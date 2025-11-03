import pypangraph as pp
from Bio import SeqIO,SeqRecord,Seq,motifs
import math
from collections import Counter
from dampa.tools.union_coverage import union_coverage
import statistics as stats
from Bio.Align import MultipleSeqAlignment,AlignInfo



def non_recursive_findcomponent(pangraph, newpaths=None, searched=None, components=None, c=0):
    if searched is None:
        searched = set()
        newpaths = []
        components = {}

    while True:
        if len(newpaths) == 0:
            allcomponents = [item for sublist in components.values() for item in sublist]
            remainingpaths = set(pangraph.paths.ids).difference(allcomponents)
            if len(remainingpaths) == 0:
                return components
            else:
                if c > 0:
                    pass  # print(f"component {c} done, {len(components[c])} paths, {len(remainingpaths)} paths remaining")
                c += 1
                searched = set()
                searchls = [remainingpaths.pop()]
                pathids = [pangraph.paths.id_to_pos[x] for x in searchls]
                ppnodedf = pangraph.nodes.df
                pathblocks = ppnodedf.loc[ppnodedf["path_id"].isin(pathids), "block_id"]
                linkedpaths = set(ppnodedf.loc[ppnodedf["block_id"].isin(
                    pathblocks), "path_id"].unique().tolist())
                pathnames = set([pangraph.paths.idx_to_name[x] for x in linkedpaths])
                newpaths = list(pathnames.difference(searched))
                continue
        else:
            searched.update(newpaths)
            if c not in components:
                components[c] = newpaths
            else:
                components[c].extend(newpaths)
            pathids = [pangraph.paths.id_to_pos[x] for x in newpaths]
            ppnodedf = pangraph.nodes.df
            pathblocks = ppnodedf.loc[ppnodedf["path_id"].isin(pathids), "block_id"]
            linkedpaths = set(ppnodedf.loc[ppnodedf["block_id"].isin(
                pathblocks), "path_id"].unique().tolist())
            pathnames = set([pangraph.paths.idx_to_name[x] for x in linkedpaths])
            newpaths = list(pathnames.difference(searched))
            continue



def find_species(components):
    species = {}
    for c in components:
        for path in components[c]:
            # s = "_".join(path.split(".")[0].split("_")[:2])
            s = path.split("_")[0]
            if c not in species:
                species[c] = set()
                species[c].add(s)
            else:
                species[c].add(s)
    return species

def shannon_entropy(seq: str) -> float:
    """Calculate Shannon entropy of a DNA sequence string."""
    counts = Counter(seq.upper())
    total = sum(counts[base] for base in "ACGT")
    if total == 0:
        return 0.0
    entropy = 0.0
    for base in "ACGT":
        p = counts[base] / total
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


def softmask_low_complexity(seqrec: SeqRecord, window: int = 120, cutoff: float = 1.6) -> Seq:
    """
    Soft-mask (lowercase) bases in low-complexity regions of a Bio.Seq object.

    seq    : Bio.Seq.Seq object
    window : sliding window size
    cutoff : entropy threshold below which sequence is considered low complexity
    """
    seq = seqrec.seq
    seq_str = str(seq).upper()
    masked = list(seq_str)  # mutable list of bases

    for i in range(0, len(seq_str) - window + 1):
        subseq = seq_str[i:i + window]
        h = shannon_entropy(subseq)
        if h < cutoff:
            for j in range(i, i + window):
                masked[j] = masked[j].lower()
    seqrec.seq = Seq.Seq("".join(masked))
    return seqrec

def rm_N_blocks(graph, nthresh):
    rmblocks = set()
    for block in graph.blocks:
        Ncount = block.consensus().count("N") + block.consensus().count("n")
        nproportion = float(Ncount) / len(block.consensus())
        if nproportion > nthresh:
            blockid = block.id
            rmblocks.add(blockid)
    return rmblocks

def filter_short_terminal_blocks(ingraph,lenthresh,rmblocks):
    for pathname,path in ingraph.paths.items():
        terminalnode = [x[1] for x in enumerate(path.nodes) if x[0] == 0 or x[0] == len(path.nodes)-1]
        for node in terminalnode:
            blockid = ingraph.nodes.df.at[str(node),"block_id"]
            if blockid not in rmblocks:
                block=ingraph.blocks[blockid]
                if len(block.consensus()) < lenthresh:
                    rmblocks.add(blockid)
                    # print(f"rm {blockid} small")

    return rmblocks

def greedy_cover(df,blocktoignore, colA="A", colB="B"):
    # Build sets of B per A
    sets = df.groupby(colA)[colB].apply(set).to_dict()
    all_B = set(df[colB].unique()) - blocktoignore

    covered = set()
    chosen = []
    strain_to_added_block = {}

    while covered != all_B:
        # pick the A with the largest number of *new* B values and return new B values in list

        best_A = max(sets, key=lambda a: len(sets[a] - covered - blocktoignore))
        new_values = sets[best_A] - covered - blocktoignore
        chosen.append(best_A)
        addset = sets[best_A] - covered - blocktoignore
        strain_to_added_block[best_A] = {
            "all_values": addset,
            "new_values": new_values
        }
        covered |= addset

        # optional: remove chosen A to speed up
        sets.pop(best_A)

    return chosen,strain_to_added_block

def get_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot unique probes and nucleotide identity across alignment positions.")
    parser.add_argument('--graphjson', type=str, required=True, help='Path to the json graph output from pangraph.')
    parser.add_argument('--nthresh', type=float, default=0.1, help='proportion of graph block consensus allowed to be N before block is ignored.')
    parser.add_argument('--lenthresh', type=int, default=120,
                        help='minimum length of terminal blocks to keep them, otherwise they are ignored.')
    parser.add_argument('--outfile', type=str, required=True, help='output fasta file of path consensus sequences covering all nodes.')
    return parser.parse_args()

def pad_short_paths(graph,usedpaths):
    """
    steps to add padding to short paths based on their position in longer paths
    1. for each path, get the block
    2. for each block, get the paths that contain it
    3. for each path that contains the block, get the length of the path
    4. if the path is shorter than some % of the longest path, add padding to the short path equivalent to the nodes covered by the longer path but missing in the shorter
    5. write out the padded paths to a new fasta file
    """
    outtargets = {s.id:s for s in usedpaths}
    pathdf = graph.nodes.df

    # usedpaths = ["LY720046.1"]
    usedpathstrim = [x.id.replace("_path_consensus","") for x in usedpaths]
    usedpathidx_to_name = {x:graph.paths.id_to_pos[x] for x in usedpathstrim}
    usedpathidx = usedpathidx_to_name.values()

    usedpathlens = {graph.paths.list[p].name: graph.paths.list[p].nuc_len for p in usedpathidx}
    # get path names that have a length < 95% of the longest path
    if len(usedpathlens) < 4:
        thresh = max(usedpathlens.values())
    else:
        thresh = stats.quantiles(usedpathlens.values(),n=4)[2]
    shortpaths = [p for p in usedpathlens if usedpathlens[p] < thresh]

    outseqs = []
    padded=[]
    for pathname in shortpaths:
        longestoverlappingpath = ("x", 0)
        pathindex = graph.paths.id_to_pos[pathname]
        shortblocks = pathdf[pathdf["path_id"] == pathindex]

        for blockid in shortblocks["block_id"]:
            blockpaths = pathdf[pathdf["block_id"] == blockid]["path_id"].tolist()
            blockpaths = [p for p in blockpaths if p in usedpathidx if p != pathindex]
            if len (blockpaths) == 0:
                continue
            pathlens = {p: graph.paths.list[p].nuc_len for p in blockpaths}
            maxlen = max(pathlens.values())
            maxlenid = max(pathlens, key=pathlens.get)

            if maxlen > longestoverlappingpath[1]:
                longestoverlappingpath = (maxlenid,maxlen)

        if longestoverlappingpath[0] == "x":
            outseqs.append(outtargets[pathname+"_path_consensus"])
            continue

        # print("newlongest",graph.paths.list[longestoverlappingpath[0]].name,longestoverlappingpath[0],longestoverlappingpath[1],pathname,pathindex,graph.paths.list[pathindex].nuc_len)
        """
        get overlapping blocks from longest path.
        """
        shortblockslist = shortblocks["block_id"].tolist()

        longpathid = longestoverlappingpath[0]
        longpathname = graph.paths.list[longestoverlappingpath[0]].name
        longpathblocks = pathdf[pathdf["path_id"] == longpathid]["block_id"].tolist()
        overlaps = [b for b in shortblockslist if b in longpathblocks]

        longpathoverlapnodes = pathdf[(pathdf["path_id"] == longpathid)&(pathdf["block_id"].isin(overlaps))]
        shortpathoverlapnodess = pathdf[(pathdf["path_id"] == pathindex)&(pathdf["block_id"].isin(overlaps))]
        """get first and last block_id in longpathoverlapnodes based on start and end positions"""
        longstartblock = longpathoverlapnodes.sort_values("start").iloc[0]
        longendblock = longpathoverlapnodes.sort_values("end").iloc[-1]
        shortstartblock = shortpathoverlapnodess.sort_values("start").iloc[0]
        shortendblock = shortpathoverlapnodess.sort_values("end").iloc[-1]
        """get consensus for each block using graph.blocks[blockid].consensus()"""
        longstartcons = graph.blocks[longstartblock["block_id"]].consensus()
        longendcons = graph.blocks[longendblock["block_id"]].consensus()
        shortstartcons = graph.blocks[shortstartblock["block_id"]].consensus()
        shortendcons = graph.blocks[shortendblock["block_id"]].consensus()
        """ get the sequences from outtargets[pathname+"_path_consensus"].seq"""
        shortseq = outtargets[pathname+"_path_consensus"].seq
        longseq = outtargets[longpathname+"_path_consensus"].seq
        """ find the position of the start and end consensus in each sequence"""
        shortoverlapstart = shortseq.find(shortstartcons)
        shortoverlapend = shortseq.find(shortendcons)
        shortlen = len(shortseq)
        longoverlapstart = longseq.find(longstartcons)
        longoverlapend = longseq.find(longendcons)
        longlen = len(longseq)
        """ calculate padding needed on left and right"""
        padleft = longoverlapstart - shortoverlapstart
        padright = (longlen - longoverlapend) - (shortlen - shortoverlapend)
        if padleft < 0:
            padleft = 0
        if padright < 0:
            padright = 0
        newseq = padleft*"N"+shortseq+padright*"N"
        newseqobj = SeqRecord.SeqRecord(Seq.Seq(newseq), id=pathname+"_path_consensus", description="")
        padded.append(pathname)
        outseqs.append(newseqobj)
        # print(f"padded {pathname} from {shortlen} to {len(newseq)} with {padleft} left and {padright} right Ns based on {longpathname} of length {longlen}")

    unchangedseqs = [outtargets[x] for x in outtargets if x.replace("_path_consensus","") not in shortpaths]
    alloutseqs = unchangedseqs + outseqs
    # SeqIO.write(alloutseqs, outpath.replace(".fa","_padded.fa").replace(".fasta","_padded.fasta"), "fasta")
    return alloutseqs,padded


def single_species_targetgen(graph,pathlist,nthresh,lenthresh,outfile,logger=None):
    pathidlist = [graph.paths.id_to_pos[x] for x in pathlist]
    filt_nodedf = graph.nodes.df[graph.nodes.df["path_id"].isin(pathidlist)]
    nodetoblock = dict(zip(filt_nodedf.index, filt_nodedf["block_id"]))
    blockno = filt_nodedf["block_id"].nunique()

    block_to_ignore = rm_N_blocks(graph, nthresh)
    # print(f"removed {len(block_to_ignore)} blocks for N proportion over {nthresh}")
    block_to_ignore = filter_short_terminal_blocks(graph, lenthresh, block_to_ignore)
    if blockno > 1:
        chosen_paths, strain_to_added_block = greedy_cover(filt_nodedf, block_to_ignore, colA="path_id", colB="block_id")
    else:
        chosen_paths = [pathidlist[0]]
    if len(chosen_paths) == 0:
        return []
    targetseqs = []
    for pathid in chosen_paths:
        pathname = graph.paths.idx_to_name[pathid]
        pathobj = graph.paths[pathname]
        outseq = ""
        blockids = []
        for node in pathobj.nodes:
            blockid = nodetoblock[str(node)]
            blockids.append(blockid)
            block = graph.blocks[blockid]
            conseq = block.consensus()
            outseq += str(conseq)


        outseqobj = SeqRecord.SeqRecord(Seq.Seq(outseq), id=f"{pathname}_path_consensus", description="")
        maskedseq = softmask_low_complexity(outseqobj, window=120, cutoff=1.6)  # TODO parameterize
        targetseqs.append(maskedseq)
        # print("masked",maskedseq.id)
    if len(targetseqs) == 0:
        return []
    if blockno > 1:
        original, final, removed, nfilt, keptidls, kepttargets = union_coverage(targetseqs, outfile, min_ident=99.0,
                                                                                min_covered_frac=0.99,
                                                                                non_n_stretch=90)  # TODO parameterize

        if len(kepttargets) == 0:
            return []
        paddedseqs, paddedids = pad_short_paths(graph, kepttargets)
    else:
        return targetseqs
    # print(f"writing {len(targetseqs)} outputs")

    return paddedseqs

def multi_species_targetgen(graph,pathlist,nthresh,lenthresh,outfile,logger=None):
    """
    For each component
        for each block
            if block is single species
                store consensus as target blockid_species
            else
                get sequences for each species
                    graph.blocks
                get consensus for each species
                save N sequences for the block if it has N species

    """

    component_nodes = graph.nodes.df[graph.nodes.df["path_name"].isin(pathlist)]
    unique_species_blocks = component_nodes.groupby('block_id')['path_species'].unique()
    block_to_species = unique_species_blocks.to_dict()
    # single_species_blocks = unique_species_blocks[unique_species_blocks.apply(len) == 1].index.tolist()
    outseq = []
    for blockid in block_to_species:
        speciesls = block_to_species[blockid]
        if len(speciesls) == 1:
            blockconsensus = graph.blocks[blockid].consensus()
            species = block_to_species[blockid]
            blockseq = SeqRecord.SeqRecord(Seq.Seq(blockconsensus),id=f"{blockid}-{species[0]}-single",name="",description="multispecies_component")
            outseq.append(blockseq)
        elif len(speciesls) > 1:
            blockalign = graph.blocks[blockid].to_biopython_alignment()
            for species in speciesls:
                speciesnodes = component_nodes[(component_nodes["block_id"]==blockid)&(component_nodes["path_species"]==species)].index.to_list()
                subset_records = [record for record in blockalign if record.id in speciesnodes]
                subset_alignment = MultipleSeqAlignment(subset_records)
                m = motifs.create([rec.seq for rec in subset_alignment])  # build Motif from the sequences
                consensus = m.consensus
                blockseq = SeqRecord.SeqRecord(Seq.Seq(consensus),id=f"{blockid}-{species}-multi",name="",description="multispecies_component")
                outseq.append(blockseq)

    outseq = [softmask_low_complexity(x, window=120, cutoff=1.6) for x in outseq]
    return outseq


def add_species_to_nodedf(graph):
    graph.nodes.df["path_name"] = graph.nodes.df["path_id"].apply(lambda x: graph.paths.idx_to_name[x])
    graph.nodes.df["path_species"] = graph.nodes.df["path_name"].apply(lambda x: str(x).split("_")[0])


def generate_targets(injson, nthresh,lenthresh,outfile,logger=None):
    graph = pp.Pangraph.from_json(injson)

    # components = recursive_findcomponent(graph)
    components = non_recursive_findcomponent(graph)
    component_species = find_species(components) # TODO separate species definition table when added to castanet
    add_species_to_nodedf(graph)

    """ 
    1 - Identify mixed species components
    2 - for single species components continue as normal
    3 - for mixed species components
        single species blocks use consensus
        multispecies blocks
        get consensus for each species
        save as blockid_species
        in castanet counts will be collapsed to blockid but initial consensus gen will generate seq for both species - allows matching
    """
    speciescount = {}
    mixedcount = 1
    outtargets = []
    for component in components:
        pathlist = components[component]
        species = list(component_species[component])
        if len(species) == 1:
            s = species[0]
            if s not in speciescount:
                speciescount[s] = 1
            else:
                speciescount[s] += 1
            paddedseqs = single_species_targetgen(graph,pathlist,nthresh,lenthresh,outfile,logger=logger)
            if not paddedseqs:
                if logger:
                    logger.warning(f"No targets generated for component {component} species {s} all nodes may have been filtered out.")
                continue
            for sequence in paddedseqs:
                sequence.id = f"{s}-component{component}-{sequence.id}"
                sequence.description = f"singlespecies_component"
            outtargets.extend(paddedseqs)
            if logger:
                logger.debug(
                    f"target number for component {component} species {s} : {len(paddedseqs)}")
        else:
            blockseqs = multi_species_targetgen(graph,pathlist,nthresh,lenthresh,outfile,logger=logger)
            outtargets.extend(blockseqs)
            if logger:
                logger.debug(
                    f"target number for component {component} species mixed: {len(blockseqs)} in {len(species)}")
            mixedcount+=1
    SeqIO.write(outtargets, outfile, "fasta")

    return outfile,len(outtargets)

def main():
    args = get_args()
    generate_targets(args.graphjson,args.nthresh,args.lenthresh,args.outfile)

if __name__ == "__main__":
    main()