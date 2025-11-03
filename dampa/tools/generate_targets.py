import pypangraph as pp
from Bio import SeqIO,SeqRecord,Seq
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

def generate_targets(injson, nthresh,lenthresh,outfile,logger=None):
    graph = pp.Pangraph.from_json(injson)
    nodetoblock = dict(zip(graph.nodes.df.index, graph.nodes.df["block_id"]))
    block_to_ignore = rm_N_blocks(graph,nthresh)
    # print(f"removed {len(block_to_ignore)} blocks for N proportion over {nthresh}")
    block_to_ignore = filter_short_terminal_blocks(graph,lenthresh,block_to_ignore,)
    chosen_paths,strain_to_added_block = greedy_cover(graph.nodes.df,block_to_ignore,colA="path_id", colB="block_id")
    targetseqs = []
    for pathid in chosen_paths:
        pathname = graph.paths.idx_to_name[pathid]
        pathobj = graph.paths[pathname]
        outseq = ""
        for node in pathobj.nodes:
            blockid = nodetoblock[str(node)]
            block = graph.blocks[blockid]
            conseq = block.consensus()
            outseq += str(conseq)

        outseqobj = SeqRecord.SeqRecord(Seq.Seq(outseq), id=pathname+"_path_consensus", description="")
        maskedseq = softmask_low_complexity(outseqobj, window=120, cutoff=1.6)#TODO parameterize
        targetseqs.append(maskedseq)
    # print(f"writing {len(targetseqs)} outputs")
    SeqIO.write(targetseqs, outfile, "fasta")

    original, final,removed, nfilt = union_coverage(outfile,outfile,min_ident=99.0,min_covered_frac=0.99,non_n_stretch=90)#TODO parameterize
    if logger:
        logger.info(f"Reduced target number: Original: {original}, Kept: {final}, Removed: {removed}, N-filtered: {nfilt}")
    return outfile,final

def main():
    args = get_args()
    generate_targets(args.graphjson,args.nthresh,args.lenthresh,args.outfile)