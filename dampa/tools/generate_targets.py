import pypangraph as pp
from Bio import SeqIO,SeqRecord,Seq,motifs
import math
from collections import Counter
from dampa.tools.union_coverage import union_coverage
import statistics as stats
import pandas as pd
import subprocess
from io import StringIO
import sys

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
                    pass
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


def softmask_low_complexity(seqrec: SeqRecord.SeqRecord, window: int = 120, cutoff: float = 1.6) -> SeqRecord.SeqRecord:
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
    rmblocks = []
    for block in graph.blocks:
        Ncount = block.consensus().count("N") + block.consensus().count("n")
        nproportion = float(Ncount) / len(block.consensus())
        conslen = len(block.consensus())
        if nproportion > nthresh or conslen < 90:
            blockid = block.id
            if blockid not in rmblocks:
                rmblocks.append(blockid)
    return rmblocks

def filter_short_terminal_blocks(ingraph,lenthresh,rmblocks):
    nodedf = ingraph.nodes.df.sort_values(["path_id", "start", "end"])

    # Rank blocks within each path
    nodedf["order"] = nodedf.groupby("path_id").cumcount()

    # Count how many blocks per path
    nodedf["n_blocks"] = nodedf.groupby("path_id")["block_id"].transform("count")

    # Assign labels
    def classify(row):
        if row["n_blocks"] == 1:
            return "first"  # or "only" if you prefer
        if row["order"] == 0:
            return "first"
        if row["order"] == row["n_blocks"] - 1:
            return "last"
        return "middle"

    nodedf["position"] = nodedf.apply(classify, axis=1)
    pos_per_block = nodedf.groupby("block_id")["position"].unique()

    # Keep block_ids that do NOT include "middle"
    only_first_or_last = pos_per_block[pos_per_block.apply(lambda x: "middle" not in x)].index
    for block in ingraph.blocks:
        if block.id in only_first_or_last:
            if block.id not in rmblocks:
                if len(block.consensus()) < lenthresh:
                    rmblocks.append(block.id)

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

def get_start_end_padding(graph,pathdf,seqs,blocktoignore):
    #TODO remove blocks that are ignored from pathdf before running
    """chosen_subclustpath"""
    outtargets = {s.id: s for s in seqs}
    #for each unique value in pathdf["chosen_subclustpath"] get pathdf["path_id"], make unique list of path ids
    sub_to_path = pathdf[pathdf["union_subclustpath"]].groupby('subcluster_path_id')['path_name'].first().to_dict()
    pathname_to_id = pathdf.groupby('path_name')['path_id'].first().to_dict()
    id_to_pathname = {v: k for k, v in pathname_to_id.items()}
    path_to_sub = {v: k for k, v in sub_to_path.items()}
    correspondingpaths = list(sub_to_path.values())
    usedpathlens = {graph.paths.list[pathname_to_id[p]].name: graph.paths.list[pathname_to_id[p]].nuc_len for p in correspondingpaths}
    strblocktoignore = list(map(str,blocktoignore))
    if len(usedpathlens) < 4:
        thresh = max(usedpathlens.values())
    else:
        thresh = stats.quantiles(usedpathlens.values(),n=4)[2]
    shortpaths = [p for p in usedpathlens if usedpathlens[p] < thresh]
    longpaths = [p for p in usedpathlens if usedpathlens[p] >= thresh]

    padding = {}
    for longpath in longpaths:
        padding[longpath] = (0,0)
    for pathname in shortpaths:
        longestoverlappingpath = ("x", 0)
        pathindex = graph.paths.id_to_pos[pathname]
        shortblocks = pathdf[pathdf["path_id"] == pathindex]

        for blockid in shortblocks["block_id"]:
            if blockid not in blocktoignore:
                blockpaths = pathdf[pathdf["block_id"].astype(str) == str(blockid)]["path_id"].tolist()
                # blockpaths = pathdf[pathdf["block_id"].astype(str) == str(blockid)]["path_name"].tolist()
                ...
                blockpaths = [p for p in blockpaths if id_to_pathname[p] in correspondingpaths if p != pathindex]
                if len (blockpaths) == 0:
                    continue
                pathlens = {p: graph.paths.list[p].nuc_len for p in blockpaths}
                maxlen = max(pathlens.values())
                maxlenid = max(pathlens, key=pathlens.get)

                if maxlen > longestoverlappingpath[1]:
                    longestoverlappingpath = (maxlenid,maxlen)

        if longestoverlappingpath[0] == "x":
            padding[pathname] = (0,0)
            continue

        """
        get overlapping blocks from longest path.
        """
        shortblockslist = shortblocks["block_id"].tolist()
        shortblockslist = [b for b in shortblockslist if str(b) not in strblocktoignore]
        longpathid = longestoverlappingpath[0]
        longpathname = graph.paths.list[longestoverlappingpath[0]].name
        longpathblocks = pathdf[pathdf["path_id"] == longpathid]["block_id"].tolist()
        longpathblocks = [b for b in longpathblocks if str(b) not in strblocktoignore]
        overlaps = [b for b in shortblockslist if b in longpathblocks]

        longpathoverlapnodes = pathdf[(pathdf["path_id"] == longpathid)&(pathdf["block_id"].isin(overlaps))]
        shortpathoverlapnodess = pathdf[(pathdf["path_id"] == pathindex)&(pathdf["block_id"].isin(overlaps))]
        """get first and last block_id in longpathoverlapnodes based on start and end positions"""
        longstartblock = longpathoverlapnodes.sort_values("start").iloc[0]
        longendblock = longpathoverlapnodes.sort_values("end").iloc[-1]
        shortstartblock = shortpathoverlapnodess.sort_values("start").iloc[0]
        shortendblock = shortpathoverlapnodess.sort_values("end").iloc[-1]
        """get consensus for each block using graph.blocks[blockid].consensus()"""
        try:
            longstartcons = graph.blocks[int(longstartblock["block_id"])].consensus()
        except:
            print("error getting longstartcons for block",longstartblock["block_id"])

            continue
        longendcons = graph.blocks[int(longendblock["block_id"])].consensus()
        shortstartcons = graph.blocks[int(shortstartblock["block_id"])].consensus()
        shortendcons = graph.blocks[int(shortendblock["block_id"])].consensus()
        """ get the sequences from outtargets[pathname+"_path_consensus"].seq"""
        shortpathseqname = path_to_sub[pathname]
        longpathseqname = path_to_sub[longpathname]
        shortseq = outtargets[shortpathseqname].seq
        longseq = outtargets[longpathseqname].seq
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
        padding[pathname] = (padleft,padright)

    paddedseqs = []
    for seq in seqs:
        pathname = sub_to_path[seq.id]
        padleft,padright = padding[pathname]
        newseq = padleft*"N"+str(seq.seq)+padright*"N"
        newseqobj = SeqRecord.SeqRecord(Seq.Seq(newseq), id=seq.id, description=seq.description)
        paddedseqs.append(newseqobj)

    return paddedseqs

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
    return alloutseqs,padded

"""
for each selected path get set of other paths that share the same blocks
run clustering on these paths and produce representative sequences for each cluster
"""

def assign_block_subcluster(component_nodes,component_blocks,block_to_ignore,graph,outfile,subnode_cluster_thresh):
    nodetosubclust = {}
    c = 0
    node_to_subcluster_rep = {}
    node_to_repid = {}
    for blockid in component_blocks:
        if blockid not in block_to_ignore:
            blockalign = graph.blocks[blockid].to_biopython_alignment()
            blockseqs = [
                SeqRecord.SeqRecord(Seq.Seq(str(record.seq).replace("-", "")), id=record.id, description=record.description)
                for record in blockalign
            ]

            reps,clusters = runvsearch_cluster(blockseqs, outfile, cutoff=subnode_cluster_thresh)


            for cno,cluster in enumerate(clusters):
                blockclustername = f"{blockid}_subcluster_{cno+1}"
                for node in cluster:
                    nodetosubclust[node.id] = blockclustername
                    node_to_subcluster_rep[node.id] = reps[cno]
                    node_to_repid[node.id] = reps[cno].id

            sys.stderr.write(f"block {blockid} had {len(reps)} subclusters from {len(blockseqs)} sequences\n")
            c+=1
            sys.stderr.write(f"completed {c} of {len(component_blocks)} blocks\n")
    component_nodes["subcluster_id"] = component_nodes.index.map(nodetosubclust)
    component_nodes["subcluster_rep"] = component_nodes.index.map(node_to_repid)
    return component_nodes,node_to_subcluster_rep

def assign_final_clusters(df):
    """
    Given a dataframe with columns:
        block_id (string or convertible to string),
        path_id,
        subcluster_id,
        path_group,
        ignoreblock (boolean)

    - Keeps block_id as string
    - Ignores rows where ignoreblock == True
    - Only includes block_ids actually present for that path_id
    - Produces a path-specific final cluster string
    """

    df = df.copy()
    originalindex = df.index

    # Ensure block_id is always a string (kept as string permanently)
    df["block_id"] = df["block_id"].astype(str)



    # Pivot: creates columns named by block_id (as strings)
    pivot = (
        df.pivot_table(
            index=["path_group", "path_id"],
            columns="block_id",               # stays string
            values="subcluster_id",
            aggfunc="first"
        )
        .reset_index()
    )

    # Block columns are the pivoted string block_id values
    block_cols = [
        col for col in pivot.columns
        if col not in ["path_group", "path_id"]
    ]

    # Sort string block_ids numerically but WITHOUT converting pivot columns
    block_cols_sorted = sorted(block_cols, key=lambda x: int(x))

    # Build signature using only blocks that exist for each path
    def build_signature(row):
        parts = []
        for b in block_cols_sorted:
            val = row[b]
            if not pd.isna(val):         # include only blocks this path actually has
                parts.append(str(val))
        return "_".join(parts) if parts else "no_blocks"

    pivot["final_cluster"] = pivot.apply(build_signature, axis=1)

    # Merge final clusters back to original df
    df_out = df.merge(
        pivot[["path_group", "path_id", "final_cluster"]],
        on=["path_group", "path_id"],
        how="left"
    )
    df_out["node_id"] = originalindex

    return df_out



def single_species_targetgen_node_subclust(graph,pathlist,nthresh,lenthresh,outfile,subnode_cluster_thresh,target_unioncov,logger=None):
    block_to_ignore = rm_N_blocks(graph, nthresh)
    block_to_ignore = list(filter_short_terminal_blocks(graph, lenthresh, block_to_ignore))
    ## before doing subnode clustering etc get set of paths that cover all non ignored blocks like in old function, then can remove all other paths from component_nodes
    nodes_to_ignore = graph.nodes.df[graph.nodes.df["block_id"].isin(block_to_ignore)].index.tolist()
    if "15205707089436575906" in block_to_ignore:
        ...

    component_nodes = graph.nodes.df[graph.nodes.df["path_name"].isin(pathlist)]

    component_nodes["ignoreblock"] = component_nodes['block_id'].isin(block_to_ignore)
    component_nodes = component_nodes[component_nodes["ignoreblock"] == False]
    component_blocks = component_nodes['block_id'].unique()
    if len(component_blocks) == 0: # if filtering removed all blocks try without filtering
        component_nodes = graph.nodes.df[graph.nodes.df["path_name"].isin(pathlist)]
        component_blocks = component_nodes['block_id'].unique()
        block_to_ignore = [x for x in block_to_ignore if x not in component_blocks]
    try:
        componentspecies = component_nodes['path_species'].unique()[0]
    except:
        component_nodes2 = graph.nodes.df[graph.nodes.df["path_name"].isin(pathlist)]
        ...
    """
    todo
    save sub block consensus sequences for each block and which paths are in those subconsensuses
    for each path type check each node to see if the path type includes multiple subconsensuses
    if they do split paths into groups of paths based on subconsensus assigment
    then for each subsequent block repeat the process until all blocks in the path type are covered.
    
    data structure required:
    df of
    block_id | path_id | subconsensus_id | subgroup
    
    then get list of blocks that make up a path type
    for each block if there is more than one subconsensus for that block in the path type
        split path type into subgroups based on subconsensus id
    """
    ...
    # Use only non-ignored blocks for cluster generation

    # for each block get subclusters at subnode_cluster_thresh identity and add as column to component_nodes
    component_nodes,node_to_subcluster_rep = assign_block_subcluster(component_nodes,component_blocks,block_to_ignore,graph,outfile,subnode_cluster_thresh)
    # for groups of paths that share same set of blocks assign path group ids
    component_nodes = assign_path_group(component_nodes)
    """for each block in a path group split the path group by subcluster_id in that block then in next block split the already split path groups again by subcluster_id in that block"""
    component_nodes = assign_final_clusters(component_nodes)

    unique_vals = component_nodes["final_cluster"].unique()
    mapping = {v: f"subcluster_path_id_{i + 1}" for i, v in enumerate(unique_vals)}
    if logger:
        logger.info(f"identified {len(unique_vals)} unique subcluster paths in species {componentspecies} covering component")
    else:
        print(f"identified {len(unique_vals)} unique subcluster paths in species {componentspecies} covering component")
    component_nodes["subcluster_path_id"] = component_nodes["final_cluster"].map(mapping)

    blockno = component_nodes["block_id"].nunique()
    if blockno > 1:
        blockset = set(block_to_ignore)
        chosen_paths, strain_to_added_block = greedy_cover(component_nodes, blockset, colA="subcluster_path_id",
                                                           colB="subcluster_id")
    else:
        chosen_paths = [component_nodes["subcluster_path_id"].tolist()[0]]
    if logger:
        logger.info(f"after greedy cover {len(chosen_paths)} subcluster paths chosen to cover all nodes in component for species {componentspecies}")
    else:
        print(f"after greedy cover {len(chosen_paths)} subcluster paths chosen to cover all nodes in component for species {componentspecies}")

    component_nodes["chosen_subclustpath"] = component_nodes["subcluster_path_id"].isin(chosen_paths)

    # get order of blocks for each path group from grapj.paths["path_name"].nodes then getting corresponding block ids from graph.nodes.df

    subpg_seqs = []
    subpgno = 1
    for subpg in chosen_paths:
        pathgroup_paths = component_nodes[component_nodes["subcluster_path_id"] == subpg]["path_id"].unique()
        example_pathid = pathgroup_paths[0]
        example_pathname = graph.paths.idx_to_name[example_pathid]
        example_pathobj = graph.paths[example_pathname]
        exampleorder = example_pathobj.nodes
        blockorder = [graph.nodes[x].block_id for x in exampleorder]
        blockorder = [x for x in blockorder if x not in block_to_ignore]
        outseq = ""
        for block in blockorder:
            rep1 = component_nodes[((component_nodes["subcluster_path_id"] == subpg)&(component_nodes["block_id"].astype(str) == str(block)))]
            rep = rep1["subcluster_rep"].unique().tolist()[0]
            if str(rep) not in node_to_subcluster_rep and logger:
                logger.warning(f"missing block to subcluster - possible vsearch issue or very short node, {rep}, {block}, {subpg}")
            elif str(rep) not in node_to_subcluster_rep and not logger:
                print("missing block to subcluster - possible vsearch issue or very short node",rep,block,subpg)
            else:
                # outseq = outseq + rep
                subclusterrep = node_to_subcluster_rep[str(rep)]
                outseq += str(subclusterrep.seq)
        outseqobj = SeqRecord.SeqRecord(Seq.Seq(outseq), id=subpg, description="singlespecies_component_subpg")
        maskedseq = softmask_low_complexity(outseqobj, window=120, cutoff=1.6)  # TODO parameterize
        subpg_seqs.append(maskedseq)
        subpgno += 1

    if blockno > 1:
        if target_unioncov:
            original, final, removed, nfilt, keptidls, kepttargets = union_coverage(subpg_seqs, outfile, min_ident=99.0,
                                                                                    min_covered_frac=0.99,
                                                                                    non_n_stretch=90)  # TODO parameterize
            if logger:
                logger.info(f"after union coverage{len(kepttargets)}sequences kept from{len(subpg_seqs)}")
            else:
                print("after union coverage",len(kepttargets),"sequences kept from",len(subpg_seqs))
            component_nodes["union_subclustpath"] = component_nodes["subcluster_path_id"].isin(list(keptidls))

            if len(kepttargets) == 0:
                return []
            paddedseqs = get_start_end_padding(graph,component_nodes,kepttargets,block_to_ignore)
        else:
            component_nodes["union_subclustpath"] =component_nodes["chosen_subclustpath"]
            paddedseqs = get_start_end_padding(graph,component_nodes,subpg_seqs,block_to_ignore)
        return paddedseqs
    else:
        return subpg_seqs

def single_species_targetgen(graph,pathlist,nthresh,lenthresh,outfile,subnode_cluster_thresh,logger=None):
    intact_pathtargets = True

    if intact_pathtargets:
        pathidlist = [graph.paths.id_to_pos[x] for x in pathlist]
        filt_nodedf = graph.nodes.df[graph.nodes.df["path_id"].isin(pathidlist)]

        # shared_paths = get_same_block_paths(filt_nodedf)
        # for i in shared_paths.values():
        #     print(len(i),i)
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
            return paddedseqs
        else:
            return targetseqs

    else:
        component_nodes = graph.nodes.df[graph.nodes.df["path_name"].isin(pathlist)]
        component_blocks = component_nodes['block_id'].unique()
        componentspecies = component_nodes['path_species'].unique()[0]
        outseq = []
        for blockid in component_blocks:
            blockalign = graph.blocks[blockid].to_biopython_alignment()
            blockseqs = [
                SeqRecord.SeqRecord(Seq.Seq(str(record.seq).replace("-", "")), id=record.id, description=record.description)
                for record in blockalign
            ]
            reps = runvsearch_consout(blockseqs, outfile, cutoff=subnode_cluster_thresh)
            if len(reps) == 1:
                blockconsensus = graph.blocks[blockid].consensus()
                blockseq = SeqRecord.SeqRecord(Seq.Seq(blockconsensus), id=f"{blockid}_{componentspecies}-single_consensus", name="",
                                               description="GRAPH_singlespecies_component")
                outseq.append(blockseq)
            else:
                r = 1
                for rep in reps:
                    rep.id = f"{blockid}_{componentspecies}-single_rep_{r}"
                    rep.description = "GRAPH_singlespecies_component"
                    outseq.append(rep)
                    r += 1
        outseq = [softmask_low_complexity(x, window=120, cutoff=1.6) for x in outseq]
        return outseq



def all_sequences_identical(alignment):
    """Return True if all sequences in a Biopython alignment are identical."""
    seqs = [str(rec.seq) for rec in alignment]
    first = seqs[0]
    return all(seq == first for seq in seqs[1:])

def runvsearch_consout(seqs,outfile,cutoff=0.95,):

    if all_sequences_identical(seqs):
        # All identical → single representative
        # print("All sequences are identical; returning one representative.")
        return [seqs[0]]

    tmpprefix = outfile.replace(".fasta","_vsearch_tmp.fasta")
    SeqIO.write(seqs, tmpprefix, "fasta")
    tmpreps = outfile.replace(".fasta","_vsearch_reps.fasta")
    try:
        result = subprocess.run([
            "vsearch",
            "--cluster_fast", "-",  # read input from stdin
            "--id", str(cutoff),
            "--consout", "-",  # write representatives to stdout
            "--threads", "4"
        ], input=open(tmpprefix).read(), text=True, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr, file=sys.stderr)
        raise
    # Parse FASTA from stdout
    reps = list(SeqIO.parse(StringIO(result.stdout), "fasta"))
    # print(f"{len(reps)} representative sequences")
    return reps


def runvsearch_cluster(seqs, outfile, cutoff=0.95):
    """
    Cluster sequences using vsearch, returning both the representatives (centroids)
    and the list of sequences in each cluster.

    Returns:
        reps: list of Bio.SeqRecord (cluster representatives)
        clusters: list of lists of Bio.SeqRecord (members per cluster)
    """
    if all_sequences_identical(seqs):
        # All identical → single representative
        return [seqs[0]], [seqs]

    tmpprefix = outfile.rsplit(".",1)[0] + "_vsearch_tmp.fasta"#.replace(".fasta", "_vsearch_tmp.fasta")
    tmpreps = outfile.rsplit(".",1)[0] + "_vsearch_reps.fasta"#.replace(".fasta", "_vsearch_reps.fasta")
    tmpuc = outfile.rsplit(".",1)[0] + "_vsearch.uc"#.replace(".fasta", "_vsearch.uc")

    SeqIO.write(seqs, tmpprefix, "fasta")

    try:
        subprocess.run([
            "vsearch",
            "--cluster_fast", tmpprefix,
            "--id", str(cutoff),
            "--centroids", tmpreps,
            "--uc", tmpuc,
            "--threads", "4",
            "--minwordmatches", "10"
        ], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print("VSEARCH failed:", e.stderr, file=sys.stderr)
        raise

    # Parse representatives (centroids)
    reps = [x for x in SeqIO.parse(tmpreps, "fasta")]

    # Parse UC file to get cluster membership
    cluster_members = {}
    seq_dict = {rec.id: rec for rec in seqs}

    with open(tmpuc) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            record_type = parts[0]  # 'S' = seed, 'H' = hit
            cluster_id = parts[1]
            seq_id = parts[8]

            cluster_members.setdefault(cluster_id, []).append(seq_dict.get(seq_id))

    # Sort clusters to match centroid order if desired
    clusters = [cluster_members[cid] for cid in sorted(cluster_members, key=int)]

    return reps, clusters

def assign_path_group(nodedf):
    # Map each path_id to the immutable set of block_ids it covers
    path_to_blocks = nodedf.groupby('path_id')['block_id'].apply(lambda x: frozenset(x))

    # Assign a unique integer id to each distinct frozenset of block IDs
    unique_sets = {s: idx for idx, s in enumerate(path_to_blocks.unique(), start=1)}

    # Map each path_id -> path_group id
    path_group_series = path_to_blocks.map(unique_sets)

    # Add column to original dataframe (map by path_id for each row)
    nodedf['path_group'] = nodedf['path_id'].map(path_group_series.to_dict())
    return nodedf

def multi_species_targetgen(graph,pathlist,nthresh,lenthresh,outfile,subnode_cluster_thresh,component,graphid,logger=None):
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

            species = block_to_species[blockid]

            blockalign = graph.blocks[blockid].to_biopython_alignment()
            blockseqs = [
                SeqRecord.SeqRecord(Seq.Seq(str(record.seq).replace("-", "")), id=record.id,
                                    description=record.description)
                for record in blockalign
            ]
            reps = runvsearch_consout(blockseqs, outfile, cutoff=subnode_cluster_thresh)
            if len(reps) == 1:
                blockconsensus = graph.blocks[blockid].consensus()
                sid = f"graph{graphid}block{blockid}_{species[0]}-single_consensus"
                blockseq = SeqRecord.SeqRecord(Seq.Seq(blockconsensus), id=sid, name="",
                                               description=f"GRAPH_multispecies_component_{component}")
                outseq.append(blockseq)
            else:
                r = 1
                for rep in reps:
                    rep.id = f"graph{graphid}block{blockid}_{species[0]}-single_rep_{r}"

                    rep.description = f"GRAPH_multispecies_component_{component}"
                    outseq.append(rep)
                    r += 1
            ...
        elif len(speciesls) > 1:

            species = block_to_species[blockid]
            blockalign = graph.blocks[blockid].to_biopython_alignment()
            blockseqs = [
                SeqRecord.SeqRecord(Seq.Seq(str(record.seq).replace("-", "")), id=record.id,
                                    description=record.description)
                for record in blockalign
            ]
            for species in speciesls:
                speciesnodes = component_nodes[(component_nodes["block_id"]==blockid)&(component_nodes["path_species"]==species)].index.to_list()
                subset_records = [record for record in blockseqs if record.id in speciesnodes]
                reps = runvsearch_consout(subset_records, outfile, cutoff=subnode_cluster_thresh)
                r=1
                for rep in reps:
                    rep.id = f"graph{graphid}block{blockid}_{species}-multi_rep_{r}"
                    rep.description = f"GRAPH_multispecies_component_{component}"
                    outseq.append(rep)
                    r+=1

    outseq = [softmask_low_complexity(x, window=120, cutoff=1.6) for x in outseq]
    return outseq


def add_species_to_nodedf(graph):
    graph.nodes.df["path_name"] = graph.nodes.df["path_id"].apply(lambda x: graph.paths.idx_to_name[x])
    graph.nodes.df["path_species"] = graph.nodes.df["path_name"].apply(lambda x: str(x).split("_")[0])


def generate_targets(injson, nthresh,lenthresh,outfile,logger=None):
    graph = pp.Pangraph.from_json(injson)
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
            paddedseqs = single_species_targetgen_node_subclust(graph,pathlist,nthresh,lenthresh,outfile,subnode_cluster_thresh,target_unioncov,logger=logger)
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
            blockseqs = multi_species_targetgen(graph, pathlist, nthresh, lenthresh, outfile, subnode_cluster_thresh,
                                                component, graphid, logger=logger)
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