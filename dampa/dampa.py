import random
import shutil
from pathlib import Path
from Bio import SeqIO,SeqRecord,Seq
from collections import defaultdict
import subprocess
import argparse
import os
import logging
import time
import sys
from collections import Counter
import platform
import subprocess
import glob
import re
import copy
import importlib.resources
# from tools.gather_probe_depth_stats import make_stats, make_propsplot, process_count_pt
from dampa.vis.plot_over_genomelen import make_genome_plots,replace_short_zeros
from dampa.tools.gather_probe_depth_stats import make_stats, make_propsplot, process_count_pt
from dampa.tools.cdhit_rep_noN import longest_or_fewest_ns_representatives
from dampa import __version__ as dampaversion


def mmseqs_subset(args,filtinput):
    """
    run the following commands
    mmseqs createdb ingenomes.fasta alltypesdb

    mmseqs cluster alltypesdb $pref""DB tmp --min-seq-id clusterident -c clustercov

    mmseqs createtsv alltypesdb $pref""DB $pref""cluster.tsv

    mmseqs createsubdb $pref""DB alltypesdb $pref""Dbreps

    mmseqs convert2fasta $pref""DBreps $pref""DB.fasta

    """
    outloc = f"{args.outputfolder}/{args.outputprefix}"
    mmseqs_log = open(outloc + "_mmseqs.log", "w")

    cmd = f"mmseqs createdb {filtinput} {outloc}_alltypesdb"
    subprocess.run(cmd, shell=True, stdout=mmseqs_log, stderr=mmseqs_log)
    cmd = f"mmseqs cluster {outloc}_alltypesdb {outloc}_DB tmp --min-seq-id {args.clusterident} -c {args.clustercov} --cov-mode 1"
    subprocess.run(cmd, shell=True, stdout=mmseqs_log, stderr=mmseqs_log)
    # cmd = f"mmseqs createtsv {outloc}_alltypesdb {outloc}_DB {outloc}_cluster.tsv"
    # subprocess.run(cmd, shell=True, stdout=mmseqs_log, stderr=mmseqs_log)
    cmd = f"mmseqs createsubdb  {outloc}_DB {outloc}_alltypesdb {outloc}_Dbreps"
    subprocess.run(cmd, shell=True, stdout=mmseqs_log, stderr=mmseqs_log)
    mmseqsreps = f"{outloc}_reps.fasta"
    cmd = f"mmseqs convert2fasta {outloc}_Dbreps {mmseqsreps}"
    subprocess.run(cmd, shell=True, stdout=mmseqs_log, stderr=mmseqs_log)

    if os.path.exists(mmseqsreps):
        logger.info(f"Mmseqs ran successfully")
        if not args.keeplogs:
            os.remove(outloc + "_mmseqs.log")
        return mmseqsreps
    else:
        logger.error(f"mmseqs output file {mmseqsreps} not present. Check for error in capture log.")


    return mmseqsreps


def cdhit_subset(args,filtinput):
    inseqs = SeqIO.to_dict(SeqIO.parse(filtinput, "fasta"))
    outloc = f"{args.outputfolder}/{args.outputprefix}"
    cdhit_log = open(outloc + "_cdhit.log", "w")

    cmd = f"cd-hit-est -M 0 -c {args.clusterident} -T {args.threads} -aS {args.clustercov} -i {filtinput} -o {outloc}_cdhit_reps -d 0"
    subprocess.run(cmd, shell=True, stdout=cdhit_log, stderr=cdhit_log)
    cdhitrepsout = f"{outloc}_cdhit_reps.fasta"
    clusters,cdhitreps = longest_or_fewest_ns_representatives(inseqs,f"{outloc}_cdhit_reps.clstr")
    # Write results
    with open(cdhitrepsout, "w") as out_f:
        SeqIO.write(cdhitreps, out_f, "fasta")


    if os.path.exists(cdhitrepsout):
        logger.info(f"cdhit ran successfully")
        if not args.keeplogs:
            os.remove(outloc + "_cdhit.log")
        return cdhitrepsout
    else:
        logger.error(f"cdhit output file {cdhitrepsout} not present. Check for error in capture log.")

    return cdhitrepsout


def cdhit_2step_clust(args,filtinput,filtered):
    """
    Runs a two-step clustering process using CD-HIT to generate representative sequences.

    Args:
        args (argparse.Namespace): The arguments passed to the script, containing various settings and file paths.
        filtinput (str): Path to the input filtered genomes file.

    Returns:
        str: Path to the output file containing representative sequences.
    """
    inseqs = SeqIO.to_dict(SeqIO.parse(filtinput, "fasta"))
    outloc = f"{args.outputfolder}/{args.outputprefix}"
    cdhit_log = open(outloc + "_cdhit.log", "w")
    cdhithighreps = f"{outloc}_cdhit_reps_high"
    cdhithighrepsout = f"{cdhithighreps}.fasta"
    cmd = f"cd-hit-est -M 0 -c {args.clusterident} -T {args.threads} -aS {args.clustercov} -i {filtinput} -o {cdhithighreps} -d 0"
    subprocess.run(cmd, shell=True, stdout=cdhit_log, stderr=cdhit_log)
    highclusters,highcdhitreps = longest_or_fewest_ns_representatives(inseqs,f"{cdhithighreps}.clstr")
    strain_to_highcluster = {strain: cluster for cluster, strains in highclusters.items() for strain in strains}
    strain_to_highcluster_size = {strain: len(strains) for strains in highclusters.values() for strain in strains}
    strain_to_highcluster_size
    logger.info(f"High cutoff clustering complete. {len(highcdhitreps)} representative sequences generated.")

    cdhitlowreps = f"{outloc}_cdhit_reps_low"
    cdhitlowrepsout = f"{cdhitlowreps}.fasta"
    cmd = f"cd-hit-est -M 0 -c {args.lowclusterident} -T {args.threads} -aS {args.clustercov} -i {cdhithighreps} -o {cdhitlowreps} -d 0"
    subprocess.run(cmd, shell=True, stdout=cdhit_log, stderr=cdhit_log)
    lowclusters,lowcdhitreps = longest_or_fewest_ns_representatives(inseqs, f"{cdhitlowreps}.clstr")

    logger.info(f"Low cutoff clustering complete. {len(lowcdhitreps)} representative sequences generated.")

    regular,outliers = split_to_reg_or_outliers(lowclusters, strain_to_highcluster_size,highcdhitreps,inseqs,outlier_size_limit=args.outliersizelimit)

    for outlier in outliers:
        seqlen = len(outlier.seq)
        newrow = {"genome id": outlier.id, "genome description": outlier.description, "filter reason": "outlier singleton genomes", "Genome length": seqlen,
                  "nonstandard proportion": ""}
        new_df = pd.DataFrame([newrow])
        filtered = pd.concat([filtered, new_df], ignore_index=True)

    if len(regular) > 1:
        regstrains = f"{outloc}_cdhit_reg.fasta"
        SeqIO.write(regular, regstrains, "fasta")
    elif len(regular) == 1:
        regstrains = f"{outloc}_cdhit_reg.fasta"
        dupe = copy.copy(regular[0])
        dupe.id = dupe.id+"_dupe"
        regular = [regular[0],dupe]
        SeqIO.write(regular, regstrains, "fasta")
    else:
        logger.error(f"\n#########\nNo regular strains found in high cutoff clusters.\n#########\n")
        sys.exit()
    if len(outliers) > 1:
        outliersstrains = f"{outloc}_cdhit_outlier.fasta"
        SeqIO.write(outliers, outliersstrains, "fasta")
    elif len(outliers) == 1:
        outliersstrains = f"{outloc}_cdhit_outlier.fasta"
        dupe = copy.copy(outliers[0])
        dupe.id = dupe.id+"_dupe"
        outliers = [outliers[0],dupe]
        SeqIO.write(outliers, outliersstrains, "fasta")
    else:
        logger.info(f"No outlier strains found in high cutoff clusters. No outlier fasta file or probes will be generated.")
        outliersstrains = False


    # TODO for singleton lowclusters check if also singletons in high clusters - if so add to alternate rep fasta file
    # use alternate rep fasta file to generate probes in separate dampa run

    cdhit_log.close()

    if os.path.exists(regstrains):
        logger.info(f"cdhit ran successfully")
        if not args.keeplogs:
            os.remove(outloc + "_cdhit.log")
        return regstrains,outliersstrains,filtered
    else:
        logger.error(f"cdhit output file {regstrains} not present. Check for error in capture log.")

def split_to_reg_or_outliers(lowclusters, strain_to_highcluster_size,highcdhitreps,inseqs,outlier_size_limit):
    altstrains = []
    regstrains = []
    altclustscount = 0
    for cluster in lowclusters:
        lowclustersize = len(lowclusters[cluster])
        if lowclustersize <= outlier_size_limit:
            # if only x sequence in low cluster, check if also x seq in high cluster
            strain = lowclusters[cluster][0]
            highclustersize = strain_to_highcluster_size[strain]
            # print(strain,highclustersize)
            if highclustersize == lowclustersize:
                # print(f"low cluster {cluster} with {lowclustersize} sequences is also a high cluster with {highclustersize} sequences. This is an outlier cluster.")
                altclustscount += 1
                for i in lowclusters[cluster]:
                    altstrains.append(inseqs[i])
                    # print(inseqs[i])

    regstraintotal = 0
    altstrainids = [altstrain.id for altstrain in altstrains]
    for strain in highcdhitreps:
        if strain.id not in altstrainids:
            regstrains.append(strain)
            if strain.id not in strain_to_highcluster_size:
                logger.error(f"strain {strain.id} not found in strain_to_highcluster_size dictionary. This should not happen.")
                for strain2 in strain_to_highcluster_size:
                    print(strain2,strain_to_highcluster_size[strain2])
            regstraintotal += strain_to_highcluster_size[strain.id]

    prop_alt_clust= float(altclustscount) / float(len(highcdhitreps))
    prop_alt_strains = float(len(altstrains)) / float(len(inseqs.keys()))
    if prop_alt_strains == 1.0:
        logger.info(f"All strains are in low cutoff clusters. No outliers will be generated.")
        regstrains = copy.deepcopy(altstrains)
        altstrains = []
    elif prop_alt_clust*100 >= 0.1:
        logger.info(f"Over 10% of high cutoff clusters ({prop_alt_clust*100:.2f}%) are also singletons at low cutoff. Outliers will not be generated.")
        regstrains = copy.deepcopy(highcdhitreps)
        altstrains = []
    if altclustscount > 0:
        logger.warning(f"{prop_alt_clust*100:.2f}% ({altclustscount}/{len(highcdhitreps)}) of high cutoff clusters are also outliers at low cutoff.\n"
                    f"{prop_alt_strains*100:.2f}% ({len(altstrains)}/{len(inseqs.keys())}) of all strains are in these clusters.\n"
                    f"These may be incorrect species or artificially modified genomes (i.e codon optimised).\n"
                    f"Please check excluded strain file: _filtered_genomes.tsv for details.\n"
                    f"These strains will be used to separately generate probes with the _outliers marker.")
    return regstrains,altstrains

def make_padded_probes(pangenomefa,probefasta,minlen,probeprefix=""):
    """
    Generates padded probes for sequences in the pangenome that are shorter than a specified length.

    Args:
        pangenomefa (str): Path to the input pangenome FASTA file.
        probefasta (str): Path to the output probes FASTA file.
        minlen (int): Minimum length of sequences to be padded.

    Returns:
        None
    """
    inpangenome = SeqIO.parse(pangenomefa, 'fasta')
    toadd = []
    for pancontig in inpangenome:
        if len(pancontig.seq) < 120:
            rmn = str(pancontig.seq).replace("N","").replace("n","")
            if len(rmn) >= minlen:
                tadd = "T"*(120-len(rmn))
                newseq = Seq.Seq(tadd + rmn)
                probeprefixm = str(probeprefix+"_")
                outpancontig = SeqRecord.SeqRecord(newseq,probeprefixm+pancontig.id+"_padded_probe",description="")
                toadd.append(outpancontig)
    with open(probefasta, "a") as fasta_file:
        SeqIO.write(toadd, fasta_file, "fasta")
    logger.info(f"Padding small sequences generated {len(toadd)} additional probes")

def nucleotide_proportions(sequences,propinit):
    """
    Calculates the nucleotide proportions in a list of sequences.

    Args:
        sequences (dict): Dict of SeqRecord objects.

    Returns:
        dict: Dictionary with nucleotide proportions.
    """
    overallseq = ""

    for i in sequences:
        try:
            sequence = sequences[i]
        except Exception as e:
            logger.error(e)
            print(i)
            print(sequences)
        overallseq += str(sequence.seq).upper()
    counts = Counter(overallseq)  # Count occurrences of each nucleotide
    total = sum(counts.values())  # Total nucleotides
    proportions = {nt: count / total for nt, count in counts.items()}# Calculate proportions
    for x in propinit:
        if x not in proportions:
            proportions[x] = 0.0
    return proportions

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

def filter_for_nonstandard_inputs(genomes,outloc,maxnonspandard,filtered):
    """
    Filters out genomes with a high proportion of non-standard nucleotides.
    Also trims genomes for trailing As or Ns, and removes genomes with excess Ns or too short sequences.
    Args:
        genomes (str): Path to the input genomes FASTA file.
        outloc (str): output path and prefix.
        maxnonspandard (float): Maximum allowed proportion of non-standard nucleotides.

    Returns:
        tuple: Path to the filtered genomes file and overall nucleotide proportions.
    """
    ingenomes = SeqIO.parse(genomes, "fasta")
    ingenomes = to_dict_remove_dups(ingenomes)
    outgenomes = []
    excluded = 0
    included = 0
    allowed = ["N", "A", "C", "G", "T", "n", "a", "c", "g", "t"]
    propinit = {x:0 for x in allowed}
    overallprops = nucleotide_proportions(ingenomes,propinit)
    added = []
    descriptions = {}
    notrimmed = 0
    for id,i in ingenomes.items():
        propinit = {x: 0 for x in allowed}
        props = nucleotide_proportions({id:i},propinit)
        allowed = ["N","A","C","G","T","n","a","c","g","t"]
        nonallowed = [x for x in i.seq if x not in allowed]
        propnonstandard = sum([props[x] for x in props if x not in allowed])
        propn = sum([props[x] for x in props if x in ["n","N"]])
        i.id = i.id.split(" ")[0]
        descriptions[i.id] = str(i.description)
        i.description = ""
        seqlen = len(i.seq)
        if propn > 0.05:
            excluded += 1
            newrow = {"genome id": i.id, "genome description": i.description, "filter reason": "excess N", "Genome length": seqlen,"nonstandard proportion":propn}
            new_df = pd.DataFrame([newrow])
            filtered = pd.concat([filtered, new_df], ignore_index=True)
        elif float(propnonstandard) >= float(maxnonspandard) and i.id not in added:
            excluded += 1
            newrow = {"genome id": i.id, "genome description": i.description, "filter reason": "excess non-standard", "Genome length": seqlen,"nonstandard proportion":propnonstandard}
            new_df = pd.DataFrame([newrow])
            filtered = pd.concat([filtered, new_df], ignore_index=True)
        elif seqlen < 100:
            excluded += 1
            newrow = {"genome id": i.id, "genome description": i.description, "filter reason": "too short (<100bp)", "Genome length": seqlen,"nonstandard proportion":""}
            new_df = pd.DataFrame([newrow])
            filtered = pd.concat([filtered, new_df], ignore_index=True)
        elif (props["C"]+props["c"]) < 0.01:
            excluded += 1
            newrow = {"genome id": i.id, "genome description": i.description, "filter reason": "very low C content, probable 3base genome from genetic signatures", "Genome length": seqlen,"nonstandard proportion":""}
            new_df = pd.DataFrame([newrow])
            filtered = pd.concat([filtered, new_df], ignore_index=True)
        else:
            oldi = copy.copy(i)
            i = strip_polyA(i)
            if str(i.seq) != str(oldi.seq):
                trimmed = len(oldi.seq) - len(i.seq)
                # logger.info(f"{trimmed}bp have been trimmed from genome {i.id}")
                newrow = {"genome id": i.id, "genome description": i.description,
                          "filter reason": "trimmed polyA (not removed)",
                          "Genome length": seqlen, "nonstandard proportion / metric": trimmed}
                new_df = pd.DataFrame([newrow])
                filtered = pd.concat([filtered, new_df], ignore_index=True)
                notrimmed +=1
            if len(i.seq) < 100:
                excluded += 1
                newrow = {"genome id": i.id, "genome description": i.description, "filter reason": "too short (<100bp)",
                          "Genome length": seqlen, "nonstandard proportion": ""}
                new_df = pd.DataFrame([newrow])
                filtered = pd.concat([filtered, new_df], ignore_index=True)
                notrimmed -= 1
            else:
                outgenomes.append(i)
                added.append(i.id)
                included += 1

    outpath = outloc+"_filtered_input.fasta"
    SeqIO.write(outgenomes, outpath, "fasta")
    logger.info(f"total genomes with trailing polyA >5bp trimmed: {notrimmed}")
    logger.info(f"total genomes excluded for due to excess non ATGCN nucleotides: {excluded}")
    logger.info(f"See *_filtered_genomes.tsv file for details")
    return outpath,overallprops,included,filtered,descriptions

class RuntimeFormatter(logging.Formatter):
    """
    Custom logging formatter to include runtime in log messages.
    """
    def __init__(self, fmt=None, datefmt=None):
        super().__init__(fmt, datefmt)
        self.start_time = time.time()  # Record the start time

    def format(self, record):
        """
        Formats the log record to include runtime.

        Args:
            record (logging.LogRecord): The log record to format.

        Returns:
            str: The formatted log message.
        """
        # Calculate elapsed time
        elapsed_time = time.time() - self.start_time
        record.runtime = f"{elapsed_time:.2f}s"  # Add runtime to the record
        return super().format(record)

def select_pangraph_binary():
    arch = platform.machine().lower()
    system = platform.system().lower()

    # Normalize architecture
    if arch in ["x86_64", "amd64"]:
        arch = "x86_64"
    elif arch in ["aarch64", "arm64"]:
        arch = "aarch64"
    else:
        logger.error(f"Unsupported architecture: {arch}")
        return None

    # Select OS and libc variant
    if system == "darwin":
        if arch == "x86_64":
            return "pangraph-x86_64-apple-darwin"
        elif arch == "aarch64":
            return "pangraph-aarch64-apple-darwin"
    elif system == "linux":
        # Check for musl or gnu
        try:
            ldd_output = subprocess.check_output("ldd --version", shell=True, stderr=subprocess.STDOUT).decode()
            if "musl" in ldd_output:
                if arch == "x86_64":
                    return "pangraph-x86_64-unknown-linux-musl"
                elif arch == "aarch64":
                    return "pangraph-aarch64-unknown-linux-musl"
            else:
                if arch == "x86_64":
                    return "pangraph-x86_64-unknown-linux-gnu"
                elif arch == "aarch64":
                    return "pangraph-aarch64-unknown-linux-gnu"
        except Exception:
            # Default to gnu if musl check fails
            if arch == "x86_64":
                return "pangraph-x86_64-unknown-linux-gnu"
            elif arch == "aarch64":
                return "pangraph-aarch64-unknown-linux-gnu"
    elif system == "windows":
        if arch == "x86_64":
            return "pangraph-x86_64-pc-windows-gnu.exe"
    else:
        logger.error(f"Unsupported system: {system}")
        return None




def linear_transitive_chain_merge_fasta(inp, fasta_file, outp):
    """
    Merge sequences in a GFA file based on linear transitive chains.

    Args:
        inp (str): Path to the GFA file.
        fasta_file (str): Path to the FASTA file containing sequences.
        outp (str): Path to the output merged FASTA file.
    """

    # === Step 0: Check input files ===
    if not os.path.exists(inp):
        raise FileNotFoundError(f"GFA file {inp} does not exist.")
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file {fasta_file} does not exist.")
    # === Step 1: Parse the GFA ===
    gfa_file = inp
    with open(gfa_file) as f:
        lines = f.readlines()

    out_edges = defaultdict(list)
    in_edges = defaultdict(list)
    links_raw = []

    for line in lines:
        if line.startswith("L"):
            parts = line.strip().split("\t")
            a, a_orient, b, b_orient = parts[1], parts[2], parts[3], parts[4]
            out_edges[(a, a_orient)].append((b, b_orient))
            in_edges[(b, b_orient)].append((a, a_orient))
            links_raw.append(parts)

    # === Step 2: Build forward-only merge chains A + → B + where no other edges touch A+ or B+
    visited = set()
    merge_chains = []

    for (node, orient) in out_edges:
        if orient != '+':
            continue
        if (node in visited) or len(out_edges[(node, '+')]) != 1:
            continue

        next_node, next_orient = out_edges[(node, '+')][0]
        if next_orient != '+' or len(in_edges[(next_node, '+')]) != 1:
            continue

        # Walk forward to find full chain
        chain = [node]
        current = next_node
        while True:
            if current in chain:
                break  # avoid cycles
            if len(in_edges[(current, '+')]) == 1 and len(out_edges.get((current, '+'), [])) == 1:
                next_target, next_target_orient = out_edges[(current, '+')][0]
                if next_target_orient != '+' or len(in_edges[(next_target, '+')]) != 1:
                    break
                chain.append(current)
                current = next_target
            else:
                chain.append(current)
                break

        # Register all nodes in the chain
        for n in chain:
            visited.add(n)
        if len(chain) > 1:
            merge_chains.append(chain)

    # === Step 3: Read FASTA and merge chains ===

    records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Build merged sequences
    merged_records = []
    used_nodes = set()

    for chain in merge_chains:
        merged_id = "_".join(chain)
        merged_seq = Seq.Seq("").join([records[n].seq for n in chain])
        merged_records.append(SeqRecord.SeqRecord(merged_seq, id=merged_id, description=""))
        used_nodes.update(chain)

    # Add unmerged records
    unmerged_records = [rec for name, rec in records.items() if name not in used_nodes]

    # Write output
    output_file = outp
    SeqIO.write(unmerged_records + merged_records, output_file, "fasta")
    return output_file

def remove_belowdepth_nodes(gfa_file, depthcut):

    # Parse the GFA file
    with open(gfa_file) as f:
        lines = f.readlines()

    # Identify nodes with depth below the threshold
    nodes_to_remove = set()
    for line in lines:
        if line.startswith("S"):
            parts = line.strip().split("\t")
            node_id = parts[1]
            totlen= float(parts[3].split(":")[-1])
            len = float(parts[4].split(":")[-1])
            depth = totlen / len
            if depth <= depthcut:
                nodes_to_remove.add(node_id)
    return list(nodes_to_remove)


def run_pangraph(args,filtinput,outloc):
    """
    Runs the PanGraph tool to generate a pangenome graph and associated files.

    Args:
        args (argparse.Namespace): The arguments passed to the script, containing various settings and file paths.

    Returns:
        None
    """
    logger.info("Starting pangenome generation with PanGraph")
    topdir = os.path.dirname(os.path.abspath(__file__))  # Get the directory of the current script
    outloc = f"{args.outputfolder}/{args.outputprefix}"  # Define the output location for PanGraph files
    pangraph_log = open(outloc + "_pangraph.log", "w")  # Open a log file for PanGraph output

    if args.pangraphstrict:
        adds= " -S"
    else:
        adds = ""
    if args.pangraphmmseqs:
        addmmseqs = " -k mmseqs"
    else:
        addmmseqs = ""

    if args.maxdiv and platform.system().lower() == "darwin":
        logger.info("--maxdiv enables strict identity threshold for Pangraph (arm macOS version only) which is still in development. Beware this may not work as expected.")
        pangraphex = "pangraph-maxdiv-aarch64-darwin"
        pangraphpath = importlib.resources.files("dampa").joinpath("tools/pangraph/"+pangraphex)
        cmd = f"""{pangraphpath} build -s {args.pangraphident} -a {args.pangraphalpha} -b {args.pangraphbeta} -l {args.pangraphlen} -j {args.threads}{adds} {filtinput} > {outloc}.json && {pangraphpath} export gfa -o {outloc}_pangenome.gfa {outloc}.json && {pangraphpath} export block-consensus -o {outloc}_pangenome.fa  {outloc}.json"""
    elif args.maxdiv and platform.system().lower() != "darwin":
        logging.error("strict identity threshold only available in arm macOS version of pangraph (change once pangraph main branch updated)")
    else:
        pangraphex = select_pangraph_binary() # TODO may be issues when conda is installed as x86 but running on arm64
        pangraphpath = importlib.resources.files("dampa").joinpath("tools/pangraph/"+pangraphex)
        cmd = f"""{pangraphpath} build -s {args.pangraphident} -a {args.pangraphalpha} -b {args.pangraphbeta} -l {args.pangraphlen} -j {args.threads} {filtinput} > {outloc}.json && {pangraphpath} export gfa -o {outloc}_pangenome.gfa {outloc}.json && {pangraphpath} export block-consensus -o {outloc}_pangenome.fa  {outloc}.json"""
    subprocess.run(cmd, shell=True, stdout=pangraph_log, stderr=pangraph_log)
    if os.path.exists(f"{outloc}_pangenome.fa") and os.path.exists(f"{outloc}_pangenome.gfa") and os.path.exists(f"{outloc}.json"):
        logger.info("Pangraph ran successfully")
        if not args.keeplogs:
            os.remove(outloc + "_pangraph.log")
        # TODO improve linear merge (appears to be missing some possible merges)
        if args.pangraphdepth > 0:
            lowdepthnodes = remove_belowdepth_nodes(f"{outloc}_pangenome.gfa", args.pangraphdepth)  # Remove nodes below a certain depth
        else:
            lowdepthnodes = []
        linear_transitive_chain_merge_fasta(f"{outloc}_pangenome.gfa",f"{outloc}_pangenome.fa",f"{outloc}_pangenome_lin.fa",lowdepthnodes)# Remove the log file if not keeping logs
        # terminal_node_additional_clustering()

        logger.info("Pangenome graph linear chain merging completed")
    else:
        logger.error(f"One or more of pangraph outputs ({args.outputprefix}_pangenome.gfa, {args.outputprefix}_pangenome.fa, {args.outputprefix}.json) in {args.outputfolder} are not present. Check for error in pangraph log")
    return

def run_finalprobetools(args, inprobes,originput):
    """
    Runs the final probe design step using the Probetools tool.

    Args:
        args (argparse.Namespace): The arguments passed to the script, containing various settings and file paths.
        inprobes (str): Path to the input probes file.

    Returns:
        tuple: A tuple containing the path to the final capture file and the total number of probes generated.
    """
    outloc = f"{args.outputfolder}/{args.outputprefix}"  # Define the output location for Probetools files
    finalpref = f"{outloc}_probetools_final"  # Define the prefix for final Probetools output files
    finalprobes = f"{finalpref}_probes.fa"  # Define the path for the final probes file
    probetools_log = open(outloc + "_probetools.log", "w")  # Open a log file for Probetools output
    # Construct the command to run Probetools with the specified parameters

    if args.nodust:
        dust=" -y N"
    else:
        dust=" -y Y"
    probetoolspath = importlib.resources.files("dampa").joinpath("tools/probetools/probetools_v_0_1_11_mod.py")
    cmd = f"python {probetoolspath} makeprobeswinput -t {originput}{dust} -b 100 -x {inprobes} -o {finalpref} -i {args.probetoolsidentity} -l {args.probetoolsalignmin} -T {args.threads} -L {args.probetools0covnmin} -c 100 -d {args.maxambig}"
    subprocess.run(cmd, shell=True, stdout=probetools_log, stderr=probetools_log)  # Execute the command

    # Check if the expected output file is present
    if os.path.exists(finalprobes):
        shutil.move(finalprobes, inprobes)  # Move the final probes file to the input probes path
        totaprobes, addedprobes = get_probeno(inprobes, "_round_")  # Get the total and added probes count
        logger.info(f"Generated {addedprobes} additional probes using probetools")
        if not args.keeplogs:
            os.remove(outloc + "_probetools.log")  # Remove the log file if not keeping logs
        return f"{finalpref}_capture.pt", totaprobes  # Return the final capture file path and total probes count
    else:
        logger.error(f"Probetools output file {finalprobes} not present. Check for error in details probetools log.")


def get_clusters(args):
    """
    Loads cluster assignments from a file.

    Args:
        args (argparse.Namespace): The arguments passed to the script, containing various settings and file paths.

    Returns:
        dict: Dictionary mapping strains to clusters.
    """
    clusterdict = {}
    cluster = ""
    if args.clustertype == "cdhit":
        clstrfile = open(args.clusterassign, "r").readlines()
        for line in clstrfile:
            if line.startswith(">"):
                cluster = line.strip().replace(">Cluster ", "")
            else:
                strain = line.split(" ")[1][1:-3]
                clusterdict[strain] = cluster
    elif args.clustertype == "tabular":
        clstrfile = open(args.clusterassign, "r").readlines()
        for line in clstrfile:
            strain = line.split("\t")[0]
            cluster = line.split("\t")[1]
            clusterdict[strain] = cluster
    return clusterdict

def load_capture_data(capture_path):
    """
    Loads capture data from a file.

    Args:
        capture_path (str): Path to the capture data file.

    Returns:
        dict: Dictionary containing capture data.
    """
    # print(f'Loading capture data from {capture_path}...')
    with open(capture_path, 'r') as input_file:
        headers, seqs, depths = [], [], []
        for line in input_file:
            if line[0] == '>':
                if not seqs == [] and seqs[-1] == '':
                    print(header[-1])
                header = line.strip().lstrip('>')
                headers.append(header)
                seqs.append('')
                depths.append('')
            elif line[0] == '$':
                seqs[-1] += line.strip()
            elif line[0] == '#':
                depths[-1] += line.strip()
    seqs = [seq.lstrip('$') for seq in seqs]
    depths = [[int(d) for d in depth.lstrip('#').split(',')] for depth in depths]
    depths = [replace_short_zeros(x, 10) for x in depths]
    if len(set(len(c) for c in [headers, seqs, depths])) != 1:
        logger.error(f'The number of headers, seqs, and probe depth lists do not match in {capture_path}.')
        exit(1)
    for header in headers:
        if headers.count(header) > 1:
            logger.error(f'Header {header} appears more than once in {capture_path}.\n')
            exit(1)
    capture_data = {header: (seq, depth) for header, seq, depth in zip(headers, seqs, depths)}
    for header in capture_data.keys():
        if len(capture_data[header][0]) != len(capture_data[header][1]):
            logger.error(f'Seq length and probe depth list length do not match for entry {header}.')
            exit(1)
    # print(f' Total targets loaded: {"{:,}".format(len(capture_data))}')
    return capture_data

def runprobetoolscapture(args,probes):
    """
    Runs the Probetools capture step.

    Args:
        args (argparse.Namespace): The arguments passed to the script, containing various settings and file paths.
        probes (str): Path to the probes file.

    Returns:
        str: Path to the capture output file.
    """
    outloc = args.outputfolder+"/"+args.outputprefix
    current_directory = os.path.dirname(os.path.abspath(__file__))
    capture_log = open(outloc + "_capture.log", "w")
    with open(os.devnull, 'w') as devnull:
        if args.nodust:
            dust = " -y Y"
        else:
            dust = " -y N"

        # outf = open("/Users/mpay0321/Dropbox/Probe_design_project/2025-01-29_integrate_probetools_probebench/stdout.txt",'w')
        cmd = f"python {current_directory}/tools/probetools/probetools_v_0_1_11_mod.py capture -t {args.input}{dust} -p {probes} -o {outloc} -i {args.probetoolsidentity} -l {args.probetoolsalignmin} -T {args.threads}"
        subprocess.run(cmd, shell=True,stdout=capture_log, stderr=capture_log)
    outf = f"{args.outputfolder}/{args.outputprefix}_capture.pt"
    if os.path.exists(outf):
        logger.info(f"Probetools capture ran successfully")
        if not args.keeplogs:
            os.remove(outloc + "_capture.log")
        return outf
    else:
        logger.error(f"Probetools capture output file {outf} not present. Check for error in capture log.")

def get_ambig_count(seq):
    """
    Counts the number of ambiguous bases in a sequence.

    Args:
        seq (str): The nucleotide sequence.

    Returns:
        int: The number of ambiguous bases.
    """
    ambig = len([x for x in seq if x not in ["A","T","G","C","a","t","c","g"]])
    return ambig


def strip_polyA(s):
    '''Strip poly-A/N/M tail, but only if tail is longer than 5 bases.'''
    seq_str = str(s.seq)
    match = re.match(r'(.*[^AaNnMm])([AaNnMm]+)$', seq_str)

    if match:
        core, tail = match.groups()
        if len(tail) > 5:
            s.seq = Seq.Seq(core)

    return s

def split_pangenome_into_probes(input_fasta, output_fasta, probe_length,probe_step,maxambig,probeprefix=""):
    """
    Splits a pangenome into probes of specified length and step size.

    Args:
        input_fasta (str): Path to the input pangenome FASTA file.
        output_fasta (str): Path to the output probes FASTA file.
        probe_length (int): Length of each probe.
        probe_step (int): Step size for generating probes.
        maxambig (int): Maximum number of ambiguous bases allowed in a probe.

    Returns:
        None
    """
    outprobes = []
    totalprobes = 0
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequence = str(record.seq)
        seq_id = record.id
        seq_length = len(sequence)
        if seq_length >= probe_length:
            probes = []  # Store generated probes

            # Generate probes with the specified step size
            for i in range(0, seq_length - probe_length + 1, probe_step):
                probes.append((i, i + probe_length))

            # Ensure the last part of the sequence is covered
            last_probe_start = seq_length - probe_length
            if probes and probes[-1][1] < seq_length:
                probes.append((last_probe_start, seq_length))
            seqno = 1
            # Print probes
            for start, end in probes:
                piece = sequence[start:end]
                numambig = get_ambig_count(piece)
                if numambig <= maxambig:
                    # Create a new SeqRecord for each piece
                    piece_id = f"{seq_id}_piece_{seqno}"
                    seqno+=1
                    probeprefixm = str(probeprefix+"_")
                    piece_record = SeqRecord.SeqRecord(Seq.Seq(piece),id=probeprefixm+piece_id,description="")
                    outprobes.append(piece_record)
                    totalprobes+=1
    # Write the probes to the output file
    SeqIO.write(outprobes, output_fasta, "fasta")
    logger.info(f"Extracted {totalprobes} probes from pangenome")

def make_summaries(args,ptcountfile,totalprobes):
    """
    Generates summary statistics and plots for probe coverage.

    Args:
        args (argparse.Namespace): The arguments passed to the script, containing various settings and file paths.
        ptcountfile (str): Path to the probe count file.
        totalprobes (int): Total number of probes.

    Returns:
        None
    """
    outloc = f"{args.outputfolder}/{args.outputprefix}"
    clusterdict = {}
    if args.clusterassign:
        clusterdict = get_clusters(args)
    depthls = []
    names = []
    meancovs = []
    propcovs = []
    maxcheck = int(args.maxdepth_describe)
    ptdata = load_capture_data(ptcountfile)
    for name, data in ptdata.items():
        stats = process_count_pt(data, maxcheck)
        propcovs.append(stats[0])
        depthls.append(stats[2])
        names.append(name)
        meancovs.append(stats[1])


    logger.info("Making summaries and figures")
    make_stats(propcovs, names, maxcheck, meancovs, outloc,totalprobes,args.report0covperc)
    make_propsplot(propcovs, names, maxcheck, outloc,clusters=clusterdict)
    make_genome_plots(ptdata, outloc)



def setup_logging():
    """
    Sets up logging with a custom formatter to include runtime.

    Returns:
        logging.Logger: The configured logger.
    """
    formatter = RuntimeFormatter("%(asctime)s - [Runtime: %(runtime)s] - %(levelname)s - %(message)s")
    formatter.default_time_format = "%Y-%m-%d %H:%M:%S"
    formatter.default_msec_format = ""
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger("runtime_logger")
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger

def get_probeno(probefile,subsetstr=""):
    """
    Counts the total number of probes and the number of added probes in a file.

    Args:
        probefile (str): Path to the probes file.
        subsetstr (str): Substring to identify added probes.

    Returns:
        tuple: Total number of probes and number of added probes.
    """
    outprobes = open(probefile, "r").read().splitlines()
    if subsetstr != "":
        addedprobes = [x for x in outprobes if x.startswith(">") and subsetstr in x]
    else:
        addedprobes = []
    totalprobes = len([x for x in outprobes if x.startswith(">")])
    return totalprobes,len(addedprobes)


def convert_ambig_nucs(props,seq):
    ambig = {"R":["A","G"],
             "K":['G','T'],
             'S':['G','C'],
             'Y':['C','T'],
             'M':['A','C'],
             'W':['A','T'],
             'B':['C','G','T'],
             'H':['A','C','T'],
             'D':['A','G','T'],
             'V':['A','C','G'],
             'N':['A','T','C','T']}
    ambigconv = {"R":["A","G"],
             "K":['G','T'],
             'S':['G','C'],
             'Y':['C','T'],
             'M':['A','C'],
             'W':['A','T'],
             'B':['C','G','T'],
             'H':['A','C','T'],
             'D':['A','G','T'],
             'V':['A','C','G'],
             'N':['A','T','C','T']}



    return seq

def subambig(probes,props):

    ambig = {"R":["A","G"],
             "K":['G','T'],
             'S':['G','C'],
             'Y':['C','T'],
             'M':['A','C'],
             'W':['A','T'],
             'B':['C','G','T'],
             'H':['A','C','T'],
             'D':['A','G','T'],
             'V':['A','C','G'],
             'N':['A','T','C','T']}

    rawprobes = SeqIO.parse(probes,"fasta")
    outprobes = []
    for probe in rawprobes:
        pseq = str(probe.seq).upper()
        outseq = ""
        for n in pseq:
            if n not in ['A','T','C','G']:
                choices = ambig[n]
                weights = [props[x] for x in choices]
                subnuc = random.choices(choices, weights=weights, k=1)[0]
                outseq += subnuc
            else:
                outseq += n
        probe.seq=Seq.Seq(outseq)
        outprobes.append(probe)
    SeqIO.write(outprobes,probes,"fasta")

def cleanup(args,filtgenomes=""):
    if args.keeptmp:
        logger.info("Keeping temporary files")
    else:
        logger.info("Cleaning up temporary files (use --keeptmp to keep pangenome graph and other temporary files)")
        outloc = f"{args.outputfolder}/{args.outputprefix}"
        tormsuffixes = ["_pangenome.gfa",".json","_pangenome.fa",
                       "_probetools_final_long_stats_report.tsv","_probetools_final_summary_stats_report.tsv",
                       "_probetools_final_capture.pt","_probetools_final_low_cov_seqs.fa","_pangenome_lin.fa"]
        for i in tormsuffixes:
            print(f"{outloc}{i}")
            if os.path.exists(f"{outloc}{i}"):
                os.remove(f"{outloc}{i}")
        if filtgenomes != "":
            filtgenomesrm=glob.glob(f"{filtgenomes}*")
            for i in filtgenomesrm:
                os.remove(i)
        mmseqsrm = glob.glob(f"{outloc}_alltypesdb*") + glob.glob(f"{outloc}_DB.*") + glob.glob(f"{outloc}_Dbreps*")
        for i in mmseqsrm:
            if os.path.exists(i):
                os.remove(i)
    logger.info(f"Cleaned up tmp files in {args.outputprefix}")

def design_probes(args,filtered_input,probeprefix,overallprops,rminp,outloc):
    probename = outloc + "_probes.fasta"
    pangenomefasta = outloc + "_pangenome_lin.fa"
    run_pangraph(args,
                 filtered_input,outloc)  # TODO possibly add check where probes are mapped onto each other in a progressive way. each time coverage of a probe by other probes is >1 across full length (at some high identity) then remove probe that is covered would remove lots of similar probes from ends of pancontigs?


    split_pangenome_into_probes(pangenomefasta, probename, args.probelen, args.probestep, args.maxambig, probeprefix)

    if not args.skip_padding:
        make_padded_probes(pangenomefasta, probename, args.minlenforpadding, probeprefix=probeprefix)
    if not args.skip_probetoolsfinal:
        finalcaptureout, totalprobes = run_finalprobetools(args, probename, filtered_input,outloc)
        if not args.skipsubambig:
            subambig(probename, overallprops)
        if not args.skip_summaries:
            make_summaries(args, finalcaptureout, totalprobes,outloc)
    else:
        totalprobes = get_probeno(probename)
        if not args.skipsubambig:
            subambig(probename, overallprops)
        if not args.skip_summaries:
            finalcaptureout = runprobetoolscapture(args, probename,outloc, filtered_input)
            make_summaries(args, finalcaptureout, totalprobes,outloc)
    if rminp:
        cleanup(args, filtered_input)
    else:
        cleanup(args)

    return totalprobes

def write_filtered_genomes(filtered, outloc,descriptions):
    """
    Writes filtered genomes to a file.

    Args:
        filtered (pandas.DataFrame): DataFrame containing filtered genome information.
        outloc (str): Output location for the filtered genomes file.

    Returns:
        str: Path to the written filtered genomes file.
    """
    outpath = f"{outloc}_filtered_genomes.tsv"

    # filtered["genome description"] = descriptions
    filtered["genome description"] = filtered["genome id"].map(descriptions)

    filtered.to_csv(outpath, sep="\t", index=False)
    # logger.info(f"Filtered genomes written to {outpath}")
    return outpath

def get_args():
    def File(MyFile):
        if not os.path.isfile(MyFile):
            raise argparse.ArgumentTypeError(MyFile + ' does not exist or is not a file.')
        return MyFile

    # Define a function to check directories exist, as a type for the argparse.
    def Dir(MyDir):
        if not os.path.isdir(MyDir):
            raise argparse.ArgumentTypeError(MyDir + \
                                             ' does not exist or is not a directory.')
        return MyDir
    def DirorFolder(Mypath):
        if not os.path.isdir(Mypath) and not os.path.isfile(Mypath):
            raise argparse.ArgumentTypeError(Mypath + \
                                             ' does not exist or is not a file or directory.')
        return Mypath



    cwd = Path.cwd()
    parser = argparse.ArgumentParser(description="DAMPA: Diversity Aware Metagenomic Panel Assignment",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    design = subparsers.add_parser("design", help="Design probes from input genomes",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    design_inputs = design.add_argument_group("Input/Output options")

    design_inputs.add_argument("-g", "--input", required=True, help="Either folder containing individual genome fasta files OR a single fasta file containing all genomes (files must end in .fna, .fa or .fasta)",type=DirorFolder)
    design_inputs.add_argument("-c", "--clusterassign", help="clstr file from cd-hit",
                        type=File)
    design_inputs.add_argument("--clustertype",
                        help="type of cluster file input cdhit (produced by cdhit) or tabular (genome and cluster tab delimited) ", choices=['cdhit','tabular'],
                        default='tabular')
    design_inputs.add_argument("--maxnonspandard",
                        help="maximum proportion of genome that can be non ATGC (0-1)",type=float,default=0.01)

    design_inputs.add_argument("-o", "--outputfolder", type=Dir,
                        help="path to output folder",default=f"{cwd}/")
    design_inputs.add_argument("-p", "--outputprefix", default="probebench_run",
                        help="prefix for all output files and folders")

    general = design.add_argument_group("General settings")

    general.add_argument("-l", "--probelen", type=int, default=120,help="length of output probes")
    general.add_argument("-s", "--probestep", type=int, default=120, help="step of probes (for no overlap set to same as probelen)")
    general.add_argument("--skipsubambig",
                        help="do NOT substitute ambiguous nucleotides (by default N or other ambiguous nucleotides are substituted for ATGC in a random selection weighted by proportions in input genomes",action='store_true')



    pangraphsettings = design.add_argument_group("Pangraph settings")

    pangraphsettings.add_argument("--pangraphident", type=int, default=20,choices=[5,10,20],help="Pangenome percentage identity setting allowable values are 5,10 or 20")
    pangraphsettings.add_argument("--pangraphalpha", type=float, default=0,help="Energy cost for splitting a block during alignment merger. Controls graph fragmentation")
    pangraphsettings.add_argument("--pangraphbeta", type=float, default=10,help="Energy cost for diversity in the alignment. A high value prevents merging of distantly-related sequences in the same block")
    pangraphsettings.add_argument("--pangraphlen", type=int, default=90,help="Minimum length of a node to allow in pangenome graph")
    pangraphsettings.add_argument("--pangraphstrict", help="enable the -S strict identity option which limits merges to 1/pangraphbeta divergence",action='store_true')
    pangraphsettings.add_argument("--pangraphmmseqs", action='store_true',
                                    help="use mmseqs2 (slower but better at low identity) in pangraph instead of minimap")
    pangraphsettings.add_argument("--pangraphdepth", type=int, default=0,help="Minimum depth of a node to allow in pangenome graph. Nodes with depth below this will be removed from the graph. Set to 0 to not remove any nodes based on depth")


    probetoolssettings = design.add_argument_group("Probetools settings")

    probetoolssettings.add_argument("--probetoolsidentity", type=int, default=85,
                                    help="Minimum identity in probe match to target to call probe binding")
    probetoolssettings.add_argument("--probetoolsalignmin", type=int, default=90,help="Minimum length (bp) of probe-target binding to allow call of binding")
    probetoolssettings.add_argument("--probetools0covnmin", type=int, default=20,
                                    help="Minimum length (bp) of 0 coverage region in input genomes to trigger design of additional probes")
    probetoolssettings.add_argument("--maxambig",help="The maximum number of ambiguous bases allowed in a probe",type=int,default=10)
    probetoolssettings.add_argument("--nodust", help="Do not run low complexity filter in BLAST (within probetools). If sample has very low GC or is very repetitive this option can be enabled to prevent low complexity regions from being removed",action='store_true')

    preclustersettings = design.add_argument_group("preclustering settings")
    preclustersettings.add_argument("--clusterer", default="cdhit",
                                    help="precluster using cdhitest or mmseqs2", choices=['cdhit','mmseqs'])
    preclustersettings.add_argument("--clustering_inputno_trigger",
                                  help="if number of input sequences exceeds this number then mmseqs will be used to deduplcate genomes above 99.9 percent identity",type=int,
                                  default=5000)
    preclustersettings.add_argument("--clusterident", type=float, default=0.999,
                                    help="Minimum identity to cluster genomes")
    preclustersettings.add_argument("--lowclusterident", type=float, default=0.85,
                                    help="Minimum identity to cluster genomes at low threshold (used in cluster2step)")
    preclustersettings.add_argument("--clustercov", type=float, default=1,
                                    help="Minimum coverage of genomes over which clusterident must apply (0-1)")
    preclustersettings.add_argument("--cluster2step",action='store_true',
                                    help="perform initial high threshold clustering followed by low threshold clustering of representatives, remove low level clusters composed of one high level cluster with <= outliersizelimit members")
    preclustersettings.add_argument("--outliersizelimit",
                                    help="In two step clustering if a cluster is <=outliersizelimit at both high and low identity then it is treated as an outlier",type=int,
                                    default=1)

    additionalsettings = design.add_argument_group("Additional settings")

    additionalsettings.add_argument("--skip_padding",
                        help="do not generate additional probes for pangenome nodes between pangraphlen and probelen in length. i.e. if padding is run 30*T would be added to the end of a 90bp pancontig",action='store_true')
    additionalsettings.add_argument("--padding_nuc",
                        help="nucleotide to use for padding probes to args.probelen", choices=['A',"T","C","G"],
                        default='T')
    additionalsettings.add_argument("--minlenforpadding",
                        help="minimum length for a pancontig for it to be padded (WARNING setting this below ~80 may result in probes that do not effectively bind, leave these small sequences for final probetools step)", type=int,
                        default='90')
    additionalsettings.add_argument("--skip_probetoolsfinal",
                        help="do NOT run final probe design step. i.e. this step uses probetools to design probes to regions that are not represented in the pangenome",action='store_true')


    additionalsettings.add_argument("-t","--threads",
                            help="number of threads",
                            type=int,
                            default=1)
    additionalsettings.add_argument("--keeplogs",
                        help="keep logs containing output from pangraph and probetools",action='store_true')
    additionalsettings.add_argument("--keeptmp",
                        help="keep intermediate files from pangraph and probetools",action='store_true')
    additionalsettings.add_argument("--skip_summaries",
                        help="do NOT run visualisation generaton of dampa probes relative to input genomes",action='store_true')
    additionalsettings.add_argument("--maxdepth_describe",
                        default=1,help="Maximum depth of probe coverage to describe separately. i.e. if 1 there will be 0,1 and >1 depth categories")
    additionalsettings.add_argument("--report0covperc",
                        help="threshold above which genomes are reported as having too much of their genome not covered by any probes",type=float,default=1)
    additionalsettings.add_argument("--version",
                                    help="print version and exit",
                                    action='store_true')
    additionalsettings.add_argument("--maxdiv",
                                    help="use new maxdiv pangraph version",
                                    action='store_true')

    evaluate = subparsers.add_parser("eval", help="Evaluate performance of a probe set against a set of genomes",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    eval_inputs = evaluate.add_argument_group("Input/Output options")

    eval_inputs.add_argument("-g", "--input", required=True, help="Genomes to check probe coverage. \n"
                                                                  "If genomes either folder containing individual genome fasta files OR a single fasta file containing all genomes (files must end in .fna, .fa or .fasta)\n"
                                                                  "If capture file then a pt file from a previous pangraph design or pangraph eval run",type=DirorFolder)
    eval_inputs.add_argument("--inputtype",
                        help="type of cluster file input cdhit (produced by cdhit) or tabular (genome and cluster tab delimited) ", choices=['genomes','capture'],
                        default='genomes')
    eval_inputs.add_argument("-q", "--probes", required=True,
                             help="Fasta file containing probes to evaluate (files must end in .fna, .fa or .fasta)",
                             type=File)
    eval_inputs.add_argument("-c", "--clusterassign", help="clstr file from cd-hit",
                        type=File)
    eval_inputs.add_argument("--clustertype",
                        help="type of cluster file input cdhit (produced by cdhit) or tabular (genome and cluster tab delimited) ", choices=['cdhit','tabular'],
                        default='tabular')
    eval_inputs.add_argument("--filtnonstandard",
                        help="remove genomes with non standard nucleotides i.e. not A,T,G,C or N",action='store_true')

    eval_inputs.add_argument("-o", "--outputfolder", type=Dir,
                        help="path to output folder",default=f"{cwd}/")
    eval_inputs.add_argument("-p", "--outputprefix", default="probebench_run",
                        help="prefix for all output files and folders")

    probetooleval = evaluate.add_argument_group("Probetools settings")

    probetooleval.add_argument("--probetoolsidentity", type=int, default=85,
                                    help="Minimum identity in probe match to target to call probe binding")
    probetooleval.add_argument("--probetoolsalignmin", type=int, default=90,help="Minimum length (bp) of probe-target binding to allow call of binding")
    probetooleval.add_argument("--nodust",
                                    help="Do not run low complexity filter in BLAST (within probetools). If sample has very low GC or is very repetitive this option can be enabled to prevent low complexity regions from being removed",
                                    action='store_true')

    additionaleval= evaluate.add_argument_group("Additional settings")

    additionaleval.add_argument("-t","--threads",
                            help="number of threads",
                            type=int,
                            default=1)
    additionaleval.add_argument("--keeplogs",
                        help="keep logs containing output from pangraph and probetools",action='store_true')

    additionaleval.add_argument("--maxdepth_describe",
                                    default=1,
                                    help="Maximum depth of probe coverage to describe separately. i.e. if 1 there will be 0,1 and >1 depth categories")
    additionaleval.add_argument("--report0covperc",
                                    help="threshold above which genomes are reported as having too much of their genome not covered by any probes",
                                    type=float, default=1)
    additionaleval.add_argument("--version",
                                    help="print version and exit",
                                    action='store_true')
    additionaleval.add_argument("--keeptmp",
                        help="keep intermediate files from pangraph and probetools",action='store_true')

    args = parser.parse_args()

    return args

def main():
    """
    get args
    setup outputs dir and tmp folder (use uuid)
    optionally run splitfasta function
    run
    """
    global logger
    logger = setup_logging()
    logger.info("Starting")

    args = get_args()

    if args.command == "design":
        if args.version:
            print(f"version {dampaversion}")
            sys.exit(0)
        logger.info("Running dampa design")
        args.filtnonstandard = True

        probeprefix = args.input.split("/")[-1].replace(".fasta","").replace(".fa","").replace(".fna","")

        logger.info("Filter genomes with too many non standard nucleotides")

        filtered = pd.DataFrame(columns=["genome id","genome description","filter reason","Genome length","nonstandard proportion"])

        filtered_input,overallprops,included,filtered,descriptions = filter_for_nonstandard_inputs(args.input, outloc,args.maxnonspandard,filtered)
        rminp = False
        if args.cluster2step:
            rminp = True
            logger.info("Using cdhit to cluster genomes in two stages")
            regstrains,outliersstrains,filtered = cdhit_2step_clust(args, filtered_input,filtered)
            logger.info("2 stage clustering finished")

            regprobes = design_probes(args, regstrains, probeprefix, overallprops,rminp,outloc)
            logger.info("regular probes designed")
            altprobes = 0
            if outliersstrains != False:
                outloc2 = f"{args.outputfolder}/{args.outputprefix}_outliers"
                altprobes = design_probes(args, outliersstrains, probeprefix, overallprops,rminp,outloc2)
                logger.info("outlier probes designed")
            logger.info(f"dampa design finished. Regular probes: {regprobes}. Outliers probes: {altprobes}")
        else:
            if included > args.clustering_inputno_trigger:
                if args.clusterer == "mmseqs":
                    rminp = True
                    logger.info("Using mmseqs to cluster genomes")
                    filtered_input = mmseqs_subset(args,filtered_input)
                elif args.clusterer == "cdhit":
                    rminp = True
                    logger.info("Using cdhit to cluster genomes")
                    filtered_input = cdhit_subset(args, filtered_input)

            totalprobes = design_probes(args, filtered_input, probeprefix, overallprops,rminp,outloc)

            logger.info(f"dampa design finished. Total probes: {totalprobes}")

        write_filtered_genomes(filtered, outloc,descriptions)
    elif args.command == "eval":
        if args.version:
            print(f"version {dampaversion}")
            sys.exit(0)
        logger.info("Running dampa eval")
        totalprobes = get_probeno(args.probes)[0]
        if args.inputtype != "capture":
            finalcaptureout = runprobetoolscapture(args, args.probes)
        else:
            finalcaptureout = args.input
        make_summaries(args, finalcaptureout,totalprobes)
        cleanup(args)
        logger.info("dampa eval finished")
if __name__ == "__main__":
    main()