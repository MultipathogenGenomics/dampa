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
import importlib.resources
# from tools.gather_probe_depth_stats import make_stats, make_propsplot, process_count_pt
from dampa.vis.plot_over_genomelen import make_genome_plots,replace_short_zeros
from dampa.tools.gather_probe_depth_stats import make_stats, make_propsplot, process_count_pt
from dampa import __version__ as dampaversion


def mmseqs_subset(args,filtinput):
    """
    run the following commands
    mmseqs createdb ingenomes.fasta alltypesdb

    mmseqs cluster alltypesdb $pref""DB tmp --min-seq-id mmident -c mmcov

    mmseqs createtsv alltypesdb $pref""DB $pref""cluster.tsv

    mmseqs createsubdb $pref""DB alltypesdb $pref""Dbreps

    mmseqs convert2fasta $pref""DBreps $pref""DB.fasta

    """
    outloc = f"{args.outputfolder}/{args.outputprefix}"
    mmseqs_log = open(outloc + "_mmseqs.log", "w")

    cmd = f"mmseqs createdb {filtinput} {outloc}_alltypesdb"
    subprocess.run(cmd, shell=True, stdout=mmseqs_log, stderr=mmseqs_log)
    cmd = f"mmseqs cluster {outloc}_alltypesdb {outloc}_DB tmp --min-seq-id {args.mmident} -c {args.mmcov} --cov-mode 1"
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

def detect_os_arch():
    """
        Detects the operating system and architecture of the current machine, and attempts to determine the type of libc used on Linux systems.

    Returns:
        dict: A dictionary containing the following keys:
            - "os": The operating system (e.g., "linux", "windows", "darwin").
            - "arch": The machine architecture (e.g., "x86_64", "arm64").
            - "libc": The type of libc used on Linux systems ("gnu", "musl", or None if undetermined).
    """
    system = platform.system().lower()  # "linux", "windows", "darwin" (macOS)
    machine = platform.machine().lower()  # "x86_64", "arm64", etc.

    libc_type = None

    if system == "linux":
        # Check for GNU libc explicitly
        try:
            result = subprocess.run(["ldd", "--version"], capture_output=True, text=True)
            if "GNU" in result.stdout:
                libc_type = "gnu"
        except FileNotFoundError:
            pass  # ldd might not exist, continue checking for musl

        # Check for musl if GNU was not detected
        if libc_type is None:
            try:
                with open("/proc/self/maps", "rb") as f:
                    if b"musl" in f.read():
                        libc_type = "musl"
            except FileNotFoundError:
                pass  # Unable to determine libc

    return {
        "os": system,
        "arch": machine,
        "libc": libc_type
    }


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

def nucleotide_proportions(sequences):
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
        overallseq += str(sequence.seq)
    counts = Counter(overallseq)  # Count occurrences of each nucleotide
    total = sum(counts.values())  # Total nucleotides
    proportions = {nt: count / total for nt, count in counts.items()}# Calculate proportions
    return proportions

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

def filter_for_nonstandard_inputs(genomes,outfolder,maxnonspandard):
    """
    Filters out genomes with a high proportion of non-standard nucleotides.

    Args:
        genomes (str): Path to the input genomes FASTA file.
        outfolder (str): Path to the output folder.
        maxnonspandard (float): Maximum allowed proportion of non-standard nucleotides.

    Returns:
        tuple: Path to the filtered genomes file and overall nucleotide proportions.
    """
    ingenomes = SeqIO.parse(genomes, "fasta")
    ingenomes = to_dict_remove_dups(ingenomes)
    outgenomes = []
    excluded = 0
    included = 0
    overallprops = nucleotide_proportions(ingenomes)
    added = []
    for id,i in ingenomes.items():
        props = nucleotide_proportions({id:i})
        allowed = ["N","A","C","G","T","n","a","c","g","t"]
        nonallowed = [x for x in i.seq if x not in allowed]
        propnonstandard = sum([props[x] for x in props if x not in allowed])
        propn = sum([props[x] for x in props if x in ["n","N"]])
        i.id = i.id.split(" ")[0]
        i.description = ""
        seqlen = len(i.seq)
        if propn > 0.05:
            logger.info(f"genome {i.id} has excess N and has been excluded")
            excluded += 1
        elif float(propnonstandard) >= float(maxnonspandard) and i.id not in added:
            notallowedstr = ",".join(list(set(nonallowed)))
            logger.info(f"genome {i.id} has non standard chraracters: {notallowedstr} and has been excluded")
            excluded += 1
        elif seqlen < 150:
            logger.info(f"genome {i.id} is too short ({seqlen}) and has been excluded")
            excluded += 1
        else:
            outgenomes.append(i)
            added.append(i.id)
            included += 1


    outpath = outfolder + "/" + genomes.split("/")[-1].replace(".fasta","").replace(".fa","").replace(".fna","")
    outpath = outpath + "_filt.fasta"
    SeqIO.write(outgenomes, outpath, "fasta")
    logger.info(f"total genomes excluded for due to excess non ATGCN nucleotides: {excluded}")
    return outpath,overallprops,included

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

def get_pangraphex(osarch):
    """
    Determines the appropriate PanGraph executable based on the operating system and architecture.

    Args:
        osarch (dict): Dictionary containing OS, architecture, and libc information.

    Returns:
        str: Path to the PanGraph executable.
    """
    if osarch["os"] == "linux":
        if osarch["arch"] == "arm64":
            if osarch["libc"] == "gnu":
                return "pangraph-aarch64-linux-gnu"
            elif osarch["libc"] == "musl":
                return "pangraph-aarch64-linux-musl"
            else:
                logger.error(f"Unsupported os/architecture/lobc combination {osarch["os"]}/{osarch['arch']}/{osarch['libc']}")
        elif osarch["arch"] == "x86_64":
            if osarch["libc"] == "gnu":
                return "pangraph-x86_64-linux-gnu"
            elif osarch["libc"] == "musl":
                return "pangraph-x86_64-linux-musl"
            else:
                logger.error(f"Unsupported os/architecture/lobc combination {osarch["os"]}/{osarch['arch']}/{osarch['libc']}")
    elif osarch["os"] == "darwin":
        if osarch["arch"] == "arm64":
            return "pangraph-aarch64-darwin"
        elif osarch["arch"] == "x86_64":
            return "pangraph-x86_64-darwin"
        else:
            logger.error(f"Unsupported os/architecture combination {osarch["os"]}/{osarch['arch']}")
    elif osarch["os"] == "windows":
        return "pangraph-x86_64-pc-windows-gnu.exe"
    else:
        logger.error(f"Unsupported os/architecture/lobc combination {osarch["os"]}/{osarch['arch']}/{osarch['libc']}")


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

def run_pangraph(args,filtinput):
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

    osarch = detect_os_arch()

    if args.pangraphstrict:
        adds= " -S"
    else:
        adds = ""

    if args.maxdiv and osarch["os"] == "darwin":
        pangraphex = "pangraph-maxdiv-aarch64-darwin"
        pangraphpath = importlib.resources.files("dampa").joinpath("tools/pangraph/"+pangraphex)
        cmd = f"""{pangraphpath} build -s {args.pangraphident} -a {args.pangraphalpha} -b {args.pangraphbeta} -l {args.pangraphlen} -j {args.threads}{adds} {filtinput} > {outloc}.json && {pangraphpath} export gfa -o {outloc}_pangenome.gfa {outloc}.json && {pangraphpath} export block-consensus -o {outloc}_pangenome.fa  {outloc}.json"""
    elif args.maxdiv and osarch["os"] != "darwin":
        logging.error("strict identity threshold only available in arm macOS version of pangraph (change once pangraph main branch updated)")
    else:
        pangraphex = get_pangraphex(osarch) # TODO may be issues when conda is installed as x86 but running on arm64
        pangraphpath = importlib.resources.files("dampa").joinpath("tools/pangraph/"+pangraphex)
        cmd = f"""{pangraphpath} build -s {args.pangraphident} -a {args.pangraphalpha} -b {args.pangraphbeta} -l {args.pangraphlen} -j {args.threads} {filtinput} > {outloc}.json && {pangraphpath} export gfa -o {outloc}_pangenome.gfa {outloc}.json && {pangraphpath} export block-consensus -o {outloc}_pangenome.fa  {outloc}.json"""
    subprocess.run(cmd, shell=True, stdout=pangraph_log, stderr=pangraph_log)
    if os.path.exists(f"{outloc}_pangenome.fa") and os.path.exists(f"{outloc}_pangenome.gfa") and os.path.exists(f"{outloc}.json"):
        logger.info("Pangraph ran successfully")
        if not args.keeplogs:
            os.remove(outloc + "_pangraph.log")
        # TODO improve linear merge (appears to be missing some possible merges)
        linear_transitive_chain_merge_fasta(f"{outloc}_pangenome.gfa",f"{outloc}_pangenome.fa",f"{outloc}_pangenome_lin.fa")# Remove the log file if not keeping logs
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
        pseq = str(probe.seq)
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
    pangraphsettings.add_argument("--pangraphalpha", type=float, default=100,help="Energy cost for splitting a block during alignment merger. Controls graph fragmentation")
    pangraphsettings.add_argument("--pangraphbeta", type=float, default=10,help="Energy cost for diversity in the alignment. A high value prevents merging of distantly-related sequences in the same block")
    pangraphsettings.add_argument("--pangraphlen", type=int, default=90,help="Minimum length of a node to allow in pangenome graph")
    pangraphsettings.add_argument("--pangraphstrict", help="enable the -S strict identity option which limits merges to 1/pangraphbeta divergence",action='store_true')

    probetoolssettings = design.add_argument_group("Probetools settings")

    probetoolssettings.add_argument("--probetoolsidentity", type=int, default=85,
                                    help="Minimum identity in probe match to target to call probe binding")
    probetoolssettings.add_argument("--probetoolsalignmin", type=int, default=90,help="Minimum length (bp) of probe-target binding to allow call of binding")
    probetoolssettings.add_argument("--probetools0covnmin", type=int, default=20,
                                    help="Minimum length (bp) of 0 coverage region in input genomes to trigger design of additional probes")
    probetoolssettings.add_argument("--maxambig",help="The maximum number of ambiguous bases allowed in a probe",type=int,default=10)
    probetoolssettings.add_argument("--nodust", help="Do not run low complexity filter in BLAST (within probetools). If sample has very low GC or is very repetitive this option can be enabled to prevent low complexity regions from being removed",action='store_true')

    mmseqssettings = design.add_argument_group("mmseqs settings")
    mmseqssettings.add_argument("--mmseqs_inputno_trigger",
                                  help="if number of input sequences exceeds this number then mmseqs will be used to deduplcate genomes above 99.9 percent identity",type=int,
                                  default=5000)
    mmseqssettings.add_argument("--mmident", type=float, default=0.999,
                                    help="Minimum identity to cluster genomes")
    mmseqssettings.add_argument("--mmcov", type=float, default=1,
                                    help="Minimum coverage of genomes over which mmident must apply (0-1)")

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

        args.input,overallprops,included = filter_for_nonstandard_inputs(args.input, args.outputfolder,args.maxnonspandard)
        originput = str(args.input)
        rminp = False
        if included > args.mmseqs_inputno_trigger:
            rminp = True
            args.input = mmseqs_subset(args,args.input)
        run_pangraph(args,args.input)#TODO possibly add check where probes are mapped onto each other in a progressive way. each time coverage of a probe by other probes is >1 across full length (at some high identity) then remove probe that is covered would remove lots of similar probes from ends of pancontigs?
        probename = args.outputfolder + "/" + args.outputprefix + "_probes.fasta"
        pangenomefasta = f"{args.outputfolder}/{args.outputprefix}_pangenome_lin.fa"

        split_pangenome_into_probes(pangenomefasta, probename,args.probelen, args.probestep,args.maxambig,probeprefix)

        if not args.skip_padding:
            make_padded_probes(pangenomefasta, probename,args.minlenforpadding,probeprefix=probeprefix)
        if not args.skip_probetoolsfinal:
            finalcaptureout,totalprobes = run_finalprobetools(args,probename,originput)
            if not args.skipsubambig:
                subambig(probename,overallprops)
            if not args.skip_summaries:
                make_summaries(args, finalcaptureout,totalprobes)
        else:
            totalprobes = get_probeno(probename)
            if not args.skipsubambig:
                subambig(probename,overallprops)
            if not args.skip_summaries:
                finalcaptureout = runprobetoolscapture(args,probename)
                make_summaries(args, finalcaptureout,totalprobes)
        if rminp:
            cleanup(args,args.input)
        else:
            cleanup(args)
        logger.info(f"dampa design finished. Total probes: {totalprobes}")
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