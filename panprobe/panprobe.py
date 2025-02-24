import random
import shutil
from pathlib import Path
from Bio import SeqIO,SeqRecord,Seq

import subprocess
import argparse
import os
import logging
import time
import sys
from collections import Counter
import platform
import subprocess
from tools.gather_probe_depth_stats import make_stats, make_propsplot, make_meanplot, process_count_pt
from vis.plot_over_genomelen import make_genome_plots,replace_short_zeros




def detect_os_arch():
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


def make_padded_probes(pangenomefa,probefasta,minlen):
    inpangenome = SeqIO.parse(pangenomefa, 'fasta')
    toadd = []
    for pancontig in inpangenome:
        if len(pancontig.seq) < 120:
            rmn = str(pancontig.seq).replace("N","").replace("n","")
            if len(rmn) >= minlen:
                tadd = "T"*(120-len(rmn))
                newseq = Seq.Seq(tadd + rmn)
                outpancontig = SeqRecord.SeqRecord(newseq,pancontig.id+"_padded_probe",description="")
                toadd.append(outpancontig)
    with open(probefasta, "a") as fasta_file:
        SeqIO.write(toadd, fasta_file, "fasta")
    logger.info(f"Padding small sequences generated {len(toadd)} additional probes")

def nucleotide_proportions(sequences):
    overallseq = ""
    for i in sequences:
        overallseq += str(i.seq)
    counts = Counter(overallseq)  # Count occurrences of each nucleotide
    total = sum(counts.values())  # Total nucleotides
    proportions = {nt: count / total for nt, count in counts.items()}# Calculate proportions
    return proportions

def filter_for_nonstandard_inputs(genomes,outfolder,maxnonspandard):
    ingenomes = list(SeqIO.parse(genomes, "fasta"))
    outgenomes = []
    excluded = 0
    overallprops = nucleotide_proportions(ingenomes)
    for i in ingenomes:
        props = nucleotide_proportions([i])
        allowed = ["N","A","C","G","T","n","a","c","g","t"]
        nonallowed = [x for x in i.seq if x not in allowed]
        propnonstandard = sum([props[x] for x in props if x not in allowed])

        i.id = i.id.split(" ")[0]
        i.description = ""
        if float(propnonstandard) < float(maxnonspandard):
            outgenomes.append(i)
        else:
            notallowedstr = ",".join(list(set(nonallowed)))
            logger.info(f"genome {i.id} has non standard chraracters: {notallowedstr} and has been excluded")
            excluded += 1
    outpath = outfolder + "/" + genomes.split("/")[-1].replace(".fasta","").replace(".fa","").replace(".fna","")
    outpath = outpath + "_filt.fasta"
    SeqIO.write(outgenomes, outpath, "fasta")
    logger.info(f"total genomes excluded for due to excess non ATGCN nucleotides: {excluded}")
    return outpath,overallprops

class RuntimeFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None):
        super().__init__(fmt, datefmt)
        self.start_time = time.time()  # Record the start time

    def format(self, record):
        # Calculate elapsed time
        elapsed_time = time.time() - self.start_time
        record.runtime = f"{elapsed_time:.2f}s"  # Add runtime to the record
        return super().format(record)

def make_prefix(args):

    if args.skip_padding:
        pad = "_nopadding"
    else:
        pad = ""
    if args.skip_probetoolsfinal:
        ptf = "_noprobetools"
    else:
        ptf = ""

    outname = args.input.split("/")[-1].replace(".fasta","").replace(".fa","").replace(".fna","")

    outprefix = f"{outname}_probelen{args.probelen}_probestep{args.probestep}_i{args.pangraphident}_a{args.pangraphalpha}_b{args.pangraphbeta}_{args.pangraphlen}min{pad}{ptf}"

    return outprefix

def get_pangraphex(osarch):
    if osarch["os"] == "linux":
        if osarch["arch"] == "arm64":
            if osarch["libc"] == "gnu":
                return "tools/pangraph/pangraph-aarch64-linux-gnu"
            elif osarch["libc"] == "musl":
                return "tools/pangraph/pangraph-aarch64-linux-musl"
            else:
                logger.error(f"Unsupported os/architecture/lobc combination {osarch["os"]}/{osarch['arch']}/{osarch['libc']}")
        elif osarch["arch"] == "x86_64":
            if osarch["libc"] == "gnu":
                return "tools/pangraph/pangraph-x86_64-linux-gnu"
            elif osarch["libc"] == "musl":
                return "tools/pangraph/pangraph-x86_64-linux-musl"
            else:
                logger.error(f"Unsupported os/architecture/lobc combination {osarch["os"]}/{osarch['arch']}/{osarch['libc']}")
    elif osarch["os"] == "darwin":
        if osarch["arch"] == "arm64":
            return "tools/pangraph/pangraph-aarch64-darwin"
        elif osarch["arch"] == "x86_64":
            return "tools/pangraph/pangraph-x86_64-darwin"
        else:
            logger.error(f"Unsupported os/architecture combination {osarch["os"]}/{osarch['arch']}")
    elif osarch["os"] == "windows":
        return "tools/pangraph/pangraph-x86_64-pc-windows-gnu.exe"
    else:
        logger.error(f"Unsupported os/architecture/lobc combination {osarch["os"]}/{osarch['arch']}/{osarch['libc']}")



def run_pangraph(args):
    logger.info("Starting pangenome generation with PanGraph")
    topdir = os.path.dirname(os.path.abspath(__file__))
    outloc = f"{args.outputfolder}/{args.outputprefix}"
    pangraph_log = open(outloc + "_pangraph.log", "w")

    osarch = detect_os_arch()
    pangraphex = get_pangraphex(osarch)
    cmd = f"""{topdir}/{pangraphex} build -s {args.pangraphident} -a {args.pangraphalpha} -b {args.pangraphbeta} -l {args.pangraphlen} -j {args.threads} {args.input} > {outloc}.json && {topdir}/tools/pangraph export gfa -o {outloc}_pangenome.gfa {outloc}.json && {topdir}/tools/pangraph export block-consensus -o {outloc}_pangenome.fa  {outloc}.json"""
    subprocess.run(cmd, shell=True, stdout=pangraph_log, stderr=pangraph_log)
    if os.path.exists(f"{outloc}_pangenome.fa") and os.path.exists(f"{outloc}_pangenome.gfa") and os.path.exists(f"{outloc}.json"):
        logger.info("Pangraph ran successfully")
        if not args.keeplogs:
            os.remove(outloc + "_pangraph.log")
    else:
        logger.error(f"One or more of pangraph outputs ({args.outputprefix}_pangenome.gfa, {args.outputprefix}_pangenome.fa, {args.outputprefix}.json) in {args.outputfolder} are not present. Check for error in pangraph log")
    return

def run_finalprobetools(args,inprobes):
    outloc = f"{args.outputfolder}/{args.outputprefix}"
    finalpref = f"{outloc}_probetools_final"
    finalprobes= f"{finalpref}_probes.fa"
    topdir = os.path.dirname(os.path.abspath(__file__))
    # cmd = f"python {current_directory}/tools/probetools/probetools_v_0_1_11.py makeprobes -t /Users/mpay0321/Dropbox/Probe_design_project/2025-02-14_add_probetools_makeprobeswinput/cluster0_pangenome_kminimap2_a200_b12.5_s20_90min_rep1_179.fasta -b 10 -o /Users/mpay0321/Dropbox/Probe_design_project/2025-02-14_add_probetools_makeprobeswinput/rep1_probetoolscomplete -c 100 -l 90 -L 20 -T 10"
    probetools_log = open(outloc + "_probetools.log","w")
    cmd = f"python {topdir}/tools/probetools/probetools_v_0_1_11.py makeprobeswinput -t {args.input} -b {args.probetoolsbatch} -x {inprobes} -o {finalpref} -i {args.probetoolsidentity} -l {args.probetoolsalignmin} -T {args.threads} -L {args.probetools0covnmin} -c 100 -d {args.maxambig}"

    subprocess.run(cmd, shell=True, stdout=probetools_log, stderr=probetools_log)


    if os.path.exists(finalprobes):
        shutil.move(finalprobes, inprobes)
        totaprobes, addedprobes = get_probeno(inprobes, "_round_")
        logger.info(f"Generated {addedprobes} additional probes using probetools")
        if not args.keeplogs:
            os.remove(outloc + "_probetools.log")
        return f"{finalpref}_capture.pt",totaprobes
    else:
        logger.error(f"Probetools output file {finalprobes} not present. Check for error in details probetools log.")

def get_clusters(args):
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
    outloc = args.outputfolder+"/"+args.outputprefix
    current_directory = os.path.dirname(os.path.abspath(__file__))
    capture_log = open(outloc + "_capture.log", "w")
    with open(os.devnull, 'w') as devnull:
        # outf = open("/Users/mpay0321/Dropbox/Probe_design_project/2025-01-29_integrate_probetools_probebench/stdout.txt",'w')
        cmd = f"python {current_directory}/tools/probetools/probetools_v_0_1_11.py capture -t {args.input} -p {probes} -o {outloc} -i {args.probetoolsidentity} -l {args.probetoolsalignmin} -T {args.threads}"
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
    ambig = len([x for x in seq if x not in ["A","T","G","C","a","t","c","g"]])
    return ambig

def split_pangenome_into_probes(input_fasta, output_fasta, probe_length,probe_step,maxambig):
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
                    piece_record = SeqRecord.SeqRecord(Seq.Seq(piece),id=piece_id,description="")
                    outprobes.append(piece_record)
                    totalprobes+=1
    # Write the probes to the output file
    SeqIO.write(outprobes, output_fasta, "fasta")
    logger.info(f"Extracted {totalprobes} probes from pangenome")

def make_summaries(args,ptcountfile,totalprobes):
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

    debug = False
    if debug:
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        args = parser.parse_args()
        ## DESIGN
        # args.input = "/Users/mpay0321/Dropbox/Probe_design_project/2025-02-17_panprobe_w_probetools_results/flu/ref/H3N2_complete_rename.fasta"
        # args.clusterassign = False
        # args.threads = 10
        # args.maxnonspandard =0.01
        # args.outputfolder="/Users/mpay0321/Dropbox/Probe_design_project/2025_02_19_probetools_dontcountN/"
        # args.outputprefix="flu"
        # args.probelen = 120
        # args.probestep = 120
        # args.skipsubambig = False
        # args.pangraphident = 20
        # args.pangraphalpha = 100
        # args.pangraphbeta = 10
        # args.pangraphlen = 90
        # args.probetoolsidentity = 85
        # args.probetoolsalignmin = 90
        # args.probetools0covnmin = 20
        # args.probetoolsbatch = 100
        # args.maxambig = 10
        # args.run_padding = True
        # args.padding_nuc = 'T'
        # args.minlenforpadding = 90
        # args.skip_probetoolsfinal = False
        # args.genprefix = False
        # args.threads = 10
        # args.keeplogs = False
        # args.skip_summaries = False
        # args.maxdepth_describe = 1
        # args.report0covperc = 1
        # args.version = False
        # args.command = "design"

        ### EVAL

        args.input = "/Users/mpay0321/Dropbox/Probe_design_project/2025-02-17_panprobe_w_probetools_results/dengue/dengue_catch_capture.pt"
        args.inputtype = "capture"
        args.probes = "/Users/mpay0321/Dropbox/Probe_design_project/2025-02-17_panprobe_w_probetools_results/dengue/dengue_catch.fasta"
        args.clusterassign = False
        args.clustertype="tabular"
        args.threads = 10
        args.maxnonspandard =0.01
        args.outputfolder="/Users/mpay0321/Dropbox/Probe_design_project/2025_02_19_probetools_dontcountN/"
        args.outputprefix="denguetest"
        args.filtnonstandard=True
        args.probetoolsidentity = 85
        args.probetoolsalignmin = 90
        args.keeplogs = False
        args.maxdepth_describe = 1
        args.report0covperc = 1
        args.version = False
        args.command = "eval"

        return args
    else:
        cwd = Path.cwd()
        parser = argparse.ArgumentParser(description="PanProbe - probe panel design using pangenome graphs",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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

        probetoolssettings = design.add_argument_group("Probetools settings")

        probetoolssettings.add_argument("--probetoolsidentity", type=int, default=85,
                                        help="Minimum identity in probe match to target to call probe binding")
        probetoolssettings.add_argument("--probetoolsalignmin", type=int, default=90,help="Minimum length (bp) of probe-target binding to allow call of binding")
        probetoolssettings.add_argument("--probetools0covnmin", type=int, default=20,
                                        help="Minimum length (bp) of 0 coverage region in input genomes to trigger design of additional probes")
        probetoolssettings.add_argument("--probetoolsbatch", type=int, default=100,
                                        help="probetools -b option, number of probes to design in each iteration")
        probetoolssettings.add_argument("--maxambig",help="The maximum number of ambiguous bases allowed in a probe",type=int,default=10)


        additionalsettings = design.add_argument_group("Additional settings")

        additionalsettings.add_argument("--run_padding",
                            help="generate additional probes for pangenome nodes between pangraphlen and probelen in length. i.e. if padding is run 30*T would be added to the end of a 90bp pancontig",action='store_true')
        additionalsettings.add_argument("--padding_nuc",
                            help="nucleotide to use for padding probes to args.probelen", choices=['A',"T","C","G"],
                            default='T')
        additionalsettings.add_argument("--minlenforpadding",
                            help="minimum length for a pancontig for it to be padded (WARNING setting this below ~80 may result in probes that do not effectively bind, leave these small sequences for final probetools step)", type=int,
                            default='90')
        additionalsettings.add_argument("--skip_probetoolsfinal",
                            help="do NOT run final probe design step. i.e. this step uses probetools to design probes to regions that are not represented in the pangenome",action='store_true')

        additionalsettings.add_argument("--genprefix",
                            help="generate output prefix from input names and settings",action='store_true')

        additionalsettings.add_argument("-t","--threads",
                                help="number of threads",
                                type=int,
                                default=1)
        additionalsettings.add_argument("--keeplogs",
                            help="keep logs containing output from pangraph and probetools",action='store_true')
        additionalsettings.add_argument("--skip_summaries",
                            help="do NOT run visualisation generaton of panprobe probes relative to input genomes",action='store_true')
        additionalsettings.add_argument("--maxdepth_describe",
                            default=1,help="Maximum depth of probe coverage to describe separately. i.e. if 1 there will be 0,1 and >1 depth categories")
        additionalsettings.add_argument("--report0covperc",
                            help="threshold above which genomes are reported as having too much of their genome not covered by any probes",type=float,default=1)
        additionalsettings.add_argument("--version",
                                        help="print version and exit",
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

        args = parser.parse_args()

        return args

def main():
    """
    get args
    setup outputs dir and tmp folder (use uuid)
    optionally run splitfasta function
    run
    """
    version = "0.1.0"
    global logger
    logger = setup_logging()
    logger.info("Starting")

    args = get_args()

    if args.command == "design":
        if args.version:
            print(f"version {version}")
            sys.exit(0)
        logger.info("Running panprobe design")
        if args.genprefix:
            args.outputprefix = make_prefix(args)
        args.filtnonstandard = True

        logger.info("Filter genomes with too many non standard nucleotides")
        args.input,overallprops = filter_for_nonstandard_inputs(args.input, args.outputfolder,args.maxnonspandard)

        run_pangraph(args)#TODO possibly add check where probes are mapped onto each other in a progressive way. each time coverage of a probe by other probes is >1 across full length (at some high identity) then remove probe that is covered would remove lots of similar probes from ends of pancontigs?
        probename = args.outputfolder + "/" + args.outputprefix + "_probes.fasta"
        pangenomefasta = f"{args.outputfolder}/{args.outputprefix}_pangenome.fa"

        split_pangenome_into_probes(pangenomefasta, probename,args.probelen, args.probestep,args.maxambig)

        if  args.run_padding:
            make_padded_probes(pangenomefasta, probename,args.minlenforpadding)
        if not args.skip_probetoolsfinal:
            finalcaptureout,totalprobes = run_finalprobetools(args,probename)
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
        logger.info(f"panprobe design finished. Total probes: {totalprobes}")
    elif args.command == "eval":
        if args.version:
            print(f"version {version}")
            sys.exit(0)
        logger.info("Running panprobe eval")
        totalprobes = get_probeno(args.probes)[0]
        if args.inputtype != "capture":
            finalcaptureout = runprobetoolscapture(args, args.probes)
        else:
            finalcaptureout = args.input
        make_summaries(args, finalcaptureout,totalprobes)
        logger.info("panprobe eval finished")
if __name__ == "__main__":
    main()