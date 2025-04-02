# Dampa
Dampa: Diversity Aware Metagenomic Panel Assignment 

## Usage

### Design Probes

```
python dampa.py design [options]
```

#### Input/Output options
- `-g, --input`: Either folder containing individual genome fasta files OR a single fasta file containing all genomes (files must end in .fna, .fa or .fasta) (required)
- `-c, --clusterassign`: clstr file from cd-hit or tab delimited file with 1st column genonme name and 2nd column cluster name
- `--clustertype`: type of cluster file input cdhit (produced by cdhit) or tabular (genome and cluster tab delimited) (default: tabular)
- `--maxnonspandard`: maximum proportion of genome that can be non ATGC (0-1) (default: 0.01)
- `-o, --outputfolder`: path to output folder (default: current working directory)
- `-p, --outputprefix`: prefix for all output files and folders (default: probebench_run)

#### General settings
- `-l, --probelen`: length of output probes (default: 120)
- `-s, --probestep`: step of probes (for no overlap set to same as probelen) (default: 120)
- `--skipsubambig`: do NOT substitute ambiguous nucleotides (by default N or other ambiguous nucleotides are substituted for ATGC in a random selection weighted by proportions in input genomes)

#### Pangraph settings
- `--pangraphident`: Pangenome percentage identity setting allowable values are 5,10 or 20 (default: 20)
- `--pangraphalpha`: Energy cost for splitting a block during alignment merger. Controls graph fragmentation (default: 100)
- `--pangraphbeta`: Energy cost for diversity in the alignment. A high value prevents merging of distantly-related sequences in the same block (default: 10)
- `--pangraphlen`: Minimum length of a node to allow in pangenome graph (default: 90)

#### Probetools settings
- `--probetoolsidentity`: Minimum identity in probe match to target to call probe binding (default: 85)
- `--probetoolsalignmin`: Minimum length (bp) of probe-target binding to allow call of binding (default: 90)
- `--probetools0covnmin`: Minimum length (bp) of 0 coverage region in input genomes to trigger design of additional probes (default: 20)
- `--probetoolsbatch`: probetools -b option, number of probes to design in each iteration (default: 100)
- `--maxambig`: The maximum number of ambiguous bases allowed in a probe (default: 10)

#### Additional settings
- `--skip_padding`: do not generate additional probes for pangenome nodes between pangraphlen and probelen in length
- `--padding_nuc`: nucleotide to use for padding probes to args.probelen (choices: A, T, C, G) (default: T)
- `--minlenforpadding`: minimum length for a pancontig for it to be padded (default: 90)
- `--skip_probetoolsfinal`: do NOT run final probe design step
- `--genprefix`: generate output prefix from input names and settings
- `-t, --threads`: number of threads (default: 1)
- `--keeplogs`: keep logs containing output from pangraph and probetools
- `--keeptmp`: keep intermediate files from pangraph and probetools
- `--skip_summaries`: do NOT run visualisation generation of dampa probes relative to input genomes
- `--maxdepth_describe`: Maximum depth of probe coverage to describe separately (default: 1)
- `--report0covperc`: threshold above which genomes are reported as having too much of their genome not covered by any probes (default: 1)
- `--version`: print version and exit

### Evaluate Probes

```
python dampa.py eval [options]
```

#### Input/Output options
- `-g, --input`: Genomes to check probe coverage. If genomes either folder containing individual genome fasta files OR a single fasta file containing all genomes (files must end in .fna, .fa or .fasta). If capture file then a pt file from a previous pangraph design or pangraph eval run (required)
- `--inputtype`: type of cluster file input cdhit (produced by cdhit) or tabular (genome and cluster tab delimited) (choices: genomes, capture) (default: genomes)
- `-q, --probes`: Fasta file containing probes to evaluate (files must end in .fna, .fa or .fasta) (required)
- `-c, --clusterassign`: clstr file from cd-hit
- `--clustertype`: type of cluster file input cdhit (produced by cdhit) or tabular (genome and cluster tab delimited) (default: tabular)
- `--filtnonstandard`: remove genomes with non standard nucleotides i.e. not A, T, G, C or N
- `-o, --outputfolder`: path to output folder (default: current working directory)
- `-p, --outputprefix`: prefix for all output files and folders (default: probebench_run)

#### Probetools settings
- `--probetoolsidentity`: Minimum identity in probe match to target to call probe binding (default: 85)
- `--probetoolsalignmin`: Minimum length (bp) of probe-target binding to allow call of binding (default: 90)

#### Additional settings
- `-t, --threads`: number of threads (default: 1)
- `--keeplogs`: keep logs containing output from pangraph and probetools
- `--maxdepth_describe`: Maximum depth of probe coverage to describe separately (default: 1)
- `--report0covperc`: threshold above which genomes are reported as having too much of their genome not covered by any probes (default: 1)
- `--version`: print version and exit
