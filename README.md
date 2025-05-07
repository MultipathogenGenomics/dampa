# Dampa
Dampa: Diversity Aware Metagenomic Panel Assignment 

Targeted metagenomics is an approach that uses sets of oligonucleotide probes (up to 100,000s) to selectively sequence targeted loci or complete genomes from hundreds of pathogen species, removing the need to identify or isolate the pathogen first.
When designing a probeset for targeted metagenomics, the optimal solution includes the minimum number of probes required to effectively enrich the desired set of target organisms or loci. DAMPA (Diversity Aware Metagenomic Panel Assignment), is a user-friendly probeset design pipeline that combines novel target selection algorithms with existing probe design tools to minimise probe number while maintaining unbiased coverage of target genomes. The key approach used by DAMPA is the generation of a pangenome graph to describe the diversity of targeted loci. A pangenome graph is a representation of all unique sequences in a species, or in this case, a dataset. Crucially, conserved regions are collapsed while diverse regions are separated. The graph therefore presents the ideal target for probe design as it will ensure that probes have sufficient identity to bind to all genomic regions regardless of their diversity and natively incorporates support for recombination. DAMPA applies existing probe design tools to the graph to design probes while also ensuring all regions of input genomes are covered.

## Dependencies
- python >=3.10
- biopython >=1.84
- pandas >=2.0
- matplotlib >= 3.10
- seaborn >=0.13
- blast >=2.16
- vsearch >=2
  
## Installation

`conda install -c bioconda -c conda-forge dampa`

[![Anaconda-Server Badge](https://anaconda.org/bioconda/dampa/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/dampa/badges/downloads.svg)](https://anaconda.org/bioconda/dampa)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/dampa/badges/version.svg)](https://anaconda.org/bioconda/dampa)

## Usage

### Design Probes
dampa design is the core module of dampa used to generate probe sets
```
dampa design [options]
```

#### Input/Output options
- `-g, --input`: Either folder containing individual genome fasta files OR a single fasta file containing all genomes (files must end in .fna, .fa or .fasta) (required)
- `-c, --clusterassign`: clstr file from cd-hit or tab delimited file with 1st column genome name and 2nd column cluster name
- `--clustertype`: type of cluster file input cdhit (produced by cdhit) or tabular (genome and cluster tab delimited) (default: tabular)
- `--maxnonspandard`: maximum proportion of genome that can be non ATGC (0-1) (default: 0.01)
- `-o, --outputfolder`: path to output folder (default: current working directory)
- `-p, --outputprefix`: prefix for all output files and folders (default: dampa_run)

#### General settings
- `-l, --probelen`: length of output probes (default: 120)
- `-s, --probestep`: step of probes (for no overlap set to same as probelen) (default: 120)
- `--skipsubambig`: do NOT substitute ambiguous nucleotides (by default N or other ambiguous nucleotides are substituted for ATGC in a random selection weighted by proportions in input genomes)

#### Pangraph settings - these modify the generation of the underlying pangenome graph used to generate probes
- `--pangraphident`: Pangenome percentage identity setting allowable values are 5,10 or 20 (default: 20)
- `--pangraphalpha`: Energy cost for splitting a block during alignment merger. Controls graph fragmentation (default: 100)
- `--pangraphbeta`: Energy cost for diversity in the alignment. A high value prevents merging of distantly-related sequences in the same block (default: 10)
- `--pangraphlen`: Minimum length of a node to allow in pangenome graph (default: 90)

#### Probetools settings - these modify the settings for probetools which is used to fill in gaps of genome coverage after pangenome probe design
- `--probetoolsidentity`: Minimum identity in probe match to target to call probe binding (default: 85)
- `--probetoolsalignmin`: Minimum length (bp) of probe-target binding to allow call of binding (default: 90)
- `--probetools0covnmin`: Minimum length (bp) of 0 coverage region in input genomes to trigger design of additional probes (default: 20)
- `--maxambig`: The maximum number of ambiguous bases allowed in a probe (default: 10)
- `--skip_probetoolsfinal`: do NOT run final probe design step (only produce probes purely from the pangenome)
- 
#### small seq padding settings
A pangraph node of length 110 will be padded with 10 Ts to make it 120bp
- `--skip_padding`: do not generate additional probes for pangenome nodes between pangraphlen and probelen in length (i.e. if pangraphlen is 100 abd probelen is 120 a pangraph node of length 110 will be padded with 10 Ts to make it 120bp)
- `--padding_nuc`: nucleotide to use for padding probes to args.probelen (choices: A, T, C, G) (default: T)
- `--minlenforpadding`: minimum length for a pancontig for it to be padded (default: 90)
- 
#### summary descriptions settings
- `--skip_summaries`: do NOT run visualisation generation of dampa probes relative to input genomes
- `--maxdepth_describe`: Maximum depth of probe coverage to describe separately (default: 1)
- `--report0covperc`: threshold above which genomes are reported as having too much of their genome not covered by any probes (default: 1)

#### Additional settings
- `-t, --threads`: number of threads (default: 1)
- `--keeplogs`: keep logs containing output from pangraph and probetools
- `--keeptmp`: keep intermediate files from pangraph and probetools
- `--version`: print version and exit

### Evaluate Probes

dampa eval will evaluate the performance of a probe set against a set of genomes and generate summary statistics and visualisations

```
dampa eval [options]
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
