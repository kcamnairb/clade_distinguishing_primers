Clade distinguishing primers
============================

This script takes a set of genomes as input and outputs primers that can distinguish 
between the clades present among the genomes. The following steps are executed in this
script:
  1. All the genomes are individually aligned to a genome designated as the reference using [C-Sibelia](https://github.com/bioinf/Sibelia).
  2. The records in the VCF files are converted to allelic primitives using [vcflib](https://github.com/vcflib/vcflib). 
  3. The SNPs and indels are separated using [bcftools](https://github.com/samtools/bcftools).https://github.com/primer3-org/primer3
  4. The SNPs from all the genomes are merged and a tree is produced using [VCF-kit](https://github.com/AndersenLab/VCF-kit).
  5. Clades present among the tree are identified using [TreeCluster](https://github.com/niemasd/TreeCluster).
  6. The lengths of the indels are binned in 500 bp intervals and summed.
  7. Primers are designed using [Primer3](https://github.com/primer3-org/primer3) to amplify regions that are variable among the genomes.
  8. Amplicon sizes are checked using [Simulate_PCR](https://sourceforge.net/projects/simulatepcr/), and primers are discarded that amplify off targets.
  9. Primers are scored by correlating amplicon sizes with the clades using Cramer's V.

## Usage
```bash
usage: clade_distinguishing_primers.py [-h] [-r REFERENCE_FASTA]
                                       [-q QUERY_FASTAS [QUERY_FASTAS ...]]
                                       [-g [GROUPS_FILE]] [-t [NUM_THREADS]]

Designs primers that will distinguish between clades or user defined groups.
If a groups file is not supplied then the SNPs identified by C-Sibelia will be
used to create a tree using VCF-kit phylo.

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE_FASTA, --reference_fasta REFERENCE_FASTA
                        Reference genome in fasta format.
  -q QUERY_FASTAS [QUERY_FASTAS ...], --query_fastas QUERY_FASTAS [QUERY_FASTAS ...]
                        Query genomes in fasta format that will be aligned to
                        the reference genome using Sibelia.
  -g [GROUPS_FILE], --groups_file [GROUPS_FILE]
                        An optional comma separated file with the first column
                        containing the sample name and the second column
                        containing the group number.
  -t [NUM_THREADS], --num_threads [NUM_THREADS]
                        Number of threads to use. The default is the total
                        number of cpus.
```
