Clade distinguishing primers
============================

This script takes a set of genomes as input and outputs primers that can distinguish 
between the clades present among the genomes. The following steps are executed in this
script:
  1. All the genomes are individually aligned to a genome designated as the reference using C-Sibelia.
  2. The records in the VCF files are converted to allelic primitives using vcflib. 
  3. The SNPs and indels are separated using bcftools.
  4. The SNPs from all the genomes are merged and a tree is produced using VCF-kit.
  5. Clades present among the tree are identified using TreeCluster.
  6. The lengths of the indels are binned in 500 bp intervals and summed.
  7. Primers are designed to amplify regions that are variable among the genomes.
  8. Amplicon sizes are checked using Simulate_PCR, and primers are discarded that amplify off targets.
  9. Primers are scored by correlating amplicon sizes with the clades using Cramer's V.
