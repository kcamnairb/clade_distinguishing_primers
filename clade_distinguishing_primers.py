from primer3 import bindings
from pprint import pprint
from Bio import SeqIO
import argparse
import os.path
import collections
import traceback
import multiprocessing
import subprocess
import pandas as pd
import re
import scipy.stats as ss
import numpy as np

parser = argparse.ArgumentParser(description='Designs primers that will distinguish \
    between clades or user defined groups. If a groups file is not supplied then \
    the SNPs identified by C-Sibelia will be used to create a tree using \
    VCF-kit phylo.')
parser.add_argument('-r', '--reference_fasta', help='Reference genome in \
    fasta format.')
parser.add_argument('-q', '--query_fastas', nargs='+', help='Query genomes \
    in fasta format that will be aligned to the reference genome using Sibelia.')
parser.add_argument('-g', '--groups_file', nargs='?', default=None, help='An \
    optional comma separated file with the first column containing \
    the sample name and the second column containing the group number.')
parser.add_argument('-t', '--num_threads', nargs='?', type=int, 
    default=multiprocessing.cpu_count(), help='Number of threads to use. \
    The default is the total number of cpus.')
args = parser.parse_args()

reference_basename = os.path.splitext(os.path.basename(args.reference_fasta))[0]

def clean_fasta(fasta_in):
    # Ambiguous DNA sequence is converted to an N, otherwise Sibelia will throw an error.
    filename = os.path.basename(fasta_in)
    fasta_out = 'fastas_clean/' + filename
    if not os.path.isfile(fasta_out):
        subprocess.call("sed '/^[^>]/ s/[^AGTC]/N/gi' " + fasta_in + " > " + fasta_out, shell=True)
    return fasta_out
    
def run_sibelia(query_fasta):
    query_basename = os.path.splitext(os.path.basename(query_fasta))[0]
    prefix = 'sibelia/' + query_basename + '_vs_' + reference_basename
    if not os.path.isfile(prefix + '.vcf'):
        command = 'C-Sibelia.py -p 1 -u' + prefix + '_unmapped.fasta ' + \
        '-v ' + prefix + '.vcf ' + reference_fasta_cleaned + ' ' + query_fasta + ' > ' + prefix + '.log'
        subprocess.call(command, shell=True)
    return prefix + '.vcf'
    
if not os.path.exists('sibelia'):
    os.makedirs('sibelia')
if not os.path.exists('fastas_clean'):
    os.makedirs('fastas_clean')
query_fastas_cleaned = [clean_fasta(f) for f in args.query_fastas]
reference_fasta_cleaned = clean_fasta(args.reference_fasta)
pool = multiprocessing.Pool(processes=args.num_threads)
vcfs_sibelia = pool.map(run_sibelia, query_fastas_cleaned)

def run_vcfallelicprimitives(vcf):
    # If multiple gaps or mismatches are present in a line, 
    # this splits the record into primitive alleles on multiple lines.
    base = vcf.split('_vs_')[0]
    name = os.path.basename(base)
    vcf_prim = base + '_prim.vcf'
    command = 'vcfallelicprimitives -kg ' + vcf + ' | vcf-sort > ' + vcf_prim
    subprocess.call(command, shell=True)
    prim_with_gt = base + '_prim_gt.vcf'
    with open(vcf_prim) as in_f, open(prim_with_gt,'w') as out_f:
        for line in in_f:
            if '#CHROM' in line:
                # Add Genotype column so vcfs can be merged in the next step
                out_f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                out_f.write(line.strip() + '\tFORMAT\t' +  name + '\n')
            elif '##' in line:
                out_f.write(line)
            else:
                out_f.write(line.rstrip('\tGT\n') + '\tGT\t1\n')
    return prim_with_gt
    
pool = multiprocessing.Pool(processes=args.num_threads)
vcfs_prim = pool.map(run_vcfallelicprimitives, vcfs_sibelia)

def get_snps(vcf):
    # Outputs just snps in zipped format
    base = os.path.splitext(vcf)[0]
    vcf_split = base + '_snps.vcf.gz'
    command = 'bcftools view -O z --types snps' + ' -o ' + vcf_split + ' ' + vcf
    subprocess.call(command, shell=True)
    subprocess.call('bcftools index ' + vcf_split, shell=True)
    return vcf_split
    
def get_indels(vcf):
    # Outputs just indels
    base = os.path.splitext(vcf)[0]
    vcf_split = base + '_indels.vcf.gz'
    command = 'bcftools view -O z --types indels' + ' -o ' + vcf_split + ' ' + vcf
    subprocess.call(command, shell=True)
    return vcf_split    
    
pool = multiprocessing.Pool(processes=args.num_threads)
vcfs_snps = pool.map(get_snps, vcfs_prim)
pool = multiprocessing.Pool(processes=args.num_threads)
vcfs_indels = pool.map(get_indels, vcfs_prim)
# Merge snp vcfs
merge_command = 'bcftools merge --missing-to-ref -o snps_only_merged.vcf ' + ' '.join(vcfs_snps)
subprocess.call(merge_command, shell=True)
with open('snps_only_merged.vcf') as in_f, open('snps_only_merged_with_ref.vcf','w') as out_f:
    for line in in_f:
        if '#CHROM' in line:
            out_f.write(line.strip() + '\t' + reference_basename + '\n')
        elif '##' in line:
            out_f.write(line)
        else:
            out_f.write(line.strip() + '\t0\n')
convert_command = 'bcftools view -O z -o snps_only_merged_with_ref.vcf.gz snps_only_merged_with_ref.vcf'            
subprocess.call(convert_command, shell=True)            
subprocess.call('vk phylo tree upgma snps_only_merged_with_ref.vcf.gz > tree_upgma.nwk', shell=True)
subprocess.call("perl -p -i -e 's/\n//g' tree_upgma.nwk", shell=True)
if args.groups_file:
    subprocess.call('TreeCluster.py -t 0.145 -tf argmax_clusters -i tree_upgma.nwk -o tree_upgma_clusters.txt', shell=True)
    clusters = pd.read_csv(args.groups_file)
else: 
    clusters = pd.read_csv('tree_upgma_clusters.txt', sep='\t')
clusters = pd.Series(clusters.ClusterNumber.tolist(), index=clusters.SequenceName.tolist())
clusters = clusters[clusters > 0]

def merge_vcfs(vcfs):
    # Performs outer merge on vcfs, and
    # returns dataframe of genotype lengths
    sibelia_merge = pd.DataFrame()
    for f in vcfs:
        isolate = f.replace('sibelia/','').replace('_prim_gt_indels.vcf.gz','')
        if isolate not in clusters:
            continue
        vcf = pd.read_csv(f, comment='#', sep='\t', header=None, compression='gzip')
        vcf.columns = 'CHROM POS ID REF ALT QUAL FILTER INFO'.split() + ['FORMAT', isolate]
        vcf[isolate] = vcf.ALT.str.len() - vcf.REF.str.len()
        vcf = vcf.loc[vcf[isolate].abs() > 0,]
        vcf = vcf.loc[:,['CHROM', 'POS', isolate]]
        if sibelia_merge.shape[0] == 0:
            sibelia_merge = vcf
        else:
            sibelia_merge = sibelia_merge.merge(vcf, "outer")
    sibelia_merge = sibelia_merge.fillna(0)
    if reference_basename in clusters:
        sibelia_merge[reference_basename] = 0
    return sibelia_merge
    
def bin_vcfs(sibelia_merge, bin_size=500):
    # Sums the genotype lengths within each bin
    def bin_group(g):
        g.insert(loc=2, column='bins', value=pd.cut(g.iloc[:,1], range(1, g.POS.max() + 1, bin_size), right=False))
        return g  
    sibelia_merge = sibelia_merge.groupby('CHROM', as_index=False).apply(bin_group)
    sibelia_binned = sibelia_merge.groupby(['CHROM', 'bins']).sum().reset_index()
    diffs = sibelia_binned.loc[:,clusters.index].apply(lambda row: row.max() - row.min(), axis=1)
    sibelia_binned = sibelia_binned.loc[(diffs <= 6000) & (diffs >= 200),:]
    sibelia_binned.to_csv('sibelia_binned.csv', index=False)
    return sibelia_binned
    
sibelia_merge = merge_vcfs(vcfs_indels)
gt_lengths = bin_vcfs(sibelia_merge)
gt_lengths.index = gt_lengths.CHROM + gt_lengths.bins.astype(str)
gt_lengths = gt_lengths.transpose().drop(['CHROM', 'bins'], axis=0)

def custom_round(x, base=5):
    return int(base * round(float(x)/base))

def cramers_v(x, y):
    # Calculates correlation for nominal variables
    #https://stackoverflow.com/questions/20892799/using-pandas-calculate-cram%C3%A9rs-coefficient-matrix
    confusion_matrix = pd.crosstab(x,y)
    chi2 = ss.chi2_contingency(confusion_matrix)[0]
    n = confusion_matrix.sum().sum()
    phi2 = chi2/n
    r,k = confusion_matrix.shape
    phi2corr = max(0, phi2-((k-1)*(r-1))/(n-1))
    rcorr = r-((r-1)**2)/(n-1)
    kcorr = k-((k-1)**2)/(n-1)
    return np.sqrt(phi2corr/min((kcorr-1),(rcorr-1)))

# Round genotype lengths to nearest 40 bp so that small differences don't skew Cramer's V statistic
gt_lengths_rounded = gt_lengths.applymap(lambda x: custom_round(x, base=40))    
gt_lengths_cor = gt_lengths_rounded.apply(lambda col: cramers_v(col, clusters), axis=0)
gt_lengths_best_cor = pd.concat((clusters,gt_lengths.loc[:, gt_lengths_cor > 0.6]),  axis=1)
gt_lengths_best_cor.rename(columns={0:'clusters'}, inplace=True)
gt_lengths_best_cor = gt_lengths_best_cor.sort_values('clusters')
contigs_dict = SeqIO.to_dict(SeqIO.parse(args.reference_fasta, 'fasta'))

def make_primers(region):
    # Designs primers to span a 500 bp region
    contig, start, stop = re.search('(^.+?)\[(\d+), (\d+)\)', region).groups()
    start = int(start)
    stop = int(stop)
    seq = contigs_dict[contig]
    index=['PRIMER_LEFT_0_SEQUENCE', 'PRIMER_RIGHT_0_SEQUENCE', 
           'PRIMER_PAIR_0_PENALTY', 'PRIMER_PAIR_0_PRODUCT_SIZE',
           'amplicon_start', 'amplicon_end']
    try:
        primers = bindings.designPrimers({
            'SEQUENCE_ID': region,
            'SEQUENCE_TEMPLATE': str(seq.seq),
            'SEQUENCE_INCLUDED_REGION':[start - 500, 2100],
            'SEQUENCE_TARGET':[[start - 200, 400]],
            },
            {'PRIMER_PRODUCT_SIZE_RANGE': [[75, 2100]],
            'PRIMER_EXPLAIN_FLAG': 1,
            'PRIMER_MAX_TM': 68,
            'PRIMER_MIN_TM': 52,
            'PRIMER_PICK_INTERNAL_OLIGO': 0})
        return pd.Series([primers['PRIMER_LEFT_0_SEQUENCE'], 
            primers['PRIMER_RIGHT_0_SEQUENCE'], primers['PRIMER_PAIR_0_PENALTY'], 
            primers['PRIMER_PAIR_0_PRODUCT_SIZE'], primers['PRIMER_LEFT_0'][0], 
            primers['PRIMER_RIGHT_0'][0] + primers['PRIMER_RIGHT_0'][1]], index=index)
    except:
        return pd.Series([np.nan] * 6, index=index)
        
pool = multiprocessing.Pool(processes=args.num_threads)
primers = pool.map(make_primers, gt_lengths_best_cor.columns[1:])
primers_df = pd.concat([gt_lengths_best_cor.iloc[:,1:].transpose(),
                     pd.DataFrame.from_records(primers, index=gt_lengths_best_cor.columns[1:])], 
                     axis=1, sort=False).dropna(axis=0)
primers_df.index = [x.replace('[','_').replace(', ','_').replace(')','') for x in primers_df.index]

def make_blast_db(fastas):
    with open('combined.fasta','w') as out_f:
        for fasta in fastas:
            # prepend genome name to each sequence name
            with open(fasta) as in_f:
                basename = os.path.splitext(os.path.basename(fasta))[0]
                name_prepended_seqs = in_f.read().replace('>', '>' + basename + '|')
                out_f.write(name_prepended_seqs)
    subprocess.call('makeblastdb -dbtype nucl -in combined.fasta -out combined_fastas', shell=True)
    
make_blast_db([args.reference_fasta] + args.query_fastas)
with open('primers.fasta', 'w') as out_f:
    for idx, row in primers_df.iterrows():
        out_f.write('>' + idx + '|F\n')
        out_f.write(row['PRIMER_LEFT_0_SEQUENCE'] + '\n')
        out_f.write('>' + idx + '|R\n')
        out_f.write(row['PRIMER_RIGHT_0_SEQUENCE'] + '\n')
# Predict amplicons using simulate_PCR
primer_blast_command = 'simulate_PCR -db combined_fastas -primers primers.fasta \
    -minlen 1 -maxlen 10000 -max_target_seqs 100000 -mux 0 -evalue 2000 \
    -num_threads ' + str(args.num_threads)
subprocess.call(primer_blast_command, shell=True)
sim_pcr = pd.read_csv('primers.fasta.pair.combined_fastas.amplicons', sep='\t') 
sim_pcr['locus'] = sim_pcr.FP_ID.str.extract('(.*)\|')
sim_pcr['isolate'] = sim_pcr.HitName.str.extract('(.*)\|')
sample_amplicon_counts = sim_pcr.groupby(
    ['locus','isolate'], as_index=False)['amplicon_len'].agg('count')
loci_with_multiple_amplicons_in_a_sample = sample_amplicon_counts.loc[
    sample_amplicon_counts.amplicon_len > 1,'locus']
# Remove primers that amplify multiple amplicons within a genome    
sim_pcr = sim_pcr.loc[
    ~sim_pcr.locus.isin(loci_with_multiple_amplicons_in_a_sample.unique()),:]
sim_pcr = sim_pcr.loc[:,['locus','isolate','amplicon_len']].merge(
    primers_df, left_on='locus', right_index=True)
amplicon_lengths = sim_pcr.loc[:,['locus','isolate','amplicon_len']].pivot(
    columns='locus', index='isolate', values='amplicon_len').fillna(0).astype('int32')
amplicon_lengths_rounded = amplicon_lengths.applymap(lambda x: custom_round(x, base=250))
# Correlate rounded amplicon lengths with clades
amplicon_lengths_cor = amplicon_lengths_rounded.apply(
    lambda col: cramers_v(col, clusters), axis=0)
amplicon_lengths_cor = amplicon_lengths_cor.sort_values(ascending=False)
amplicon_lengths = amplicon_lengths.transpose()
amplicon_lengths['cramers_v'] = amplicon_lengths_cor
amplicon_lengths = amplicon_lengths.sort_values('cramers_v', ascending=False)
amplicon_lengths_w_primers = pd.concat(
    [pd.DataFrame(clusters, columns=['clusters']).transpose(),
     amplicon_lengths], axis=0)
amplicon_lengths_w_primers = amplicon_lengths_w_primers.merge(
    primers_df.drop(clusters.index, axis=1), how='left', left_index=True, right_index=True)
amplicon_lengths_w_primers = amplicon_lengths_w_primers.loc[:,
    ['PRIMER_LEFT_0_SEQUENCE','PRIMER_RIGHT_0_SEQUENCE','PRIMER_PAIR_0_PENALTY',
    'PRIMER_PAIR_0_PRODUCT_SIZE','amplicon_start','amplicon_end', 'cramers_v'] +
    list(clusters.sort_values().index)]
amplicon_lengths_w_primers.to_csv('distinguishing_primers.csv')