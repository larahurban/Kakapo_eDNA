#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 19:16:09 2022

@author: lurban
"""


## 00. porechop v0.2.4 and Nanofilt per fastq_pass folder in nanopore folder 
# -> fuses reads into one fastq file
porechop -i ./fastq_pass -b ./porechop0.fastq

# run NanoFilt v2.6: :
NanoFilt porechop0.fastq -q 7 > porechop.fastq

# raw data at this stage is provided via NCBI SRA (NCBI BioProject ID PRJNA810861)


# 0. align reads using minimap2 v2.17 and process with samtools v1.13
minimap2 -ax map-ont NCBI_Strigops_habroptilus_Reference.fasta porechop.fastq > porechop.sam
samtools view -S -b porechop.sam > porechop.bam
samtools sort -o porechop_sorted.bam porechop.bam
samtools view -b -F 4 porechop_sorted.bam > porechop_mapped.bam
samtools index porechop_mapped.bam


# 1. call variants using medaka v1.2.5
# NCBI_Strigops_habroptilus_Reference.fasta can be obtained from NCBI: https://www.ncbi.nlm.nih.gov/genome/?term=txid2489341[orgn]
medaka_variant -f NCBI_Strigops_habroptilus_Reference.fasta -i porechop_mapped.bam -P 0 -o porechop_mapped
medaka snp NCBI_Strigops_habroptilus_Reference.fasta ./porechop_mapped/round_1_hap_1_probs.hdf ./porechop_mapped/round_1_hap_2_probs.hdf ./porechop_mapped/round_1_gvcf.vcf

# index vcf and gvcf files using tabix v0.2.6

bgzip -c round_1.vcf > round_1.vcf.gz
tabix -f -p vcf round_1.vcf.gz

bgzip -c round_1_gvcf.vcf > round_1_gvcf.vcf.gz
tabix -f -p vcf round_1_gvcf.vcf.gz



# 3. analyse variant callsets (vcf.gz and vcf.gz.tbi files) in python and compare
# with existing population genomic variant callset for individual identification
import pysam
from pysam import VariantFile
from tqdm import tqdm
import h5py
import pandas as pd

# this population-wide genomic variant callset has been produced by applying
# DeepVariant to ~30x coverage Illumina sequencing data of all extant kakapo
# on the sampled island, Whenua Hou
bcf_in = VariantFile("Kakapo_population_callset.vcf.gz")

# load the medaka vcf and gvcf variant callset (apply the following procedure to both)
bcf_in1 = VariantFile("porechop_mapped/round_1.vcf.gz")
bcf_in2 = VariantFile("porechop_mapped/round_1_gvcf.vcf.gz") 

# a. first process bcf_in1
df0 = pd.DataFrame(columns=['contig','pos','ref','alt'])

# overlap variant positions between population genomic callset and medaka variants
ct = 0
for rec in bcf_in1.fetch():
    df0.loc[ct] = np.nan
    df0['pos'][ct] = rec.pos
    df0['contig'][ct] = rec.contig
    df0['ref'][ct] = rec.alleles[0]
    df0['alt'][ct] = rec.alleles[1]
    ct += 1
    
df1 = pd.DataFrame(columns=['contig','pos'])
contigsu = np.unique(df0.contig)

# overlap alleles between population genomic callset and medaka variants
ct = 0
for chromosome in contigsu:
    print(ct)
    positions = df0[df0.contig==chromosome].pos.values
    for rec in bcf_in.fetch(chromosome):
        if rec.pos in positions:
            if len(rec.alleles) == 2:
                if (rec.alleles[1] == df0[(df0.contig==chromosome)&(df0.pos==rec.pos)]['alt'].values[0]):      
                    df1.loc[ct] = np.nan
                    df1['pos'][ct] = rec.pos
                    df1['contig'][ct] = rec.contig
                    print(rec.pos)
                    print(rec.contig)
                    ct += 1

df1 = df1.drop_duplicates() 


# b. first process bcf_in2 in the same way
df0 = pd.DataFrame(columns=['contig','pos','ref','alt'])

ct = 0
for rec in bcf_in2.fetch():
    df0.loc[ct] = np.nan
    df0['pos'][ct] = rec.pos
    df0['contig'][ct] = rec.contig
    df0['ref'][ct] = rec.alleles[0]
    df0['alt'][ct] = rec.alleles[1]
    ct += 1
    
df2 = pd.DataFrame(columns=['contig','pos'])
contigsu = np.unique(df0.contig)

# overlap alleles between population genomic callset and medaka variants
ct = 0
for chromosome in contigsu:
    print(ct)
    positions = df0[df0.contig==chromosome].pos.values
    for rec in bcf_in.fetch(chromosome):
        if rec.pos in positions:
            if len(rec.alleles) == 2:
                if (rec.alleles[1] == df0[(df0.contig==chromosome)&(df0.pos==rec.pos)]['alt'].values[0]):      
                    df2.loc[ct] = np.nan
                    df2['pos'][ct] = rec.pos
                    df2['contig'][ct] = rec.contig
                    print(rec.pos)
                    print(rec.contig)
                    ct += 1

df2 = df2.drop_duplicates() 


# c. fuse df1 and df2 
# these two data frames will contain highly overlapping variants already, but 
# overlapping them will allow us to include all variation and all indels
df = pd.concat([df1,df2])
df.index = range(0,df.shape[0])
df = df.drop_duplicates()


# d. overlap genotypes

summ = pd.DataFrame(index=rec.samples.keys())

for index, row in df.iterrows():
    col = row['contig']+'_'+str(row['pos'])
    for rec in bcf_in.fetch(row['contig'], row['pos']-1, row['pos']):
        gts = [s['GT'] for s in rec.samples.values()]
    summ[col] = np.nan
    for i in range(0, len(gts)):
        if not((gts[i][0]==None) | (gts[i][1]==None)):
            summ[col][i] = np.sum(gts[i])
        
summ.loc['soil'] = np.nan

for index, row in df.iterrows():
    col = row['contig']+'_'+str(row['pos'])
    #if exists()
    try:
        for rec in bcf_in1.fetch(row['contig'], row['pos']-1, row['pos']):
            gts = [s['GT'] for s in rec.samples.values()]
    except:
        for rec in bcf_in2.fetch(row['contig'], row['pos']-1, row['pos']):
            gts = [s['GT'] for s in rec.samples.values()]
    if not((gts[0][0]==None) | (gts[0][1]==None)):
        summ[col].loc['soil'] = np.sum(gts)
        
        
# look at genotypes from soil
for index, row in df.iterrows():
    col = row['contig']+'_'+str(row['pos'])
    #if exists()
    try:
        for rec in bcf_in1.fetch(row['contig'], row['pos']-1, row['pos']):
            gts = [s['GT'] for s in rec.samples.values()]
    except:
        for rec in bcf_in2.fetch(row['contig'], row['pos']-1, row['pos']):
            gts = [s['GT'] for s in rec.samples.values()]
    print(gts)


# e. get haplotype information 

haps = h5py.File('porechop_mapped/round_1_hap_1_probs.hdf', 'r')

reads = pd.DataFrame(list(haps['samples/data'].keys()))

reads['chr'] = [x.split(':')[0] for x in reads['0'].values]
reads['start'] = [x.split(':')[1].split('-')[0].split('.')[0] for x in reads['0'].values]
reads['end'] = [x.split(':')[1].split('-')[1].split('.')[0] for x in reads['0'].values]
reads['start'] = reads['start'].astype('int')
reads['end'] = reads['end'].astype('int')
del reads['0']

df['reads'] = 0
for i in range(0,df.shape[0]):
    for j in range(0,reads.shape[0]):
        if reads.chr.values[j] == df.contig.values[i]:
            if ((df.pos.values[i] >= reads.start.values[j]) & (df.pos.values[i] <= reads.end.values[j])):
                df['reads'][i] = j

readsi = np.unique(df.reads)

for index, row in df.iterrows():
    col = row['contig']+'_'+str(row['pos'])
    for rec in bcf_in.fetch(row['contig'], row['pos']-1, row['pos']):
        gts = [s['GT'] for s in rec.samples.values()]
    summ[col] = np.nan
    for i in range(0, len(gts)):
        if not((gts[i][0]==None) | (gts[i][1]==None)):
            summ[col][i] = np.sum(gts[i])
       
summ_reads = pd.DataFrame(index=summ_r.index.values)
ct = 0

for reads in readsi:
    dfi = df[df.reads==reads]
    summ_reads[str(reads)] = np.nan
    summ_ph = pd.DataFrame(index=summ_r.index.values)
    for i in dfi.pos.values: 
        for rec in bcf_in.fetch(dfi['contig'].values[0], i-1, i):
            gts = [s['GT'] for s in rec.samples.values()]
            summ_ph[str(i)] = gts
    for j in range(0,summ_ph.shape[0]):
        for i in range(0,summ_ph.shape[1]):  
            if not((summ_ph.iloc[j,i][0]==None) | (summ_ph.iloc[j,i][1]==None)):
                if np.sum(summ_ph.iloc[j,i]) == 0:
                    summ_reads[str(reads)][j] = 0
                    break
                else:
                    summ_reads[str(reads)][j] = 1

                               
# f. haplotype agreement score
                               
# calculate haplotype agreement scores
summ_reads_r = pd.DataFrame(index=summ_reads.index, columns=['ov']) 
for i in summ_reads.index:
    summi = summ_reads.loc[i]
    summi = summi.dropna()
    summ_reads_r.loc[i] = np.sum(summi==0)/summi.shape[0]

# plot haplotype agreement scores
fig, ax = plt.subplots(figsize=(8, 4), facecolor='white')
plt.hist(1-summ_reads_r.ov.values,bins=30, color='teal',histtype='bar', ec='powderblue')
ax.set_xlabel('haplotype agreement [%]'); ax.set_ylabel('number of individuals')
ax.spines['right'].set_visible(False); ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left'); ax.xaxis.set_ticks_position('bottom')
                        
