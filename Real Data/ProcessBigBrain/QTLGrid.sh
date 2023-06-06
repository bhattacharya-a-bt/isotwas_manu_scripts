#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=30G,h_rt=12:00:00
#$ -o /u/scratch/a/abtbhatt/renormalize_0429/outFiles
#$ -t 100-300:25

module add qtltools/1.3.1 
cd /u/scratch/a/abtbhatt/renormalize_0429

k=${SGE_TASK_ID}

QTLtools cis --vcf BB_hmp_forQTL.vcf.gz \
--bed gene_exp_bigbrain_042923_nohcp.bed.gz \
--cov cov_hcp${k}_gene.txt --nominal 0.0001 --normal \
--std-err \
--out hcp${k}_cis_gene.txt 

QTLtools cis --vcf BB_hmp_forQTL.vcf.gz \
--bed tx_exp_bigbrain_042923_nohcp.bed.gz \
--cov cov_hcp${k}_isoform.txt --nominal 0.0001 --normal \
--std-err \
--out hcp${k}_cis_isoform.txt 



module add qtltools/1.3.1 
cd /u/scratch/a/abtbhatt/renormalize_0429

k=175

QTLtools cis --vcf BB_hmp_forQTL_allSNPs.vcf.gz \
--bed gene_exp_bigbrain_042923_nohcp.bed.gz \
--cov cov_hcp${k}_gene.txt --nominal 1.0 --normal \
--std-err \
--out /u/project/pasaniuc/abtbhatt/hcp${k}_cis_gene_all.txt 

k=150

QTLtools cis --vcf BB_hmp_forQTL_allSNPs.vcf.gz \
--bed tx_exp_bigbrain_042923_nohcp.bed.gz \
--cov cov_hcp${k}_isoform.txt --nominal 1.0 --normal \
--std-err \
--out /u/project/pasaniuc/abtbhatt/hcp${k}_cis_isoform_all.txt 