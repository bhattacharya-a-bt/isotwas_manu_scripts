require(data.table)
require(bigsnpr)

setwd('/u/project/gandalm/shared/GenomicDatasets-processed/PEC_AMP_AD_Arjun/Metadata')

rnaseq = fread('rnaseq_20230501.csv')
genotype = fread('genotype.csv')
individual = fread('individual_20230421.csv')

bed = fread('/u/scratch/a/abtbhatt/renormalize_0429/gene_exp_bigbrain_042923_nohcp.bed',
            nrow = 10)


snps = snp_attach(snp_readBed2('/u/scratch/a/abtbhatt/renormalize_0429/bigbrain_hmp.bed',
                               backingfile = tempfile()))

individual = subset(individual,genotypingID %in% snps$fam$sample.ID)
individual = subset(individual,individualID %in% colnames(bed))
individual = individual[!duplicated(individual$individualID),]
snps = snp_attach(subset(snps,
                         ind.row = which(snps$fam$sample.ID %in% 
                                           individual$genotypingID)))
individual = individual[match(snps$fam$sample.ID,
                              individual$genotypingID),]
snps$fam$family.ID = snps$fam$sample.ID = individual$individualID

snp_writeBed(snps,'/u/scratch/a/abtbhatt/renormalize_0429/BB_hmp_forQTL.bed')

### IN BASH
# cd /u/scratch/a/abtbhatt/renormalize_0429
# plink --bfile BB_hmp_forQTL_allSNPs --recode vcf-fid --out BB_hmp_forQTL_allSNPs
# module add htslib
# bgzip -c BB_hmp_forQTL_allSNPs.vcf > BB_hmp_forQTL_allSNPs.vcf.gz
# tabix -p vcf BB_hmp_forQTL_allSNPs.vcf.gz 
# bgzip gene_exp_bigbrain_042923_nohcp.bed && tabix -p bed gene_exp_bigbrain_042923_nohcp.bed.gz
# bgzip tx_exp_bigbrain_042923_nohcp.bed && tabix -p bed tx_exp_bigbrain_042923_nohcp.bed.gz

