### This script takes in an isoTWAS model file,
### a GWAS summary statistics file, and calculates gene-level and
### isoform-specific trait associations using step-wise association testing.
## This script is built to be used
### in a job array through SGE or SLURM.

library("optparse")
option_list <- list(
  make_option(c("-f", "--file"), action="store_true", default=TRUE,
              help="file for model",type='character')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

filename = opt$file[1]

require(data.table)
temp_gwas = '/pine/scr/a/b/abhattac/GTEx_isotwas/summStats/temp'
dir.create(temp_gwas)
twas_results = '/pine/scr/a/b/abhattac/PE_fetal/SCZ_results'
dir.create(twas_results)
model = fread(filename)
sumStat = '/pine/scr/a/b/abhattac/GTEx_isotwas/summStats/PGC3_SCZ_wave3_public.v2.tsv'
tissue = 'Fetal_Brain'
gene = unique(model$Gene)[1]

chr = unique(model$Chromosome)[1]
min = min(model$Position)
max = max(model$Position)
system(paste0("awk '$1==",chr,"' ",sumStat,' > ',
              paste0(temp_gwas,'/',gene,'_gwas.tsv')))

require(data.table)
gwas = fread(paste0(temp_gwas,'/',gene,'_gwas.tsv'))
colnames(gwas) = colnames(fread(sumStat,nrow=2))
gwas = gwas[,c('CHR','SNP','BP','A1','A2','OR','SE','P')]
gwas$BETA = log(gwas$OR)
gwas$ChrPos = paste(gwas$CHR,gwas$BP,sep=':')

require(bigsnpr)
system(paste0('plink --bfile ',paste0('/pine/scr/a/b/abhattac/1000GP_Phase3/plink_eur/chr',chr),
              ' --chr ',chr,
              ' --from-bp ',max(0,min(min)-5e6),
              ' --to-bp ',max(max)+5e6,
              ' --make-bed --out ',file.path(temp_gwas,gene)))
ld_bed = paste0(file.path(temp_gwas,gene),'.bed')
ld_snp = snp_attach(snp_readBed2(ld_bed,
                                 backingfile = tempfile()))
ld_snp$map$marker.ID = paste(ld_snp$map$chromosome,
                             ld_snp$map$physical.pos,
                             sep=':')
isoform_df = as.data.frame(matrix(ncol = 11,
                                  nrow = 0))
colnames(isoform_df) = c('Gene',
                         'Transcript',
                         'h2',
                         'h2.se',
                         'h2.p',
                         'R2',
                         'TWAS.Effect',
                         'TWAS.SE',
                         'TWAS.Z',
                         'TWAS.P',
                         'Permutation.P')

for (i in 1:length(unique(model$Transcript))){
  this_model = subset(model,Transcript == unique(model$Transcript)[i])
  this_model$ChrPos = paste(this_model$Chromosome,
                            this_model$Position,
                            sep=':')
  int_snp = intersect(gwas$ChrPos,this_model$ChrPos)
  if (length(int_snp) > 0){
    gwas_cur = gwas[match(int_snp,
                          gwas$ChrPos),]
    model_df = this_model[match(int_snp,
                                this_model$ChrPos),]
    ld_snp_cur = snp_attach(subset(ld_snp,
                                   ind.col = which(ld_snp$map$marker.ID %in%
                                                     int_snp),
                                   backingfile = tempfile()))

    model_df = model_df[match(ld_snp_cur$map$marker.ID,
                              model_df$ChrPos),]
    gwas_cur = gwas_cur[match(ld_snp_cur$map$marker.ID,
                              gwas_cur$ChrPos),]

    gwas_cur$EFFECT = ifelse(gwas_cur$A1 == model_df$ALT,
                             gwas_cur$BETA,
                             -1*gwas_cur$BETA)

    gwas_cur$Z = gwas_cur$EFFECT/gwas_cur$SE

    effect = model_df$Weight %*% gwas_cur$Z
    se = sqrt(as.numeric((model_df$Weight %*%
                            snp_cor(ld_snp_cur$genotypes) %*%
                            model_df$Weight)))

    calculateTWAS <- function(effects,
                              Z,
                              LD,
                              indices){
      effects = effects[indices]
      twasZ = as.numeric(effects %*% Z)
      twasr2pred = as.numeric(effects %*% LD %*% effects)
      if (twasr2pred > 0){
        twas = as.numeric(twasZ/sqrt(twasr2pred))
      } else {
        twas = 0
      }
      return(twas)
    }

    if (abs(effect/se) > 2.2){
      permutationLD = boot::boot(data = model_df$Weight,
                                 statistic = calculateTWAS,
                                 R = 10000,
                                 sim = 'permutation',
                                 Z = gwas_cur$Z,
                                 LD = snp_cor(ld_snp_cur$genotypes))
      permute.p = (10000 * mean(abs(permutationLD$t) >
                                  abs(permutationLD$t0)) + 1)/(10000+1)
    } else {
      permute.p = 1
    }

    df = data.frame(Gene = gene,
                    Transcript = this_model$Transcript[1],
                    h2 = this_model$h2[1],
                    h2.se = this_model$h2.se[1],
                    h2.p = this_model$h2.p[1],
                    R2 = this_model$R2[1],
                    TWAS.Effect = effect,
                    TWAS.SE = se,
                    TWAS.Z = effect/se,
                    TWAS.P = 2*pnorm(-1*abs(effect/se)),
                    Permutation.P = permute.p)
    isoform_df = rbind(isoform_df,df)}


}
isoform_df = subset(isoform_df,
                    h2.p < .1 &
                      R2 > 0.01)

if (nrow(isoform_df) > 0){
  isoform_df$Tissue = tissue

  require(ACAT)
  isoform_gene_df = data.frame(Gene = gene,
                               Tissue = tissue,
                               ACAT.P = ACAT(isoform_df$TWAS.P),
                               stageR.P = min(isoform_df$TWAS.P),
                               chi2.P = 1-pchisq(sum(isoform_df$TWAS.Z^2),
                                                 df = nrow(isoform_df)))

  fwrite(isoform_df,file.path(twas_results,'SCZ_isoform_isoTWAS_assoc.tsv'),
         sep='\t',append=T,row.names=F,quote=F)
  fwrite(isoform_gene_df,file.path(twas_results,'SCZ_gene_isoTWAS_assoc.tsv'),
         sep='\t',append=T,row.names=F,quote=F)
}

setwd(temp_gwas)
fff= list.files()
if (any(grepl(gene,fff))){
  fff = fff[grepl(gene,fff)]
  file.remove(fff)
}



require(data.table)
temp_gwas = '/pine/scr/a/b/abhattac/GTEx_isotwas/summStats/temp'
dir.create(temp_gwas)
twas_results = '/pine/scr/a/b/abhattac/PE_fetal/BD_results'
dir.create(twas_results)
model = fread(filename)
sumStat = '/pine/scr/a/b/abhattac/GTEx_isotwas/summStats/BD.PGC.2018.gwas_sumstats.txt'
tissue = 'Fetal_Brain'
gene = unique(model$Gene)[1]

chr = unique(model$Chromosome)[1]
min = min(model$Position)
max = max(model$Position)
system(paste0("awk '$1==",chr,"' ",sumStat,' > ',
              paste0(temp_gwas,'/',gene,'_gwas.tsv')))

require(data.table)
gwas = fread(paste0(temp_gwas,'/',gene,'_gwas.tsv'))
colnames(gwas) = colnames(fread(sumStat,nrow=2))
gwas = gwas[,c('CHR','SNP','BP','A1','A2','OR','SE','P')]
gwas$BETA = log(gwas$OR)
gwas$ChrPos = paste(gwas$CHR,gwas$BP,sep=':')

require(bigsnpr)
system(paste0('plink --bfile ',paste0('/pine/scr/a/b/abhattac/1000GP_Phase3/plink_eur/chr',chr),
              ' --chr ',chr,
              ' --from-bp ',max(0,min(min)-5e6),
              ' --to-bp ',max(max)+5e6,
              ' --make-bed --out ',file.path(temp_gwas,gene)))
ld_bed = paste0(file.path(temp_gwas,gene),'.bed')
ld_snp = snp_attach(snp_readBed2(ld_bed,
                                 backingfile = tempfile()))
ld_snp$map$marker.ID = paste(ld_snp$map$chromosome,
                             ld_snp$map$physical.pos,
                             sep=':')
isoform_df = as.data.frame(matrix(ncol = 11,
                                  nrow = 0))
colnames(isoform_df) = c('Gene',
                         'Transcript',
                         'h2',
                         'h2.se',
                         'h2.p',
                         'R2',
                         'TWAS.Effect',
                         'TWAS.SE',
                         'TWAS.Z',
                         'TWAS.P',
                         'Permutation.P')

for (i in 1:length(unique(model$Transcript))){
  this_model = subset(model,Transcript == unique(model$Transcript)[i])
  this_model$ChrPos = paste(this_model$Chromosome,
                            this_model$Position,
                            sep=':')
  int_snp = intersect(gwas$ChrPos,this_model$ChrPos)
  if (length(int_snp) > 0){
    gwas_cur = gwas[match(int_snp,
                          gwas$ChrPos),]
    model_df = this_model[match(int_snp,
                                this_model$ChrPos),]
    ld_snp_cur = snp_attach(subset(ld_snp,
                                   ind.col = which(ld_snp$map$marker.ID %in%
                                                     int_snp),
                                   backingfile = tempfile()))

    model_df = model_df[match(ld_snp_cur$map$marker.ID,
                              model_df$ChrPos),]
    gwas_cur = gwas_cur[match(ld_snp_cur$map$marker.ID,
                              gwas_cur$ChrPos),]

    gwas_cur$EFFECT = ifelse(gwas_cur$A1 == model_df$ALT,
                             gwas_cur$BETA,
                             -1*gwas_cur$BETA)

    gwas_cur$Z = gwas_cur$EFFECT/gwas_cur$SE

    effect = model_df$Weight %*% gwas_cur$Z
    se = sqrt(as.numeric((model_df$Weight %*%
                            snp_cor(ld_snp_cur$genotypes) %*%
                            model_df$Weight)))

    calculateTWAS <- function(effects,
                              Z,
                              LD,
                              indices){
      effects = effects[indices]
      twasZ = as.numeric(effects %*% Z)
      twasr2pred = as.numeric(effects %*% LD %*% effects)
      if (twasr2pred > 0){
        twas = as.numeric(twasZ/sqrt(twasr2pred))
      } else {
        twas = 0
      }
      return(twas)
    }

    if (abs(effect/se) > 2.2){
      permutationLD = boot::boot(data = model_df$Weight,
                                 statistic = calculateTWAS,
                                 R = 10000,
                                 sim = 'permutation',
                                 Z = gwas_cur$Z,
                                 LD = snp_cor(ld_snp_cur$genotypes))
      permute.p = (10000 * mean(abs(permutationLD$t) >
                                  abs(permutationLD$t0)) + 1)/(10000+1)
    } else {
      permute.p = 1
    }

    df = data.frame(Gene = gene,
                    Transcript = this_model$Transcript[1],
                    h2 = this_model$h2[1],
                    h2.se = this_model$h2.se[1],
                    h2.p = this_model$h2.p[1],
                    R2 = this_model$R2[1],
                    TWAS.Effect = effect,
                    TWAS.SE = se,
                    TWAS.Z = effect/se,
                    TWAS.P = 2*pnorm(-1*abs(effect/se)),
                    Permutation.P = permute.p)
    isoform_df = rbind(isoform_df,df)}


}
isoform_df = subset(isoform_df,
                    h2.p < .1 &
                      R2 > 0.01)

if (nrow(isoform_df) > 0){
  isoform_df$Tissue = tissue

  require(ACAT)
  isoform_gene_df = data.frame(Gene = gene,
                               Tissue = tissue,
                               ACAT.P = ACAT(isoform_df$TWAS.P),
                               stageR.P = min(isoform_df$TWAS.P),
                               chi2.P = 1-pchisq(sum(isoform_df$TWAS.Z^2),
                                                 df = nrow(isoform_df)))

  fwrite(isoform_df,file.path(twas_results,'BD_isoform_isoTWAS_assoc.tsv'),
         sep='\t',append=T,row.names=F,quote=F)
  fwrite(isoform_gene_df,file.path(twas_results,'BD_gene_isoTWAS_assoc.tsv'),
         sep='\t',append=T,row.names=F,quote=F)
}

setwd(temp_gwas)
fff= list.files()
if (any(grepl(gene,fff))){
  fff = fff[grepl(gene,fff)]
  file.remove(fff)
}





require(data.table)
temp_gwas = '/pine/scr/a/b/abhattac/GTEx_isotwas/summStats/tempASD'
dir.create(temp_gwas)
twas_results = '/pine/scr/a/b/abhattac/PE_fetal/ASD_results'
dir.create(twas_results)
model = fread(filename)
sumStat = '/pine/scr/a/b/abhattac/GTEx_isotwas/summStats/ASD.iPSYCHPGC.2018.sumstats.gz'
tissue = 'Fetal_Brain'
gene = unique(model$Gene)[1]


chr = unique(model$Chromosome)[1]
min = min(model$Position)
max = max(model$Position)

require(data.table)
gwas = fread(sumStat)
gwas = subset(gwas,SNP %in% model$SNP)

require(bigsnpr)
system(paste0('plink --bfile ',paste0('/pine/scr/a/b/abhattac/1000GP_Phase3/plink_eur/chr',chr),
              ' --chr ',chr,
              ' --from-bp ',max(0,min(min)-5e6),
              ' --to-bp ',max(max)+5e6,
              ' --make-bed --out ',file.path(temp_gwas,gene)))
ld_bed = paste0(file.path(temp_gwas,gene),'.bed')
ld_snp = snp_attach(snp_readBed2(ld_bed,
                                 backingfile = tempfile()))
isoform_df = as.data.frame(matrix(ncol = 11,
                                  nrow = 0))
colnames(isoform_df) = c('Gene',
                         'Transcript',
                         'h2',
                         'h2.se',
                         'h2.p',
                         'R2',
                         'TWAS.Effect',
                         'TWAS.SE',
                         'TWAS.Z',
                         'TWAS.P',
                         'Permutation.P')

for (i in 1:length(unique(model$Transcript))){
  this_model = subset(model,Transcript == unique(model$Transcript)[i])
  this_model$ChrPos = paste(this_model$Chromosome,
                            this_model$Position,
                            sep=':')
  int_snp = intersect(gwas$SNP,this_model$SNP)
  if (length(int_snp) > 0){
    gwas_cur = gwas[match(int_snp,
                          gwas$SNP),]
    model_df = this_model[match(int_snp,
                                this_model$SNP),]
    ld_snp_cur = snp_attach(subset(ld_snp,
                                   ind.col = which(ld_snp$map$marker.ID %in%
                                                     int_snp),
                                   backingfile = tempfile()))

    model_df = model_df[match(ld_snp_cur$map$marker.ID,
                              model_df$SNP),]
    gwas_cur = gwas_cur[match(ld_snp_cur$map$marker.ID,
                              gwas_cur$SNP),]

    gwas_cur$Z = ifelse(gwas_cur$A1 == model_df$ALT,
                        gwas_cur$Z,
                        -1*gwas_cur$Z)

    effect = model_df$Weight %*% gwas_cur$Z
    se = sqrt(as.numeric((model_df$Weight %*%
                            snp_cor(ld_snp_cur$genotypes) %*%
                            model_df$Weight)))

    calculateTWAS <- function(effects,
                              Z,
                              LD,
                              indices){
      effects = effects[indices]
      twasZ = as.numeric(effects %*% Z)
      twasr2pred = as.numeric(effects %*% LD %*% effects)
      if (twasr2pred > 0){
        twas = as.numeric(twasZ/sqrt(twasr2pred))
      } else {
        twas = 0
      }
      return(twas)
    }

    if (abs(effect/se) > 2.2){
      permutationLD = boot::boot(data = model_df$Weight,
                                 statistic = calculateTWAS,
                                 R = 10000,
                                 sim = 'permutation',
                                 Z = gwas_cur$Z,
                                 LD = snp_cor(ld_snp_cur$genotypes))
      permute.p = (10000 * mean(abs(permutationLD$t) >
                                  abs(permutationLD$t0)) + 1)/(10000+1)
    } else {
      permute.p = 1
    }

    df = data.frame(Gene = gene,
                    Transcript = this_model$Transcript[1],
                    h2 = this_model$h2[1],
                    h2.se = this_model$h2.se[1],
                    h2.p = this_model$h2.p[1],
                    R2 = this_model$R2[1],
                    TWAS.Effect = effect,
                    TWAS.SE = se,
                    TWAS.Z = effect/se,
                    TWAS.P = 2*pnorm(-1*abs(effect/se)),
                    Permutation.P = permute.p)
    isoform_df = rbind(isoform_df,df)}


}
isoform_df = subset(isoform_df,
                    h2.p < .1 &
                      R2 > 0.01)

if (nrow(isoform_df) > 0){
  isoform_df$Tissue = tissue

  require(ACAT)
  isoform_gene_df = data.frame(Gene = gene,
                               Tissue = tissue,
                               ACAT.P = ACAT(isoform_df$TWAS.P),
                               stageR.P = min(isoform_df$TWAS.P),
                               chi2.P = 1-pchisq(sum(isoform_df$TWAS.Z^2),
                                                 df = nrow(isoform_df)))

  fwrite(isoform_df,file.path(twas_results,'ASD_isoform_isoTWAS_assoc.tsv'),
         sep='\t',append=T,row.names=F,quote=F)
  fwrite(isoform_gene_df,file.path(twas_results,'ASD_gene_isoTWAS_assoc.tsv'),
         sep='\t',append=T,row.names=F,quote=F)
}

setwd(temp_gwas)
fff= list.files()
if (any(grepl(gene,fff))){
  fff = fff[grepl(gene,fff)]
  file.remove(fff)
}






require(data.table)
temp_gwas = '/pine/scr/a/b/abhattac/GTEx_isotwas/summStats/tempALZ'
dir.create(temp_gwas)
twas_results = '/pine/scr/a/b/abhattac/PE_fetal/ALZ_results'
dir.create(twas_results)
model = fread(filename)
sumStat = '/pine/scr/a/b/abhattac/GTEx_isotwas/summStats/PGC_ALZ2_sumstats_hg38.txt'
tissue = 'Fetal_Brain'
gene = unique(model$Gene)[1]
model$ChrPos = paste(model$Chromosome,
                     model$Position,
                     sep=':')


chr = unique(model$Chromosome)[1]
min = min(model$Position)
max = max(model$Position)

system(paste0("awk '$1==",chr,"' ",sumStat,' > ',
              paste0(temp_gwas,'/',gene,'_gwas.tsv')))

require(data.table)
gwas = fread(paste0(temp_gwas,'/',gene,'_gwas.tsv'))
colnames(gwas) = colnames(fread(sumStat,nrow=2))
### LIFTOVER
library(rtracklayer)
library(liftOver)
path = '/pine/scr/a/b/abhattac/hg38ToHg19.over.chain'
ch = import.chain(path)
gwas$SNP = paste0('SNP',1:nrow(gwas))

snp_df = data.frame(SNP = gwas$SNP,
                    Position = gwas$POS)
snp_df = snp_df[!duplicated(gwas$SNP),]

df_out <- data.frame(chr=paste0('chr',gwas$CHR[1]),
                     start=snp_df$Position,
                     end=snp_df$Position + 1,
                     strand='+',
                     score=1:nrow(snp_df),
                     names = as.character(snp_df$SNP))

grl <- makeGRangesListFromDataFrame(df_out,names.field = 'names')
names(grl) = df_out$names
cur19 = liftOver(grl, ch)
cur19 = as.data.frame(cur19)
liftover = data.frame(SNP = cur19$group_name,
                      Position_hg19 = cur19$start)

gwas = merge(gwas,liftover,by='SNP')
gwas$ChrPos = paste(gwas$CHR,gwas$Position_hg19,sep=':')
gwas = subset(gwas,ChrPos %in% model$ChrPos)

require(bigsnpr)
system(paste0('plink --bfile ',paste0('/pine/scr/a/b/abhattac/1000GP_Phase3/plink_eur/chr',chr),
              ' --chr ',chr,
              ' --from-bp ',max(0,min(min)-5e6),
              ' --to-bp ',max(max)+5e6,
              ' --make-bed --out ',file.path(temp_gwas,gene)))
ld_bed = paste0(file.path(temp_gwas,gene),'.bed')
ld_snp = snp_attach(snp_readBed2(ld_bed,
                                 backingfile = tempfile()))
ld_snp$map$marker.ID = paste(ld_snp$map$chromosome,
                             ld_snp$map$physical.pos,
                             sep=':')
isoform_df = as.data.frame(matrix(ncol = 11,
                                  nrow = 0))
colnames(isoform_df) = c('Gene',
                         'Transcript',
                         'h2',
                         'h2.se',
                         'h2.p',
                         'R2',
                         'TWAS.Effect',
                         'TWAS.SE',
                         'TWAS.Z',
                         'TWAS.P',
                         'Permutation.P')

for (i in 1:length(unique(model$Transcript))){
  this_model = subset(model,Transcript == unique(model$Transcript)[i])
  this_model$ChrPos = paste(this_model$Chromosome,
                            this_model$Position,
                            sep=':')
  int_snp = intersect(gwas$ChrPos,this_model$ChrPos)
  if (length(int_snp) > 0){
    gwas_cur = gwas[match(int_snp,
                          gwas$ChrPos),]
    model_df = this_model[match(int_snp,
                                this_model$ChrPos),]
    ld_snp_cur = snp_attach(subset(ld_snp,
                                   ind.col = which(ld_snp$map$marker.ID %in%
                                                     int_snp),
                                   backingfile = tempfile()))

    model_df = model_df[match(ld_snp_cur$map$marker.ID,
                              model_df$ChrPos),]
    gwas_cur = gwas_cur[match(ld_snp_cur$map$marker.ID,
                              gwas_cur$ChrPos),]

    gwas_cur$Z = ifelse(gwas_cur$A1 == model_df$ALT,
                        gwas_cur$Z,
                        -1*gwas_cur$Z)

    effect = model_df$Weight %*% gwas_cur$Z
    se = sqrt(as.numeric((model_df$Weight %*%
                            snp_cor(ld_snp_cur$genotypes) %*%
                            model_df$Weight)))

    calculateTWAS <- function(effects,
                              Z,
                              LD,
                              indices){
      effects = effects[indices]
      twasZ = as.numeric(effects %*% Z)
      twasr2pred = as.numeric(effects %*% LD %*% effects)
      if (twasr2pred > 0){
        twas = as.numeric(twasZ/sqrt(twasr2pred))
      } else {
        twas = 0
      }
      return(twas)
    }

    if (abs(effect/se) > 2.2){
      permutationLD = boot::boot(data = model_df$Weight,
                                 statistic = calculateTWAS,
                                 R = 10000,
                                 sim = 'permutation',
                                 Z = gwas_cur$Z,
                                 LD = snp_cor(ld_snp_cur$genotypes))
      permute.p = (10000 * mean(abs(permutationLD$t) >
                                  abs(permutationLD$t0)) + 1)/(10000+1)
    } else {
      permute.p = 1
    }

    df = data.frame(Gene = gene,
                    Transcript = this_model$Transcript[1],
                    h2 = this_model$h2[1],
                    h2.se = this_model$h2.se[1],
                    h2.p = this_model$h2.p[1],
                    R2 = this_model$R2[1],
                    TWAS.Effect = effect,
                    TWAS.SE = se,
                    TWAS.Z = effect/se,
                    TWAS.P = 2*pnorm(-1*abs(effect/se)),
                    Permutation.P = permute.p)
    isoform_df = rbind(isoform_df,df)}


}
isoform_df = subset(isoform_df,
                    h2.p < .1 &
                      R2 > 0.01)

if (nrow(isoform_df) > 0){
  isoform_df$Tissue = tissue

  require(ACAT)
  isoform_gene_df = data.frame(Gene = gene,
                               Tissue = tissue,
                               ACAT.P = ACAT(isoform_df$TWAS.P),
                               stageR.P = min(isoform_df$TWAS.P),
                               chi2.P = 1-pchisq(sum(isoform_df$TWAS.Z^2),
                                                 df = nrow(isoform_df)))

  fwrite(isoform_df,file.path(twas_results,'ALZ_isoform_isoTWAS_assoc.tsv'),
         sep='\t',append=T,row.names=F,quote=F)
  fwrite(isoform_gene_df,file.path(twas_results,'ALZ_gene_isoTWAS_assoc.tsv'),
         sep='\t',append=T,row.names=F,quote=F)
}

setwd(temp_gwas)
fff= list.files()
if (any(grepl(gene,fff))){
  fff = fff[grepl(gene,fff)]
  file.remove(fff)
}








require(data.table)
temp_gwas = '/pine/scr/a/b/abhattac/GTEx_isotwas/summStats/tempADHD'
dir.create(temp_gwas)
twas_results = '/pine/scr/a/b/abhattac/PE_fetal/ADHD_results'
dir.create(twas_results)
model = fread(filename)
sumStat = '/pine/scr/a/b/abhattac/GTEx_isotwas/summStats/adhd_jul2017'
tissue = 'Fetal_Brain'
gene = unique(model$Gene)[1]
model$ChrPos = paste(model$Chromosome,
                     model$Position,
                     sep=':')


chr = unique(model$Chromosome)[1]
min = min(model$Position)
max = max(model$Position)

system(paste0("awk '$1==",chr,"' ",sumStat,' > ',
              paste0(temp_gwas,'/',gene,'_gwas.tsv')))

require(data.table)
gwas = fread(paste0(temp_gwas,'/',gene,'_gwas.tsv'))
colnames(gwas) = colnames(fread(sumStat,nrow=2))
gwas$ChrPos = paste(gwas$CHR,gwas$BP,sep=':')
gwas = subset(gwas,ChrPos %in% model$ChrPos)
gwas$BETA = log(gwas$OR)

require(bigsnpr)
system(paste0('plink --bfile ',paste0('/pine/scr/a/b/abhattac/1000GP_Phase3/plink_eur/chr',chr),
              ' --chr ',chr,
              ' --from-bp ',max(0,min(min)-5e6),
              ' --to-bp ',max(max)+5e6,
              ' --make-bed --out ',file.path(temp_gwas,gene)))
ld_bed = paste0(file.path(temp_gwas,gene),'.bed')
ld_snp = snp_attach(snp_readBed2(ld_bed,
                                 backingfile = tempfile()))
ld_snp$map$marker.ID = paste(ld_snp$map$chromosome,
                             ld_snp$map$physical.pos,
                             sep=':')
isoform_df = as.data.frame(matrix(ncol = 11,
                                  nrow = 0))
colnames(isoform_df) = c('Gene',
                         'Transcript',
                         'h2',
                         'h2.se',
                         'h2.p',
                         'R2',
                         'TWAS.Effect',
                         'TWAS.SE',
                         'TWAS.Z',
                         'TWAS.P',
                         'Permutation.P')

for (i in 1:length(unique(model$Transcript))){
  this_model = subset(model,Transcript == unique(model$Transcript)[i])
  this_model$ChrPos = paste(this_model$Chromosome,
                            this_model$Position,
                            sep=':')
  int_snp = intersect(gwas$ChrPos,this_model$ChrPos)
  if (length(int_snp) > 0){
    gwas_cur = gwas[match(int_snp,
                          gwas$ChrPos),]
    model_df = this_model[match(int_snp,
                                this_model$ChrPos),]
    ld_snp_cur = snp_attach(subset(ld_snp,
                                   ind.col = which(ld_snp$map$marker.ID %in%
                                                     int_snp),
                                   backingfile = tempfile()))

    model_df = model_df[match(ld_snp_cur$map$marker.ID,
                              model_df$ChrPos),]
    gwas_cur = gwas_cur[match(ld_snp_cur$map$marker.ID,
                              gwas_cur$ChrPos),]

    gwas_cur$EFFECT = ifelse(gwas_cur$A1 == model_df$ALT,
                             gwas_cur$BETA,
                             -1*gwas_cur$BETA)
    gwas_cur$Z = gwas_cur$EFFECT/gwas_cur$SE

    effect = model_df$Weight %*% gwas_cur$Z
    se = sqrt(as.numeric((model_df$Weight %*%
                            snp_cor(ld_snp_cur$genotypes) %*%
                            model_df$Weight)))

    calculateTWAS <- function(effects,
                              Z,
                              LD,
                              indices){
      effects = effects[indices]
      twasZ = as.numeric(effects %*% Z)
      twasr2pred = as.numeric(effects %*% LD %*% effects)
      if (twasr2pred > 0){
        twas = as.numeric(twasZ/sqrt(twasr2pred))
      } else {
        twas = 0
      }
      return(twas)
    }

    if (abs(effect/se) > 2.2){
      permutationLD = boot::boot(data = model_df$Weight,
                                 statistic = calculateTWAS,
                                 R = 10000,
                                 sim = 'permutation',
                                 Z = gwas_cur$Z,
                                 LD = snp_cor(ld_snp_cur$genotypes))
      permute.p = (10000 * mean(abs(permutationLD$t) >
                                  abs(permutationLD$t0)) + 1)/(10000+1)
    } else {
      permute.p = 1
    }

    df = data.frame(Gene = gene,
                    Transcript = this_model$Transcript[1],
                    h2 = this_model$h2[1],
                    h2.se = this_model$h2.se[1],
                    h2.p = this_model$h2.p[1],
                    R2 = this_model$R2[1],
                    TWAS.Effect = effect,
                    TWAS.SE = se,
                    TWAS.Z = effect/se,
                    TWAS.P = 2*pnorm(-1*abs(effect/se)),
                    Permutation.P = permute.p)
    isoform_df = rbind(isoform_df,df)}


}
isoform_df = subset(isoform_df,
                    h2.p < .1 &
                      R2 > 0.01)

if (nrow(isoform_df) > 0){
  isoform_df$Tissue = tissue

  require(ACAT)
  isoform_gene_df = data.frame(Gene = gene,
                               Tissue = tissue,
                               ACAT.P = ACAT(isoform_df$TWAS.P),
                               stageR.P = min(isoform_df$TWAS.P),
                               chi2.P = 1-pchisq(sum(isoform_df$TWAS.Z^2),
                                                 df = nrow(isoform_df)))

  fwrite(isoform_df,file.path(twas_results,'ADHD_isoform_isoTWAS_assoc.tsv'),
         sep='\t',append=T,row.names=F,quote=F)
  fwrite(isoform_gene_df,file.path(twas_results,'ADHD_gene_isoTWAS_assoc.tsv'),
         sep='\t',append=T,row.names=F,quote=F)
}

setwd(temp_gwas)
fff= list.files()
if (any(grepl(gene,fff))){
  fff = fff[grepl(gene,fff)]
  file.remove(fff)
}
