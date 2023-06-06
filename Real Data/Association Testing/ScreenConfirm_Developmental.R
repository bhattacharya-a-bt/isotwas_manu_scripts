require(data.table)
require(vroom)
setwd('Walker')

fff = list.files(recursive = T)

traits = c('ADHD',
           'ALZ',
           'ASD',
           'BP',
           'BV',
           'CortTH',
           'ICV',
           'MDD',
           'NTSM',
           'PANIC',
           'SCZ',
           'PTSD',
           'OCD',
           'AN',
           'CDG')

require(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

list_intersect = list()
allGenes = c()
require(rlist)

for (tr in traits){
  
  iso_f = paste0(tr,'/',tr,'_isoTWAS.tsv')
  twas_f = paste0(tr,'/',tr,'_TWAS.tsv')
  isotwas_res = vroom::vroom(iso_f)
  isotwas_res = isotwas_res[complete.cases(isotwas_res),]
  isotwas_res = isotwas_res[order(abs(isotwas_res$Z),decreasing = T),]
  isotwas_res = isotwas_res[!duplicated(isotwas_res$Transcript) & 
                              abs(isotwas_res$Z) < Inf,]
  colnames(isotwas_res)[8] = 'R2'
  
  twas_res = vroom::vroom(twas_f)
  twas_res = twas_res[complete.cases(twas_res) & abs(twas_res$Z) < Inf,]
  twas_res = twas_res[!duplicated(twas_res$Gene),]
  
  intersectGenes = intersect(twas_res$Gene,isotwas_res$Gene)
  list_intersect = list.append(list_intersect,intersectGenes)
  allGenes = unique(c(allGenes,
                      twas_res$Gene,
                      isotwas_res$Gene))
  
}

names(list_intersect) = traits
saveRDS(list_intersect,'intersect_gene_tested_developmental.RDS')



bm = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol',
                          'chromosome_name','start_position','end_position',
                          'gene_biotype'),
           filters = 'ensembl_gene_id',
           values = allGenes, 
           mart = ensembl)
colnames(bm) = c('Gene','HGNC','Chromosome','Start','End','Biotype')



for (i in 1:15){
  tr = traits[i]
  iso_f = paste0(tr,'/',tr,'_isoTWAS.tsv')
  twas_f = paste0(tr,'/',tr,'_TWAS.tsv')
  isotwas_res = vroom::vroom(iso_f)
  isotwas_res = isotwas_res[complete.cases(isotwas_res),]
  isotwas_res = isotwas_res[order(abs(isotwas_res$Z),decreasing = T),]
  isotwas_res = isotwas_res[!duplicated(isotwas_res$Transcript) & 
                              abs(isotwas_res$Z) < Inf,]
  colnames(isotwas_res)[8] = 'R2'
  
  twas_res = vroom::vroom(twas_f)
  twas_res = twas_res[complete.cases(twas_res) & abs(twas_res$Z) < Inf,]
  twas_res = twas_res[!duplicated(twas_res$Gene),]
  
  require(bacon)
  bc_twas = bacon(twas_res$Z)
  bc_isotwas = bacon(isotwas_res$Z)
  outFile_bacon = 'bacon_chisquare_compAdult.tsv'
  df_measures = data.frame(Trait = tr,
                           Method = c('TWAS','isoTWAS'),
                           Median_Inflation = c(median(bc_twas@traces[,7,1]),
                                                median(bc_isotwas@traces[,7,1])),
                           Lower_Inflation = c(quantile(bc_twas@traces[,7,1],.025),
                                               quantile(bc_isotwas@traces[,7,1],.025)),
                           Upper_Inflation = c(quantile(bc_twas@traces[,7,1],.975),
                                               quantile(bc_isotwas@traces[,7,1],.975)))
  
  isotwas_res = merge(bm,isotwas_res,by='Gene')
  twas_res = merge(bm,twas_res,by='Gene')
  
  www = which(isotwas_res$Chromosome == 6 &
                isotwas_res$Start < 35e6 &
                isotwas_res$End > 27e6)
  isotwas_res = isotwas_res[-www,]
  
  www = which(twas_res$Chromosome == 6 &
                twas_res$Start < 35e6 &
                twas_res$End > 27e6)
  twas_res = twas_res[-www,]
  
  
  gene = data.frame(Gene = unique(isotwas_res$Gene),
                    HGNC = NA,
                    Chromosome = NA,
                    Start = NA,
                    End = NA,
                    Biotype = NA,
                    Screen.P = NA)
  
  require(tidyverse)
  gene = isotwas_res %>%
    group_by(Gene) %>%
    summarise(HGNC = unique(HGNC),
              Chromosome = unique(Chromosome),
              Start = unique(Start),
              End = unique(End),
              Biotype = unique(Biotype),
              Screen.P = isotwas::p_screen(P))
  
  alpha1=.05
  G = nrow(gene)
  gene$Screen.P.Adjusted = p.adjust(gene$Screen.P,method = 'fdr')
  R = length(unique(gene$Gene[gene$Screen.P.Adjusted < alpha1]))
  alpha2 = (R*alpha1)/G
  isoform_new = as.data.frame(matrix(nrow = 0,
                                     ncol = ncol(isotwas_res)+2))
  colnames(isoform_new) = c(colnames(isotwas_res),'Screen.P','Confirmation.P')
  gene = gene[order(gene$Screen.P),]
  ttt = merge(isotwas_res,
              gene[,c('Gene','Screen.P',
                      'Screen.P.Adjusted')])
  isoform_new = ttt %>%
    group_by(Gene) %>%
    summarise(Transcript = Transcript,
              Confirmation.P = isotwas::p_confirm(P,alpha = alpha2))
  isoform_new = merge(isoform_new,ttt,by=c('Gene','Transcript'))
  isoform_new$Confirmation.P = ifelse(isoform_new$Screen.P.Adjusted < 0.05,
                                      isoform_new$Confirmation.P,
                                      1)
  
  intersect_gene = intersect(isoform_new$Gene,twas_res$Gene)
  isoform_new = isoform_new[order(isoform_new$Confirmation.P),]
  twas_res = subset(twas_res,Gene %in% intersect_gene)
  isoform_new_intersect = subset(isoform_new,Gene %in% intersect_gene)
  twas_res$FDR = p.adjust(twas_res$P,method = 'fdr')
  
  df_comp = merge(twas_res[,c('Gene','P')],
                  gene[,c('Gene','Screen.P')],
                  by='Gene')
  df_comp$TWAS_ChiSquare = qchisq(df_comp$P,df=1)
  df_comp$isoTWAS_ChiSquare = qchisq(df_comp$Screen.P,df=1)
  df_comp = subset(df_comp,TWAS_ChiSquare > 1 & isoTWAS_ChiSquare > 1 &
                     TWAS_ChiSquare != isoTWAS_ChiSquare)
  
  percdiff = 
    (df_comp$isoTWAS_ChiSquare - 
       df_comp$TWAS_ChiSquare)/(df_comp$TWAS_ChiSquare) * 100
  df_measures$Mean_ChiSquare = mean(percdiff)
  df_measures$SD_ChiSquare = Monte.Carlo.se::jack.se(percdiff,mean)
  
  fwrite(df_measures,outFile_bacon,sep='\t',
         row.names=F,quote=F,append=T)
  
  
  aaa = data.frame(Threshold = 'Adjusted P < 0.05',
                   Number = c(length(unique(isoform_new$Gene[isoform_new$Screen.P.Adjusted < alpha1 &
                                                               isoform_new$Confirmation.P < alpha2 &
                                                               isoform_new$permute.P < 0.05])),
                              length(unique(isoform_new$Transcript[isoform_new$Screen.P.Adjusted < alpha1 &
                                                                     isoform_new$Confirmation.P < alpha2 &
                                                                     isoform_new$permute.P < 0.05])),
                              sum(twas_res$FDR < .05 & twas_res$permute.P < .05)),
                   Trait = tr,
                   Feature = c('Gene','Isoform','Gene'),
                   Method = c('isoTWAS','isoTWAS','TWAS'))
  aaa$Tissue = 'Developing brain (PEC)'
  fwrite(aaa,'NumberDevPECWalker_allstudied.tsv',sep='\t',
         append=T,row.names=F,quote=F)
  
  aaa = data.frame(Threshold = 'Adjusted P < 0.05',
                   Number = c(length(unique(isoform_new_intersect$Gene[isoform_new_intersect$Screen.P.Adjusted < alpha1 &
                                                                         isoform_new_intersect$Confirmation.P < alpha2 &
                                                                         isoform_new_intersect$permute.P < 0.05])),
                              length(unique(isoform_new_intersect$Transcript[isoform_new_intersect$Screen.P.Adjusted < alpha1 &
                                                                               isoform_new_intersect$Confirmation.P < alpha2 &
                                                                               isoform_new_intersect$permute.P < 0.05])),
                              sum(twas_res$FDR < .05 & twas_res$permute.P < .05)),
                   Trait = tr,
                   Feature = c('Gene','Isoform','Gene'),
                   Method = c('isoTWAS','isoTWAS','TWAS'))
  aaa$Tissue = 'Developing brain (PEC)'
  fwrite(aaa,'NumberDevPECWalker_intersect.tsv',sep='\t',
         append=T,row.names=F,quote=F)
  
  gnomad = 
    fread('/u/scratch/a/abtbhatt/gnomad_plof.txt')
  gnomad$HGNC = gnomad$gene
  
  isoform_sig = subset(isoform_new, 
                       Screen.P < alpha1 &
                         Confirmation.P < alpha2 &
                         permute.P < 0.05)
  
  if (nrow(isoform_sig) > 1){
  isoform_sig$Chromosome = as.numeric(isoform_sig$Chromosome)
  
  isoform_sig = isoform_sig[order(isoform_sig$Chromosome,
                                  isoform_sig$Start),]
  isoform_sig$Trait = tr
  isoform_sig = isoform_sig[,c('Trait','Gene','HGNC','Chromosome',
                               'Start','End','Biotype','Transcript',
                               'Z','P','permute.P','topSNP','topSNP.P',
                               'Screen.P.Adjusted','Confirmation.P','R2')]
  colnames(isoform_sig) = c('Trait','Gene','HGNC','Chromosome',
                            'Start','End','Biotype','Transcript',
                            'Z','P','Permutation P','Top GWAS SNP','Top GWAS P',
                            'Screening Adjusted P','Confirmation P','R2')
  
  isoform_sig = merge(isoform_sig,
                      gnomad[,c('HGNC','pLI')],
                      by = 'HGNC',
                      all.x=T)
  
  xlsx::write.xlsx(isoform_sig,
                   file = 'Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx',
                   append = T,
                   sheetName = tr,
                   col.names = T,
                   row.names = F)
  }
  
  twas_sig = subset(twas_res,FDR < .05 & permute.P < .05)
  
  if (nrow(twas_sig) > 1){
  twas_sig = twas_sig[order(as.numeric(twas_sig$Chromosome),
                            twas_sig$Start),]
  twas_sig = twas_sig[,c('Gene','HGNC','Chromosome','Start','End',
                         'Biotype','Transcript','Z','P','permute.P',
                         'topSNP','topSNP.P','FDR','R2')]
  colnames(twas_sig) = c('Gene','HGNC','Chromosome','Start',
                         'End','Biotype','Transcript','Z','P','Permutation P',
                         'Top GWAS SNP','Top GWAS P','Adjusted P','R2')
  
  
  twas_sig = merge(twas_sig,
                   gnomad[,c('HGNC','pLI')],
                   by='HGNC',
                   all.x=T)
  xlsx::write.xlsx(twas_sig,
                   file = 'Developmental_Associations_BeforeFOCUS_TWAS.xlsx',
                   append = T,
                   sheetName = tr,
                   col.names = T,
                   row.names = F)
  }
  
  
  
  
}

library(readxl)
path <- "Developmental_Associations_BeforeFOCUS_TWAS.xlsx"
excel_sheets(path = path)


library(readxl)
path <- "Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx"
excel_sheets(path = path)

df = data.frame(Trait = traits,
                TWAS = 0,
                TWAS_total = 0,
                isoTWAS = 0,
                isoTWAS_total = 0)
for (t in 1:length(traits)){
  print(t)
  tr = traits[t]
  if (tr %in% excel_sheets(path = 'Developmental_Associations_BeforeFOCUS_TWAS.xlsx')){
    
    twas = read_xlsx('Developmental_Associations_BeforeFOCUS_TWAS.xlsx',
                     sheet = tr)
    df$TWAS[t] = sum(twas$pLI > 0.9 & !is.na(twas$pLI))
    df$TWAS_total[t] = length(unique(twas$Gene))
    
  }
  
  if (tr %in% excel_sheets('Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx')){
  isotwas = read_xlsx('Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx',
                      sheet = tr)
  df$isoTWAS[t] = length(unique(isotwas$Gene[isotwas$pLI > 0.9 &
                                               !is.na(isotwas$pLI)]))
  
  df$isoTWAS_total[t] = length(unique(isotwas$Gene))
  }
  
  
}

test.this = matrix(c(sum(df$isoTWAS),sum(df$TWAS),
                     sum(df$isoTWAS_total - df$isoTWAS),
                     sum(df$TWAS_total - df$TWAS)),
                   ncol = 2)
p.value = fisher.test(test.this)$p.value
df$FisherP = p.value
fwrite(df,'pLI_count_test_allStudied.tsv',sep='\t',
       col.names=T,row.names=F,quote=F)




df = data.frame(Trait = traits,
                TWAS = 0,
                TWAS_total = 0,
                isoTWAS = 0,
                isoTWAS_total = 0)
for (t in 1:length(traits)){
  print(t)
  
  intGenes = list_intersect[[t]]
  
  tr = traits[t]
  if (tr %in% excel_sheets(path = 'Developmental_Associations_BeforeFOCUS_TWAS.xlsx')){
    
    twas = read_xlsx('Developmental_Associations_BeforeFOCUS_TWAS.xlsx',
                     sheet = tr)
    df$TWAS[t] = sum(twas$pLI > 0.9 & !is.na(twas$pLI) & 
                       twas$Gene %in% intGenes)
    df$TWAS_total[t] = length(unique(twas$Gene[twas$Gene %in% intGenes]))
    
  }
  
  if (tr %in% excel_sheets('Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx')){
    isotwas = read_xlsx('Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx',
                        sheet = tr)
    df$isoTWAS[t] = length(unique(isotwas$Gene[isotwas$pLI > 0.9 &
                                                 !is.na(isotwas$pLI) &
                                                 isotwas$Gene %in% intGenes]))
    
    df$isoTWAS_total[t] = length(unique(isotwas$Gene[isotwas$Gene %in% 
                                                       intGenes]))
  }
  
  
}

df_reg = data.frame(Method = rep(c('TWAS','isoTWAS'),each = 15),
                    Number = c(df$TWAS,df$isoTWAS))
reg = lm(Number ~ Method,data = df_reg)

sapply(list_intersect,length)

test.this = matrix(c(sum(df$isoTWAS),sum(df$TWAS),
                     sum(df$isoTWAS_total - df$isoTWAS),
                     sum(df$TWAS_total - df$TWAS)),
                   ncol = 2)
p.value = fisher.test(test.this)$p.value
df$FisherP = p.value
fwrite(df,'pLI_count_test_intersectGenes.tsv',sep='\t',
       col.names=T,row.names=F,quote=F)










require(data.table)
require(vroom)
setwd('Walker')
library(readxl)
path <- "Developmental_Associations_BeforeFOCUS_TWAS.xlsx"
excel_sheets(path = path)

traits = c('ADHD',
           'ALZ',
           'ASD',
           'BP',
           'BV',
           'CortTH',
           'ICV',
           'MDD',
           'NTSM',
           'PANIC',
           'SCZ',
           'PTSD',
           'OCD',
           'AN',
           'CDG')


library(readxl)
require(LDlinkR)
require(bigsnpr)

snps = snp_attach(snp_readBed2('/u/scratch/a/abtbhatt/Walker/PEC_fetal_626_hmp.bed',
                               backingfile = tempfile()))

path <- "Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx"
excel_sheets(path = path)


df_out = data.frame(Trait = NA,
                    Method = NA,
                    Tissue = NA,
                    Number = NA)
for (t in 1:length(traits)){
  print(t)
  tr = traits[t]
  twas_GWAS_SNP = c()
  if (tr %in% excel_sheets(path = 'Developmental_Associations_BeforeFOCUS_TWAS.xlsx')){
    
    twas = read_xlsx('Developmental_Associations_BeforeFOCUS_TWAS.xlsx',
                     sheet = tr)
    twas$Method = 'TWAS'
    
  } else {
    twas = data.frame(Gene = 'Dummy',
                      Chromosome = 1,
                      Start = 1,
                      End = 1,
                      Z = 1,
                      'Top GWAS SNP' = 'Dummy',
                      'Top GWAS P' = 1,
                      Method = 'TWAS')
    colnames(twas)[6:7] = c('Top GWAS SNP',
                            'Top GWAS P')
  }
  
  isotwas_GWAS_SNP = c()
  if (tr %in% excel_sheets('Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx')){
    isotwas = read_xlsx('Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx',
                        sheet = tr)
    isotwas$Method = 'isoTWAS'
  } else {
    
    isotwas = data.frame(Gene = 'Dummy',
                         Chromosome = 1,
                         Start = 1,
                         End = 1,
                         Z = 1,
                         'Top GWAS SNP' = 'Dummy',
                         'Top GWAS P' = 1,
                         Method = 'isoTWAS')
    colnames(isotwas)[6:7] = c('Top GWAS SNP',
                               'Top GWAS P')
  }
  
  total = rbind(twas[,c('Gene','Chromosome','Start','End',
                        'Z','Top GWAS SNP','Top GWAS P','Method')],
                isotwas[,c('Gene','Chromosome','Start','End',
                           'Z','Top GWAS SNP','Top GWAS P','Method')])
  
  total = as.data.frame(total[order(total$Chromosome,
                                    total$Start),])
  total = subset(total,`Top GWAS P` < 5e-8)
  if (nrow(total) > 0){
    total$Group = 1
    ggg = 1
    for (i in 1:(nrow(total)-1)){
      if (total$End[i] <= total$Start[i+1] - 1e6){
        ggg = ggg+1
        total$Group[(i+1):nrow(total)] = ggg
      }
    }
    
    both_n = twas_n = isotwas_n = 0
    for (ggg in unique(total$Group)){
      
      ttt = subset(total,Group == ggg)
      if (all(c('TWAS','isoTWAS') %in% ttt$Method)){
        both_n = both_n + 1
      }
      
      if ('TWAS' %in% ttt$Method &
          !'isoTWAS' %in% ttt$Method){
        twas_n = twas_n + 1
      }
      
      if (!'TWAS' %in% ttt$Method &
          'isoTWAS' %in% ttt$Method){
        isotwas_n = isotwas_n + 1
      }
      
    }
  } else {
    both_n = 0
    twas_n = 0
    isotwas_n = 0
  }
  
  df_out = rbind(df_out,
                 data.frame(Trait = tr,
                            Method = c('Both','TWAS','isoTWAS'),
                            Tissue = 'Walker',
                            Number = c(both_n,
                                       twas_n,
                                       isotwas_n)))
  
  rm(twas,isotwas,total)
  
}
df_out = df_out[complete.cases(df_out),]

require(tidyverse)
df_total = df_out %>%
  group_by(Method,Tissue) %>%
  summarise(Number = sum(Number))
data.table::fwrite(df_total,
                   'Walker/numGWASTagged.tsv',
                   sep='\t',
                   col.names=T,
                   row.names=F)


data.table::fwrite(df_out,
                   'Walker/numGWASTagged_ByTrait.tsv',
                   sep='\t',
                   col.names=T,
                   row.names=F)

















require(data.table)
require(vroom)
setwd('Walker')

fff = list.files(recursive = T)
path <- "Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx"
traits = excel_sheets(path = path)


allGenes = c()
for (tr in traits){
  setwd(tr)
  
  fff = list.files()
  splitwas = fread(fff[grepl('spli',fff)])
  allGenes = unique(c(allGenes,splitwas$Gene))
  setwd('Walker')
}


require(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

bm = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol',
                          'chromosome_name','start_position','end_position',
                          'gene_biotype'),
           filters = 'ensembl_gene_id',
           values = allGenes, 
           mart = ensembl)
colnames(bm) = c('Gene','HGNC','Chromosome','Start','End','Biotype')



for (i in 1:length(traits)){
  tr = traits[i]
  
  setwd(tr)
  fff = list.files()
  
  full_isoform_res = readxl::read_excel('Walker/Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx',
                                        sheet = tr)
  
  isotwas_res = fread(fff[grepl('spli',fff)])
  isotwas_res = isotwas_res[complete.cases(isotwas_res),]
  isotwas_res = isotwas_res[order(abs(isotwas_res$Z),decreasing = T),]
  isotwas_res = isotwas_res[!duplicated(isotwas_res$Transcript) & 
                              abs(isotwas_res$Z) < Inf,]
  isotwas_res = merge(bm,isotwas_res,by='Gene')
  full_isoform_res = subset(full_isoform_res,
                            Gene %in% isotwas_res$Gene)
  
  www = which(isotwas_res$Chromosome == 6 &
                isotwas_res$Start < 35e6 &
                isotwas_res$End > 27e6)
  isotwas_res = isotwas_res[-www,]
  
  
  gene = data.frame(Gene = unique(isotwas_res$Gene),
                    HGNC = NA,
                    Chromosome = NA,
                    Start = NA,
                    End = NA,
                    Biotype = NA,
                    Screen.P = NA)
  
  require(tidyverse)
  gene = isotwas_res %>%
    group_by(Gene) %>%
    summarise(HGNC = unique(HGNC),
              Chromosome = unique(Chromosome),
              Start = unique(Start),
              End = unique(End),
              Biotype = unique(Biotype),
              Screen.P = isotwas::p_screen(P))
  
  alpha1=.05
  G = nrow(gene)
  gene$Screen.P.Adjusted = p.adjust(gene$Screen.P,method = 'bonf')
  R = length(unique(gene$Gene[gene$Screen.P.Adjusted < alpha1]))
  alpha2 = (R*alpha1)/G
  isoform_new = as.data.frame(matrix(nrow = 0,
                                     ncol = ncol(isotwas_res)+2))
  colnames(isoform_new) = c(colnames(isotwas_res),'Screen.P','Confirmation.P')
  gene = gene[order(gene$Screen.P),]
  ttt = merge(isotwas_res,
              gene[,c('Gene','Screen.P',
                      'Screen.P.Adjusted')])
  isoform_new = ttt %>%
    group_by(Gene) %>%
    summarise(Transcript = Transcript,
              Confirmation.P = isotwas::p_confirm(P,alpha = alpha2))
  isoform_new = merge(isoform_new,ttt,by=c('Gene','Transcript'))
  isoform_new$Confirmation.P = ifelse(isoform_new$Screen.P.Adjusted < 0.05,
                                      isoform_new$Confirmation.P,
                                      1)
  
  
  
  
  gnomad = 
    fread('/u/scratch/a/abtbhatt/gnomad_plof.txt')
  gnomad$HGNC = gnomad$gene
  
  isoform_sig = subset(isoform_new, 
                       Screen.P < alpha1 &
                         Confirmation.P < alpha2 &
                         permute.P < 0.05)
  
  xlsx::write.xlsx(isoform_sig,
                     file = 
                       'Walker/Developmental_Associations_BeforeFOCUS_spliceTWAS.xlsx',
                     append = T,
                     sheetName = tr,
                     col.names = T,
                     row.names = F)
  
  aaa = data.frame(Threshold = 'Adjusted P < 0.05',
                   Number = c(length(unique(isoform_sig$Gene)),
                              length(unique(full_isoform_res$Gene))),
                   Trait = tr,
                   Feature = c('Gene','Gene'),
                   Method = c('splice-TWAS','isoTWAS'))
  aaa$Tissue = 'Developing brain (PEC)'
  fwrite(aaa,'/u/scratch/a/abtbhatt/Walker/Associations/NumberDevPECWalker_allstudied_spliceIso.tsv',sep='\t',
         append=T,row.names=F,quote=F)
  
  
  setwd('/u/scratch/a/abtbhatt/Walker/Associations')
  
  
}






require(data.table)
require(vroom)
setwd('Walker')
library(readxl)
path <- "Developmental_Associations_BeforeFOCUS_spliceTWAS.xlsx"




library(readxl)
require(LDlinkR)
require(bigsnpr)

snps = snp_attach(snp_readBed2('/u/scratch/a/abtbhatt/PEC_fetal_626_hmp.bed',
                               backingfile = tempfile()))

path <- "Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx"
excel_sheets(path = path)
traits = excel_sheets(path = path)

df_out = data.frame(Trait = NA,
                    Method = NA,
                    Tissue = NA,
                    Number = NA)
for (t in 1:length(traits)){
  print(t)
  tr = traits[t]
  twas_GWAS_SNP = c()
  if (tr %in% excel_sheets(path = 'Developmental_Associations_BeforeFOCUS_spliceTWAS.xlsx')){
    
    twas = read_xlsx('Developmental_Associations_BeforeFOCUS_spliceTWAS.xlsx',
                     sheet = tr)
    twas$Method = 'splice-TWAS'
    colnames(twas)[12:13] = c('Top GWAS SNP',
                              'Top GWAS P')
    
  } else {
    twas = data.frame(Gene = 'Dummy',
                      Chromosome = 1,
                      Start = 1,
                      End = 1,
                      Z = 1,
                      'Top GWAS SNP' = 'Dummy',
                      'Top GWAS P' = 1,
                      Method = 'splice-TWAS')
    colnames(twas)[6:7] = c('Top GWAS SNP',
                            'Top GWAS P')
  }
  
  isotwas_GWAS_SNP = c()
  if (tr %in% excel_sheets('Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx')){
    isotwas = read_xlsx('Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx',
                        sheet = tr)
    isotwas$Method = 'isoTWAS'
  } else {
    
    isotwas = data.frame(Gene = 'Dummy',
                         Chromosome = 1,
                         Start = 1,
                         End = 1,
                         Z = 1,
                         'Top GWAS SNP' = 'Dummy',
                         'Top GWAS P' = 1,
                         Method = 'isoTWAS')
    colnames(isotwas)[6:7] = c('Top GWAS SNP',
                               'Top GWAS P')
  }
  
  total = rbind(twas[,c('Gene','Chromosome','Start','End',
                        'Z','Top GWAS SNP','Top GWAS P','Method')],
                isotwas[,c('Gene','Chromosome','Start','End',
                           'Z','Top GWAS SNP','Top GWAS P','Method')])
  
  total = as.data.frame(total[order(total$Chromosome,
                                    total$Start),])
  total = subset(total,`Top GWAS P` < 5e-8)
  if (nrow(total) > 0){
    total$Group = 1
    ggg = 1
    for (i in 1:(nrow(total)-1)){
      if (total$End[i] <= total$Start[i+1] - 1e6){
        ggg = ggg+1
        total$Group[(i+1):nrow(total)] = ggg
      }
    }
    
    both_n = twas_n = isotwas_n = 0
    for (ggg in unique(total$Group)){
      
      ttt = subset(total,Group == ggg)
      if (all(c('splice-TWAS','isoTWAS') %in% ttt$Method)){
        both_n = both_n + 1
      }
      
      if ('splice-TWAS' %in% ttt$Method &
          !'isoTWAS' %in% ttt$Method){
        twas_n = twas_n + 1
      }
      
      if (!'splice-TWAS' %in% ttt$Method &
          'isoTWAS' %in% ttt$Method){
        isotwas_n = isotwas_n + 1
      }
      
    }
  } else {
    both_n = 0
    twas_n = 0
    isotwas_n = 0
  }
  
  df_out = rbind(df_out,
                 data.frame(Trait = tr,
                            Method = c('Both','splice-TWAS','isoTWAS'),
                            Tissue = 'Walker',
                            Number = c(both_n,
                                       twas_n,
                                       isotwas_n)))
  
  rm(twas,isotwas,total)
  
}
df_out = df_out[complete.cases(df_out),]

require(tidyverse)
df_total = df_out %>%
  group_by(Method,Tissue) %>%
  summarise(Number = sum(Number))
data.table::fwrite(df_total,
                   'Walker/numGWASTagged_splice.tsv',
                   sep='\t',
                   col.names=T,
                   row.names=F)


data.table::fwrite(df_out,
                   'Walker/numGWASTagged_splice_byTrait.tsv',
                   sep='\t',
                   col.names=T,
                   row.names=F)


