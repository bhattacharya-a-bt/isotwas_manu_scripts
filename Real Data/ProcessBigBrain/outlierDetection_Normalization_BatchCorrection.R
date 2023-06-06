require(sva)
require(DESeq2)
require(SummarizedExperiment)
require(tximeta)
outFolder = '/u/scratch/a/abtbhatt/renormalize_0429'
setwd(outFolder)

tx_raw = readRDS('AdultBigBrain_tx_exp_raw_042923.RDS')

### REMOVE LOW EXPRESSED TRANSCRIPTS
www = which(rowMeans(assays(tx_raw)[[2]] > 0.1) > 0.25)
tx_raw = tx_raw[www,]

### AGGREGATE TO GENE LEVEL
gene_raw = summarizeToGene(tx_raw)
saveRDS(gene_raw,
        'AdultBigBrain_gene_exp_raw_042923.RDS')

www = which(rowMeans(assays(gene_raw)[[2]] > 0.1) > 0.25)
gene_raw = gene_raw[www,]

print(dim(tx_raw))
print(dim(gene_raw))

### VST NORMALIZATION
tx_cts = DESeq2::vst(round(assays(tx_raw)[[1]]))
gene_cts = DESeq2::vst(round(assays(gene_raw)[[1]]))


### REMOVE OUTLIERS
require(WGCNA)
outliers_gene = c()

### FIRST WITH GENE EXPRESSION
for (b in unique(gene_raw$batch)){
  
  print(b)
  this_batch = which(gene_raw$batch == b)
  gene_cts_batch = gene_cts[,this_batch]
  normadj <- adjacency(gene_cts_batch,type = 'signed',corFnc = 'bicor')   #Calculate network adjacency
  netsummary <- fundamentalNetworkConcepts(normadj)
  C <- netsummary$Connectivity   #Extract connectivity of each sample
  Z.C <- (C-mean(C))/sqrt(var(C))   #Covert to Z-score
  outliers <- (Z.C < -3)
  outliers_gene = c(outliers_gene,
                    names(which(outliers == T)))
  
}

### NOW WITH TX EXPRESSION
for (b in unique(gene_raw$batch)){
  
  print(b)
  this_batch = which(tx_raw$batch == b)
  tx_cts_batch = tx_cts[,this_batch]
  normadj <- adjacency(tx_cts_batch,type = 'signed',corFnc = 'bicor')   #Calculate network adjacency
  netsummary <- fundamentalNetworkConcepts(normadj)
  C <- netsummary$Connectivity   #Extract connectivity of each sample
  Z.C <- (C-mean(C))/sqrt(var(C))   #Covert to Z-score
  outliers <- (Z.C < -3)
  outliers_gene = c(outliers_gene,
                    names(which(outliers == T)))
  
}

length(outliers_gene)
tx_cts = tx_cts[,!colnames(tx_cts) %in% outliers_gene]
gene_cts = gene_cts[,!colnames(gene_cts) %in% outliers_gene]


tx_raw = tx_raw[,!colnames(tx_raw) %in% outliers_gene]
gene_raw = gene_raw[,!colnames(gene_raw) %in% outliers_gene]



### REMOVE GENES WITH 0 VARIANCE IN A BATCH
outdf = as.data.frame(matrix(nrow = nrow(gene_cts),
                             ncol = 10))
colnames(outdf) = c('Gene',sort(
  unique(unlist(colData(gene_raw)$batch))
))
outdf$Gene = rownames(gene_cts)

require(tidyverse)
for (i in 1:nrow(gene_cts)){
  print(i)
  x = gene_cts[i,]
  ddd = data.frame(Gene = x,
                   Study = unlist(colData(gene_raw)$batch))
  tot = ddd %>%
    group_by(Study) %>%
    summarise(Variance = var(Gene))
  outdf[i,-1] = tot$Variance
}

rrr = which(rowSums(outdf[,-1] != 0) == 9)
gene_cts = gene_cts[rrr,]
gene_raw = gene_raw[rrr,]



outdf = as.data.frame(matrix(nrow = nrow(tx_cts),
                             ncol = 10))
colnames(outdf) = c('Gene',sort(
  unique(unlist(colData(tx_raw)$batch))
))
outdf$Gene = rownames(tx_cts)

for (i in 1:nrow(outdf)){

  print(i)
  ddd = data.frame(Gene = tx_cts[i,],
                   Study = unlist(colData(tx_raw)$batch))
  tot = ddd %>%
    group_by(Study) %>%
    summarise(Variance = var(Gene))
  outdf[i,-1] = tot$Variance
}

rrr = which(rowSums(outdf[,-1] != 0) == 9)
tx_cts = tx_cts[rrr,]
tx_raw = tx_raw[rrr,]

### BATCH CORRECTION WITH COMBAT
tx_cbt = ComBat(
  tx_cts,
  batch = colData(tx_raw)$batch
)
gene_cbt = ComBat(
  gene_cts,
  batch = colData(gene_raw)$batch
)


### WRITE OUT NORMALIZED, BATCH CORRECTED EXPRESSION
require(biomaRt)
tx_bed = as.data.frame(cbind(pid = rep(NA,nrow(tx_cbt)),
                             gid = rep(NA,nrow(tx_cbt)),
                             strand = '+'))
tx_bed = cbind(tx_bed,tx_cbt)
tx_bed$gid = as.character(pbapply::pbsapply(strsplit(as.character(unlist(rowData(tx_raw)$gene_id)),
                                                         '[.]'),
                                                function(x) x[1]))
tx_bed$pid = pbapply::pbsapply(strsplit(as.character(unlist(rowData(tx_raw)$tx_name)),
                                          '[.]'),
                                 function(x) x[1])


gene_bed = as.data.frame(cbind(pid = rep(NA,nrow(gene_cts)),
                               gid = rep(NA,nrow(gene_cts)),
                               strand = '+'))
gene_bed = cbind(gene_bed,gene_cbt)
gene_bed$pid = gene_bed$gid = sapply(strsplit(unlist(rowData(gene_raw)$gene_id),
                                              '[.]'),
                                     function(x) x[1])


ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl")
bm_gene = getBM(attributes = c('ensembl_gene_id',
                               'chromosome_name',
                               'start_position',
                               'end_position'),
                filters = c('ensembl_gene_id'),
                mart = ensembl,
                values = unique(c(tx_bed$gid,
                                  gene_bed$gid)))
colnames(bm_gene) = c('gid','#Chr','start','end')
gene_bed = merge(bm_gene,gene_bed,by='gid')
tx_bed = merge(bm_gene,tx_bed,by='gid')
gene_bed = subset(gene_bed,`#Chr` %in% 1:22)
tx_bed = subset(tx_bed,`#Chr` %in% 1:22)

gene_bed = gene_bed[,c('#Chr','start','end','pid','gid','strand',
                       colnames(gene_bed)[7:ncol(gene_bed)])]
gene_bed = gene_bed[order(gene_bed$`#Chr`,
                          gene_bed$start),]
require(data.table)
fwrite(gene_bed,'/u/scratch/a/abtbhatt/renormalize_0429/gene_exp_bigbrain_042923_nohcp.bed',
       sep='\t',col.names=T,row.names=F)

tx_bed = tx_bed[,c('#Chr','start','end','pid','gid','strand',
                       colnames(tx_bed)[7:ncol(tx_bed)])]
tx_bed = tx_bed[order(tx_bed$`#Chr`,
                          tx_bed$start),]

fwrite(tx_bed,'/u/scratch/a/abtbhatt/renormalize_0429/tx_exp_bigbrain_042923_nohcp.bed',
       sep='\t',col.names=T,row.names=F)



