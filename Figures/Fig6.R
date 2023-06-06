require(data.table)
require(ggplot2)
require(bigsnpr)

dir.create('isoqtlPlots')
setwd('/u/scratch/a/abtbhatt/isoqtlPlots')
eqtl_file = 'hcp175_cis_gene.txt'
isoqtl_file = 'hcp150_cis_isoform.txt'
gen_file = 'BB_hmp_forQTL.bed'
tx_exp_file = 'tx_exp_bigbrain_0503_rint.bed.gz'
gene_exp_file = 'gene_exp_bigbrain_0503_rint.bed.gz'

examples = fread('tryTheseisoQTLPlots_InGene.tsv')

examples = subset(examples,(InGene == 1 &
                              `Top GWAS P` < 5e-8 &
                              Dataset == 'Adult') | 
                    HGNC %in% c('AKT3',
                                'KMT5A',
                                'SFMBT1'))
examples = subset(examples,chr_name %in% 1:22)
examples = examples[!duplicated(examples),]


getIsoQTLPlot = function(i,isoform = NULL){
  
  gene = examples$Gene[i]
  if (is.null(isoform)){
    isoform = examples$Transcript[i]
  }
  trait = examples$Trait[i]
  hgnc = examples$HGNC[i]
  print(hgnc)
  print(trait)
  
  outFile = paste0('index',i,'_',
                   trait,'_Gene',hgnc,'_Isoform',
                   isoform,'_QTL')
  
  gwasfile = file.path('./01-reheadered-sumstats',
                       list.files('./01-reheadered-sumstats'))
  gwasfile = gwasfile[grepl(trait,gwasfile)]
  gwasfile = gwasfile[!grepl('/z_',gwasfile)]
  
  if (trait == 'SCZ'){
    gwasfile = gwasfile[grepl('2021',gwasfile)]
  }
  
  gwasSummStats_hmp = 
    gwasfile[1]
  
  ### extract gene and isoform qtls
  system(paste('grep',gene,eqtl_file,' > eqtl_this.txt'))
  system(paste('grep',isoform,isoqtl_file,' > isoqtl_this.txt'))
  
  eqtl = fread('eqtl_this.txt')
  isoqtl = fread('isoqtl_this.txt')
  gwas = fread(gwasSummStats_hmp)
  
  if (trait == 'NTSM'){
    gwas$FRQ = gwas$EAF_UKB
  }
  
  if (trait %in% c('BV','NTSM')){
    gwas$BETA = gwas$Z/(sqrt(2*gwas$FRQ*(1-gwas$FRQ)*(gwas$N + gwas$Z)))
    gwas$SE = 1/(sqrt(2*gwas$FRQ*(1-gwas$FRQ)*(gwas$N + gwas$Z)))
  }
  
  colnames(eqtl) = colnames(isoqtl) = 
    c('Phenotype',
      'Chromosome',
      'Start',
      'End',
      'Strand',
      'TotalN',
      'Distance',
      'SNP',
      'Chromosome_SNP',
      'Position',
      'Position2',
      'P',
      'R2',
      'Beta',
      'SE',
      'Top')
  eqtl$Type = 'Gene'
  isoqtl$Type = 'Isoform'
  
  require(biomaRt)
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                    host="grch37.ensembl.org", 
                    path="/biomart/martservice" ,
                    dataset="hsapiens_gene_ensembl")
  bm = getBM(attributes = c('ensembl_gene_id','chromosome_name',
                            'start_position','end_position','hgnc_symbol'),
             filters = c('ensembl_gene_id'),
             values = gene,
             mart = ensembl)
  
  if ('POS' %in% colnames(gwas) & !'BP' %in% colnames(gwas)){
    gwas$BP = gwas$POS
  }
  gwas = subset(gwas,CHR == bm$chromosome_name[1] &
                  BP < bm$start_position[1] + 5e5 &
                  BP > max(1,bm$end_position[1] - 5e5))
  
  if (!'OR' %in% colnames(gwas) &
      'BETA' %in% colnames(gwas)){
    gwas$OR = exp(gwas$BETA)
  }
  gwas = gwas[,c('SNP','BP','A1','A2','OR','SE','P','N')]
  
  eqtl = eqtl[,c('Phenotype',
                 'Chromosome','Start','End',
                 'SNP','Position',
                 'P','Beta','SE')]
  colnames(eqtl)[-5] = paste0(colnames(eqtl)[-5],'_Gene')
  isoqtl = isoqtl[,c('Phenotype',
                     'Chromosome','Start','End',
                     'SNP','Position',
                     'P','Beta','SE')]
  colnames(isoqtl)[-5] = paste0(colnames(isoqtl)[-5],'_Isoform')
  
  total_all = merge(merge(eqtl,isoqtl,
                          by = 'SNP'),gwas,by='SNP',all.y = T)
  
  
  ### Make coloc plot
  total_all_melt = reshape2::melt(
    total_all,
    measure.vars = c('P','P_Gene','P_Isoform')
  )
  
  
  colnames(total_all_melt)[22:23] = c('Phenotype','P')
  total_all_melt$Phenotype = as.factor(total_all_melt$Phenotype)
  levels(total_all_melt$Phenotype) = c(trait,
                                       hgnc,
                                       isoform)
  
  total_all_melt = total_all_melt[,c('SNP','A1','A2','BP','Chromosome_Gene',
                                     'Start_Gene','End_Gene','Position_Gene',
                                     'Phenotype','P')]
  total_all_melt$Sig = ifelse(total_all_melt$Phenotype == trait,
                              -log10(5e-8),
                              -log10(1e-6))
  
  chr = unique(total_all_melt$Chromosome_Gene)[!is.na(unique(total_all_melt$Chromosome_Gene))]
  snps = snp_attach(snp_readBed2(paste0('./plink_1KG_bypop/EUR/chr',
                                        chr,'.bed'),
                                 backingfile = tempfile()))
  snps = snp_attach(subset(snps,
                           ind.col = which(snps$map$marker.ID %in%
                                             total_all_melt$SNP)))
  
  ld_all = snp_cor(snps$genotypes)
  colnames(ld_all) = rownames(ld_all) = snps$map$marker.ID
  total_all_melt = subset(total_all_melt,SNP %in% colnames(ld_all))
  ld_all = ld_all[total_all_melt$SNP,
                  total_all_melt$SNP]
  
  best_isoform_snp = 
    subset(total_all_melt,
           Phenotype == isoform)$SNP[which.min(subset(total_all_melt,
                                                      Phenotype == isoform)$P)]
  total_all_melt$Label = ifelse(total_all_melt$SNP == 
                                  best_isoform_snp,
                                total_all_melt$SNP,'')
  
  ld_df = data.frame(marker.ID = colnames(ld_all),
                     LD = as.numeric(ld_all[isoqtl$SNP[which.min(isoqtl$P_Isoform)],]))
  
  ld_df = merge(ld_df,
                snps$map,
                by='marker.ID')
  colnames(ld_df)[c(1,2,5)] = c('SNP','LD','FinalPosition')
  
  total_all_melt = merge(total_all_melt,
                         ld_df,
                         by = 'SNP')
  total_all_melt$LDBlock = cut(total_all_melt$LD, 
                               breaks = c(0,.2,.4,.6,.8,1),
                               include.lowest=TRUE)
  
  all_out = as.data.frame(matrix(ncol = ncol(total_all_melt),
                                 nrow = 0))
  colnames(all_out) = colnames(total_all_melt)
  
  for (q in unique(total_all_melt$Phenotype)){
    
    this_out = subset(total_all_melt,Phenotype == q)
    this_out = this_out[!duplicated(this_out),]
    this_out$Label[duplicated(this_out$Label)] = ''
    all_out = rbind(all_out,this_out)
    
  }
  
  require(liftOver)
  ch = import.chain('hg19ToHg38.over.chain')
  df.snp <- data.frame(chr=paste0('chr',all_out$chromosome), 
                       start=all_out$FinalPosition, 
                       end=all_out$FinalPosition, 
                       score=1:nrow(all_out),
                       id = all_out$SNP)
  
  gr  = makeGRangesFromDataFrame(df.snp, 
                                 ignore.strand=TRUE,
                                 keep.extra.columns=TRUE)
  cur38 = as.data.frame(unlist(liftOver(gr, ch)))
  
  genelocs = data.frame(SNP = cur38$id,
                        chr = sapply(strsplit(as.character(cur38$seqnames),'r'),
                                     function(x) x[2]),
                        FinalPosition_38 = cur38$start,
                        right = cur38$end)
  genelocs = genelocs[!duplicated(genelocs),]
  all_out = merge(all_out,
                  genelocs[,c('SNP','FinalPosition_38')],
                  by = 'SNP')
  all_out$Trait = trait
  all_out$Gene = hgnc
  all_out$Isoform = isoform
  fwrite(all_out,
         '/u/scratch/a/abtbhatt/isoqtlPlots/finalSummStats.tsv',
         append=T,row.names=F,sep='\t')
  
  
  
  require(ggrepel)
  locuszoom = ggplot(data = subset(all_out,!is.na(LDBlock)),
                     aes(x = FinalPosition_38,
                         y = -log10(P),
                         label = Label,
                         color = LDBlock)) +
    facet_wrap(~Phenotype,ncol = 1,scales = 'free_y') +
    geom_hline(aes(yintercept = Sig),
               linetype = 2,
               color = 'grey') +
    geom_vline(xintercept = c(eqtl$Start_Gene[1],
                              eqtl$End_Gene[1]),
               linetype = 2,
               color = 'grey') +
    geom_point(aes(size = ifelse(Label != '',
                                 'not lead','lead'),
                   shape = ifelse(Label != '',
                                  'not lead','lead'))) +
    theme_minimal() +
    theme(axis.text=element_text(size=7.5),
          axis.title=element_text(size=8),
          plot.title = element_text(size = 8),
          legend.title=element_text(size=8),
          legend.text=element_text(size=7),
          strip.text = element_text(size=8),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = NA,
                                      fill = NA, size = .1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.line = element_line(colour = "black")) +
    geom_label_repel(size = 2,
                     color = 'black',
                     max.overlaps = 1000000,
                     force = 10) +
    scale_size_discrete("level", range=c(.8,2)) +
    guides(size = 'none',
           shape = 'none') +
    xlab(paste0('Position on Chromosome ',
                unique(total_all_melt$Chromosome_Gene)[!is.na(unique(total_all_melt$Chromosome_Gene))])) +
    ylab(expression(-log[10]~"P-value")) +
    theme(legend.position = 'bottom',
          legend.direction = 'horizontal') +
    labs(color = 'LD') +
    scale_color_manual(values = c('navyblue',
                                  'lightblue',
                                  'darkgreen',
                                  'orange',
                                  'red'))
  
  require(cowplot)
  ld_legend = get_legend(locuszoom)
  
  locuszoom_tot = plot_grid(locuszoom + guides(color = 'none'),
                            ld_legend,
                            ncol=1,
                            rel_heights = c(1,.05))
  
  
  ## BOXPLOT
  isoform_exp = fread(tx_exp_file)
  isoform_exp = subset(isoform_exp,
                       id == isoform)
  gene_exp = fread(gene_exp_file)
  gene_exp = subset(gene_exp,
                    id == gene)
  snps = snp_attach(snp_readBed2(gen_file,
                                 backingfile = tempfile()))
  snps = snp_attach(subset(snps,
                           ind.col = which(snps$map$marker.ID %in%
                                             best_isoform_snp)))
  
  
  boxplot_df = data.frame(ID = snps$fam$family.ID,
                          SNP = as.numeric(snps$genotypes[]))
  bp_gene_df = data.frame(ID = colnames(gene_exp)[7:ncol(gene_exp)],
                          Gene =
                            as.numeric(gene_exp[1,
                                                7:ncol(gene_exp)]))
  bp_tx_df = data.frame(ID = colnames(isoform_exp)[7:ncol(isoform_exp)],
                        Isoform =
                          as.numeric(isoform_exp[1,
                                                 7:ncol(isoform_exp)]))
  boxplot_df = merge(merge(boxplot_df,bp_gene_df,by='ID'),
                     bp_tx_df,by='ID')
  boxplot_df_melt = reshape2::melt(boxplot_df,
                                   measure.vars = c('Gene','Isoform'))
  
  colnames(boxplot_df_melt) = c('ID','SNP','Phenotype','Expression')
  boxplot_df_melt$Phenotype = factor(boxplot_df_melt$Phenotype,
                                     c('Gene','Isoform'))
  levels(boxplot_df_melt$Phenotype) = c(paste0(bm$hgnc_symbol[1]),
                                        paste0(isoform))
  boxplot_df_melt$Phenotype = relevel(boxplot_df_melt$Phenotype,
                                      ref = paste0(bm$hgnc_symbol[1]))
  boxplot_df_melt$Genotype = as.factor(boxplot_df_melt$SNP)
  ref = all_out[all_out$SNP == best_isoform_snp,]$A1[1]
  alt = all_out[all_out$SNP == best_isoform_snp,]$A2[1]
  levels(boxplot_df_melt$Genotype) = c(paste0(ref,'/',ref),
                                       paste0(ref,'/',alt),
                                       paste0(alt,'/',alt))
  reg_df = data.frame(Phenotype = c(paste0(bm$hgnc_symbol[1]),
                                    paste0(isoform)),
                      Slope = c(coef(lm(Gene ~ SNP,
                                        data = boxplot_df))[2],
                                coef(lm(Isoform ~ SNP,
                                        data = boxplot_df))[2]),
                      Intercept = c(coef(lm(Gene ~ SNP,
                                            data = boxplot_df))[1],
                                    coef(lm(Isoform ~ SNP,
                                            data = boxplot_df))[1]))
  reg_df$Phenotype = factor(reg_df$Phenotype,
                            levels = c(paste0(bm$hgnc_symbol[1]),
                                       paste0(isoform)))
  
  
  boxplot_full = ggplot(data = boxplot_df_melt,
                        aes(x = Genotype,
                            y = Expression,
                            fill = Phenotype)) +
    geom_boxplot(outlier.shape = NA,
                 alpha = .3) +
    geom_jitter(width = 0.2,
                size = .5,
                alpha = .2) +
    facet_wrap(~Phenotype,
               ncol = 2) +
    theme_minimal() +
    theme(axis.text=element_text(size=7),
          axis.title=element_text(size=8),
          plot.title = element_text(size = 8),
          legend.title=element_text(size=8),
          legend.text=element_text(size=7),
          strip.text = element_text(size=7),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = NA,
                                      fill = NA, size = .1),
          legend.position="bottom",
          legend.box = "horizontal",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(n=3,'Set1')) +
    guides(fill = 'none') +
    geom_abline(data = reg_df,
                aes(slope = Slope,
                    intercept = Intercept),
                color = 'black',
                linetype = 2,
                size = 1)
  
  
  ### make forest plot
  forest_df = data.frame(Phenotype = c(trait,hgnc,isoform),
                         Beta = c(log(gwas$OR[gwas$SNP == best_isoform_snp]),
                                  eqtl$Beta_Gene[eqtl$SNP == best_isoform_snp],
                                  isoqtl$Beta_Isoform[isoqtl$SNP == best_isoform_snp]),
                         SE = c(gwas$SE[gwas$SNP == best_isoform_snp],
                                eqtl$SE_Gene[eqtl$SNP == best_isoform_snp],
                                isoqtl$SE_Isoform[isoqtl$SNP == best_isoform_snp]),
                         P = c(gwas$P[gwas$SNP == best_isoform_snp],
                               eqtl$P_Gene[eqtl$SNP == best_isoform_snp],
                               isoqtl$P_Isoform[isoqtl$SNP == best_isoform_snp]))
  
  forest_df$Lower = forest_df$Beta - abs(qnorm(1-1e-6)) * forest_df$SE
  forest_df$Upper = forest_df$Beta + abs(qnorm(1-1e-6)) * forest_df$SE
  forest_df$Label = signif(forest_df$P,2)
  forest_df$Phenotype = factor(forest_df$Phenotype,
                               c(trait,hgnc,isoform))
  
  forest = ggplot(data = forest_df,
                  aes(x = Phenotype,
                      y = Beta,
                      color = Phenotype)) +
    geom_point() +
    theme_minimal() +
    theme(axis.text=element_text(size=7),
          axis.title=element_text(size=8),
          plot.title = element_text(size = 8),
          legend.title=element_text(size=8),
          legend.text=element_text(size=7),
          strip.text = element_text(size=8),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = NA,
                                      fill = NA, size = .1),
          legend.position="bottom",
          legend.box = "horizontal",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(angle = 45)) +
    guides(color = 'none') +
    xlab('') +
    ylab('Effect') +
    geom_errorbar(aes(ymin = Lower,
                      ymax = Upper),
                  width = .1) +
    scale_color_manual(values = c('black',
                                  RColorBrewer::brewer.pal(n = 3,
                                                           'Set1')[1:2])) +
    geom_hline(yintercept = 0,
               color = 'grey',
               linetype = 2) +
    coord_flip() +
    geom_label_repel(aes(label = prettyNum(P,format='e',digits=2)),
                     size = 2,
                     force = 2,
                     nudge_x = .1)
  
  total_all_plot = plot_grid(locuszoom_tot,
                             plot_grid(boxplot_full,forest,
                                       ncol = 2,
                                       rel_widths = c(1.5,1)),
                             ncol = 1,
                             rel_heights = c(2,1))
  
  
  ggsave(plot = total_all_plot,
         filename = paste0(outFile,'.pdf'),
         height = 5,width = 5)
  
  
  ggsave(plot = total_all_plot,
         filename = paste0(outFile,'.png'),
         height = 5,width = 5)
  
  return(total = total_all_plot)
  
}

akt3_scz = getIsoQTLPlot(19)
cul3_scz = getIsoQTLPlot(1)
hspd1_scz = getIsoQTLPlot(14)
pclo_cdg = getIsoQTLPlot(20)
cnnm_cdg = getIsoQTLPlot(4)
kmt5a_scz = getIsoQTLPlot(8)
sfmbt1_scz = getIsoQTLPlot(11)

akt3_bv = getIsoQTLPlot(3)
akt3_bv_alt = getIsoQTLPlot(3,isoform = 'ENST00000492957')

a = fread('/u/scratch/a/abtbhatt/isoqtlPlots/finalSummStats.tsv')
writexl::write_xlsx(a,
                    '/u/scratch/a/abtbhatt/isoqtlPlots/SupplementalData14.xlsx')

require(cowplot)
total = plot_grid(akt3_scz,
                  cul3_scz,
                  hspd1_scz,
                  pclo_cdg,
                  ncol=2,
                  labels = c('a',
                             'b',
                             'c',
                             'd'))
ggsave(plot = total,
       filename = 'mainFig_Examples.png',
       height = 10,width = 10)


require(cowplot)
akt_bv_total = plot_grid(akt3_bv,akt3_bv_alt,
                         ncol=2,
                         labels = c('a',
                                    'b'))
ggsave(plot = akt_bv_total,
       filename = 'akt3_bv_total.png',
       height = 5,width = 10)


library(magrittr)
library(dplyr)
library(ggplot2)
library(ggtranscript)
gtf <- rtracklayer::import('gencode.v38.annotation.gtf')
gtf$transcript_id = sapply(strsplit(gtf$transcript_id,'[.]'),
                           function(x) x[1])
gtf = subset(gtf, 
             transcript_id %in% c(examples$Transcript[19],
                                  examples$Transcript[3]))
gtf <- gtf %>% dplyr::as_tibble()
akt3_exons <- gtf %>% dplyr::filter(type == "exon")


akt3_tx_plot = akt3_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    aes(fill = transcript_id)
  ) +
  geom_intron(
    data = to_intron(akt3_exons, "transcript_id"),
    aes(strand = strand)
  ) +
  theme_minimal() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=11),
        strip.text = element_text(size=10),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position="bottom",
        legend.box = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(angle = 45)) +
  xlab('Coordinates') +
  ylab('Transcript') +
  guides(fill = 'none')

require(cowplot)
akt_bv_total = plot_grid(plot_grid(akt3_bv,akt3_bv_alt,
                                   ncol=2,
                                   labels = c('a',
                                              'b')),
                         akt3_tx_plot,ncol=1,rel_heights = c(1,.2))
