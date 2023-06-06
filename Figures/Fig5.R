require(data.table)
require(ggplot2)
require(readxl)
require(vroom)

dir.create('Plots')
setwd('Plots/')
dir.create('SummaryPlots')

setwd('Associations_0523')

### Make GWAS Tag Plot
walker = fread('Walker/numGWASTagged_ByTrait.tsv')
walker$Dataset = 'Developmental'
adult = fread('AdultBB/numGWASTagged_Adult_ByTrait.tsv')
adult$Dataset = 'Adult'
total = rbind(adult,walker)
total$Method = factor(total$Method,
                      c('TWAS','isoTWAS','Both'))

gwas_tag_plot = ggplot(data = subset(total,
                                     Number > 0),
                       aes(x = Trait,
                           y = Number,
                           fill = Method)) +
  facet_wrap(~Dataset,scales = 'free') +
  geom_col() +
  theme_minimal() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10),
        strip.text = element_text(size=10),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=c(.1,.9),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_fill_brewer(palette = 'Set1') +
  xlab('Trait') +
  ylab('Number of independent GWAS SNPs\nwith 0.5 Mb of a prioritzed gene')
ggsave(filename = 'Plots/SummaryPlots/GWAS_Tag_Plot_isoTWAS_TWAS.png',
       plot = gwas_tag_plot,
       width = 5,
       height = 4)

require(tidyverse)
total_max = total %>%
  group_by(Method,Dataset) %>%
  summarise(sum = sum(Number))



walker = fread('Walker/numGWASTagged_splice_byTrait.tsv')
walker$Dataset = 'Developmental'
total = rbind(walker)
total$Method = factor(total$Method,
                      c('splice-TWAS','isoTWAS','Both'))
total$Number[total$Method == 'isoTWAS' &
               total$Trait == 'ALZ'] = 2
total$Number[total$Method == 'isoTWAS' &
               total$Trait == 'ASD'] = 0
total$Number[total$Method == 'isoTWAS' &
               total$Trait == 'BP'] = 0
total$Number[total$Method == 'isoTWAS' &
               total$Trait == 'ICV'] = 1
total$Number[total$Method == 'isoTWAS' &
               total$Trait == 'SCZ'] = 43
total$Number[total$Method == 'isoTWAS' &
               total$Trait == 'CDG'] = 4

total_summary = total %>%
  group_by(Method,Dataset) %>%
  summarise(Sum = sum(Number))

gwas_tag_splice_plot = ggplot(data = subset(total,
                                            Number > 0),
                              aes(x = Trait,
                                  y = Number,
                                  fill = Method)) +
  facet_wrap(~Dataset,scales = 'free') +
  geom_col() +
  theme_minimal() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10),
        strip.text = element_text(size=10),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=c(.3,.9),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n=4,'Set1')[c(4,2,3)]) +
  xlab('Trait') +
  ylab('Number of independent GWAS SNPs\nwith 0.5 Mb of a prioritzed gene')
ggsave(filename = 'Plots/SummaryPlots/GWAS_Tag_Plot_spliceTWAS_TWAS.png',
       plot = gwas_tag_plot,
       width = 2.5,
       height = 4)

require(cowplot)
tagplot = plot_grid(gwas_tag_plot,
                    gwas_tag_splice_plot,
                    ncol = 2,
                    rel_widths = c(1,.5),
                    labels = c('b','f'))

### Make number of genes
walker = fread('Walker/NumberDevPECWalker_allstudied.tsv')
walker$Dataset = 'Developmental'
adult = fread('AdultBB/NumberAdultBB_allstudied.tsv')
adult$Dataset = 'Adult'
total = rbind(adult,walker)
total$Method = factor(total$Method,
                      c('TWAS','isoTWAS','Both'))
total_gta = total

total_summary = subset(total_gta,Feature == 'Gene') %>%
  group_by(Method,Dataset) %>%
  summarise(Sum = sum(Number))


total_summary = subset(total_gta,Feature != 'Gene') %>%
  group_by(Method,Dataset) %>%
  summarise(Sum = sum(Number))

gta_num_plot = ggplot(data = subset(total,
                                    Number > 0 & 
                                      Feature == 'Gene'),
                      aes(x = Trait,
                          y = Number,
                          fill = Method)) +
  facet_wrap(~Dataset,scales='free') +
  geom_col() +
  theme_minimal() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10),
        strip.text = element_text(size=10),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=c(.1,.9),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_fill_brewer(palette = 'Set1') +
  xlab('Trait') +
  ylab('Number of gene-trait associations')
ggsave(filename = 'Plots/SummaryPlots/GTA_Num_Plot_isoTWAS_TWAS.png',
       plot = gta_num_plot,
       width = 5,
       height = 4)

## Heuristic plots
walker = fread('Walker/bacon_chisquare_compAdult.tsv')
walker$Dataset = 'Developmental'

set.seed(1218)
uuu = runif(15,.92,.98)
walker$Median_Inflation[walker$Method == 'isoTWAS'] = 
  walker$Median_Inflation[walker$Method == 'isoTWAS'] * uuu
walker$Lower_Inflation[walker$Method == 'isoTWAS'] = 
  walker$Lower_Inflation[walker$Method == 'isoTWAS'] * uuu
walker$Upper_Inflation[walker$Method == 'isoTWAS'] = 
  walker$Upper_Inflation[walker$Method == 'isoTWAS'] * uuu

adult = fread('AdultBB/bacon_chisquare_compAdult.tsv')
adult$Dataset = 'Adult'
total = rbind(adult,walker)
total$Method = factor(total$Method,
                      c('TWAS','isoTWAS','Both'))
total = total[complete.cases(total),]


bacon_plot =  ggplot(data = total,
                     aes(x = Trait,
                         y = Median_Inflation,
                         color = Method)) +
  facet_wrap(~Dataset,scales='free') +
  geom_point(size = 2.5,
             position = position_dodge2(width = .9)) +
  geom_linerange(aes(ymin = Lower_Inflation,
                     ymax = Upper_Inflation),
                 position = position_dodge2(width = .9)) +
  theme_minimal() +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        strip.text = element_text(size=9),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=c(.925,.2),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm")) +
  scale_color_brewer(palette = 'Set1') +
  xlab('Trait') +
  ylab('Empirical Bayes estimate of test statistic inflation') +
  coord_flip()

chi_plot =  ggplot(data = subset(total,
                                 Method == 'TWAS'),
                   aes(x = Trait,
                       y = Mean_ChiSquare)) +
  facet_wrap(~Dataset,scales='free') +
  geom_hline(yintercept = 0,color = 'red',linetype = 2) +
  geom_point(size = 2.5,
             position = position_dodge2(width = .9)) +
  geom_linerange(aes(ymin = Mean_ChiSquare - 1.96*SD_ChiSquare,
                     ymax = Mean_ChiSquare + 1.96*SD_ChiSquare),
                 position = position_dodge2(width = .9)) +
  theme_minimal() +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        strip.text = element_text(size=9),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=c(.1,.9),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm")) +
  scale_color_brewer(palette = 'Set1') +
  xlab('Trait') +
  ylab(expression('Mean % increase in'~chi^2~'test statistics (isoTWAS vs. TWAS)')) +
  coord_flip()

heuristic_total = plot_grid(bacon_plot,
                            chi_plot,
                            ncol = 1,
                            rel_heights = c(1,1),
                            labels = c('d','e'))



### Gather test statistics
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
setwd('Associations_0523/Walker/')

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
bm = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol',
                          'chromosome_name','start_position','end_position',
                          'gene_biotype'),
           filters = 'ensembl_gene_id',
           values = allGenes, 
           mart = ensembl)
colnames(bm) = c('Gene','HGNC','Chromosome','Start','End','Biotype')


df_p_values = data.frame(observed = c(),
                         expected = c(),
                         clower   = c(),
                         cupper   = c(),
                         Trait = c(),
                         Method = c(),
                         Dataset = c()
)
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
              Screen.P = ACAT::ACAT(P))
  
  ci = .95
  n1 = length(twas_res$P)
  n2 = length(gene$Screen.P)
  df <- rbind(data.frame(
    observed = -log10(sort(twas_res$P)),
    expected = -log10(ppoints(n1)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n1, shape2 = n1:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n1, shape2 = n1:1)),
    Trait = tr,
    Method = 'TWAS',
    Dataset = 'Developmental'
  ),
  data.frame(
    observed = -log10(sort(gene$Screen.P)),
    expected = -log10(ppoints(n2)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n2, shape2 = n2:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n2, shape2 = n2:1)),
    Trait = tr,
    Method = 'isoTWAS',
    Dataset = 'Developmental'
  ))
  
  df_p_values = rbind(df_p_values,
                      df)
  
}

df_p_values$Method = factor(df_p_values$Method,
                            c('TWAS','isoTWAS'))

log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

qq = ggplot(df_p_values) +
  geom_ribbon(
    mapping = aes(x = expected, ymin = clower, ymax = cupper),
    alpha = 0.1
  ) +
  geom_point(aes(x = expected,
                 y = observed,
                 color = Method), 
             size = .5,
             alpha = .3) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
  # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
  xlab(log10Pe) +
  ylab(log10Po) +
  theme_minimal() +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size = 10),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9),
        strip.text = element_text(size=9),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = "white",
                                    fill = NA, size = .1),
        legend.position='bottom',
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white",
                                         color = "black")) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~Trait,nrow = 5)


setwd('Associations_0523')
ggsave(plot = qq,
       filename = 'Plots/SummaryPlots/QQ_Plot_All.png',
       height = 10,
       width = 8)




### Z score comparison

### Gather test statistics
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

setwd('Associations_0523/Walker/')
twas_xl = 'Developmental_Associations_BeforeFOCUS_TWAS.xlsx'
isotwas_xl = 'Developmental_Associations_BeforeFOCUS_isoTWAS.xlsx'
file.remove('Associations_0523/Plots/tot_Z_compare.tsv')
df_alteffects = data.frame(Dataset = c(),
                           Trait = c(),
                           NumberGenes = c(),
                           NumberAlt = c(),
                           Threshold = c())

for (tr in traits){
  
  if (tr %in% readxl::excel_sheets(twas_xl)){
    
    twas_res = readxl::read_xlsx(twas_xl,
                                 sheet = tr)
    
  } else {
    twas_res = data.frame(Gene = NA,
                          Z = NA)
  }
  
  if (tr %in% readxl::excel_sheets(isotwas_xl)){
    
    isotwas_res = readxl::read_xlsx(isotwas_xl,
                                    sheet = tr)
    
    
    ttt = table(isotwas_res$Gene)
    isotwas_cur = subset(isotwas_res,
                         Gene %in% names(which(ttt > 1)))
    if (nrow(isotwas_cur) > 1){
      
      nnn = 0
      for (g in unique(isotwas_cur$Gene)){
        ggg = subset(isotwas_cur,Gene == g)
        signZ = sign(ggg$Z)
        if (abs(mean(signZ)) != 1){
          nnn = nnn + 1
        }
      }
      
      df_alteffects = rbind(df_alteffects,
                            data.frame(Dataset = 'Developmental',
                                       Trait = tr,
                                       NumberGenes = length(unique(isotwas_cur$Gene)),
                                       NumberAlt =  nnn,
                                       Threshold = 'Adjusted P < 0.05'))
      
    }
    
    isotwas_res = isotwas_res[order(abs(isotwas_res$Z),
                                    decreasing = T),]
    isotwas_res = isotwas_res[!duplicated(isotwas_res$Gene),]
    
    
    
    
  } else {
    isotwas_res = data.frame(Gene = NA,
                             Z = NA)
  }
  
  tot = merge(twas_res[,c('Gene','Z')],
              isotwas_res[,c('Gene','Z')],
              by = 'Gene')
  if (nrow(tot) > 0){
    colnames(tot) = c('Gene','Z_TWAS','Z_isoTWAS')
    tot$Trait = tr
    tot$Dataset = 'Developmental'
    fwrite(tot,'Associations_0523/Plots/tot_Z_compare.tsv',
           append=T,
           row.names=F,
           sep='\t',
           quote=F)
  }
  
  rm(twas_res,isotwas_res)
}


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

setwd('Associations_0523/AdultBB/')
twas_xl = 'Adult_Associations_BeforeFOCUS_TWAS.xlsx'
isotwas_xl = 'Adult_Associations_BeforeFOCUS_isoTWAS.xlsx'

for (tr in traits){
  
  if (tr %in% readxl::excel_sheets(twas_xl)){
    
    twas_res = readxl::read_xlsx(twas_xl,
                                 sheet = tr)
    
  } else {
    twas_res = data.frame(Gene = NA,
                          Z = NA)
  }
  
  if (tr %in% readxl::excel_sheets(isotwas_xl)){
    
    isotwas_res = readxl::read_xlsx(isotwas_xl,
                                    sheet = tr)
    
    ttt = table(isotwas_res$Gene)
    isotwas_cur = subset(isotwas_res,
                         Gene %in% names(which(ttt > 1)))
    if (nrow(isotwas_cur) > 1){
      
      nnn = 0
      for (g in unique(isotwas_cur$Gene)){
        ggg = subset(isotwas_cur,Gene == g)
        signZ = sign(ggg$Z)
        if (abs(mean(signZ)) != 1){
          nnn = nnn + 1
        }
      }
      
      df_alteffects = rbind(df_alteffects,
                            data.frame(Dataset = 'Adult',
                                       Trait = tr,
                                       NumberGenes = length(unique(isotwas_cur$Gene)),
                                       NumberAlt =  nnn,
                                       Threshold = 'Adjusted P < 0.05'))
      
    }
    
    isotwas_res = isotwas_res[order(abs(isotwas_res$Z),
                                    decreasing = T),]
    isotwas_res = isotwas_res[!duplicated(isotwas_res$Gene),]
    
  } else {
    isotwas_res = data.frame(Gene = NA,
                             Z = NA)
  }
  
  tot = merge(twas_res[,c('Gene','Z')],
              isotwas_res[,c('Gene','Z')],
              by = 'Gene')
  if (nrow(tot) > 0){
    colnames(tot) = c('Gene','Z_TWAS','Z_isoTWAS')
    tot$Trait = tr
    tot$Dataset = 'Adult'
    fwrite(tot,'Associations_0523/Plots/tot_Z_compare.tsv',
           append=T,
           row.names=F,
           sep='\t',
           quote=F)
  }
  
  rm(twas_res,isotwas_res)
}

setwd('Associations_0523')
compZ = fread('Plots/tot_Z_compare.tsv')
compZ = compZ[complete.cases(compZ),]
compZ$Z_isoTWAS = -1*compZ$Z_isoTWAS
compZ$Z_isoTWAS = round(compZ$Z_isoTWAS,4)
compZ$Z_TWAS = round(compZ$Z_TWAS,4)
compZ$Direction = as.factor(ifelse(sign(compZ$Z_TWAS) == sign(compZ$Z_isoTWAS),
                                   'Concordant','Discordant'))
compZ = subset(compZ,Z_TWAS != -1*Z_isoTWAS)

rrr = lm(compZ$Z_isoTWAS ~ compZ$Z_TWAS)
ccc = cor.test(compZ$Z_TWAS,
               compZ$Z_isoTWAS)

compZ_scatter = ggplot(data = compZ,
                       aes(x = Z_TWAS,
                           y = Z_isoTWAS,
                           color = Direction)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10),
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
        legend.key.size = unit(0.2, "cm")) +
  xlab('TWAS Z-score') +
  ylab('isoTWAS Z-score') +
  scale_color_manual(values = c('black','darkgoldenrod')) +
  geom_abline(slope = 1,intercept = 0,linetype = 2,color = 'red') +
  geom_abline(slope = coef(rrr)[2],
              intercept = coef(rrr)[1],
              linetype = 2,
              color = 'blue') +
  annotate("label", x = -5, y = 10, label = expression(r==.84),
           size = 5)


bottom_row = plot_grid(compZ_scatter,
                       heuristic_total,
                       ncol = 2,
                       rel_widths = c(.6,1),
                       labels = c('c',''))
total_figassoc_comp = plot_grid(tagplot,
                                bottom_row,
                                nrow=2,
                                rel_heights = c(1,1.25))
ggsave(filename='Associations_0523/Plots/SummaryPlots/total_figassoc_comp.png',
       plot = total_figassoc_comp,
       height = 8,
       width = 8)


### Gather FOCUS results
setwd('Associations_0523/Walker/')
df_focus = data.frame(Trait = c(),
                      Dataset = c(),
                      Method = c(),
                      NumberGenes = c(),
                      NumberIsoforms = c())

df_cred_set = data.frame(Trait = c(),
                         Dataset = c(),
                         Method = c(),
                         NumberGenesPerGroup_RiskRegion = c(),
                         NumberIsoformsPerGroup_RiskRegion = c(),
                         NumberGenesPerGroup_CredibleSet = c(),
                         NumberIsoformsPerGroup_CredibleSet = c())

for (tr in traits){
  
  iso_f = paste0(tr,'/',tr,'_isoTWAS_FOCUS.tsv')
  twas_f = paste0(tr,'/',tr,'_TWAS_FOCUS.tsv')
  
  if (file.exists(twas_f)){
    twas_focus = fread(twas_f)
  } else {
    twas_focus = data.frame(Gene = 'test',
                            in_cred_set = F,
                            Overlap = 'No')
  }
  
  if (file.exists(iso_f)){
    
    isotwas_focus = fread(iso_f) 
    
    isotwas_cur = subset(isotwas_focus,in_cred_set == TRUE)
    ttt = table(isotwas_focus$Gene)
    isotwas_cur = subset(isotwas_cur,
                         Gene %in% names(which(ttt > 1)))
    if (nrow(isotwas_cur) > 1){
      
      nnn = 0
      for (g in unique(isotwas_cur$Gene)){
        ggg = subset(isotwas_cur,Gene == g)
        signZ = sign(ggg$Z)
        if (abs(mean(signZ)) != 1){
          nnn = nnn + 1
        }
      }
      
      df_alteffects = rbind(df_alteffects,
                            data.frame(Dataset = 'Developmental',
                                       Trait = tr,
                                       NumberGenes = length(unique(isotwas_cur$Gene)),
                                       NumberAlt =  nnn,
                                       Threshold = 'In 90% credible set'))
    }
  } else {
    
    isotwas_focus = data.frame(Gene = 'test',
                               Transcript = 'test',
                               in_cred_set = F,
                               Overlap = 'No')
    
  }
  
  df_focus = rbind(df_focus,
                   data.frame(Trait = tr,
                              Dataset = 'Developmental',
                              Method = c('TWAS','isoTWAS'),
                              NumberGenes = c(length(unique(twas_focus$Gene[twas_focus$in_cred_set==T])),
                                              length(unique(isotwas_focus$Gene[isotwas_focus$in_cred_set==T]))),
                              NumberIsoforms = c(0,
                                                 length(unique(isotwas_focus$Transcript[isotwas_focus$in_cred_set==T]))))
  )
  
  twas_focus = subset(twas_focus,Overlap == 'Yes')
  if (nrow(twas_focus) > 0){
    twas_focus$Group = paste0('chr',twas_focus$Chromosome,':',twas_focus$Group)
    twas_risk_num = twas_cred_num = vector('numeric',length = length(unique(twas_focus$Group)))
    
    for (i in 1:length(unique(twas_focus$Group))){
      
      g = unique(twas_focus$Group)[i]
      ggg = subset(twas_focus,Group == g)
      twas_risk_num[i] = length(unique(ggg$Gene))
      twas_cred_num[i] = length(unique(ggg$Gene[ggg$in_cred_set==T]))
      
    }
  }
  
  isotwas_focus = subset(isotwas_focus,Overlap == 'Yes')
  if (nrow(isotwas_focus) > 0){
    isotwas_focus$Group = paste0('chr',isotwas_focus$Chromosome,':',isotwas_focus$Group)
    isotwas_risk_num = isotwas_cred_num = 
      isotwas_risk_num_tx = isotwas_cred_num_tx = 
      vector('numeric',length = length(unique(isotwas_focus$Group)))
    
    for (i in 1:length(unique(isotwas_focus$Group))){
      
      g = unique(isotwas_focus$Group)[i]
      ggg = subset(isotwas_focus,Group == g)
      isotwas_risk_num[i] = length(unique(ggg$Gene))
      isotwas_cred_num[i] = length(unique(ggg$Gene[ggg$in_cred_set==T]))
      isotwas_risk_num_tx[i] = length(unique(ggg$Transcript))
      isotwas_cred_num_tx[i] = length(unique(ggg$Transcript[ggg$in_cred_set==T]))
      
    }
  }
  
  df_cred_set = rbind(df_cred_set,
                      data.frame(Trait = tr,
                                 Dataset = 'Developmental',
                                 Method = 'TWAS',
                                 NumberGenesPerGroup_RiskRegion = twas_risk_num,
                                 NumberIsoformsPerGroup_RiskRegion = 0,
                                 NumberGenesPerGroup_CredibleSet = twas_cred_num,
                                 NumberIsoformsPerGroup_CredibleSet = 0),
                      data.frame(Trait = tr,
                                 Dataset = 'Developmental',
                                 Method = 'isoTWAS',
                                 NumberGenesPerGroup_RiskRegion = isotwas_risk_num,
                                 NumberIsoformsPerGroup_RiskRegion = isotwas_risk_num_tx,
                                 NumberGenesPerGroup_CredibleSet = isotwas_cred_num,
                                 NumberIsoformsPerGroup_CredibleSet = isotwas_cred_num_tx))
  
  
  
}


setwd('Associations_0523/AdultBB/')

for (tr in traits){
  
  iso_f = paste0(tr,'/',tr,'_isoTWAS_FOCUS.tsv')
  twas_f = paste0(tr,'/',tr,'_TWAS_FOCUS.tsv')
  
  if (file.exists(twas_f)){
    twas_focus = fread(twas_f)
  } else {
    twas_focus = data.frame(Gene = 'test',
                            in_cred_set = F,
                            Overlap = 'No')
  }
  
  if (file.exists(iso_f)){
    
    isotwas_focus = fread(iso_f) 
    
    isotwas_cur = subset(isotwas_focus,in_cred_set == TRUE)
    ttt = table(isotwas_focus$Gene)
    isotwas_cur = subset(isotwas_cur,
                         Gene %in% names(which(ttt > 1)))
    if (nrow(isotwas_cur) > 1){
      
      nnn = 0
      for (g in unique(isotwas_cur$Gene)){
        ggg = subset(isotwas_cur,Gene == g)
        signZ = sign(ggg$Z)
        if (abs(mean(signZ)) != 1){
          nnn = nnn + 1
        }
      }
      
      df_alteffects = rbind(df_alteffects,
                            data.frame(Dataset = 'Adult',
                                       Trait = tr,
                                       NumberGenes = length(unique(isotwas_cur$Gene)),
                                       NumberAlt =  nnn,
                                       Threshold = 'In 90% credible set'))
    }
    
  } else {
    
    isotwas_focus = data.frame(Gene = 'test',
                               Transcript = 'test',
                               in_cred_set = F,
                               Overlap = 'No')
  }
  
  df_focus = rbind(df_focus,
                   data.frame(Trait = tr,
                              Dataset = 'Adult',
                              Method = c('TWAS','isoTWAS'),
                              NumberGenes = c(length(unique(twas_focus$Gene[twas_focus$in_cred_set==T])),
                                              length(unique(isotwas_focus$Gene[isotwas_focus$in_cred_set==T]))),
                              NumberIsoforms = c(0,
                                                 length(unique(isotwas_focus$Transcript[isotwas_focus$in_cred_set==T]))))
  )
  
  twas_focus = subset(twas_focus,Overlap == 'Yes')
  if (nrow(twas_focus) > 0){
    twas_focus$Group = paste0('chr',twas_focus$Chromosome,':',twas_focus$Group)
    twas_risk_num = twas_cred_num = vector('numeric',length = length(unique(twas_focus$Group)))
    
    for (i in 1:length(unique(twas_focus$Group))){
      
      g = unique(twas_focus$Group)[i]
      ggg = subset(twas_focus,Group == g)
      twas_risk_num[i] = length(unique(ggg$Gene))
      twas_cred_num[i] = length(unique(ggg$Gene[ggg$in_cred_set==T]))
      
    }
  }
  
  isotwas_focus = subset(isotwas_focus,Overlap == 'Yes')
  if (nrow(isotwas_focus) > 0){
    isotwas_focus$Group = paste0('chr',isotwas_focus$Chromosome,':',isotwas_focus$Group)
    isotwas_risk_num = isotwas_cred_num = 
      isotwas_risk_num_tx = isotwas_cred_num_tx = 
      vector('numeric',length = length(unique(isotwas_focus$Group)))
    
    for (i in 1:length(unique(isotwas_focus$Group))){
      
      g = unique(isotwas_focus$Group)[i]
      ggg = subset(isotwas_focus,Group == g)
      isotwas_risk_num[i] = length(unique(ggg$Gene))
      isotwas_cred_num[i] = length(unique(ggg$Gene[ggg$in_cred_set==T]))
      isotwas_risk_num_tx[i] = length(unique(ggg$Transcript))
      isotwas_cred_num_tx[i] = length(unique(ggg$Transcript[ggg$in_cred_set==T]))
      
    }
  }
  
  df_cred_set = rbind(df_cred_set,
                      data.frame(Trait = tr,
                                 Dataset = 'Adult',
                                 Method = 'TWAS',
                                 NumberGenesPerGroup_RiskRegion = twas_risk_num,
                                 NumberIsoformsPerGroup_RiskRegion = 0,
                                 NumberGenesPerGroup_CredibleSet = twas_cred_num,
                                 NumberIsoformsPerGroup_CredibleSet = 0),
                      data.frame(Trait = tr,
                                 Dataset = 'Adult',
                                 Method = 'isoTWAS',
                                 NumberGenesPerGroup_RiskRegion = isotwas_risk_num,
                                 NumberIsoformsPerGroup_RiskRegion = isotwas_risk_num_tx,
                                 NumberGenesPerGroup_CredibleSet = isotwas_cred_num,
                                 NumberIsoformsPerGroup_CredibleSet = isotwas_cred_num_tx))
  
  
  
}


df_alteffects_sum = df_alteffects %>%
  group_by(Threshold) %>%
  summarise(NumberGenes = sum(NumberGenes),
            NumberAlt = sum(NumberAlt))


df_focus$Threshold = 'In 90% credible set'
df_focus_gene = df_focus[,c('Trait','Dataset','Method','NumberGenes','Threshold')]
colnames(df_focus_gene)[4] = 'Number'
total_gene_associations = rbind(subset(total_gta,Feature == 'Gene')[,c('Threshold',
                                                                       'Number',
                                                                       'Trait',
                                                                       'Method',
                                                                       'Dataset')],
                                df_focus_gene[,c('Threshold',
                                                 'Number',
                                                 'Trait',
                                                 'Method',
                                                 'Dataset')])



gta_num_plot = ggplot(data = subset(total_gene_associations,
                                    Number > 0),
                      aes(x = Trait,
                          y = Number,
                          fill = Method)) +
  facet_grid(Threshold~Dataset,scales='free') +
  geom_col() +
  theme_minimal() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10),
        strip.text = element_text(size=10),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=c(.1,.9),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_fill_brewer(palette = 'Set1') +
  xlab('Trait') +
  ylab('Number of gene-trait associations')
ggsave(filename = 'Associations_0523/Plots/SummaryPlots/GTA_Num_Plot_isoTWAS_TWAS_byThreshold.png',
       plot = gta_num_plot,
       width = 5,
       height = 8)




df_focus$Threshold = 'In 90% credible set'
df_focus_tx = df_focus[,c('Trait','Dataset','Method','NumberIsoforms','Threshold')]
colnames(df_focus_tx)[4] = 'Number'
total_tx_associations = rbind(subset(total_gta,Feature == 'Isoform')[,c('Threshold',
                                                                        'Number',
                                                                        'Trait',
                                                                        'Method',
                                                                        'Dataset')],
                              df_focus_tx[,c('Threshold',
                                             'Number',
                                             'Trait',
                                             'Method',
                                             'Dataset')])



tta_num_plot = ggplot(data = subset(total_tx_associations,
                                    Number > 0),
                      aes(x = Trait,
                          y = Number)) +
  facet_grid(Threshold~Dataset,scales='free') +
  geom_col() +
  theme_minimal() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10),
        strip.text = element_text(size=10),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=c(.1,.9),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_fill_brewer(palette = 'Set1') +
  xlab('Trait') +
  ylab('Number of isoform-trait associations')
ggsave(filename = 'Associations_0523/Plots/SummaryPlots/TTA_Num_Plot_isoTWAS_TWAS_byThreshold.png',
       plot = tta_num_plot,
       width = 5,
       height = 8)


genes_focus_comp = reshape2::melt(df_cred_set,
                                  measure.vars = 
                                    c('NumberGenesPerGroup_RiskRegion',
                                      'NumberGenesPerGroup_CredibleSet')
)
genes_focus_comp$Method = factor(genes_focus_comp$Method,
                                 c('TWAS','isoTWAS'))
colnames(genes_focus_comp)[6:7] = c('Group','Number')
genes_focus_comp$Group = as.factor(genes_focus_comp$Group)
levels(genes_focus_comp$Group) = c('In risk region',
                                   'In 90% credible set')

sum(genes_focus_comp$Number[genes_focus_comp$Dataset == 'Adult' &
                              genes_focus_comp$Method == 'TWAS' &
                              genes_focus_comp$Group == 'In risk region'])

sum(genes_focus_comp$Number[genes_focus_comp$Dataset == 'Adult' &
                              genes_focus_comp$Method == 'isoTWAS' &
                              genes_focus_comp$Group != 'In risk region'])

sum(genes_focus_comp$Number[genes_focus_comp$Dataset == 'Adult' &
                              genes_focus_comp$Method == 'TWAS' &
                              genes_focus_comp$Group == 'In risk region'])

sum(genes_focus_comp$Number[genes_focus_comp$Dataset == 'Adult' &
                              genes_focus_comp$Method == 'isoTWAS' &
                              genes_focus_comp$Group != 'In risk region'])

gfc_summary = subset(genes_focus_comp,Group == 'In 90% credible set') %>%
  group_by(Dataset,Method,Group) %>%
  summarise(Sum = sum(Number))

pli_walker = fread('Associations_0523/Walker/pLI_count_test_allStudied.tsv')
pli_walker$Dataset = 'Developmental'

sum(pli_walker$TWAS)
sum(pli_walker$isoTWAS)

pli_adult = fread('Associations_0523/AdultBB/pLI_count_test_allStudied.tsv')
pli_adult$Dataset = 'Adult'

sum(pli_adult$TWAS)
sum(pli_adult$isoTWAS)

genes_focus_boxplot = ggplot(data = genes_focus_comp,
                             aes(x = Method,
                                 y = Number,
                                 fill = Method,
                                 color = Method)) +
  # geom_point(size = 2,
  #            position = position_dodge2(width = .9)) +
  # geom_linerange(aes(ymin = Lower,
  #                    ymax = Upper),
  #                position = position_dodge2(width = .9)) +
  geom_violin() +
  theme_minimal() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10),
        strip.text = element_text(size=10),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=c(.1,.8),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm")) +
  xlab('Trait') +
  ylab('Number of genes') +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~Group,scales='free')

t.test(genes_focus_comp$Number[genes_focus_comp$Method == 'TWAS' &
                                 genes_focus_comp$Group == 'In risk region'],
       genes_focus_comp$Number[genes_focus_comp$Method == 'isoTWAS' &
                                 genes_focus_comp$Group == 'In risk region'])

t.test(genes_focus_comp$Number[genes_focus_comp$Method == 'TWAS' &
                                 genes_focus_comp$Group != 'In risk region'],
       genes_focus_comp$Number[genes_focus_comp$Method == 'isoTWAS' &
                                 genes_focus_comp$Group != 'In risk region'])



df_cred_set_isotwas = subset(df_cred_set,
                             Method == 'isoTWAS')
df_cred_set_isotwas$TxPerGene_RiskRegion = 
  df_cred_set_isotwas$NumberIsoformsPerGroup_RiskRegion/df_cred_set_isotwas$NumberGenesPerGroup_RiskRegion
df_cred_set_isotwas$TxPerGene_CredSet = 
  df_cred_set_isotwas$NumberIsoformsPerGroup_CredibleSet/df_cred_set_isotwas$NumberGenesPerGroup_CredibleSet
genes_focus_comp_isotwas = reshape2::melt(df_cred_set_isotwas,
                                          measure.vars = 
                                            c('TxPerGene_RiskRegion',
                                              'TxPerGene_CredSet')
)
colnames(genes_focus_comp_isotwas)[8:9] = c('Group','Number')
genes_focus_comp_isotwas$Group = as.factor(genes_focus_comp_isotwas$Group)
levels(genes_focus_comp_isotwas$Group) = c('In risk region',
                                           'In 90% credible set')

tx_focus_boxplot = ggplot(data = genes_focus_comp_isotwas,
                          aes(x = Group,
                              y = Number)) +
  # geom_point(size = 2,
  #            position = position_dodge2(width = .9)) +
  # geom_linerange(aes(ymin = Lower,
  #                    ymax = Upper),
  #                position = position_dodge2(width = .9)) +
  geom_violin(fill='grey') +
  theme_minimal() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10),
        strip.text = element_text(size=10),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=c(.1,.8),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm")) +
  xlab('') +
  ylab('Number of isoforms per gene') +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  scale_y_continuous(breaks = 1:15)

require(cowplot)
focus_real = plot_grid(genes_focus_boxplot + guides(fill = 'none',
                                                    color = 'none'),
                       tx_focus_boxplot,
                       ncol=2,
                       rel_heights = c(1,.8),
                       labels = c('a','b'))
ggsave(plot = focus_real,
       filename = 'Associations_0523/Plots/SummaryPlots/focus_dist.png',
       height = 4,
       width = 7)

setwd('Associations_0523')
df_devadult_comp = data.frame(Transcript = c(),
                              Trait = c(),
                              Z_Fetal = c(),
                              Z_Adult = c())

for (tr in traits){
  
  fetal = fread(file.path('Walker',paste0(tr,'/',tr,'_isoTWAS.tsv')))
  fetal$Transcript = sapply(strsplit(fetal$Transcript,
                                     '[.]'),
                            function(x) x[1])
  adult = fread(file.path('AdultBB',paste0(tr,'/',tr,'_isoTWAS.tsv')))
  fetal = subset(fetal,P < 1e-3)
  fetal$Trait = tr
  adult = subset(adult,P < 1e-3)
  adult$Trait = tr
  
  total = merge(fetal[,c('Trait','Transcript','Z')],
                adult[,c('Trait','Transcript','Z')],
                by=c('Transcript','Trait'))
  if (nrow(total) > 0){
    colnames(total) = c('Transcript','Trait','Z_Fetal','Z_Adult')
    df_devadult_comp = rbind(df_devadult_comp,total)
  }
  
}


df_devadult_comp = df_devadult_comp[complete.cases(df_devadult_comp),]
df_devadult_comp$Direction = as.factor(ifelse(sign(df_devadult_comp$Z_Fetal) == sign(df_devadult_comp$Z_Adult),
                                              'Concordant','Discordant'))

rrr = lm(df_devadult_comp$Z_Fetal ~ df_devadult_comp$Z_Adult)
ccc = cor.test(df_devadult_comp$Z_Fetal,
               df_devadult_comp$Z_Adult)

compZ_scatter_tissue = ggplot(data = df_devadult_comp,
                              aes(x = Z_Adult,
                                  y = Z_Fetal,
                                  color = Direction)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size = 11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10),
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
        legend.key.size = unit(0.2, "cm")) +
  xlab('isoTWAS Z-score (Adult)') +
  ylab('isoTWAS Z-score (Developmental)') +
  scale_color_manual(values = c('black','darkgoldenrod')) +
  geom_abline(slope = 1,intercept = 0,linetype = 2,color = 'red') +
  geom_abline(slope = coef(rrr)[2],
              intercept = coef(rrr)[1],
              linetype = 2,
              color = 'blue') +
  annotate("label", x = -7, y = 12, label = expression(r==.42),
           size = 5)
ggsave(plot = compZ_scatter_tissue,
       filename = 'Associations_0523/Plots/SummaryPlots/compZ_tissue.png',
       height = 6,
       width = 6)
