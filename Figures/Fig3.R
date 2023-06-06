setwd('GTEx_Models_Raw_042523')
tissues = list.files()
tissues = tissues[!grepl('[.]',
                         tissues)]
tissues = tissues[tissues != 'DoneRaw']

require(data.table)
require(vroom)
require(tidyverse)
require(ggplot2)
require(cowplot)
require(wesanderson)


hash = rbind(fread('GTExScripts/file_location_brain.tsv'),
             fread('GTExScripts/file_location_other.tsv'))
hash = hash[,c('Tissue','SampleSize')]
hash = hash[!duplicated(hash),]
hash$R2Cutoff = 0.01

for (i in 1:nrow(hash)){
  
  print(i)
  samplesize = hash$SampleSize[i]
  r2_seq = seq(0.0001,.05,by=0.00001)
  t = sqrt(abs(r2_seq)/(1-abs(r2_seq)) * (samplesize-2))
  p_val = 1 - pt(t,samplesize-2)
  r2_cutoff = max(0.01,r2_seq[min(which(p_val < 0.05))])
  hash$R2Cutoff[i] = r2_cutoff
  
}

outFile = 'GTEx_Models_Raw_042523/isoformComparison_GTEx.tsv'
file.remove(outFile)

pbmcapply::pbmclapply(tissues,function(f){
  setwd(paste0(f,'/Comparison'))
  if (file.exists('all_tissue_res.tsv')){
    file.remove('all_tissue_res.tsv')
  }
  print(f)
  system(paste('zcat *_isoformcomparison.tsv.gz >> all_tissue_res.tsv'))
  this = vroom::vroom('all_tissue_res.tsv',show_col_types = F)
  this$Tissue = f
  this = subset(this,Gene != 'Gene' &
                  Gene != 'Transcript' &
                  Tissue != 'Tissue')
  data.table::fwrite(this,outFile,sep='\t',append=T,row.names=F,quote=F)
  
  
  setwd('GTEx_Models_Raw_042523')
},mc.cores = 4)


setwd('GTEx_Models_Raw_042523')
outFile = 'GTEx_Models_Raw_042523/geneComparison_GTEx.tsv'
file.remove(outFile)

pbmcapply::pbmclapply(tissues,function(f){
  setwd(paste0(f,'/Comparison'))
  if (file.exists('all_tissue_res.tsv')){
    file.remove('all_tissue_res.tsv')
  }
  print(f)
  system(paste('zcat *_comparison.tsv.gz >> all_tissue_res.tsv'))
  this = vroom::vroom('all_tissue_res.tsv',show_col_types = F)
  this$Tissue = f
  this = subset(this,Gene != 'Gene' &
                  Gene != 'Transcript' &
                  Tissue != 'Tissue')
  data.table::fwrite(this,outFile,sep='\t',append=T,row.names=F,quote=F)
  
  
  setwd('GTEx_Models_Raw_042523')
},mc.cores = 4)


### OVERALL PLOTS
outFolder = '/GTExPlotsTables'
outModel = 'out_models_041923'
dir.create(outFolder)

geneComp = fread('GTEx_Models_Raw_042523/geneComparison_GTEx.tsv')
isoformComp = fread('GTEx_Models_Raw_042523/isoformComparison_GTEx.tsv')
geneComp = merge(geneComp,hash,by='Tissue')
isoformComp = merge(isoformComp,hash,by='Tissue')
geneComp$TissueOut = as.factor(geneComp$Tissue)
isoformComp$TissueOut = as.factor(isoformComp$Tissue)

levels(geneComp$TissueOut) = levels(isoformComp$TissueOut) = 
  c('Adipose (subcutaneous)',
    'Adipose (visc. omentum)',
    'Adrenal gland',
    'Artery (aorta)',
    'Artery (coronary)',
    'Artery (tibial)',
    'Amygdala',
    'Ant. cingulate cortex',
    'Caudate basal ganglia',
    'Cerebellar hemisphere',
    'Cerebellum',
    'Cortex',
    'Frontal cortex',
    'Hippocampus',
    'Hypothalamus',
    'Nucleus accumbens basal ganglia',
    'Putamen basal ganglia',
    'Cervical spinal cord',
    'Substantia nigra',
    'Breast',
    'Fibroblast',
    'Lymphocytes',
    'Colon (sigmoid)',
    'Colon (transverse)',
    'Esophagus (GE junction)',
    'Esophagus (mucosa)',
    'Esophagus (muscularis)',
    'Heart (atrial appendage)',
    'Heart (left ventricle)',
    'Liver',
    'Lung',
    'Minor salivary gland',
    'Muscle (skeletal)',
    'Nerve (tibial)',
    'Ovary',
    'Pancreas',
    'Pituitary',
    'Prostate',
    'Skin (suprapubic)',
    'Skin (lower leg)',
    'Small intestine',
    'Spleen',
    'Stomach',
    'Testis',
    'Thyroid',
    'Uterus',
    'Vagina',
    'Whole blood')

all_tissues = c('Adipose (subcutaneous)',
                'Adipose (visc. omentum)',
                'Adrenal gland',
                'Artery (aorta)',
                'Artery (coronary)',
                'Artery (tibial)',
                'Amygdala',
                'Ant. cingulate cortex',
                'Caudate basal ganglia',
                'Cerebellar hemisphere',
                'Cerebellum',
                'Cortex',
                'Frontal cortex',
                'Hippocampus',
                'Hypothalamus',
                'Nucleus accumbens basal ganglia',
                'Putamen basal ganglia',
                'Cervical spinal cord',
                'Substantia nigra',
                'Breast',
                'Fibroblast',
                'Lymphocytes',
                'Colon (sigmoid)',
                'Colon (transverse)',
                'Esophagus (GE junction)',
                'Esophagus (mucosa)',
                'Esophagus (muscularis)',
                'Heart (atrial appendage)',
                'Heart (left ventricle)',
                'Liver',
                'Lung',
                'Minor salivary gland',
                'Muscle (skeletal)',
                'Nerve (tibial)',
                'Ovary',
                'Pancreas',
                'Pituitary',
                'Prostate',
                'Skin (suprapubic)',
                'Skin (lower leg)',
                'Small intestine',
                'Spleen',
                'Stomach',
                'Testis',
                'Thyroid',
                'Uterus',
                'Vagina',
                'Whole blood')

brain_tissues = c('Amygdala',
                  'Ant. cingulate cortex',
                  'Caudate basal ganglia',
                  'Cerebellar hemisphere',
                  'Cerebellum',
                  'Cortex',
                  'Frontal cortex',
                  'Hippocampus',
                  'Hypothalamus',
                  'Nucleus accumbens basal ganglia',
                  'Putamen basal ganglia',
                  'Cervical spinal cord',
                  'Substantia nigra')

### ISOFORM COMPARISON
### Percent difference isoforms
isoformComp$mrce_lasso = as.numeric(isoformComp$mrce_lasso)
isoformComp$multi_enet = as.numeric(isoformComp$multi_enet)
isoformComp$joinet = as.numeric(isoformComp$joinet)
isoformComp$spls = as.numeric(isoformComp$spls)
isoformComp$univariate = as.numeric(isoformComp$univariate)
isoformComp$multivariate = as.numeric(isoformComp$multivariate)

isoformComp$PD = (isoformComp$multivariate -
                    isoformComp$univariate)/abs(isoformComp$univariate) * 100
isoformCompTissue = isoformComp %>%
  group_by(TissueOut) %>%
  summarise(Min = quantile(PD,.1),
            Lower = quantile(PD,.25),
            Median = quantile(PD,.5),
            Upper = quantile(PD,.75),
            Max = quantile(PD,.9),
            Min_r2 = quantile(multivariate - univariate,.1),
            Lower_r2 = quantile(multivariate - univariate,.25),
            Median_r2 = quantile(multivariate - univariate,.5),
            Upper_r2 = quantile(multivariate - univariate,.75),
            Max_r2 = quantile(multivariate - univariate,.9),
            Multivariate_Models = sum(multivariate > R2Cutoff),
            Univariate_Models = sum(univariate > R2Cutoff),
            NumIsoform = sum(multivariate > R2Cutoff |
                               univariate > R2Cutoff),
            Greater = 100*mean(multivariate > univariate),
            MRCE = sum(mrce_lasso > R2Cutoff),
            `Multivariate elastic net` = sum(multi_enet > R2Cutoff),
            `joinet` = sum(joinet > R2Cutoff),
            `SPLS` = sum(spls > R2Cutoff))
isoformCompTissue$Label_Median = paste0(round(isoformCompTissue$Greater,2),'%')
isoformCompTissue = isoformCompTissue[order(isoformCompTissue$Median,
                                            decreasing = T),]
isoformCompTissue$TissueOut = factor(as.character(isoformCompTissue$TissueOut),
                                     levels = 
                                       as.character(isoformCompTissue$TissueOut))

boxplot_overall_pd = ggplot(data = isoformCompTissue,
                            aes(x = TissueOut)) +
  geom_hline(yintercept = 0,color = 'red',linetype=2) +
  geom_boxplot(aes(middle = Median,
                   lower = Lower,
                   upper = Upper,
                   ymin = Min,
                   ymax = Max),
               stat = 'identity',
               width = .5,
               fill = 'grey') +
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
  scale_fill_brewer(palette = 'Set1') +
  xlab('Tissue') +
  ylab(expression("% difference in "~R^2~"for isoforms")) +
  guides(fill = 'none') +
  geom_label(aes(y = 2000,label = Label_Median),size = 2,fill='white') +
  coord_flip()

boxplot_brain_pd = ggplot(data = subset(isoformCompTissue,
                                        TissueOut %in% brain_tissues),
                          aes(x = TissueOut)) +
  geom_hline(yintercept = 0,color = 'red',linetype=2) +
  geom_boxplot(aes(middle = Median,
                   lower = Lower,
                   upper = Upper,
                   ymin = Min,
                   ymax = Max),
               stat = 'identity',
               width = .5,
               fill = 'grey') +
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
  scale_fill_brewer(palette = 'Set1') +
  xlab('Tissue') +
  ylab(expression("% difference in "~R^2~"for isoforms")) +
  guides(fill = 'none') +
  geom_label(aes(y = 2000,label = Label_Median),size = 3,fill='white') +
  coord_flip()
ggsave(plot = boxplot_overall_pd,
       filename = file.path(outFolder,
                            'Boxplot_IsoformPD_AllTissues.png'),
       height = 10,
       width = 6)
ggsave(plot = boxplot_brain_pd,
       filename = file.path(outFolder,
                            'Boxplot_IsoformPD_BrainTissues.png'),
       height = 4,
       width = 6)

## Number of isoforms over .01
Above01_isoform = reshape2::melt(
  isoformCompTissue,
  measure.vars = c('Multivariate_Models',
                   'Univariate_Models')
)
colnames(Above01_isoform)[(ncol(Above01_isoform) - 1):ncol(Above01_isoform)] = 
  c('Method','Number')
Above01_isoform$Method = as.factor(Above01_isoform$Method)
levels(Above01_isoform$Method) = c('Multivariate','Univariate')
Above01_isoform$TissueOut = factor(as.character(Above01_isoform$TissueOut),
                                   levels = all_tissues)

Above01_isoform_ratio = Above01_isoform %>%
  group_by(TissueOut) %>%
  summarize(Ratio = Number[Method == 'Multivariate']/Number[Method == 'Univariate'])

AT_plot = ggplot(data = Above01_isoform,
                 aes(x = TissueOut,
                     y = Number,
                     fill = Method)) +
  geom_col(position = position_dodge2(width = .9)) +
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
  scale_fill_manual(values = wes_palette("Darjeeling2")) +
  xlab('Tissue') +
  ylab(expression("Number of isoforms predicted at"~R^2 >= 0.01)) +
  coord_flip()

BT_plot =  ggplot(data = subset(Above01_isoform,
                                TissueOut %in% brain_tissues),
                  aes(x = TissueOut,
                      y = Number,
                      fill = Method)) +
  geom_col(position = position_dodge2(width = .9)) +
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
  scale_fill_manual(values = wes_palette("Darjeeling2")) +
  xlab('Tissue') +
  ylab(expression("Number of isoforms predicted at"~R^2 >= 0.01)) +
  coord_flip()

num01_isoform_bt = BT_plot


ggsave(plot = AT_plot,
       file = file.path(outFolder,
                        'Barplot_IsoformAbove01_AllTissues.png'),
       height = 10,
       width = 6)
ggsave(plot = BT_plot,
       file = file.path(outFolder,
                        'Barplot_IsoformAbove01_BrainTissues.png'),
       height = 4,
       width = 6)

### Best method barplot
isoformComp$BestMethod = pbapply::pbapply(isoformComp,1,
                                          function(x){
                                            names(which.max(x[c(3:7)]))
                                          })
bestMethodOverall = data.frame(Method = unique(isoformComp$BestMethod))
bestMethodOverall$Number = NA
for (i in 1:nrow(bestMethodOverall)){
  
  print(i)
  m = bestMethodOverall$Method[i]
  bestMethodOverall$Number[i] = 
    sum(isoformComp$BestMethod[isoformComp$univariate > isoformComp$R2Cutoff |
                                 isoformComp$multivariate > isoformComp$R2Cutoff] ==
          m)
  
}
bestMethodOverall$Method = c('Multivariate elastic net',
                             'SPLS','MRCE',
                             'Best univariate','Joinet')
bestMethodOverall$Method = factor(bestMethodOverall$Method,
                                  levels = c('Best univariate',
                                             'SPLS','MRCE','Joinet',
                                             'Multivariate elastic net'))

AT_plot = ggplot(data = bestMethodOverall,
                 aes(x = Method,
                     y = Number,
                     fill = Method)) +
  geom_col(position = position_dodge2(width = .9)) +
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
  scale_fill_brewer(palette = 'Dark2') +
  xlab('Tissue') +
  ylab(expression("Number of isoforms predicted at"~R^2 >= 0.01)) +
  guides(fill = 'none')

ggsave(plot = AT_plot,
       file = file.path(outFolder,
                        'Barplot_BestMethod.png'),
       height = 5,
       width = 7)

writexl::write_xlsx(as.data.frame(isoformCompTissue),
                    path = file.path(outFolder,
                                     'IsoformComparison_ByTissue.xlsx'),
                    col_names = T)


## TWAS vs isoTWAS inclusion criteria
setwd('outModels/GTEx')
tissues = list.files()
tissues = tissues[!grepl('[.]',
                         tissues)]
tissues = tissues[tissues != 'DoneRaw']
inclusionCriteria = data.frame(Tissue = tissues,
                               TWAS = NA,
                               isoTWAS = NA)
for (i in 1:nrow(inclusionCriteria)){
  t = inclusionCriteria$Tissue[i]
  print(t)
  setwd(t)
  
  fff = list.files('TWAS/')
  inclusionCriteria$TWAS[i] = length(fff)
  fff = list.files('isoTWAS/')
  inclusionCriteria$isoTWAS[i] = length(fff)
  
  
  setwd('outModels/GTEx')
}

inclusionCriteriaMelt = reshape2::melt(inclusionCriteria,
                                       measure.vars = c('TWAS','isoTWAS'))
colnames(inclusionCriteriaMelt)[2:3] = c('Method','Number')
inclusionCriteriaMelt$TissueOut = inclusionCriteriaMelt$Tissue
inclusionCriteriaMelt$TissueOut = as.factor(inclusionCriteriaMelt$TissueOut)
levels(inclusionCriteriaMelt$TissueOut) = all_tissues
inclusionCriteriaMelt$Method = relevel(inclusionCriteriaMelt$Method,
                                       ref = 'TWAS')

AT_plot = ggplot(data = inclusionCriteriaMelt,
                 aes(x = TissueOut,
                     y = Number,
                     fill = Method)) +
  geom_col(position = position_dodge2(width = .9)) +
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
  scale_fill_brewer(palette = 'Set1') +
  xlab('Tissue') +
  ylab(expression("Number of genes passing CV thresholds")) +
  coord_flip()

BT_plot =  ggplot(data = subset(inclusionCriteriaMelt,
                                TissueOut %in% brain_tissues),
                  aes(x = TissueOut,
                      y = Number,
                      fill = Method)) +
  geom_col(position = position_dodge2(width = .9)) +
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
  scale_fill_brewer(palette = 'Set1') +
  xlab('Tissue') +
  ylab(expression("Number of genes passing CV thresholds")) +
  coord_flip()

numInclusion_isoTWAS_TWAS_bt = BT_plot

ggsave(plot = AT_plot,
       file = file.path(outFolder,
                        'Barplot_GenesPassingCV_AllTissues.png'),
       height = 10,
       width = 6)
ggsave(plot = BT_plot,
       file = file.path(outFolder,
                        'Barplot_GenesPassingCV_BrainTissues.png'),
       height = 4,
       width = 6)

## Gene expression comparison
geneComp$R2_isoTWAS = as.numeric(geneComp$R2_isoTWAS)
geneComp$R2_TWAS = as.numeric(geneComp$R2_TWAS)

geneComp$PD = (geneComp$R2_isoTWAS -
                 geneComp$R2_TWAS)/abs(geneComp$R2_TWAS) * 100
geneCompTissue = geneComp %>%
  group_by(TissueOut) %>%
  summarise(Min = quantile(PD,.1),
            Lower = quantile(PD,.25),
            Median = quantile(PD,.5),
            Upper = quantile(PD,.75),
            Max = quantile(PD,.9),
            Min_r2 = quantile(R2_isoTWAS - R2_TWAS,.1),
            Lower_r2 = quantile(R2_isoTWAS - R2_TWAS,.25),
            Median_r2 = quantile(R2_isoTWAS - R2_TWAS,.5),
            Upper_r2 = quantile(R2_isoTWAS - R2_TWAS,.75),
            Max_r2 = quantile(R2_isoTWAS - R2_TWAS,.9),
            TWAS = sum(R2_TWAS > R2Cutoff),
            isoTWAS = sum(R2_isoTWAS > R2Cutoff),
            Greater = 100*mean(R2_isoTWAS > R2_TWAS))
geneCompTissue$Label_Median = paste0(round(geneCompTissue$Greater,2),'%')
geneCompTissue = geneCompTissue[order(geneCompTissue$Median,
                                      decreasing = T),]
geneCompTissue$TissueOut = factor(as.character(geneCompTissue$TissueOut),
                                  levels = 
                                    as.character(geneCompTissue$TissueOut))

boxplot_overall_pd = ggplot(data = geneCompTissue,
                            aes(x = TissueOut)) +
  geom_hline(yintercept = 0,color = 'red',linetype=2) +
  geom_boxplot(aes(middle = Median,
                   lower = Lower,
                   upper = Upper,
                   ymin = Min,
                   ymax = Max),
               stat = 'identity',
               width = .5,
               fill = 'grey') +
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
  scale_fill_brewer(palette = 'Set1') +
  xlab('Tissue') +
  ylab(expression("% difference in "~R^2~"for genes")) +
  guides(fill = 'none') +
  geom_label(aes(y = 2000,label = Label_Median),size = 2,fill='white') +
  coord_flip()

ggsave(plot = boxplot_overall_pd,
       filename = file.path(outFolder,
                            'Boxplot_GenePD_AllTissues.png'),
       height = 10,
       width = 7)

## Number of isoforms over .01
Above01_gene = reshape2::melt(
  geneCompTissue,
  measure.vars = c('TWAS',
                   'isoTWAS')
)
colnames(Above01_gene)[(ncol(Above01_gene) - 1):ncol(Above01_gene)] = 
  c('Method','Number')
Above01_gene$Method = as.factor(Above01_gene$Method)
Above01_gene$Method = relevel(Above01_gene$Method,
                              ref = 'TWAS')
Above01_gene$TissueOut = factor(as.character(Above01_gene$TissueOut),
                                levels = all_tissues)



Above01_gene_ratio = Above01_gene %>%
  group_by(TissueOut) %>%
  summarize(Ratio = Number[Method == 'isoTWAS']/Number[Method == 'TWAS'])

df_summary_isotwas = data.frame(Tissue = rep(unique(Above01_gene$TissueOut),2),
                                Feature = rep(c('Gene','Isoform'),
                                              each =
                                                length(unique(Above01_gene$TissueOut))))
df_summary_isotwas$Number = 0
for (i in 1:nrow(df_summary_isotwas)){
  
  if (df_summary_isotwas$Feature[i] == 'Gene'){
    df_summary_isotwas$Number[i] = 
      inclusionCriteriaMelt$Number[inclusionCriteriaMelt$Method == 'isoTWAS' &
                                     inclusionCriteriaMelt$TissueOut == df_summary_isotwas$Tissue[i]]
  } else {
    df_summary_isotwas$Number[i] = 
      isoformCompTissue$NumIsoform[isoformCompTissue$TissueOut == 
                                     df_summary_isotwas$Tissue[i]]
  }
  
}


sumIsotwas = ggplot(data = 
                      df_summary_isotwas,
                    aes(x = Tissue,
                        y = Number)) +
  geom_col(position = position_dodge2(width = .9)) +
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
  scale_fill_brewer(palette = 'Set1') +
  xlab('Tissue') +
  ylab(expression("Number of features predicted at"~R^2 >= 0.01)) +
  coord_flip() +
  facet_wrap(~Feature,scales='free_x')
ggsave(plot = sumIsotwas,
       filename = '/GTExPlotsTables/NumFeatures_isoTWAS.png',
       height = 8,
       width = 6)

AT_plot = ggplot(data = Above01_gene,
                 aes(x = TissueOut,
                     y = Number,
                     fill = Method)) +
  geom_col(position = position_dodge2(width = .9)) +
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
  scale_fill_brewer(palette = 'Set1') +
  xlab('Tissue') +
  ylab(expression("Number of genes predicted at"~R^2 >= 0.01)) +
  coord_flip()

BT_plot =  ggplot(data = subset(Above01_gene,
                                TissueOut %in% brain_tissues),
                  aes(x = TissueOut,
                      y = Number,
                      fill = Method)) +
  geom_col(position = position_dodge2(width = .9)) +
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
  scale_fill_brewer(palette = 'Set1') +
  xlab('Tissue') +
  ylab(expression("Number of genes predicted at"~R^2 >= 0.01)) +
  coord_flip()

num01_gene_bt = BT_plot

ggsave(plot = AT_plot,
       file = file.path(outFolder,
                        'Barplot_GeneAbove01_AllTissues.png'),
       height = 10,
       width = 6)
ggsave(plot = BT_plot,
       file = file.path(outFolder,
                        'Barplot_GeneAbove01_BrainTissues.png'),
       height = 4,
       width = 6)


writexl::write_xlsx(as.data.frame(geneCompTissue),
                    path = file.path(outFolder,
                                     'geneComparison_ByTissue.xlsx'),
                    col_names = T)


### External
external = fread('GTEx_Models_Raw_042523/ExternalR2_Table.tsv')
external$Label = paste0(round(external$Above),
                        '%')

external$LabelAlt = paste0('Median:\n',round(external$Median,1),
                           '%')
external$Dataset = c('PEC+AMP-AD/GTEx',
                     'GTEx/PEC+AMP-AD')

percdiff_gene = ggplot(data = external,
                       aes(x = Dataset)) +
  geom_hline(color = 'red',
             linetype = 2,
             yintercept = 0) +
  geom_boxplot(aes(middle = Median,
                   lower = Lower,
                   upper = Upper,
                   ymin = Min,
                   ymax = Max),
               stat = 'identity',
               width = .2,
               fill = 'grey') +
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
  scale_fill_brewer(palette = 'Set1') +
  xlab('Training/imputation set') +
  ylab(expression("% difference in "~R^2~"for genes")) +
  geom_label(aes(y = 500,label = LabelAlt),size = 4,fill='white')


### Summary stats
setwd('GTEx_Models_Raw_042523')
gc = geneComp[,c('TissueOut','SampleSize')]
gc = gc[!duplicated(gc),]
Above01_isoform = merge(Above01_isoform,gc,
                        by = 'TissueOut')
Above01_isoform$Brain = ifelse(Above01_isoform$TissueOut %in% brain_tissues,
                               'Brain','Other')
summary(lm(Above01_isoform$Median ~ Above01_isoform$Brain + Above01_isoform$SampleSize))
brainOther = ggplot(data = subset(Above01_isoform,Method == 'Multivariate'),
                    aes(x = as.factor(Brain),
                        y = Median,
                        fill = as.factor(Brain))) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  guides(fill = 'none') +
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
  xlab('Tissue') +
  ylab(expression("Median % difference in "~R^2~"(multivariate - univariate)")) +
  scale_fill_manual(values = c('darkorchid2','grey')) +
  annotate("label", x = 1, y = 180, label = expression(P==0.011),
           size = 6)


require(cowplot)
isoform_leg = get_legend(num01_isoform_bt)
num01_isoform_bt = 
  num01_isoform_bt + ylab(expression('# isoforms with'~R^2>=0.01))
brainOther
gene_leg = get_legend(numInclusion_isoTWAS_TWAS_bt)
numInclusion_isoTWAS_TWAS_bt = numInclusion_isoTWAS_TWAS_bt +
  ylab('# genes passing CV thresholds')
num01_gene_bt =
  num01_gene_bt + ylab(expression('# genes with'~R^2>=0.01))
percdiff_gene

top = plot_grid(plot_grid(num01_isoform_bt + guides(fill = 'none'),
                          isoform_leg,
                          ncol=1,
                          rel_heights = c(1,.05)),
                plot_grid(num01_gene_bt + guides(fill = 'none'),
                          gene_leg,
                          ncol=1,
                          rel_heights = c(1,.05)),rel_widths = c(1,1),
                labels = c('a','d'))
bottom = plot_grid(brainOther,
                   plot_grid(numInclusion_isoTWAS_TWAS_bt + guides(fill = 'none'),
                             gene_leg,
                             ncol=1,
                             rel_heights = c(1,.05)),
                   percdiff_gene,
                   rel_widths = c(.7,1.2,.95),
                   nrow=1,
                   labels = c('b','c','e'))
tot_all = plot_grid(top,bottom,
                    ncol = 1,
                    rel_heights = c(1,1))
ggsave(plot = tot_all,
       filename = '/GTExPlotsTables/AllGTExPredPlot.png',
       height = 10,
       width = 12)


require(data.table)
require(ggplot2)
df_out = fread('GTEx_Models_Raw_042523/SummaryStats.tsv')
