require(data.table)
res = fread('sims_res_isoTWAS_prediction.tsv',fill=T)
res = res[complete.cases(res),]
res$n.qtl = as.numeric(res$n.qtl)
res = subset(res,n.qtl %in% c(200,500,1000))

require(tidyverse)
params = res %>%
  group_by(n.qtl, p.causal, 
           prop.shared, eqtl.h2, 
           eqtl_nongenetic_corr, n.tx) %>%
  summarise(Min = quantile(PercDiff,.1),
            Low = quantile(PercDiff,.25),
            Median = quantile(PercDiff,.5),
            High = quantile(PercDiff,.75),
            Max = quantile(PercDiff,.9),
            Min_R2 = quantile(isotwas_r2 - twas_r2,.1),
            Low_R2 = quantile(isotwas_r2 - twas_r2,.25),
            Median_R2 = quantile(isotwas_r2 - twas_r2,.5),
            High_R2 = quantile(isotwas_r2 - twas_r2,.75),
            Max_R2 = quantile(isotwas_r2 - twas_r2,.9))
params_r2 = res %>%
  group_by(n.qtl, p.causal, 
           prop.shared, eqtl.h2, 
           eqtl_nongenetic_corr, n.tx) %>%
  summarise(Min = quantile(isotwas_r2 - twas_r2,.1),
            Low = quantile(isotwas_r2 - twas_r2,.25),
            Median = quantile(isotwas_r2 - twas_r2,.5),
            High = quantile(isotwas_r2 - twas_r2,.75),
            Max = quantile(isotwas_r2 - twas_r2,.9))

head(params)
params$prop.shared = paste0('p[shared] == ',
                            params$prop.shared)
params$eqtl_nongenetic_corr = paste0('sigma[h] == ',
                                     params$eqtl_nongenetic_corr)
params$p.causal = as.factor(params$p.causal)
params$n.tx = as.factor(params$n.tx)


params_r2$prop.shared = paste0('p[shared] == ',
                               params_r2$prop.shared)
params_r2$eqtl_nongenetic_corr = paste0('sigma[h] == ',
                                        params_r2$eqtl_nongenetic_corr)
params_r2$p.causal = as.factor(params_r2$p.causal)
params_r2$n.tx = as.factor(params_r2$n.tx)

require(ggplot2)
require(ggallin)

for (ss in as.numeric(unique(params$n.qtl))){
  for (h2 in as.numeric(unique(params$eqtl.h2))){
    
    print(ss)
    params_forplots = subset(params,n.qtl == ss &
                               eqtl.h2 == h2)
    params_r2_forplots = subset(params_r2,n.qtl == ss &
                                  eqtl.h2 == h2)
    perc = ggplot(data = params_forplots,
           aes(x = p.causal,
               fill = n.tx)) +
      geom_hline(yintercept = 0,linetype = 2,color = 'red') +
      facet_grid(eqtl_nongenetic_corr ~ prop.shared, labeller = label_parsed) +
      geom_boxplot(aes(middle = Median,
                       lower = Low,
                       upper = High,
                       ymin = Min,
                       ymax = Max),
                   stat = 'identity',
                   width = .5) +
      theme_minimal() +
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=11),
            plot.title = element_text(size = 12),
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
      scale_fill_manual(values = wesanderson::wes_palette('Moonrise2')) +
      ylab(expression("% difference in adjusted"~R^2~"(log-scale)")) +
      xlab(expression(p[c])) +
      labs(fill = 'Number of transcripts') +
      scale_y_continuous(trans = pseudolog10_trans,
                         breaks = c(-100,-10,0,10,100,1000)) +
      ggtitle(bquote(N==.(ss)~"and"~h[g]^2==.(h2)))
    
    r2plot = ggplot(data = params_r2_forplots,
           aes(x = as.factor(p.causal),
               fill = as.factor(n.tx))) +
      geom_hline(yintercept = 0,linetype = 2,color = 'red') +
      facet_grid(as.factor(eqtl_nongenetic_corr) ~ as.factor(prop.shared),
                 labeller = label_parsed) +
      geom_boxplot(aes(middle = Median,
                       lower = Low,
                       upper = High,
                       ymin = Min,
                       ymax = Max),
                   stat = 'identity',
                   width = .5) +
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
      scale_fill_manual(values = wesanderson::wes_palette('Moonrise2')) +
      ylab(expression("Difference in adjusted"~R^2~"(log-scale)")) +
      xlab(expression(p[c])) +
      labs(fill = 'Number of transcripts') +
      ggtitle(bquote(N==.(ss)~"and"~h[g]^2==.(h2)))
    
    dir.create('PredictionPlots')
    ggsave(filename = file.path('PredictionPlots',
                                paste0('SampleSize',ss,'_h2',h2,'_Percent.png')),
           plot = perc,
           height = 4,
           width = 6)
    ggsave(filename = file.path('PredictionPlots',
                                paste0('SampleSize',ss,'_h2',h2,'_r2.png')),
           plot = r2plot,
           height = 5,
           width = 6)
    
    
  }
}





require(ggplot2)
require(data.table)
require(vroom)
method_diff = vroom('sims_res_isoTWAS_prediction_050123.tsv')
method_diff = as.data.frame(method_diff)
method_diff = method_diff[complete.cases(method_diff),]

method_diff_melt = reshape2::melt(method_diff,
                                  measure.vars = colnames(method_diff)[7:11])
colnames(method_diff_melt)[12:13] = c('Method','R2')
method_diff_melt$Method = as.factor(method_diff_melt$Method)
levels(method_diff_melt$Method) =  c('MRCE','Multivariate elastic net',
                                     'Best univariate','SPLS','joinet')
method_diff_melt$Method = relevel(method_diff_melt$Method,
                                  ref = 'Best univariate')



require(tidyverse)
params_r2 = method_diff_melt %>%
  group_by(p.causal, 
           prop.shared, 
           eqtl_nongenetic_corr,Method) %>%
  summarise(Min = quantile(R2,.1),
            Low = quantile(R2,.25),
            Median = quantile(R2,.5),
            High = quantile(R2,.75),
            Max = quantile(R2,.9))

params_r2$eqtl_nongenetic_corr = paste0('sigma[h] == ',
                                        params_r2$eqtl_nongenetic_corr)
params_r2$p.causal = paste0('p[c] == ',
                            params_r2$p.causal)
params_r2$prop.shared = as.factor(params_r2$prop.shared)

overallboxr2 = ggplot(data = params_r2,
                aes(x = prop.shared,
                    fill = Method)) +
  geom_hline(yintercept = 0.01,linetype = 2,color = 'red') +
  facet_grid(as.factor(eqtl_nongenetic_corr) ~ as.factor(p.causal),
             labeller = label_parsed) +
  geom_boxplot(aes(middle = Median,
                   lower = Low,
                   upper = High,
                   ymin = Min,
                   ymax = Max),
               stat = 'identity',
               width = .5) +
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
  ylab(expression('CV'~R^2)) +
  xlab(expression(p[shared])) +
  scale_color_brewer(palette = 'Dark2')






require(ggplot2)
require(data.table)
require(vroom)
method_diff = vroom('sims_res_isoTWAS_prediction_050123.tsv')
method_diff = as.data.frame(method_diff)
method_diff = method_diff[complete.cases(method_diff),]


method_diff$BestMethod = apply(method_diff,1,function(x){
  return(names(which.max(x[c(7,8,9,10,11)])))
})

params_r2 = method_diff %>%
  group_by(p.causal, 
           prop.shared, 
           eqtl_nongenetic_corr) %>%
  summarise(mrce_lasso = mean(BestMethod == 'mrce_lasso'),
            joinet = mean(BestMethod == 'joinet'),
            univariate = mean(BestMethod == 'univariate'),
            spls = mean(BestMethod == 'spls'),
            multi_enet = mean(BestMethod == 'multi_enet'))
method_diff_melt = reshape2::melt(
  params_r2,
  measure.vars = unique(method_diff$BestMethod)
)
colnames(method_diff_melt)[4:5] = c('Method','Proportion')
method_diff_melt$Method = as.factor(method_diff_melt$Method)
levels(method_diff_melt$Method) =  c('Best univariate','Multivariate elastic net',
                                     'SPLS','MRCE','joinet')
method_diff_melt$Method = relevel(method_diff_melt$Method,
                                  ref = 'Best univariate')


method_diff_melt$eqtl_nongenetic_corr = paste0('sigma[h] == ',
                                        method_diff_melt$eqtl_nongenetic_corr)
method_diff_melt$p.causal = paste0('p[c] == ',
                                   method_diff_melt$p.causal)
trackplot = ggplot(data = method_diff_melt,
                      aes(x = as.factor(prop.shared),
                          y = Proportion,
                          fill = Method,
                          color = Method,
                          group = Method)) +
  geom_hline(yintercept = 0,linetype = 2,color = 'red') +
  facet_grid(as.factor(eqtl_nongenetic_corr) ~ as.factor(p.causal),
             labeller = label_parsed) +
  geom_point(size = 3) +
  geom_line() +
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
  ylab(expression('% of simulations with best model')) +
  xlab(expression(p[shared])) +
  scale_color_brewer(palette = 'Dark2')


require(cowplot)
legend = get_legend(trackplot)
tot_plot = plot_grid(plot_grid(overallboxr2 + guides(fill = 'none'),
                               trackplot + guides(color = 'none',
                                                  fill = 'none',
                                                  shape = 'none'),
                               ncol=2,
                               rel_widths = c(1,1),
                               labels = c('a','b')),
                     legend,ncol = 1,rel_heights = c(1,.05))
ggsave(plot = tot_plot,
       filename = 'BestMethod_Sims.png',
       height = 5,
       width = 8)
