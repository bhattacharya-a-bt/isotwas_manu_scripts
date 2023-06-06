require(data.table)
params = fread('ST - GeneNotIsoPowerSim.tsv')
params$prop_shared = paste0('p[shared] == ',
                            params$prop_shared)
params$eqtl_nongenetic_corr = paste0('sigma[h] == ',
                                     params$eqtl_nongenetic_corr)
params$p.causal = as.factor(params$p.causal)
params$n.tx = as.factor(params$n.tx)


require(ggplot2)
require(ggallin)

for (gv in unique(params$gwas_var)){
  for (h2 in unique(params$eqtl_h2)){
    
    params_forplots = subset(params,gwas_var == gv &
                               eqtl_h2 == h2)
    
    params_melt = rbind(params_forplots[,1:7],
                        params_forplots[,1:7])
    params_melt$Power = c(params_forplots$Gene_Power,
                          params_forplots$Tx_Power)
    params_melt$Lower = c(params_forplots$Gene_Lower,
                          params_forplots$Tx_Lower)
    params_melt$Higher = c(params_forplots$Gene_Higher,
                           params_forplots$Tx_Higher)
    params_melt$Method = as.factor(rep(c('TWAS','isoTWAS'),
                                       each = nrow(params_forplots)))
    params_melt$Method = relevel(params_melt$Method,
                                 ref = 'TWAS')
    params_melt$n.tx = factor(as.character(params_melt$n.tx),
                              levels = c(2,5,10))
    
    perc = ggplot(data = params_melt,
                  aes(x = n.tx,
                      y = Power,
                      shape = p.causal,
                      color = Method,
                      group = interaction(p.causal,Method))) +
      geom_hline(yintercept = 0.8,linetype = 2,color = 'grey') +
      geom_point(size = 2,position = position_dodge2(width=.9)) +
      geom_line(position = position_dodge2(width = .9)) +
      facet_grid(eqtl_nongenetic_corr ~ prop_shared, labeller = label_parsed) +
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
      scale_color_brewer(palette = 'Set1') +
      ylab(expression("Power")) +
      xlab("Number of isoforms") +
      labs(color = 'Method',
           shape = expression(p[c])) +
      ggtitle(bquote(h^2==.(gv)~"and"~h[g]^2==.(h2))) +
      ylim(c(0,1))
    
    ggsave(filename = file.path('./',
                                paste0('GWASVar',gv,'_h2',h2,
                                       '_GeneNotIso.png')),
           plot = perc,
           height = 5,
           width = (6*(5/4)))
    
  }
}



require(data.table)
params = 
  fread('ST - OnlyOnePowerSim.tsv')


params$prop_shared = paste0('p[shared] == ',
                            params$prop_shared)
params$eqtl_nongenetic_corr = paste0('sigma[h] == ',
                                     params$eqtl_nongenetic_corr)
params$p.causal = as.factor(params$p.causal)
params$n.tx = as.factor(params$n.tx)


require(ggplot2)
require(ggallin)
for (pc in unique(params$p.causal)){
  for (gv in unique(params$gwas_var)){
    for (h2 in unique(params$eqtl_h2)){
      
      params_forplots = subset(params,gwas_var == gv &
                                 eqtl_h2 == h2 &
                                 p.causal == pc)
      
      params_melt = rbind(params_forplots[,1:8],
                          params_forplots[,1:8])
      params_melt$Power = c(params_forplots$Power_TWAS,
                            params_forplots$Power_isoTWAS)
      params_melt$Lower = c(params_forplots$Lower_TWAS,
                            params_forplots$Lower_isoTWAS)
      params_melt$Higher = c(params_forplots$Higher_TWAS,
                             params_forplots$Higher_isoTWAS)
      params_melt$Method = as.factor(rep(c('TWAS','isoTWAS'),
                                         each = nrow(params_forplots)))
      params_melt$Method = relevel(params_melt$Method,
                                   ref = 'TWAS')
      pc = as.numeric(as.character(pc))
      
      
      perc = ggplot(data = params_melt,
                    aes(x = usage,
                        y = Power,
                        shape = as.factor(n.tx),
                        color = Method,
                        group=interaction(n.tx, Method))) +
        geom_hline(yintercept = 0.8,linetype = 2,color = 'grey') +
        geom_point(size = 2,position = position_dodge2(width=.9)) +
        geom_line(position = position_dodge2(width=.9)) +
        facet_grid(eqtl_nongenetic_corr ~ prop_shared, labeller = label_parsed) +
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
        scale_color_brewer(palette = 'Set1') +
        scale_x_continuous(breaks = c(.1,.3,.5,.7,.9)) +
        ylab(expression("Power")) +
        xlab(expression("Usage of effect isoform")) +
        labs(color = 'Method',
             shape = 'Number of isoforms') +
        ggtitle(bquote(h^2==.(gv)~","~h[g]^2==.(h2)~"and"~p[c]==.(pc)))+
        ylim(c(0,1))
      
    ggsave(filename = file.path('./',
                                  paste0('GWASVar',gv,'_h2',h2,'_pc',
                                         pc,'_OnlyOne.png')),
             plot = perc,
             height = 5,
             width = (6*(5/4)))
      
    }
  }
}




require(data.table)
params = 
  fread('ST - TwoOffset.tsv')


params$prop_shared = paste0('p[shared] == ',
                            params$prop_shared)
params$eqtl_nongenetic_corr = paste0('sigma[h] == ',
                                     params$eqtl_nongenetic_corr)
params$p.causal = as.factor(params$p.causal)
params$n.tx = as.factor(params$n.tx)
params$ratio_effect = as.factor(params$ratio_effect)


require(ggplot2)
require(ggallin)

for (pc in unique(params$p.causal)){
  for (gv in unique(params$gwas_var)){
    for (h2 in unique(params$eqtl_h2)){
      
      params_forplots = subset(params,gwas_var == gv &
                                 eqtl_h2 == h2 &
                                 p.causal == pc)
      
      params_melt = rbind(params_forplots[,1:8],
                          params_forplots[,1:8])
      params_melt$Power = c(params_forplots$Power_TWAS,
                            params_forplots$Power_isoTWAS)
      params_melt$Lower = c(params_forplots$Lower_TWAS,
                            params_forplots$Lower_isoTWAS)
      params_melt$Higher = c(params_forplots$Higher_TWAS,
                             params_forplots$Higher_isoTWAS)
      params_melt$Method = as.factor(rep(c('TWAS','isoTWAS'),
                                         each = nrow(params_forplots)))
      params_melt$Method = relevel(params_melt$Method,
                                   ref = 'TWAS')
      pc = as.numeric(as.character(pc))
      
      
      perc = ggplot(data = params_melt,
                    aes(x = ratio_effect,
                        y = Power,
                        shape = as.factor(n.tx),
                        color = Method,
                        group=interaction(n.tx, Method))) +
        geom_hline(yintercept = 0.8,linetype = 2,color = 'grey') +
        geom_point(size = 2,position = position_dodge2(width=.9)) +
        geom_line(position = position_dodge2(width=.9)) +
        facet_grid(eqtl_nongenetic_corr ~ prop_shared, labeller = label_parsed) +
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
        scale_color_brewer(palette = 'Set1') +
        ylab(expression("Power")) +
        xlab(expression("Ratio of effect sizes")) +
        labs(color = 'Method',
             shape = 'Number of isoforms') +
        ggtitle(bquote(h^2==.(gv)~","~h[g]^2==.(h2)~"and"~p[c]==.(pc))) +
        ylim(c(0,1))
      
    ggsave(filename = file.path('./',
                                  paste0('GWASVar',gv,'_h2',h2,'_pc',
                                         pc,'_TwoOffset.png')),
             plot = perc,
             height = 5,
             width = 6*(5/4))
      
    }
  }
}

