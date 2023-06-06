library("optparse")
option_list <- list(
  make_option(c("-i", "--index"), action="store_true", default=TRUE,
              help="index",type='numeric')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

##ceiling = 648
##number = 1000

index = opt$index
source('/u/scratch/a/abtbhatt/sims_isotwas_power/GeneNotIso_Function.R')
source('/u/scratch/a/abtbhatt/sims_isotwas_power/TwoOffset_Function.R')
source('/u/scratch/a/abtbhatt/sims_isotwas_power/OnlyOne_Function.R')


#### create sims
require(bigsnpr)
require(MOSTWAS)
require(data.table)
require(isotwas)
n.qtl = c(500)
p.causal = c(0.01,.05)
prop.shared = c(0,.5,1)
eqtl_h2.loc = c(0.05,.1)
eqtl_nongenetic_corr = c(.1,.25)
gwas_var = c(.1,.25,.05)
n.tx = c(2,5,10)[3:1]
scenarios = c('gene_not_iso',
              'only_one',
              'two_offset')
params = expand.grid(n.qtl,
                     p.causal,
                     prop.shared,
                     eqtl_h2.loc,
                     eqtl_nongenetic_corr,
                     n.tx,
                     gwas_var,
                     scenarios)
colnames(params) = c('n.qtl',
                     'p.causal',
                     'prop.shared',
                     'eqtl.h2',
                     'heterogeneity',
                     'n.tx',
                     'gwas_var',
                     'scenario')

snp_folder = '/u/scratch/a/abtbhatt/sims_isotwas_prediction/snp_files'
snp_files = file.path(snp_folder,
                      list.files(snp_folder))

for (chr in 1:22){
  
  dir.create('/u/scratch/a/abtbhatt/sims_isotwas_power')
  if (params$scenario[index] == 'gene_not_iso'){
    done = data.table::fread('/u/scratch/a/abtbhatt/sims_isotwas_power/sims_power_done_genenotiso.tsv')
    done = subset(done,
                  Chromosome == chr &
                    n.qtl == params$n.qtl[index] &
                    p.causal == params$p.causal[index] &
                    prop_shared == params$prop.shared[index] &
                    eqtl_h2 == params$eqtl.h2[index] &
                    eqtl_nongenetic_corr == params$heterogeneity[index] &
                    n.tx == params$n.tx[index] &
                    gwas_var == params$gwas_var[index])
    
    
    if (nrow(done) == 0){
      
      bed_files = snp_files[grepl(paste0('_chr',chr,'_'),
                                  snp_files)]
      bed_files = bed_files[grepl('.bed',
                                  bed_files)]
      snps = snp_attach(snp_readBed2(bed_files[1],
                                     backingfile = tempfile()))
      ttt = pbapply::pbreplicate(20,
                                 power.sim_gene_not_iso(snps,
                                                        n.qtl = params$n.qtl[index],
                                                        p.causal = params$p.causal[index],
                                                        prop_shared = params$prop.shared[index],
                                                        eqtl_h2 = params$eqtl.h2[index],
                                                        eqtl_nongenetic_corr = params$heterogeneity[index],
                                                        n.tx = params$n.tx[index],
                                                        gwas_var = params$gwas_var[index],
                                                        scenario = params$scenario[index],
                                                        chr = chr,
                                                        gene = strsplit(strsplit(bed_files,
                                                                                 '/')[[1]][8],
                                                                        '_')[[1]][1],
                                                        outFile = file.path('/u/scratch/a/abtbhatt/sims_isotwas_power',
                                                                            'sims_res_isoTWAS_genenotiso.tsv')
                                                        )
                                   )
        
        data.table::fwrite(data.frame(Chromosome = chr,
                                      n.qtl = params$n.qtl[index],
                                      p.causal = params$p.causal[index],
                                      prop_shared = params$prop.shared[index],
                                      eqtl_h2 = params$eqtl.h2[index],
                                      eqtl_nongenetic_corr = params$heterogeneity[index],
                                      n.tx = params$n.tx[index],
                                      gwas_var = params$gwas_var[index]),
                           '/u/scratch/a/abtbhatt/sims_isotwas_power/sims_power_done_genenotiso.tsv',
                           append=T,
                           row.names=F,
                           quote=F,
                           sep='\t')
        }
  }
  
  
  if (params$scenario[index] == 'only_one'){
    done = data.table::fread('/u/scratch/a/abtbhatt/sims_isotwas_power/sims_power_done_onlyone.tsv')
    done = subset(done,
                  Chromosome == chr &
                    n.qtl == params$n.qtl[index] &
                    p.causal == params$p.causal[index] &
                    prop_shared == params$prop.shared[index] &
                    eqtl_h2 == params$eqtl.h2[index] &
                    eqtl_nongenetic_corr == params$heterogeneity[index] &
                    n.tx == params$n.tx[index] &
                    gwas_var == params$gwas_var[index])
    
    
    if (nrow(done) == 0){
      
      bed_files = snp_files[grepl(paste0('_chr',chr,'_'),
                                  snp_files)]
      bed_files = bed_files[grepl('.bed',
                                  bed_files)]
      snps = snp_attach(snp_readBed2(bed_files[1],
                                     backingfile = tempfile()))
      ttt = pbapply::pbreplicate(20,
                                 power.sim_only_one(snps,
                                                        n.qtl = params$n.qtl[index],
                                                        p.causal = params$p.causal[index],
                                                        prop_shared = params$prop.shared[index],
                                                        eqtl_h2 = params$eqtl.h2[index],
                                                        eqtl_nongenetic_corr = params$heterogeneity[index],
                                                        n.tx = params$n.tx[index],
                                                        gwas_var = params$gwas_var[index],
                                                        scenario = params$scenario[index],
                                                        chr = chr,
                                                        gene = strsplit(strsplit(bed_files,
                                                                                 '/')[[1]][8],
                                                                        '_')[[1]][1],
                                                        outFile = file.path('/u/scratch/a/abtbhatt/sims_isotwas_power',
                                                                            'sims_res_isoTWAS_onlyone.tsv')
                                 )
      )
      
      data.table::fwrite(data.frame(Chromosome = chr,
                                    n.qtl = params$n.qtl[index],
                                    p.causal = params$p.causal[index],
                                    prop_shared = params$prop.shared[index],
                                    eqtl_h2 = params$eqtl.h2[index],
                                    eqtl_nongenetic_corr = params$heterogeneity[index],
                                    n.tx = params$n.tx[index],
                                    gwas_var = params$gwas_var[index]),
                         '/u/scratch/a/abtbhatt/sims_isotwas_power/sims_power_done_onlyone.tsv',
                         append=T,
                         row.names=F,
                         quote=F,
                         sep='\t')
    }
  }
  
  
  if (params$scenario[index] == 'two_offset'){
    done = data.table::fread('/u/scratch/a/abtbhatt/sims_isotwas_power/sims_power_done_twooffset.tsv')
    done = subset(done,
                  Chromosome == chr &
                    n.qtl == params$n.qtl[index] &
                    p.causal == params$p.causal[index] &
                    prop_shared == params$prop.shared[index] &
                    eqtl_h2 == params$eqtl.h2[index] &
                    eqtl_nongenetic_corr == params$heterogeneity[index] &
                    n.tx == params$n.tx[index] &
                    gwas_var == params$gwas_var[index])
    
    
    if (nrow(done) == 0){
      
      bed_files = snp_files[grepl(paste0('_chr',chr,'_'),
                                  snp_files)]
      bed_files = bed_files[grepl('.bed',
                                  bed_files)]
      snps = snp_attach(snp_readBed2(bed_files[1],
                                     backingfile = tempfile()))
      ttt = pbapply::pbreplicate(20,
                                 power.sim_two_offset(snps,
                                                        n.qtl = params$n.qtl[index],
                                                        p.causal = params$p.causal[index],
                                                        prop_shared = params$prop.shared[index],
                                                        eqtl_h2 = params$eqtl.h2[index],
                                                        eqtl_nongenetic_corr = params$heterogeneity[index],
                                                        n.tx = params$n.tx[index],
                                                        gwas_var = params$gwas_var[index],
                                                        scenario = params$scenario[index],
                                                        chr = chr,
                                                        gene = strsplit(strsplit(bed_files,
                                                                                 '/')[[1]][8],
                                                                        '_')[[1]][1],
                                                        outFile = file.path('/u/scratch/a/abtbhatt/sims_isotwas_power',
                                                                            'sims_res_isoTWAS_twooffset.tsv')
                                 )
      )
      
      data.table::fwrite(data.frame(Chromosome = chr,
                                    n.qtl = params$n.qtl[index],
                                    p.causal = params$p.causal[index],
                                    prop_shared = params$prop.shared[index],
                                    eqtl_h2 = params$eqtl.h2[index],
                                    eqtl_nongenetic_corr = params$heterogeneity[index],
                                    n.tx = params$n.tx[index],
                                    gwas_var = params$gwas_var[index]),
                         '/u/scratch/a/abtbhatt/sims_isotwas_power/sims_power_done_twooffset.tsv',
                         append=T,
                         row.names=F,
                         quote=F,
                         sep='\t')
    }
  }
  
}
