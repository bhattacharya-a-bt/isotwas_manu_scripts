library("optparse")
option_list <- list(
  make_option(c("-i", "--index"), action="store_true", default=TRUE,
              help="index",type='numeric'),
  make_option(c("-n", "--number"), action="store_true", default=TRUE,
              help="number of genes",type='numeric')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

##ceiling = 43
##number = 1000

index = opt$index
number = opt$number

twas_folder="out_models_041923/AdultBB/TWAS"
isotwas_folder="out_models_041923/AdultBB/isoTWAS"
ld_folder = "out_models_041923/AdultBB/LD"
gwas_folder="GWAS"

trait = c('ADHD',
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



tr = trait[index]
fff = list.files(gwas_folder)
fff = fff[grepl(tr,fff)]
fff = fff[grepl('sumstats',fff)]
fff = fff[!grepl('.no_',fff)]
sumstats_file = file.path(gwas_folder,
                          fff[1])

library(isotwas)
library(data.table)
library(bigsnpr)


print(tr)
sumstats = fread(sumstats_file)
dir.create(file.path('AdultBB/Associations',
                     tr),recursive = T)

outFile_twas = file.path('AdultBB/Associations',
                         tr,
                         paste0(tr,'_TWAS.tsv'))
outFile_isotwas = file.path('AdultBB/Associations',
                            tr,
                            paste0(tr,'_isoTWAS.tsv'))



fff = list.files(isotwas_folder)
toDoGenes = sapply(strsplit(fff,'_iso'),
                   function(x) x[1])


if (file.exists(outFile_isotwas) &
    file.exists(outFile_twas)){
  ttt = fread(outFile_twas)
  iii = fread(outFile_isotwas)
  toDoGenes = toDoGenes[!toDoGenes %in% c(ttt$Gene,
                                          iii$Gene)]
}



doneFile = 
  paste0('AdultBB/',
         tr,'_done_dev.tsv')

file.remove(doneFile)

if (!exists('ttt') & !exists('iii')){
  
  
  fwrite(data.frame(Trait = tr,
                    Gene = 'test',
                    Index = 1),
         doneFile,
         sep='\t',
         col.names=T,
         row.names=F,
         quote=F)
  
}

if (!file.exists(doneFile)){
  
  fwrite(data.frame(Trait = tr,
                    Gene = unique(c(ttt$Gene,
                                    iii$Gene)),
                    Index = 1),
         doneFile,
         sep='\t',
         col.names=T,
         row.names=F,
         quote=F)
  
}


runAssociation = function(gene){
  
  print(gene)
  
  if (nrow(subset(vroom::vroom(doneFile,show_col_types=F),
                  Gene == gene & Trait == tr)) == 0){
    
    ### RUN ISOTWAS
    isotwas_model = readRDS(file.path(isotwas_folder,
                                      paste0(gene,'_isoTWAS.RDS')))
    isotwas_model$R2 = unlist(isotwas_model$R2)
    
    if (file.exists(file.path(twas_folder,
                              paste0(gene,'_TWAS.RDS')))){
      twas_model = readRDS(file.path(twas_folder,
                                     paste0(gene,'_TWAS.RDS')))
    } else {
      twas_model = as.data.frame(isotwas_model[1:2,])
      twas_model$R2 = 0
    }
    
    if (any(twas_model$R2 > 0.01) |
        any(unlist(isotwas_model$R2) > 0.01)){
      
      LD_isotwas = readRDS(file.path(ld_folder,
                                     paste0(gene,'_LDMatrix.RDS')))
      
      LD_twas = readRDS(file.path(ld_folder,
                                  paste0(gene,'_LDMatrix.RDS')))
      sumstats.cur = subset(sumstats,SNP %in% unique(c(twas_model$SNP,
                                                       isotwas_model$SNP)))
      tot = rbind(twas_model[,c('SNP','Chromosome','Position')],
                  isotwas_model[,c('SNP','Chromosome','Position')])
      tot = tot[!duplicated(tot$SNP),]
      sumstats.cur = merge(sumstats.cur,tot,by = 'SNP')
      sumstats.cur$Beta = sumstats.cur$Z
      sumstats.cur$SE = 1
      twas_model = subset(twas_model,SNP %in% sumstats$SNP)
      isotwas_model = subset(isotwas_model,
                             SNP %in% sumstats$SNP)
      twas_model$A1 = twas_model$REF
      twas_model$A2 = twas_model$ALT
      twas_model$Transcript = twas_model$Feature
      if (nrow(twas_model) > 1){twas_model = as.data.frame(apply(twas_model,2,unlist))}
      twas_model$Chromosome = as.numeric(twas_model$Chromosome)
      twas_model$Weight = as.numeric(twas_model$Weight)
      twas_model$R2 = as.numeric(twas_model$R2)
      twas_model$Position = as.numeric(twas_model$Position)
      isotwas_model$A1 = isotwas_model$REF
      isotwas_model$A2 = isotwas_model$ALT
      isotwas_model$Transcript = isotwas_model$Feature
      print(max(twas_model$R2))
      print(max(isotwas_model$R2))
      
      ### RUN TWAS ASSOCIATION
      if (nrow(twas_model) > 0){
        gene_df_twas = burdenTest(mod = twas_model,
                                  ld = LD_twas,
                                  gene = gene,
                                  sumStats = sumstats.cur,
                                  chr = 'Chromosome',
                                  pos = 'Position',
                                  a1 = 'A1',
                                  a2 = 'A2',
                                  Z = 'Z',
                                  beta = 'Beta',
                                  se = 'SE',
                                  R2cutoff = .01,
                                  alpha = 1e-3,
                                  nperms = 1e3,
                                  usePos = F)
        if (class(gene_df_twas) == 'data.frame'){
          gene_df_twas$R2 = unlist(twas_model$R2)[1]
          fwrite(gene_df_twas,outFile_twas,
                 append = T, sep = '\t',
                 quote = F, row.names=F)
        }
      }
      
      
      isotwas_model = subset(isotwas_model,R2 >= 0.01)
      for (tx in unique(isotwas_model$Transcript)){
        
        if (nrow(isotwas_model) > 0){
          #print(tx)
          tx_df_isotwas = burdenTest(mod = subset(isotwas_model,
                                                  Transcript == tx),
                                     ld = LD_isotwas,
                                     gene = gene,
                                     sumStats = sumstats.cur,
                                     chr = 'Chromosome',
                                     pos = 'Position',
                                     a1 = 'A1',
                                     a2 = 'A2',
                                     Z = 'Z',
                                     beta = 'Beta',
                                     se = 'SE',
                                     R2cutoff = 0.01,
                                     alpha = 1e-3,
                                     nperms = 1e3,
                                     usePos = F)
          if (class(tx_df_isotwas) == 'data.frame'){
            
            tx_df_isotwas$R2 = subset(isotwas_model,
                                      Transcript == tx)$R2[1]
            colnames(tx_df_isotwas) = c('Gene','Transcript',
                                        'Z','P','permute.P','topSNP',
                                        'topSNP.P','R2')
            
            fwrite(tx_df_isotwas,outFile_isotwas,
                   append = T, sep = '\t',
                   quote = F, row.names=F)
            
          }
        }
        
      }
      
    }
    
    
    
    fwrite(data.frame(Trait = tr,
                      Gene = gene,
                      Index = 1),
           doneFile,
           sep='\t',
           row.names=F,
           quote=F,append=T)
  }
}


require(pbmcapply)
tc_assoc = function(s){
  tryCatch(runAssociation(s),
           error = function(e) {
             print(paste0('error with ',s))
           }
  )
}

lapply(toDoGenes,
       tc_assoc)
