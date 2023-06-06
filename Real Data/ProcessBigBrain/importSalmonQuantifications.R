require(data.table)
require(rjson)
setwd('/u/project/gandalm/shared/GenomicDatasets-processed/PEC_AMP_AD_Arjun/Metadata')

rnaseq = fread('rnaseq_20230501.csv')

require(tximeta)
genotype = fread('genotype.csv')
individual = fread('individual_20230421.csv')
individual = subset(individual,ageDeath > 0)
genotype = merge(genotype,
                 individual[,c('individualID',
                               'genotypingID')],
                 by = 'genotypingID')

## ROSMAP is sample id mismatch
### Mayo is Temporal_Cortex
### BrainGVEX is sample id mismatch
### MSBB is frontal pole
### UCLA-ASD is dorsolateral prefrontal cortex
### LIBD_szcontrol/CMC_HBCC have ones in commone

dupIDs = rnaseq$individualID[duplicated(rnaseq$individualID)]
dups = subset(rnaseq,individualID %in%
                dupIDs)
rnaseq_nondup = subset(rnaseq,!individualID %in%
                         rnaseq$individualID[duplicated(rnaseq$individualID)])
sum(rnaseq_nondup$individual %in% genotype$individualID)

rnaseq_nondup = rbind(rnaseq_nondup,
                      subset(dups,
                             study == 'Mayo' &
                               tissue != 'cerebellum'))
sum(rnaseq_nondup$individual %in% genotype$individualID)

rnaseq_nondup = rbind(rnaseq_nondup,
                      subset(dups,
                             study == 'MSBB' &
                               tissue == 'frontal pole'))
sum(rnaseq_nondup$individual %in% genotype$individualID)

gvex_dups = subset(dups,
                   study == 'BrainGVEX')
www = vector('numeric',length=nrow(gvex_dups))
for (i in 1:nrow(gvex_dups)){
  if (grepl(gvex_dups$individualID[i],gvex_dups$path[i])){
    www[i] = 1
  }
}
rnaseq_nondup = rbind(rnaseq_nondup,
                      subset(gvex_dups,www==1))
sum(rnaseq_nondup$individual %in% genotype$individualID)

rnaseq_nondup = rbind(rnaseq_nondup,
                      subset(dups,study=='UCLA-ASD' &
                               tissue == 'dorsolateral prefrontal cortex'))
sum(rnaseq_nondup$individual %in% genotype$individualID)

dup_remaining = subset(dups,!study %in% c('Mayo','MSBB','BrainGVEX','UCLA-ASD') &
                         individualID != '')
dup_remaining = dup_remaining[order(dup_remaining$individualID),]
rnaseq_nondup = rbind(rnaseq_nondup,
                      subset(dup_remaining,grepl('PITT',individualID) &
                               !grepl('_BP_',specimenID)))
sum(rnaseq_nondup$individual %in% genotype$individualID)
dup_remaining = subset(dup_remaining,!grepl('PITT',individualID))
dup_remaining = dup_remaining[order(dup_remaining$RIN,decreasing = T),]
dup_remaining = dup_remaining[!duplicated(dup_remaining$individualID)]
rnaseq_nondup = rbind(rnaseq_nondup,
                      dup_remaining)
sum(rnaseq_nondup$individualID %in% genotype$individualID)

still_remaining = rnaseq_nondup$individualID[duplicated(rnaseq_nondup$individualID)]

dups = subset(rnaseq_nondup,
              individualID %in% still_remaining)
rnaseq_nondup = subset(rnaseq_nondup,
                       !individualID %in% still_remaining)

dups = dups[order(dups$RIN,decreasing = T),]
dups = dups[!duplicated(dups$individualID),]
rnaseq_nondup = rbind(rnaseq_nondup,dups)
sum(rnaseq_nondup$individualID %in% genotype$individualID)
rnaseq_nondup = subset(rnaseq_nondup, RIN > 5.5)

rnaseq = rnaseq_nondup

total_id = unique(intersect(genotype$individualID,
                            rnaseq$individualID))

rnaseq = subset(rnaseq,individualID %in% total_id)
rnaseq = rnaseq[!duplicated(rnaseq$individualID),]
genotype = subset(genotype,individualID %in% total_id)
genotype = genotype[!duplicated(genotype$individualID),]
fwrite(individual,'/u/scratch/a/abtbhatt/renormalize_0429/individual_inrnaseq.tsv',
       sep = '\t',col.names=T,row.names=F)
fwrite(genotype,'/u/scratch/a/abtbhatt/renormalize_0429/individual_hash.tsv',
       sep='\t',col.names=T,row.names=F)


fwrite(rnaseq,
       '/u/scratch/a/abtbhatt/renormalize_0429/rnaseq_050123_nodups.tsv',
       col.names=T,
       row.names=F,
       quote=F,
       sep='\t')


ddd = rnaseq[,c('individualID','bamFilePath')]
fwrite(ddd,'/u/scratch/a/abtbhatt/bigbrain_rnaseqqc/bamfilepaths.txt',
       col.names=F,row.names=F,quote=F,sep='\t')


outFolder = '/u/scratch/a/abtbhatt/renormalize_0429'
dir.create(outFolder)
setwd(outFolder)
dir.create('chunks')
require(SummarizedExperiment)

rnaseq_read = merge(rnaseq,genotype[,c('individualID',
                                       'genotypingID')],
                    by = 'individualID')
cd = data.frame(files = file.path(rnaseq_read$path,
                                  'quant.sf'),
                names = rnaseq_read$individualID,
                batch = rnaseq_read$study,
                genotypingID = rnaseq_read$genesDetect,
                bamPath = rnaseq_read$bamFilePath)
se = tximeta(cd,dropInfReps = T,
             countsFromAbundance = 'lengthScaledTPM')
saveRDS(se,
        file.path(outFolder,paste0('AdultBigBrain_tx_exp_raw_042923.RDS')))
