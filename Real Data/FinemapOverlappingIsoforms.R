### This script takes in isoTWAS results and runs FOCUS on overlapping
### transcripts

require(vroom)
require(data.table)
require(GenomicRanges)
require(Matrix)
require(bigsnpr)
require(xlsx)
require(isotwas)

traits = c('SCZ','BD','ALZ','ASD','ADHD')
setwd('/path/to/results')

print(tr)
setwd(paste0(tr,'_results'))
twas = fread('/path/to/isotwas/results')
res = twas[order(twas$Chromosome,twas$Start),]
res$GTA = paste(res$Gene,res$Tissue,sep=':')
res$TTA = paste(res$Transcript,res$Tissue,sep=':')

chr.table = table(res$Chromosome)
chr.un = unique(res$Chromosome)
keep.gta = c()
keep.tta = c()
for (c in chr.un){
  res.cur = subset(res,Chromosome == c)
  res.cur = res.cur[order(res.cur$Start),]
  if (nrow(res.cur) > 1){
    for (i in 1:(nrow(res.cur)-1)){
      if (res.cur$End[i] > res.cur$Start[i+1] - 1e6){
        keep.gta = unique(c(keep.gta,
                            c(res.cur$GTA[c(i,i+1)])))
        keep.tta = unique(c(keep.tta,
                            c(res.cur$TTA[c(i,i+1)])))
      }
    }
  }
}

new_res = subset(res,!TTA %in% keep.tta)
new_res$pip = 1
new_res$in_cred_set = T


res = subset(res,TTA %in% keep.tta)

chr.un = unique(res$Chromosome)
if (length(chr.un) >= 1){
  for (c in chr.un[1:length(chr.un)]){
    print(c)
    all.snps = c()
    omega = c()
    pos = c()
    gene = c()
    snp.chr = c()
    this.res.tot = subset(res,Chromosome == c)
    this.res.tot$Group = 1
    ggg = 1
    for (i in 1:(nrow(this.res.tot)-1)){
      if (this.res.tot$End[i] <= this.res.tot$Start[i+1] - 1e6){
        ggg = ggg+1
        this.res.tot$Group[(i+1):nrow(this.res.tot)] = ggg
      }
    }
    for (g in unique(this.res.tot$Group)){
      this.res = subset(this.res.tot,
                        Group == g)
      for (i in 1:nrow(this.res)){
        aaa=fread('/path/to/isotwas/model/tsv')
        aaa = subset(aaa,Transcript == this.res$Transcript[i])
        Model = data.frame(SNP = aaa$SNP,
                           Chromosome = aaa$Chromosome,
                           Position = aaa$Position,
                           Effect = aaa$Weight,
                           A1 = aaa$ALT,
                           A2 = aaa$REF)
        Model = subset(Model,Effect!=0)
        Model = Model[!duplicated(Model$SNP),]
        all.snps = c(all.snps,
                     as.character(Model$SNP))
        omega = c(omega,
                  as.numeric(Model$Effect))
        gene = c(gene,
                 rep(this.res$TTA[i],nrow(Model)))
        snp.chr = c(snp.chr,
                    as.numeric(Model$Chromosome))
        pos = c(pos,as.numeric(Model$Position))
      }
      tot.df = data.frame(SNP = all.snps,
                          Gene = gene,
                          Effect = omega,
                          Chromosome = snp.chr)
      model.df = as.data.frame(matrix(nrow = length(unique(all.snps)),
                                      ncol = nrow(this.res)+1))
      colnames(model.df) = c('SNP',this.res$TTA)
      model.df$SNP = as.character(unique(all.snps))
      for (g in 1:nrow(this.res)){
        print(this.res$TTA[g])
        cur.tot.df = subset(tot.df,Gene == this.res$TTA[g])
        cur.tot.df$SNP = as.character(cur.tot.df$SNP)
        for (i in 1:nrow(model.df)){
          print(i)
          w = which(cur.tot.df$SNP == model.df$SNP[i])
          model.df[i,g+1] = ifelse(length(w) != 0,
                                   cur.tot.df$Effect[w],
                                   0)
        }
      }
      model.df$Chromosome = c
      for (i in 1:nrow(model.df)){
        rrr = subset(tot.df,SNP == model.df$SNP[i])
        model.df$Chromosome[i] = rrr$Chromosome[1]
      }

      min = max(1,min(pos)-1e6)
      max = max(pos)+1e6

      system(paste0('plink ',
                    ' --bfile /path/to/plink/prefix/chr',c,
                    ' --chr ',c,
                    ' --from-bp ',min,
                    ' --to-bp ',max,
                    ' --make-bed --out temp'))

      snps = snp_attach(snp_readBed2(paste0('temp.bed'),
                                     backingfile = tempfile()))
      snp.set = subset(snps$map,
                       marker.ID %in% model.df$SNP)
      model.df = model.df[match(snp.set$marker.ID,
                                model.df$SNP),]
      snpMat = as.matrix(snp_attach(subset(snps,
                                           ind.col =
                                             which(snps$map$marker.ID %in%
                                                     model.df$SNP),
                                           backingfile=tempfile()))$genotypes[])
      rm(snp.set)
      snpMat = Matrix(t(snpMat), sparse = T)
      V = snp_cor(snp_attach(subset(snps,
                                    ind.col =
                                      which(snps$map$marker.ID %in%
                                              model.df$SNP),
                                    backingfile=tempfile()))$genotypes)

      Omega = Matrix(as.matrix(model.df[,-c(1,ncol(model.df))]))
      zscores = this.res$TWAS.Z
      m = length(zscores)
      wcor = estimate_cor(as.matrix(Omega),
                          as.matrix(V),intercept=T)[[1]]
      diag(wcor) = 1
      wcor[is.na(wcor)] = 0
      swld = estimate_cor(as.matrix(Omega),
                          as.matrix(V),
                          intercept=T)[[2]]
      null_res = m * log(1 - 1e-3)
      marginal = m * log(1 - 1e-3)
      comb_list = list()
      for (n in 1:min(3,length(zscores))){
        comb_list = c(comb_list,
                      combn(1:length(zscores),n,simplify=F))
      }
      pips = rep(0,length(zscores))
      zscores = get_resid(zscores,as.matrix(swld),as.matrix(wcor))[[1]]
      for (j in 1:length(comb_list)){
        subset = comb_list[[j]]
        local = bayes_factor(zscores,
                             idx_set = subset,
                             wcor = wcor)
        marginal = log(exp(local) + exp(marginal))
        for (idx in subset){
          if (pips[idx] == 0){
            pips[idx] = local
          } else {
            pips[idx] = log(exp(pips[idx]) + exp(local))
          }
        }
        print(pips)
        print(marginal)
      }

      iter=2
      while (any(pips == Inf) & iter >= 1){

        zscores = this.res$TWAS.Z
        m = length(zscores)
        null_res = m * log(1 - 1e-3)
        marginal = m * log(1 - 1e-3)
        comb_list = list()
        for (n in 1:min(iter,
                        length(zscores))){
          comb_list = c(comb_list,
                        combn(1:length(zscores),n,simplify=F))
        }

        pips = rep(0,length(zscores))
        for (j in 1:length(comb_list)){
          subset = comb_list[[j]]
          local = isotwas::bayes_factor(zscores,
                                        idx_set = subset,
                                        wcor = wcor)
          marginal = log(exp(local) + exp(marginal))
          for (idx in subset){
            if (pips[idx] == 0){
              pips[idx] = local
            } else {
              pips[idx] = log(exp(pips[idx]) + exp(local))
            }
          }
          print(pips)
          print(marginal)
        }
        iter = iter-1
      }

      pips = exp(pips - marginal)
      null_res = exp(null_res - marginal)
      this.res$pip = pips
      this.res = this.res[order(this.res$pip,decreasing = T),]
      npost = this.res$pip/sum(this.res$pip)
      csum = cumsum(npost)
      this.res$in_cred_set = F
      for (i in 1:nrow(this.res)){
        this.res$in_cred_set[i] = T
        if (csum[i] > .9){
          this.res$in_cred_set[i] = F
        }

      }
      if ('Group' %in% colnames(this.res)){
        w = which(colnames(this.res) == 'Group')
        this.res = as.data.frame(this.res)
        this.res = this.res[,-w]
      }
      this.res = this.res[,c('Gene','HGNC',
                             'Transcript','Chromosome',
                             'Start','End','Tissue','TWAS.Z',
                             'Screen.P',"Confirmation.P",
                             "Permutation.P","h2",
                             "h2.se","h2.p","R2",
                             "pip",'in_cred_set')]
      fwrite(this.res,'/path/to/focus/results',append=T,
             row.names= F,quote=F,sep='\t')
      }
    }
}
