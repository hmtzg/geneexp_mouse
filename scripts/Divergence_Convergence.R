library(tidyverse)
library(clusterProfiler)
library(cluster)

expr = readRDS('./data/processed/raw/expression.rds')
ind_id = readRDS('./data/processed/raw/individual_ids.rds')
age = readRDS('./data/processed/raw/ages.rds')
samp_tissue = readRDS('./data/processed/raw/tissue_ids.rds')
colnames(expr) = ind_id

############################## calculate cov without cortex (3 tissues):

exp2 = expr[, !samp_tissue == "Cortex"]
genecov2 = sapply(unique(colnames(exp2)), function(x){
  sapply(rownames(exp2),function(y){
    sd(exp2[y, colnames(exp2) == x]) / mean(exp2[y, colnames(exp2) == x])
  })
})

genecovtidy2 = reshape2::melt(genecov2) %>%
  set_names('gene_id','ind_id','CoV')

saveRDS(genecovtidy2, file="./data/processed/tidy/CoV_wo_cortex.rds")

#################### remove one individual from other tissues that lacks expression in cortex
expr = expr[, !colnames(expr) =='465']
names(age) = ind_id
samp_tissue = samp_tissue[!names(age)=="465"]
age = age[!names(age)=="465"]
ageu = age[unique(names(age))]

#################### calculate CoV among tissues for each individual
genecov = sapply(unique(colnames(expr)),function(x){
  sapply(rownames(expr),function(y){
    sd(expr[y, colnames(expr) == x]) / mean(expr[y, colnames(expr) == x])
  })
})
genecovtidy = reshape2::melt(genecov) %>%
  set_names('gene_id','ind_id','CoV')
saveRDS(genecovtidy, file="./data/processed/tidy/CoV.rds")

#################### calculate CoV-age correlations in development and ageing

covch_D = t(sapply(rownames(genecov),function(x){
  x_dev = cor.test(genecov[x, ageu<90], ageu[ageu<90], m ="s")
  c(x_dev$est, x_dev$p.val)
}))
covch_D = cbind(covch_D, p.adjust(covch_D[,2], method = "BH"))
colnames(covch_D)[c(2,3)] = c("pval"," BH")

covch_A = t(sapply(rownames(genecov),function(x){
  x_aging = cor.test(genecov[x, ageu>90], ageu[ageu>90], m ="s")
  c(x_aging$est, x_aging$p.val)
}))
covch_A = cbind(covch_A, p.adjust(covch_A[,2], method = "BH"))
colnames(covch_A)[c(2,3)] = c("pval"," BH")

covch_dev = data.frame(covch_D, period = 'development', gene_id = rownames(covch_D), row.names = NULL) %>%
  rename(CoV_change = rho, 
         FDR = X.BH)

covch_aging = data.frame(covch_A, period = 'aging', gene_id = rownames(covch_A), row.names = NULL) %>%
  rename(CoV_change = rho, 
         FDR = X.BH)

covch = rbind(covch_dev, covch_aging)
saveRDS(covch, './data/processed/tidy/CoV_change.rds')


#######################
#######################
#######################
#######################

#######################
#######################
####################### test conv/div difference with jacknife: (not dc but cd ratio is calculated)
#######################
#######################


genecov.dev = genecov[, ageu < 90]
genecov.aging = genecov[, ageu > 90]
### test development:
divconv_ratio = sapply(colnames(genecov.dev), function(exc){
  excgenecov = genecov.dev[, !colnames(genecov.dev)%in%exc]
  excage = ageu[!names(ageu)%in%exc]
  excage = excage[excage<90]
  sub_cov_dev = t(sapply(rownames(excgenecov), function(x){
    a = cor.test(excgenecov[x,], excage, m='s' )
    c(a$est, a$p.val)
  }))
  sub_cov_dev = cbind(sub_cov_dev, p.adjust(sub_cov_dev[,2], m='BH'))
  colnames(sub_cov_dev) = c('rho', 'p', 'BH')
  sub_cov_dev_sig = sub_cov_dev[sub_cov_dev[,3]<0.1, ]
  div_conv_ratio = sum(sub_cov_dev_sig[,1]<0) / sum(sub_cov_dev_sig[,1] > 0)
  return(div_conv_ratio)
})

## jackknife standard error:
dc_sig = cbind(Dev = c( sum(covch_D[covch_D[,3]< 0.1, 1]>0),
                        sum(covch_D[covch_D[,3]< 0.1, 1]<0)),
               Ageing = c( sum(covch_A[covch_A[,3]< 0.1, 1]>0),
                        sum(covch_A[covch_A[,3]< 0.1, 1]<0)))
rownames(dc_sig) = c('Divergent','Convergent')
  
obs = dc_sig[2,'Dev'] / dc_sig[1,'Dev']
emprical_average = mean(divconv_ratio, na.rm=T)
jk_se = round(( sum((divconv_ratio - emprical_average)^2, na.rm=T)*(4/5) )^0.5 , 2)
## jackknife standard error:
Q_dot = mean(divconv_ratio, na.rm=T)
jk_bias = (5-1)*(Q_dot - obs)

## http://www.comparingpartitions.info/?link=Tut9
vard = sum((divconv_ratio - emprical_average)^2, na.rm=T) / 4
CI95_dev = c(emprical_average - 2*(vard/5)^0.5,
             emprical_average + 2*(vard/5)^0.5)

#### ageing:
divconv_ratio_aging = sapply(colnames(genecov.aging), function(exc){
  excgenesd = genecov.aging[, !colnames(genecov.aging)%in%exc]
  excage = ageu[!names(ageu)%in%exc]
  excage = excage[excage>90]
  sub_cov = t(sapply(rownames(excgenesd), function(x){
    a = cor.test(excgenesd[x,], excage, m='s' )
    c(a$est, a$p.val)
  }))
  sub_cov = cbind(sub_cov, p.adjust(sub_cov[,2], m='BH'))
  colnames(sub_cov) = c('rho', 'p', 'BH')
  sub_cov_sig = sub_cov[sub_cov[,3]<0.1, ]
  div_conv_ratio = sum(sub_cov_sig[,1]<0) / sum(sub_cov_sig[,1] > 0)
  return(div_conv_ratio)
})

## jackknife standard error:
obs = dc_sig[2,'Ageing'] / dc_sig[1,'Ageing']
emprical_average = mean(divconv_ratio_aging, na.rm=T)
jk_se = round(( sum((divconv_ratio - emprical_average)^2, na.rm=T)*(5/6) )^0.5 , 2)
## jackknife standard error:
Q_dot = mean(divconv_ratio, na.rm=T)
jk_bias = (6-1)*(Q_dot - obs)

vara = sum((divconv_ratio_aging - emprical_average)^2, na.rm=T) / 6
CI95_aging = c(emprical_average - 2*(vara/6)^0.5,
               emprical_average + 2*(vara/6)^0.5)

saveRDS(cbind(CI95_dev, CI95_aging),
        './data/processed/raw/cov_ratio_jk_CI.rds')

#################### 
#################### 
#################### 
#################### 


#################### 
#################### 
#################### Divergent-Convergent genes GSEA among development divergent genes:
#################### 
div_dev = covch_D[covch_D[,1] > 0, 1]
ddc_genes = div_dev*covch_A[names(div_dev),1]

saveRDS(ddc_genes, './data/processed/raw/dev_divergent_genes_dc_rank.rds')

dc_gse = gseGO(geneList = sort(ddc_genes, decreasing = T), OrgDb = "org.Mm.eg.db", ont = "BP", pvalueCutoff = 1,
                keyType = "ENSEMBL", nPerm = 1000, minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BY',
                verbose = F)
saveRDS(dc_gse, file="./data/processed/raw/dc_gse.rds")

# dc_gse = readRDS(file="./data/processed/raw/dc_gse.rds")
dc_gse_genelist =  strsplit(dc_gse@result[,11], split = '/')
names(dc_gse_genelist) = dc_gse@result[,'ID']
dc_gse_genelist =reshape2::melt(dc_gse_genelist) %>% 
  set_names(c('gene_id','GO_ID'))
dc_gse_table = list(enrichment = dc_gse@result[,1:10],
                    genelist = dc_gse_genelist)
write.xlsx(dc_gse_table, './data/Table_S10.xlsx')

#################### 
####################  Test enrichment of reversal genes among DC genes with fisher test
#################### 
####################


#################### shared reversal genes vs dc genes (among development divergent genes)
rev.shared = readRDS('./data/processed/tidy/revgenes.shared.rds')

OR1 = rev.shared %>%
  right_join(data.frame(gene_id=names(ddc_genes), DCpattern = ddc_genes, row.names=NULL) ) %>% 
  mutate(DCpattern = ifelse(DCpattern<0,'DC', 'DD')) %>%
  mutate(Revness = ifelse(Pattern == 'UpDown'| Pattern == 'DownUp','Rev','NonRev')) %>%
  mutate(Revness = replace_na(Revness, 'NonRev')) %>%
  select(-Pattern) %>%
  mutate(Revness = replace_na(Revness, 'NonRev')) %>%
  group_by(Revness, DCpattern) %>%
  summarise(n=n()) %>% 
  spread(Revness, n) %>%
  summarise(OR = fisher.test(cbind(Rev,NonRev))$est,
            pval = fisher.test(cbind(Rev,NonRev))$p.val)

#################### reversal genes (for each tissue) vs dc genes (among development divergent genes)
revgenes = readRDS('./data/processed/tidy/revgenes.tissues.rds')

OR2 = data.frame(gene_id=rep(names(ddc_genes),4), DCpattern = rep(ddc_genes,4) ,
           tissue= rep(c('Cortex','Lung', 'Liver', 'Muscle'),each = length(ddc_genes)), row.names=NULL) %>%
  mutate(DCpattern = ifelse(DCpattern<0,'DC', 'DD')) %>%
  left_join(revgenes) %>%
  mutate(Revness = ifelse(direction == 'UpDown'| direction == 'DownUp','Rev','NonRev')) %>% 
  select(-direction) %>% 
  mutate(Revness = replace_na(Revness, 'NonRev')) %>%
  group_by(Revness, tissue, DCpattern) %>%
  summarise(n=n()) %>%
  spread(Revness, n) %>%
  group_by(tissue) %>%
  summarise(OR = fisher.test(cbind(Rev, NonRev))$est,
            p = fisher.test(cbind(Rev, NonRev))$p.val)
  
####################
#################### reversal genes (for each tissue) vs dc genes (among development divergent genes)
#################### test updown and downup genes separately
####################


####################
####################


