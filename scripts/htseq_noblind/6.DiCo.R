library(tidyverse)
library(clusterProfiler)
library(cluster)

expr = readRDS('./data/htseq/expr_mat_no_blind.rds')
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

# saveRDS(genecovtidy2, file="./data/processed/tidy/CoV_wo_cortex.rds")
saveRDS(genecovtidy2, file="./data/htseq/no_blind/CoV_wo_cortex.rds")

#################### remove one individual from other tissues that lacks expression in cortex
expr = expr[, !colnames(expr) =='465']
names(age) = ind_id
samp_tissue = samp_tissue[!names(age)=="465"]
age = age[!names(age)=="465"]
ageu = age[unique(names(age))]


############################## calculate cov without each tissue (3 tissues, using 15 individuals):

ts = unique(samp_tissue)
cov3ts = lapply(unique(samp_tissue), function(i){
  print(paste('Calculating CoV without:',i))
  cov3tsX = sapply(unique(colnames(expr)), function(x){
    expX = expr[, !samp_tissue%in%i]
    sapply(rownames(expX), function(y){
      sd(expX[y, colnames(expX)==x]) / mean( expX[y, colnames(expX) == x] )
    })
  })
  cov3tsX
})

names(cov3ts) = paste0('wo_', ts)

cov3tstidy  = reshape2::melt(cov3ts) %>% 
  set_names(c('gene_id', 'ind_id', 'CoV', 'Excluded')) %>%
  mutate(ind_id = as.character(ind_id)) %>%
  left_join(data.frame(age=ageu, ind_id=names(ageu))) %>%
  mutate(period = ifelse(age<90,'Development', 'Ageing') )

# saveRDS(cov3tstidy, file='./data/processed/tidy/CoV_wo_eachtissue.rds')
saveRDS(cov3tstidy, file='./data/htseq/no_blind/CoV_wo_eachtissue.rds')

#################### calculate CoV among 4 tissues for each individual

# CoV_f = function(mat){
#   # this function calculates coeffient of variation (CoV) using the expression levels from different tissues
#   # of the same individual
#   
#   sapply(unique(colnames(mat)), function(x){
#     
#   })
#   
# }


genecov = sapply(unique(colnames(expr)), function(x){
  sapply(rownames(expr),function(y){
    sd(expr[y, colnames(expr) == x]) / mean(expr[y, colnames(expr) == x])
  })
})
genecovtidy = reshape2::melt(genecov) %>%
  set_names('gene_id','ind_id','CoV')
# saveRDS(genecovtidy, file="./data/processed/tidy/CoV.rds")
saveRDS(genecovtidy, file="./data/htseq/no_blind/CoV.rds")

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
# saveRDS(covch, './data/processed/tidy/CoV_change.rds')
saveRDS(covch, './data/htseq/no_blind/CoV_change.rds')

#################### CoV-age correlations for 3 tissues;

# head(cov3ts$wo_Cortex)
# head(cov3tstidy)
# 
# covch3ts = cov3tstidy %>%
#   group_by(period, Excluded, gene_id) %>%
#   summarise(rho  = cor.test(CoV, age, m='s')$est,
#             pval = cor.test(CoV, age, m='s')$p.val)
# 
# covch3ts = covch3ts %>%
#   group_by(period) %>%
#   mutate(BH = p.adjust(pval, method = 'BH'))

# saveRDS(covch3ts, './data/processed/tidy/CoV_change_wo_eachtissue.rds')


#######################
#######################
#######################
#######################

#######################
#######################
####################### test conv/div ratio with jacknife: (not dc but cd ratio is calculated)
####################### (for significant genes)
#######################


# jk_CoDi_ratio_test = function(mat, age){
#   # exclude each colu
#   
# }

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

# saveRDS(cbind(CI95_dev, CI95_aging),
#         './data/processed/raw/cov_ratio_jk_CI.rds')

saveRDS(cbind(CI95_dev, CI95_aging),
                 './data/htseq/no_blind/cov_ratio_jk_CI.rds')
#################### 
#################### 
#################### test for all genes
#################### 

### test development:
# divconv_ratio_all = sapply(colnames(genecov.dev), function(exc){
#   excgenecov = genecov.dev[, !colnames(genecov.dev)%in%exc]
#   excage = ageu[!names(ageu)%in%exc]
#   excage = excage[excage<90]
#   sub_cov_dev = t(sapply(rownames(excgenecov), function(x){
#     a = cor.test(excgenecov[x,], excage, m='s' )
#     c(a$est, a$p.val)
#   }))
#   colnames(sub_cov_dev) = c('rho', 'p')
#   #sub_cov_dev_sig = sub_cov_dev[sub_cov_dev[,3]<0.1, ]
#   div_conv_ratio = sum(sub_cov_dev[,1]<0) / sum(sub_cov_dev[,1] > 0)
#   return(div_conv_ratio)
# })

## jackknife standard error:
# dc_all = cbind(Dev = c( sum(covch_D[, 1]>0),
#                         sum(covch_D[, 1]<0)),
#                Ageing = c( sum(covch_A[, 1]>0),
#                            sum(covch_A[, 1]<0)))
# rownames(dc_all) = c('Divergent','Convergent')
# 
# obs_all = dc_all[2,'Dev'] / dc_all[1,'Dev']
# emprical_average_all = mean(divconv_ratio_all)
# jk_se_all = round(( sum((divconv_ratio_all - emprical_average_all)^2)*(6/7) )^0.5 , 2)
# ## jackknife standard error:
# Q_dot_all = mean(divconv_ratio_all)
# jk_bias_all = (7-1)*(Q_dot_all - obs_all)
# 
# ## http://www.comparingpartitions.info/?link=Tut9
# vard_all = sum((divconv_ratio_all - emprical_average_all)^2) / 6
# CI95_dev_all = c(emprical_average_all - 2*(vard_all/7)^0.5,
#                  emprical_average_all + 2*(vard_all/7)^0.5)

## aging
# divconv_ratio_aging_all = sapply(colnames(genecov.aging), function(exc){
#   excgenesd = genecov.aging[, !colnames(genecov.aging)%in%exc]
#   excage = ageu[!names(ageu)%in%exc]
#   excage = excage[excage>90]
#   sub_cov = t(sapply(rownames(excgenesd), function(x){
#     a = cor.test(excgenesd[x,], excage, m='s' )
#     c(a$est, a$p.val)
#   }))
#   #sub_cov = cbind(sub_cov, p.adjust(sub_cov[,2], m='BH'))
#   colnames(sub_cov) = c('rho', 'p')
#   div_conv_ratio = sum(sub_cov[,1]<0) / sum(sub_cov[,1] > 0)
#   return(div_conv_ratio)
# })

## jackknife standard error:
# obs_all = dc_all[2,'Ageing'] / dc_all[1,'Ageing']
# emprical_average_all = mean(divconv_ratio_aging_all)
# jk_se_all = round(( sum((divconv_ratio_aging_all - emprical_average_all)^2, na.rm=T)*(6/7) )^0.5 , 2)
## jackknife standard error:
# Q_dot_all = mean(divconv_ratio_aging_all)
# jk_bias_all = (7-1)*(Q_dot_all - obs_all)
# 
# vara_all = sum((divconv_ratio_aging_all - emprical_average_all)^2, na.rm=T) / 6
# CI95_aging_all = c(emprical_average_all - 2*(vara_all/7)^0.5,
#                    emprical_average_all + 2*(vara_all/7)^0.5)

# saveRDS(cbind(CI95_dev_all, CI95_aging_all),
#        './data/processed/raw/cov_ratio_jk_CI_nosig.rds')

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

# saveRDS(ddc_genes, './data/processed/raw/dev_divergent_genes_dc_rank.rds')
saveRDS(ddc_genes, './data/htseq/no_blind/dev_divergent_genes_dc_rank.rds')
ddc_genes_qn = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
length(intersect(names(ddc_genes), names(ddc_genes_qn) ) )
length(intersect( names(which(ddc_genes_qn<0)), names(which(ddc_genes<0)) ))

dico_old = names(which(ddc_genes_qn<0))
dico = names(which(ddc<0))

sg = intersect(names(ddc), names(ddc_genes_qn) ) 
sg = intersect(dico, dico_old ) 
head(ddc)
head(ddc_genes_qn)
cor(ddc[sg], ddc_genes_qn[sg], m='s')

#


# run here:
dc_gse = gseGO(geneList = sort(ddc_genes, decreasing = T), OrgDb = "org.Mm.eg.db", ont = "BP", pvalueCutoff = 1,
               keyType = "ENSEMBL", nPerm = 1000, minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BY',
               verbose = F)
# saveRDS(dc_gse, file="./data/processed/raw/dc_gse.rds")
saveRDS(dc_gse, file="./data/htseq/no_blind/dc_gse.rds")

# dc_gse = readRDS(file="./data/processed/raw/dc_gse.rds")

dc_gse_genelist =  strsplit(dc_gse@result[,11], split = '/')
names(dc_gse_genelist) = dc_gse@result[,'ID']
dc_gse_genelist =reshape2::melt(dc_gse_genelist) %>% 
  set_names(c('gene_id','GO_ID'))
dc_gse_table = list(enrichment = dc_gse@result[,1:10],
                    genelist = dc_gse_genelist)
# write.xlsx(dc_gse_table, './results/SI_tables/TableS11.xlsx')

#################### 
####################  Test enrichment of reversal genes among DC genes with fisher test
#################### 
####################


#################### shared reversal genes vs dc genes (among development divergent genes)

rev.shared = readRDS('./data/htseq/no_blind/revgenes.shared.rds')

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

# OR1: 0.5, p< 1e-10

#################### reversal genes (for each tissue) vs dc genes (among development divergent genes)
revgenes = readRDS('./data/htseq/no_blind/revgenes.tissues_tidy.rds')

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

OR2 

# tissue    OR        p
# <chr>  <dbl>    <dbl>
#   1 Cortex  1.07 7.14e- 2
# 2 Liver   1.71 7.18e-43
# 3 Lung    1.28 2.47e-10
# 4 Muscle  1.56 1.80e-30
####################
#################### reversal genes (for each tissue) vs dc genes (among development divergent genes)
#################### test updown and downup genes separately
####################

OR_UD = data.frame(gene_id=rep(names(ddc_genes),4), DCpattern = rep(ddc_genes,4) ,
                   tissue= rep(c('Cortex','Lung', 'Liver', 'Muscle'),each = length(ddc_genes)), row.names=NULL) %>%
  mutate(DCpattern = ifelse(DCpattern<0,'DC', 'DD')) %>%
  left_join(revgenes) %>%
  mutate(Revness = ifelse(direction == 'UpDown', 'Rev','NonRev')) %>% 
  select(-direction) %>% 
  mutate(Revness = replace_na(Revness, 'NonRev')) %>%
  group_by(Revness, tissue, DCpattern) %>%
  summarise(n=n()) %>%
  spread(Revness, n) %>%
  group_by(tissue) %>%
  summarise(OR = fisher.test(cbind(Rev, NonRev))$est,
            p = fisher.test(cbind(Rev, NonRev))$p.val)
OR_UD
# tissue    OR           p
# <chr>  <dbl>       <dbl>
#   1 Cortex 0.939 0.174      
# 2 Liver  1.06  0.221      
# 3 Lung   1.20  0.0000501  
# 4 Muscle 1.25  0.000000409

OR_DU = data.frame(gene_id=rep(names(ddc_genes),4), DCpattern = rep(ddc_genes,4) ,
                   tissue= rep(c('Cortex','Lung', 'Liver', 'Muscle'),each = length(ddc_genes)), row.names=NULL) %>%
  mutate(DCpattern = ifelse(DCpattern<0,'DC', 'DD')) %>%
  left_join(revgenes) %>%
  mutate(Revness = ifelse(direction == 'DownUp', 'Rev','NonRev')) %>% 
  select(-direction) %>% 
  mutate(Revness = replace_na(Revness, 'NonRev')) %>%
  group_by(Revness, tissue, DCpattern) %>%
  summarise(n=n()) %>%
  spread(Revness, n) %>%
  group_by(tissue) %>%
  summarise(OR = fisher.test(cbind(Rev, NonRev))$est,
            p = fisher.test(cbind(Rev, NonRev))$p.val)
OR_DU
# tissue    OR        p
# <chr>  <dbl>    <dbl>
#   1 Cortex  1.19 2.92e- 4
# 2 Liver   1.73 2.64e-40
# 3 Lung    1.15 1.78e- 3
# 4 Muscle  1.37 3.95e-14

####################
####################
#################### Test DiCo genes with permutations among Dev. Divergent genes
####################
####################

# covA = readRDS("/home/black/Dropbox/projects/ageing/proc_data/perm/aging.cov.perm.rds")
# head(div_dev)

# DiCo/DiDi: 
# div_dev = covch_D[covch_D[,1] > 0, 1]
# ddc_genes = div_dev*covch_A[names(div_dev),1]
# co.t = summary(covch_A[covch_A[,1]<0,1])['3rd Qu.']
# co = covch_A[covch_A[,1]<co.t,1]

# ddc_genes2 = div_dev[names(div_dev)%in%names(co)]*covch_A[rownames(covch_A)%in%intersect(names(div_dev), names(co)), 1]
# length(ddc_genes)

# -:dico, +:didi
# obs_dico = sum(ddc_genes<0) / (sum(ddc_genes<0)+sum(ddc_genes>0))
# 
# DiCodist = sapply(1:ncol(covA), function(perm){
#   pxDiCo = sum(covA[names(div_dev),perm]<0)
#   pxDiDi = sum(covA[names(div_dev),perm]>0)
#   return(pxDiCo / (pxDiCo+pxDiDi))
# })
# names(DiCodist) = paste0('p',1:length(DiCodist))
# dico_p = round(mean(DiCodist >= obs_dico),3)
# dico_fpr = round(median(DiCodist) / obs_dico,3)
#saveRDS(data.frame(Perms= names(DiCodist), rho = DiCodist, row.names=NULL),
#        './data/processed/tidy/DiCo_sig_perm.rds')

# dico_prop_perm = data.frame(Perms= names(DiCodist), rho = DiCodist, row.names=NULL) %>% 
#   ggplot(aes(x=rho)) +
#   geom_histogram(bins=40) +
#   geom_vline(xintercept = obs_dico, linetype='dashed', col='darkred') +
#   annotate('text', x= 0.55, y=65, label=paste('Obs = ', round(obs_dico,2)), hjust = 0, size=2) + 
#   annotate('text', x= 0.55, y=63, label=paste('FPR = ',dico_fpr), hjust = 0, size=2) + 
#   annotate('text', x= 0.55, y=61, label=paste('p.val = ',dico_p), hjust = 0, size=2) +
#   xlab('DiCo Proportion') +
#   ylab('Frequency') +
#   theme_bw()

#ggsave('./results/SI_figures/Figure_S30.pdf',units='cm', width=12, height=12, useDingbats=F)
#ggsave('./results/SI_figures/Figure_S30.png',units='cm', width=12, height=12)

#############
############# DiCo with rho cutoff
#############
# -:dico, +:didi
# obs_dico = sum(ddc_genes<0) / (sum(ddc_genes<0)+sum(ddc_genes>0))
# 
# DiCodist = sapply(1:ncol(covA), function(perm){
#   pxDiCo = sum(covA[names(div_dev),perm]<0)
#   pxDiDi = sum(covA[names(div_dev),perm]>0)
#   return(pxDiCo / (pxDiCo+pxDiDi))
# })
# names(DiCodist) = paste0('p',1:length(DiCodist))
# dico_p = round(mean(DiCodist >= obs_dico),3)
# dico_fpr = round(median(DiCodist) / obs_dico,3)
# 
# dico_prop_perm = data.frame(Perms= names(DiCodist), rho = DiCodist, row.names=NULL) %>% 
#   ggplot(aes(x=rho)) +
#   geom_histogram(bins=40) +
#   geom_vline(xintercept = obs_dico, linetype='dashed', col='darkred') +
#   annotate('text', x= 0.55, y=65, label=paste('Obs = ', round(obs_dico,2)), hjust = 0, size=2) + 
#   annotate('text', x= 0.55, y=63, label=paste('FPR = ',dico_fpr), hjust = 0, size=2) + 
#   annotate('text', x= 0.55, y=61, label=paste('p.val = ',dico_p), hjust = 0, size=2) +
#   xlab('DiCo Proportion') +
#   ylab('Frequency') +
#   theme_bw()
# 
# ggsave('./results/SI_figures/Figure_S36.pdf',units='cm', width=12, height=12, useDingbats=F)
# ggsave('./results/SI_figures/Figure_S36.png',units='cm', width=12, height=12)


############
############ development convergent genes
# conv_dev = covch_D[covch_D[,1]<0,1] # 5939
# dcc_genes = conv_dev*covch_A[names(conv_dev),1] # +:coco, -:codi
# 
# obs_coco = sum(dcc_genes>0) / (sum(dcc_genes>0) + sum(dcc_genes<0))
# 
# CoCodist = sapply(1:ncol(covA), function(perm){
#   pxCoCo = sum(covA[names(conv_dev),perm]>0)
#   pxCoDi = sum(covA[names(conv_dev),perm]<0)
#   return(pxCoCo / (pxCoCo+pxCoDi))
# })
# names(CoCodist) = paste0('p',1:length(CoCodist))
# coco_p = round(mean(CoCodist >= obs_coco),3)
# coco_fpr = round(median(CoCodist) / obs_coco,3)
# 
# data.frame(Perms= names(CoCodist), rho = CoCodist, row.names=NULL) %>% 
#   ggplot(aes(x=rho)) +
#   geom_histogram(bins=40) +
#   geom_vline(xintercept = obs_coco, linetype='dashed', col='darkred') +
#   annotate('text', x= 0.55, y=65, label=paste('Obs = ', round(obs_coco,2)), hjust = 0, size=2) + 
#   annotate('text', x= 0.55, y=63, label=paste('FPR = ',coco_fpr), hjust = 0, size=2) + 
#   annotate('text', x= 0.55, y=61, label=paste('p.val = ',coco_p), hjust = 0, size=2) +
#   xlab('DiCo Proportion') +
#   ylab('Frequency') +
#   theme_bw()
# 
# 
# 












