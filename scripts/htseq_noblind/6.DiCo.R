library(tidyverse)
library(ggpubr)
source('./scripts/functions.R')

expr = readRDS('./data/htseq/expr_mat_no_blind.rds')
ind_id = readRDS('./data/processed/raw/individual_ids.rds')
age = readRDS('./data/processed/raw/ages.rds')
samp_tissue = readRDS('./data/processed/raw/tissue_ids.rds')
colnames(expr) = ind_id

############################## calculate cov without cortex (3 tissues):
# 
# exp2 = expr[, !samp_tissue == "Cortex"]
# genecov2 = sapply(unique(colnames(exp2)), function(x){
#   sapply(rownames(exp2),function(y){
#     sd(exp2[y, colnames(exp2) == x]) / mean(exp2[y, colnames(exp2) == x])
#   })
# })
# 
# genecovtidy2 = reshape2::melt(genecov2) %>%
#   set_names('gene_id','ind_id','CoV')
# 
# saveRDS(genecovtidy2, file="./data/htseq/no_blind/CoV_wo_cortex.rds")

#################### remove one individual from other tissues that lacks expression in cortex
expr = expr[, !colnames(expr) =='465']
names(age) = ind_id
samp_tissue = samp_tissue[!names(age)=="465"]
age = age[!names(age)=="465"]
ageu = age[unique(names(age))]


############################## calculate cov without each tissue (3 tissues, using 15 individuals):
# 
# ts = unique(samp_tissue)
# cov3ts = lapply(unique(samp_tissue), function(i){
#   print(paste('Calculating CoV without:',i))
#   cov3tsX = sapply(unique(colnames(expr)), function(x){
#     expX = expr[, !samp_tissue%in%i]
#     sapply(rownames(expX), function(y){
#       sd(expX[y, colnames(expX)==x]) / mean( expX[y, colnames(expX) == x] )
#     })
#   })
#   cov3tsX
# })
# 
# names(cov3ts) = paste0('wo_', ts)
# 
# cov3tstidy  = reshape2::melt(cov3ts) %>% 
#   set_names(c('gene_id', 'ind_id', 'CoV', 'Excluded')) %>%
#   mutate(ind_id = as.character(ind_id)) %>%
#   left_join(data.frame(age=ageu, ind_id=names(ageu))) %>%
#   mutate(period = ifelse(age<90,'Development', 'Ageing') )
# 
#saveRDS(cov3tstidy, file='./data/htseq/no_blind/CoV_wo_eachtissue.rds')

#################### calculate CoV among 4 tissues for each individual

genecov = sapply(unique(colnames(expr)), function(x){
  sapply(rownames(expr),function(y){
    sd(expr[y, colnames(expr) == x]) / mean(expr[y, colnames(expr) == x])
  })
})
genecovtidy = reshape2::melt(genecov) %>%
  set_names('gene_id','ind_id','CoV')

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
  rename(CoV_change = rho, FDR = X.BH)

covch_aging = data.frame(covch_A, period = 'aging', gene_id = rownames(covch_A), row.names = NULL) %>%
  rename(CoV_change = rho, FDR = X.BH)

covch = rbind(covch_dev, covch_aging)

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

genecov.dev = genecov[, ageu < 90]
genecov.aging = genecov[, ageu > 90]

### calculate leave-one-out ratios:
codi_dev = LOO_CoDiR(genecov.dev, ageu[ageu<90] )
codi_aging = LOO_CoDiR(genecov.aging, ageu[ageu>90] )

##### ## http://www.comparingpartitions.info/?link=Tut9

loo = list(dev = codi_dev, aging = codi_aging )

saveRDS(loo, file='./data/htseq/no_blind/codi_jk_pseudovalues.rds')

LOO_JK_CI(codi_dev)
LOO_JK_CI(codi_aging)

#################### 
#################### 
#################### test for all genes
#################### 

### calculate leave-one-out ratios:
codi_dev_all = LOO_CoDiR(genecov.dev, ageu[ageu<90], qval = 1.1 )
codi_aging_all = LOO_CoDiR(genecov.aging, ageu[ageu>90], qval = 1.1 )

loo_allgenes = list(dev = codi_dev_all, aging = codi_aging_all)

saveRDS(loo_allgenes, file='./data/htseq/no_blind/codi_allgenes_jk_pseudovalues.rds')

##### DiCo genes:

covch_vst = readRDS('./data/htseq/no_blind/CoV_change.rds')
covch_qn = readRDS('./data/processed/tidy/CoV_change.rds')
# covch_qn  = covch_qn %>%
#   mutate(period = gsub('development', 'Development', period)) %>%
#   mutate(period = gsub('aging', 'Ageing', period))

covch_vst %>%
  rename(CoV_change_vst = CoV_change) %>%
  inner_join(covch_qn, by=c('gene_id', 'period') ) %>%
  group_by(period) %>%
  summarise(rho = cor(CoV_change_vst, CoV_change, m='s') )
# period        rho
# <chr>       <dbl>
# 1 aging       0.690
# 2 development 0.693

covchplot = covch_vst %>%
  rename(CoV_change_vst = CoV_change) %>%
  inner_join(covch_qn, by=c('gene_id', 'period') ) %>%
  mutate(period = factor(period, levels=c('development', 'ageing'))) %>%
  ggplot(aes(x=CoV_change, y=CoV_change_vst)) +
  facet_grid(~period) +
  geom_point(size=0.5, color=adjustcolor('tomato4', alpha.f = 0.4)) +
  geom_smooth(method='lm', color='midnightblue') +
  stat_cor(r.digits = 2, method='spearman', cor.coef.name='rho') +
  xlab('Quantile normalized CoV change') +
  ylab('VST normalized data CoV change')

ggsave('./results/htseq/no_blind/CoV_ch.pdf', covchplot, units='cm', height = 10, width = 16, useDingbats=F)
ggsave('./results/htseq/no_blind/CoV_ch.png', covchplot, units='cm', height = 10, width = 16)

#################### 
####################  Test enrichment of reversal genes among DC genes with fisher test
#################### 
####################

ddc_genes = readRDS('./data/htseq/no_blind/dev_divergent_genes_dc_rank.rds')

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
OR1
# OR     pval
# <dbl>    <dbl>
#   1 0.527 1.45e-10
# shared reversal genes (updown+downup) are depleted among DiCo genes

#################### reversal genes (for each tissue) vs dc genes (among development divergent genes)
revgenes = readRDS('./data/htseq/no_blind/revgenes.tissues_tidy.rds')

OR2 = data.frame(gene_id=rep(names(ddc_genes),4), DCpattern = rep(ddc_genes,4) ,
                 tissue= rep(c('Cortex','Lung', 'Liver', 'Muscle'),
                             each = length(ddc_genes)), row.names=NULL) %>%
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
# reversal genes in each tissue (updown+downup) are enrichmed among DiCo genes

####################
#################### reversal genes (for each tissue) vs dc genes (among development divergent genes)
#################### test updown and downup genes separately
####################

# actx =  revgenes[revgenes$tissue=='Cortex',]
# actx = actx[actx$gene_id%in%names(ddc_genes),]
# ifelse(actx =='UpDown'| actx =='DownUp', T, F)
# a = data.frame(DC = ddc_genes[actx$gene_id]<0, 
#                Rev = ifelse(actx$direction =='UpDown', T, F) )
# table(a)
# fisher.test(table(a))


OR_UD = data.frame(gene_id=rep(names(ddc_genes),4),
                   DCpattern = rep(ddc_genes,4) ,
                   tissue= rep(c('Cortex','Lung', 'Liver', 'Muscle'), each = length(ddc_genes)),
                   row.names=NULL) %>%
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

####################
OR_DU = data.frame(gene_id=rep(names(ddc_genes),4), DCpattern = rep(ddc_genes,4) ,
                   tissue= rep(c('Cortex','Lung', 'Liver', 'Muscle'),each = length(ddc_genes)), 
                   row.names=NULL) %>%
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
# 1 Cortex  1.19 2.92e- 4
# 2 Liver   1.73 2.64e-40
# 3 Lung    1.15 1.78e- 3
# 4 Muscle  1.37 3.95e-14

# DownUp genes in each tissue are enriched among DiCo genes.

Rev_DiCo_enrichments = list(shdRev = OR1, Revs = OR2, UD = OR_UD, DU = OR_DU)

saveRDS(Rev_DiCo_enrichments, './data/htseq/no_blind/rev_dico_enrichments_OR.rds')

save(list=ls(), 
     file='./data/htseq/no_blind/6.DiCo.rdata')


