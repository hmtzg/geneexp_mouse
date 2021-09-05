library(tidyverse)
library(openxlsx)
source('./scripts/functions.R')

exp = readRDS("./data/htseq/expr_mat_no_blind.rds")
age = readRDS("./data/preprocess/ages.rds")
samp_tissue = readRDS("./data/preprocess/tissue_ids.rds")
sample_info = readRDS('./data/processed/tidy/sample_info.rds')

dev = sapply(unique(samp_tissue), function(b){
  ts = t(apply(exp, 1, function(c){
    chs = age < 90 & samp_tissue == b
    k = cor.test(c[chs], age[chs], m = "s")
    c(k$est, k$p.val)
  }))
  ts = cbind(ts, p.adjust(ts[,2], m = "BH"))
  colnames(ts) = c("rho", "p", "BH")
  return(ts)
},simplify = F)

saveRDS(dev, file="./data/htseq/no_blind/development_expression_change.rds")

#dev2 = lapply(dev, function(x) x[complete.cases(x),])

aging = sapply(unique(samp_tissue), function(b){
  ts = t(apply(exp, 1, function(c){
    chs = age >= 93 & samp_tissue == b
    k = cor.test(c[chs], age[chs], m = "s")
    c(k$est, k$p.val)
  }))
  ts = cbind(ts, p.adjust(ts[,2], m = "BH"))
  colnames(ts) = c("rho", "p", "BH")
  return(ts)
},simplify = F)

# saveRDS(aging, file="./data/processed/raw/ageing_expression_change.rds")
saveRDS(aging, file="./data/htseq/no_blind/ageing_expression_change.rds")

# aging = lapply(aging, function(x) x[complete.cases(x),])

expch_dev = reshape2::melt(dev) %>%
  spread(key= `Var2`, value = `value`) %>%
  mutate(period = 'development')
#expch_dev = expch_dev[complete.cases(expch_dev),]

expch_age = reshape2::melt(aging) %>%
  spread(key= `Var2`, value = `value`) %>%
  mutate(period = 'aging')

expch = rbind(expch_dev,expch_age) %>%
  set_names(c('gene_id','tissue','Expression Change','p','FDR','period'))

# saveRDS(expch,'./data/processed/tidy/expression_change.rds')
saveRDS(expch,'./data/htseq/no_blind/expression_change.rds')

age_related_genes = expch %>%
  set_names(c('gene_id','tissue','Expression Change (rho)','p','FDR','period')) %>%
  filter( FDR < 0.1) %>%
  arrange(tissue, group_by=period)

####################
####################
####################
####################

expch_qn = readRDS('./data/processed/tidy/expression_change.rds')

colnames(expch_qn)[3] = 'Exp. Ch. (QN)'

expch = expch[complete.cases(expch),]

expch %>%
  select(-p, -FDR) %>% 
  inner_join(select(expch_qn, -p, -FDR) ) %>% 
  group_by(tissue, period) %>% 
  summarise(cor = cor(`Expression Change`, `Exp. Ch. (QN)`))
  
# tissue period        cor
# <chr>  <chr>       <dbl>
#   1 Cortex aging       0.825
# 2 Cortex development 0.937
# 3 Liver  aging       0.870
# 4 Liver  development 0.810
# 5 Lung   aging       0.926
# 6 Lung   development 0.897
# 7 Muscle aging       0.845
# 8 Muscle development 0.811

expr = readRDS('./data/htseq/expr_no_blind.rds')
expr_qn = readRDS('./data/processed/tidy/expression.rds')
colnames(expr_qn)[3] = 'expression (qn)'

cdat = expr %>% 
  inner_join(expr_qn, by=c('gene_id','sample_id')) %>% 
  left_join(select(sample_info, -ind_id, -log2age), by='sample_id' ) 

## not finished

#################### buna gerek yok (simdilik)
#################### GORA for each tissue :
#################### up (or down) significant genes against all genes in that period and tissue
####################

# tissue_up_go = sapply(unique(expch$tissue), function(i){
#   sapply(unique(expch$period), function(pr){
#     a = expch %>%
#       filter(period == pr & tissue == i) %>%
#       mutate(genelist = ifelse(FDR < 0.1 & `Expression Change` > 0, 1, 0)) %>%
#       pull(genelist, name=gene_id)
#     if(sum(a)>=10) return(go_bp_enrich.test.Mm(genelist = a, selection = 1))
#   },simplify = F)
# },simplify = F)
# 
# tissue_down_go = sapply(unique(expch$tissue), function(i){
#   sapply(unique(expch$period), function(pr){
#     a = expch %>%
#       filter(period == pr & tissue == i) %>%
#       mutate(genelist = ifelse(FDR < 0.1 & `Expression Change` < 0, 1, 0)) %>%
#       pull(genelist, name=gene_id)
#     if(sum(a)>=10) return(go_bp_enrich.test.Mm(genelist = a, selection = 1))
#   },simplify = F)
# },simplify = F)
# 
# #save(tissue_up_go, tissue_down_go, file='./data/tissue_go.rdata')
# 
# tissue_up_go = unlist(tissue_up_go, recursive = F)
# tissue_down_go = unlist(tissue_down_go, recursive = F)
# 
# table_s1 = c(list('age_related_change' = age_related_genes), Upgenes = tissue_up_go, Downgenes = tissue_down_go)

####################
####################
####################
####################

