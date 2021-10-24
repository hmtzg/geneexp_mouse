library(tidyverse)
library(openxlsx)
library(ggpubr)

source('./scripts/functions.R')
exp = readRDS("./data/htseq/expr_mat_blinded.rds")
age = readRDS("./data/preprocess/ages.rds")
samp_tissue = readRDS("./data/preprocess/tissue_ids.rds")
sample_info = readRDS('./data/processed/tidy/sample_info.rds')

theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)

tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('Up','Down'))
revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))
periodcol = c(Development = "#FE6100", Ageing ="#648FFF")

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

saveRDS(dev, file="./data/htseq/blinded/development_expression_change.rds")

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

saveRDS(aging, file="./data/htseq/blinded/ageing_expression_change.rds")

# aging = lapply(aging, function(x) x[complete.cases(x),])

expch_dev = reshape2::melt(dev) %>%
  spread(key= `Var2`, value = `value`) %>%
  mutate(period = 'Development')
#expch_dev = expch_dev[complete.cases(expch_dev),]

expch_age = reshape2::melt(aging) %>%
  spread(key= `Var2`, value = `value`) %>%
  mutate(period = 'Ageing')

expch = rbind(expch_dev,expch_age) %>%
  set_names(c('gene_id','tissue','Expression Change','p','FDR','period'))

saveRDS(expch,'./data/htseq/blinded/expression_change.rds')

age_related_genes = expch %>%
  set_names(c('gene_id','tissue','Expression Change (rho)','p','FDR','period')) %>%
  filter( FDR < 0.1) %>%
  arrange(tissue, group_by=period)

####################
####################
####################
####################



