library(tidyverse)
library(openxlsx)
source('./scripts/functions.R')

exp = readRDS("./data/processed/raw/expression.rds")
age = readRDS("./data/preprocess/ages.rds")
samp_tissue = readRDS("./data/preprocess/tissue_ids.rds")

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

saveRDS(dev, file="./data/processed/raw/development_expression_change.rds")

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

saveRDS(aging, file="./data/processed/raw/ageing_expression_change.rds")

aging = lapply(aging, function(x) x[complete.cases(x),])

expch_dev = reshape2::melt(dev) %>%
  spread(key= `Var2`, value = `value`) %>%
  mutate(period = 'development')

expch_age = reshape2::melt(aging) %>%
  spread(key= `Var2`, value = `value`) %>%
  mutate(period = 'aging')

expch = rbind(expch_dev,expch_age) %>%
  set_names(c('gene_id','tissue','Expression Change','p','FDR','period'))

saveRDS(expch,'./data/processed/tidy/expression_change.rds')

age_related_genes = expch %>%
  set_names(c('gene_id','tissue','Expression Change (rho)','p','FDR','period')) %>%
  filter( FDR < 0.1) %>%
  arrange(tissue, group_by=period)


## number of sig genes in each tissue and period:
age_related_genes %>%
  group_by(period, tissue) %>%
  summarise(n = n(), npercent = n()/ 15063 * 100)
# period      tissue     n npercent
# <chr>       <chr>  <int>    <dbl>
# 1 aging       Cortex    68   0.451 
# 2 aging       Liver    156   1.04  
# 3 aging       Lung    2319  15.4   
# 4 aging       Muscle     2   0.0133
# 5 development Cortex  6151  40.8   
# 6 development Liver   4471  29.7   
# 7 development Lung    3396  22.5   
# 8 development Muscle  1941  12.9  

####################
####################
####################
#################### expr change with pearson:

expr = readRDS('data/processed/tidy/expression.rds')
sinfo = readRDS('data/processed/tidy/sample_info.rds')
agecor = expr %>%
  left_join(sinfo) %>%
  mutate(period = ifelse(age<90,'development','ageing')) %>%
  group_by(tissue,gene_id, period) %>%
  summarise(cor = cor.test(expression, age, m='pearson')$est,
            p.val = cor.test(expression, age, m='pearson')$p.val)
agecor
saveRDS(agecor, 'data/processed/tidy/expression_change_pearson.rds')
####################
#################### GORA for each tissue :
#################### up (or down) significant genes against all genes in that period and tissue
####################

# i = unique(expch$tissue)[4]
# pr = 'development'
tissue_up_go = sapply(unique(expch$tissue), function(i){
  sapply(unique(expch$period), function(pr){
    print(i)
    print(pr)
    a = expch %>%
      filter(period == pr & tissue == i) %>% 
      mutate(genelist = ifelse(FDR < 0.1 & `Expression Change` > 0, 1, 0)) %>%
      pull(genelist, name=gene_id)
    if(sum(a)>=10) return(go_bp_enrich.test.Mm(genelist = a, selection = 1,padj = 'BH'))
  },simplify = F)
},simplify = F)

tissue_down_go = sapply(unique(expch$tissue), function(i){
  sapply(unique(expch$period), function(pr){
    print(i)
    print(pr)
    a = expch %>%
      filter(period == pr & tissue == i) %>%
      mutate(genelist = ifelse(FDR < 0.1 & `Expression Change` < 0, 1, 0)) %>%
      pull(genelist, name=gene_id)
    if(sum(a)>=10) return(go_bp_enrich.test.Mm(genelist = a, selection = 1, padj = 'BH'))
  },simplify = F)
},simplify = F)

save(tissue_up_go, tissue_down_go, file='./data/processed/raw/tissue_go.rdata')

tissue_up_go = unlist(tissue_up_go, recursive = F)
tissue_down_go = unlist(tissue_down_go, recursive = F)

table_sX = c(list('age_related_change' = age_related_genes), Upgenes = tissue_up_go, 
             Downgenes = tissue_down_go)

write.xlsx(table_sX,
           file='./results/supplementary_files/Supplementary_File_1.xlsx', row.names=F)
