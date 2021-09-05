##############
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

####################
####################
####################
####################

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

write.xlsx(table_sX, file='./results/SI_tables/TableS2.xlsx', row.names=F)



### GSEA for each tissue and each period separately:
# library(clusterProfiler)
# tissue_gsea = sapply(unique(expch$tissue), function(i){
#   sapply(unique(expch$period), function(pr){
#     a = expch %>%
#       filter(period == 'development' & tissue == 'Cortex') %>%
#       pull(rho, name = gene_id) %>%
#       sort(., decreasing = T) %>%
#       gseGO(geneList = ., OrgDb = "org.Mm.eg.db", ont = "BP", keyType = "ENSEMBL", nPerm = 1000, 
#              minGSSize = 10,
#             maxGSSize = 500, pvalueCutoff = 1, verbose = F, pAdjustMethod = 'BY')
#     return(a)
#   },simplify = F)
# },simplify = F)
# 
# library(openxlsx)
# table_s1_gsea = list('age_related_change' = age_related_genes,
#                 'Cortex_development' = cortex_dev_gsea, 'Cortex_ageing' = cortex_aging_gsea,
#                 'Lung_development' = lung_dev_gsea, 'Lung_ageing' = lung_aging_gsea,
#                 'Liver_development' = liver_dev_gsea, 'Liver_ageing' = liver_aging_gsea,
#                 'Muscle_development' = muscle_dev_gsea, 'Muscle_ageing' = muscle_aging_gsea)
# write.xlsx(table_s1_gsea, file='./data/Table_S1_gsea.xlsx', row.names=T)
