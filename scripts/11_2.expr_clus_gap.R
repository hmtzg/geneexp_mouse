library(tidyverse)
library(reshape2)
library(openxlsx)
library(ggpubr)
library(gridExtra)
library(grid)
library(cluster)
library(factoextra)

setwd('/home/hmt/Dropbox/projects/repos/geneexp_mouse/')
source('./scripts/functions.R')

#cov = readRDS('./data/processed/tidy/CoV.rds')
#covch = readRDS('./data/processed/tidy/CoV_change.rds')
expr = readRDS('./data/processed/tidy/expression.rds')
sample_info = readRDS('./data/processed/tidy/sample_info.rds')
#dgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
#dcgenes = names(which(dgenes<0))

############### scale expression values:
scexpr = expr %>%
  left_join(select(sample_info, -log2age) ) %>%
  group_by(gene_id, tissue) %>%
  mutate(scExpr = scale(expression)[,1]) %>% ungroup()
print('deneme')
############### get dc gene matrix

scdat = split(scexpr,f = scexpr$tissue)
scdat = lapply(scdat, function(x){
  x %>%
    acast(gene_id~sample_id, value.var = 'scExpr')
})

### gap statistics to find optimum K for kmeans:
k=0
gapstat = lapply(scdat, function(x){
  k<<- k+1
  print(paste0(names(scdat)[k],'...'))
  gapstatX = clusGap(x, FUNcluster = kmeans , K.max = 30, B = 100, d.power = 2, spaceH0 = 'original')
})
print('bitti')
save(list=ls(), file='/home/hmt/Dropbox/projects/repos/geneexp_mouse/results/expr_cluster_gap.rdata')
