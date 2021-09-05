library(tidyverse)
library(ggpubr)
source('./scripts/functions.R')
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
attr = readRDS('./data/mirna/attr.rds')
ddgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
covch = readRDS('./data/processed/tidy/CoV_change.rds')
ts.spec = readRDS('./data/processed/raw/ts.specQ3.genes.rds')
ddgenes = reshape2::melt(ddgenes) %>% 
  rownames_to_column(var='gene_id') %>%
  rename('DiCoRank' = value) %>%
  mutate(pattern = ifelse(DiCoRank<0, 'DiCo', 'DiDi')) %>% tibble()

attrM = attr %>% 
  select(-GeneSym, -GeneID, -`NA`, -weight, -ENS_hs) %>%
  rename('gene_id' = ENS_mm)

### Cortex

bg.ctx = attrM %>%
  filter(gene_id%in%ts.spec$Cortex) %>%
  inner_join(ddgenes)
bg.ctx

length(unique(bg.ctx$miRNA)) # 172 mirna to be tested.
length(unique(bg.ctx$gene_id)) # 343 DiCo+DiDi genes to be tested.
table(sapply(unique(bg.ctx$miRNA), function(x){ sum(is.na(bg.ctx[bg.ctx$miRNA%in%x,'pattern'])) })) 
# background check.

# number of miRNA targets for each gene
bg.ctx %>%  
  group_by(gene_id) %>%
  summarise(n=unique(length(miRNA))) %>%
  select(n) %>% table()

# number of genes each miRNA targets
bg.ctx %>%  
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  select(n) %>% table() 
# 69 miRNAs target only one gene
# 24 miRNAs target two gene
# 17 miRNAs target  three gene

# cutoff for at least 5 genes targeted by each miRNA
bgX.ctx = bg.ctx %>%
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>4) %>%
  left_join(bg.ctx) %>%
  mutate(pattern = factor(pattern))

mirTest = function(mat, miRcol = 1, genecol=3, patcol=5){
  colnames(mat)[c(miRcol, genecol, patcol)] = c('miRNA', 'gene_id', 'pattern')
  # data frame with columns;
  # miRNA: miRNA names
  # gene_id: miRNA target gene names
  # pattern: factor with two levels to be tested, first level will be used in fisher test, OR>1
  
  sapply(unique(mat$miRNA), function(x){
    mirX = mat %>% 
      filter(miRNA==x) %>%
      #count(pattern, .drop=F)
      group_by(pattern) %>%
      summarise(n = length(unique(gene_id)))
    
    mirOther = mat %>% 
      filter(!miRNA==x) %>%
      group_by(pattern) %>%
      summarise(n = length(unique(gene_id)))
    
    fisherX = mirX %>% full_join(mirOther, by='pattern') %>%
      rename( 'mirX'= n.x, 'mirOthers' = n.y) %>%
      arrange(pattern) %>% 
      column_to_rownames('pattern') %>% 
      t() 
    if(sum(is.na(fisherX))>0) fisherX = fisherX %>% replace_na(replace = 0)
    
    fisherXtest = fisher.test(fisherX)
    list(table= fisherX, Fistest = fisherXtest)
  }, USE.NAMES=T, simplify=F)
}

miRtestMm.ctx = mirTest(bgX.ctx, miRcol = 1, genecol = 4, patcol = 6)
miRtestMmresult.ctx = t(sapply(miRtestMm.ctx, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(miRtestMmresult.ctx) = c('OR', 'pvalue')
miRtestMmresult.ctx = cbind(miRtestMmresult.ctx, 'BY' = p.adjust(miRtestMmresult.ctx[,'pvalue'], method='BY'))
miRtestMmresult.ctx = as.data.frame(miRtestMmresult.ctx)

miRtestMmresult.ctx %>% arrange(OR)
miRtestMmresult.ctx %>% filter(BY < 0.1)  # non significant

### Lung
bg.lng = attrM %>%
  filter(gene_id%in%ts.spec$Lung) %>%
  inner_join(ddgenes)
bg.lng
length(unique(bg.lng$miRNA)) # 214 mirna to be tested.
length(unique(bg.lng$gene_id)) # 327 DiCo+DiDi genes to be tested.
table(sapply(unique(bg.lng$miRNA), function(x){ sum(is.na(bg.lng[bg.lng$miRNA%in%x,'pattern'])) })) 
# background check.
# number of miRNA targets for each gene
bg.lng %>%  
  group_by(gene_id) %>%
  summarise(n=unique(length(miRNA))) %>%
  select(n) %>% table()
# number of genes each miRNA targets
bg.lng %>%  
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  select(n) %>% table() 
# 86 miRNAs target only one gene
# 36 miRNAs target two gene
# 17 miRNAs target  three gene
# cutoff for at least 5 genes targeted by each miRNA
bgX.lng = bg.lng %>%
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>4) %>%
  left_join(bg.lng) %>%
  mutate(pattern = factor(pattern))

miRtestMm.lng = mirTest(bgX.lng, miRcol = 1, genecol = 4, patcol = 6)
miRtestMmresult.lng = t(sapply(miRtestMm.lng, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(miRtestMmresult.lng) = c('OR', 'pvalue')
miRtestMmresult.lng = cbind(miRtestMmresult.lng, 'BY' = p.adjust(miRtestMmresult.lng[,'pvalue'], method='BY'))
miRtestMmresult.lng = as.data.frame(miRtestMmresult.lng)

miRtestMmresult.lng %>% arrange(OR)
miRtestMmresult.lng %>% filter(BY < 0.1)  # non significant

### Liver
bg.lvr = attrM %>%
  filter(gene_id%in%ts.spec$Liver) %>%
  inner_join(ddgenes)
bg.lvr
length(unique(bg.lvr$miRNA)) # 146 mirna to be tested.
length(unique(bg.lvr$gene_id)) # 216 DiCo+DiDi genes to be tested.
table(sapply(unique(bg.lvr$miRNA), function(x){ sum(is.na(bg.lvr[bg.lvr$miRNA%in%x,'pattern'])) })) 
# background check.
# number of miRNA targets for each gene
bg.lvr %>%  
  group_by(gene_id) %>%
  summarise(n=unique(length(miRNA))) %>%
  select(n) %>% table()
# number of genes each miRNA targets
bg.lvr %>%  
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  select(n) %>% table() 
# 68 miRNAs target only one gene
# 25 miRNAs target two gene
# 11 miRNAs target  three gene
# cutoff for at least 5 genes targeted by each miRNA
bgX.lvr = bg.lvr %>%
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>4) %>%
  left_join(bg.lvr) %>%
  mutate(pattern = factor(pattern))

miRtestMm.lvr = mirTest(bgX.lvr, miRcol = 1, genecol = 4, patcol = 6)
miRtestMmresult.lvr = t(sapply(miRtestMm.lvr, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(miRtestMmresult.lvr) = c('OR', 'pvalue')
miRtestMmresult.lvr = cbind(miRtestMmresult.lvr, 'BY' = p.adjust(miRtestMmresult.lvr[,'pvalue'], method='BY'))
miRtestMmresult.lvr = as.data.frame(miRtestMmresult.lvr)

miRtestMmresult.lvr %>% arrange(OR)
miRtestMmresult.lvr %>% filter(BY < 0.1)  # non significant
#######

### Muscle
bg.ms = attrM %>%
  filter(gene_id%in%ts.spec$Muscle) %>%
  inner_join(ddgenes)
bg.ms
length(unique(bg.ms$miRNA)) # 150 mirna to be tested.
length(unique(bg.ms$gene_id)) # 284 DiCo+DiDi genes to be tested.
table(sapply(unique(bg.ms$miRNA), function(x){ sum(is.na(bg.ms[bg.ms$miRNA%in%x,'pattern'])) })) 
# background check.
# number of miRNA targets for each gene
bg.ms %>%  
  group_by(gene_id) %>%
  summarise(n=unique(length(miRNA))) %>%
  select(n) %>% table()
# number of genes each miRNA targets
bg.ms %>%  
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  select(n) %>% table() 
# 52 miRNAs target only one gene
# 25 miRNAs target two gene
# 14 miRNAs target  three gene
# cutoff for at least 5 genes targeted by each miRNA
bgX.ms = bg.ms %>%
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>4) %>%
  left_join(bg.ms) %>%
  mutate(pattern = factor(pattern))

miRtestMm.ms = mirTest(bgX.ms, miRcol = 1, genecol = 4, patcol = 6)
miRtestMmresult.ms = t(sapply(miRtestMm.ms, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(miRtestMmresult.ms) = c('OR', 'pvalue')
miRtestMmresult.ms = cbind(miRtestMmresult.ms, 'BY' = p.adjust(miRtestMmresult.ms[,'pvalue'], method='BY'))
miRtestMmresult.ms = as.data.frame(miRtestMmresult.ms)

miRtestMmresult.ms %>% arrange(OR)
miRtestMmresult.ms %>% filter(BY < 0.1)  # non significant

