library(tidyverse)
library(ggpubr)
source('./scripts/functions.R')
attr = readRDS('./data/tfbs/attr.rds')
colnames(attr)[c(3,4)] = c('tf_sym', 'tf_id')
ddgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
ts.spec = readRDS('./data/processed/raw/ts.specQ3.genes.rds')

ddgenes = reshape2::melt(ddgenes) %>% 
  rownames_to_column(var='gene_id') %>%
  rename('DiCoRank' = value) %>%
  mutate(pattern = ifelse(DiCoRank<0, 'DiCo', 'DiDi')) %>% tibble()

attrM = attr %>% 
  select(-GeneSym, -GeneID, -weight, -ENS_hs) %>%
  rename('gene_id' = ENS_mm)

### Cortex
bg.ctx = attrM %>%
  filter(gene_id%in%ts.spec$Cortex) %>%
  inner_join(ddgenes)
bg.ctx

length(unique(bg.ctx$tf_sym)) # 15 TF to be tested.
length(unique(bg.ctx$gene_id)) # 591 DiCo+DiDi genes to be tested.
table(sapply(unique(bg.ctx$tf_sym), function(x){ sum(is.na(bg.ctx[bg.ctx$tf_sym%in%x,'pattern'])) })) 
# background check.

## number of TF targets for each gene
bg.ctx %>%  
  group_by(gene_id) %>%
  summarise(n=unique(length(tf_sym))) %>%
  select(n) %>% table()

# number of genes each TF targets
bg.ctx %>%  
  group_by(tf_sym) %>%
  summarise(n=unique(length(gene_id))) %>%
  select(n) %>% table() 
# 6 TF targets only 1 genes
# 7 TF targets only 2 gene
# 4 TF targets only 3 gene
# 2 TFs target only 4 genes

# cutoff for at least 5 genes targeted by each TF
bgX.ctx = bg.ctx %>%  
  group_by(tf_sym) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>4) %>%
  left_join(bg.ctx) %>%
  mutate(pattern = factor(pattern))

mirTest = function(mat, miRcol = 1, genecol=3, patcol=5){
  colnames(mat)[c(miRcol, genecol, patcol)] = c('tfbs', 'gene_id', 'pattern')
  # data frame with columns;
  # tfbs: tfbs names
  # gene_id: tfbs target gene names
  # pattern: factor with two levels to be tested, first level will be used in fisher test, OR>1
  
  sapply(unique(mat$tfbs), function(x){
    mirX = mat %>% 
      filter(tfbs==x) %>%
      #count(pattern, .drop=F)
      group_by(pattern) %>%
      summarise(n = length(unique(gene_id)))
    
    mirOther = mat %>% 
      filter(!tfbs==x) %>%
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

tfMm.ctx = mirTest(bgX.ctx, miRcol = 1, genecol = 4, patcol = 6)
tfMm2.ctx = t(sapply(tfMm.ctx, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(tfMm2.ctx) = c('OR', 'pvalue')
tfMm2.ctx = cbind(tfMm2.ctx, 'BY' = p.adjust(tfMm2.ctx[,'pvalue'], method='BY'))
tfMm2.ctx = as.data.frame(tfMm2.ctx)

tfMm2.ctx %>% arrange(OR)
tfMm2.ctx %>% filter(BY < 0.1)  # non significant


### Lung
bg.lng = attrM %>%
  filter(gene_id%in%ts.spec$Lung) %>%
  inner_join(ddgenes)
bg.lng
length(unique(bg.lng$tf_sym)) # 151 mirna to be tested.
length(unique(bg.lng$gene_id)) # 456 DiCo+DiDi genes to be tested.
table(sapply(unique(bg.lng$tf_sym), function(x){ sum(is.na(bg.lng[bg.lng$tf_sym%in%x,'pattern'])) })) 
# background check.
# number of miRNA targets for each gene
bg.lng %>%  
  group_by(gene_id) %>%
  summarise(n=unique(length(tf_sym))) %>%
  select(n) %>% table()
# number of genes each miRNA targets
bg.lng %>%  
  group_by(tf_sym) %>%
  summarise(n=unique(length(gene_id))) %>%
  select(n) %>% table() 
# 4 tfs target only one gene
# 9 tfs target two gene
# 3 tfs target  three gene
# cutoff for at least 5 genes targeted by each miRNA
bgX.lng = bg.lng %>%
  group_by(tf_sym) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>4) %>%
  left_join(bg.lng) %>%
  mutate(pattern = factor(pattern))

tfMm.lng = mirTest(bgX.lng, miRcol = 1, genecol = 4, patcol = 6)
tfMm2.lng = t(sapply(tfMm.lng, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(tfMm2.lng) = c('OR', 'pvalue')
tfMm2.lng = cbind(tfMm2.lng, 'BY' = p.adjust(tfMm2.lng[,'pvalue'], method='BY'))
tfMm2.lng = as.data.frame(tfMm2.lng)

tfMm2.lng %>% arrange(OR)
tfMm2.lng %>% filter(BY < 0.1)  # non significant

### Liver
bg.lvr = attrM %>%
  filter(gene_id%in%ts.spec$Liver) %>%
  inner_join(ddgenes)
bg.lvr
length(unique(bg.lvr$tf_sym)) # 146 mirna to be tested.
length(unique(bg.lvr$gene_id)) # 302 DiCo+DiDi genes to be tested.
table(sapply(unique(bg.lvr$tf_sym), function(x){ sum(is.na(bg.lvr[bg.lvr$tf_sym%in%x,'pattern'])) })) 
# background check.
# number of miRNA targets for each gene
bg.lvr %>%  
  group_by(gene_id) %>%
  summarise(n=unique(length(tf_sym))) %>%
  select(n) %>% table()
# number of genes each miRNA targets
bg.lvr %>%  
  group_by(tf_sym) %>%
  summarise(n=unique(length(gene_id))) %>%
  select(n) %>% table() 
# 4 tfs target only one gene
# 7 tfs target two gene
# 4 tfs target  three gene
# cutoff for at least 5 genes targeted by each miRNA
bgX.lvr = bg.lvr %>%
  group_by(tf_sym) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>4) %>%
  left_join(bg.lvr) %>%
  mutate(pattern = factor(pattern))

tfMm.lvr = mirTest(bgX.lvr, miRcol = 1, genecol = 4, patcol = 6)
tfMm2.lvr = t(sapply(tfMm.lvr, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(tfMm2.lvr) = c('OR', 'pvalue')
tfMm2.lvr = cbind(tfMm2.lvr, 'BY' = p.adjust(tfMm2.lvr[,'pvalue'], method='BY'))
tfMm2.lvr = as.data.frame(tfMm2.lvr)

tfMm2.lvr %>% arrange(OR)
tfMm2.lvr %>% filter(BY < 0.1)  # non significant

### Muscle
bg.ms = attrM %>%
  filter(gene_id%in%ts.spec$Muscle) %>%
  inner_join(ddgenes)
bg.ms
length(unique(bg.ms$tf_sym)) # 150 mirna to be tested.
length(unique(bg.ms$gene_id)) # 413 DiCo+DiDi genes to be tested.
table(sapply(unique(bg.ms$tf_sym), function(x){ sum(is.na(bg.ms[bg.ms$tf_sym%in%x,'pattern'])) })) 
# background check.
# number of miRNA targets for each gene
bg.ms %>%  
  group_by(gene_id) %>%
  summarise(n=unique(length(tf_sym))) %>%
  select(n) %>% table()
# number of genes each miRNA targets
bg.ms %>%  
  group_by(tf_sym) %>%
  summarise(n=unique(length(gene_id))) %>%
  select(n) %>% table() 
# 8 tfs target only one gene
# 7 tfs target two gene
# 4 tfs target  three gene
# cutoff for at least 5 genes targeted by each miRNA
bgX.ms = bg.ms %>%
  group_by(tf_sym) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>4) %>%
  left_join(bg.ms) %>%
  mutate(pattern = factor(pattern))

tfMm.ms = mirTest(bgX.ms, miRcol = 1, genecol = 4, patcol = 6)
tfMm2.ms = t(sapply(tfMm.ms, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(tfMm2.ms) = c('OR', 'pvalue')
tfMm2.ms = cbind(tfMm2.ms, 'BY' = p.adjust(tfMm2.ms[,'pvalue'], method='BY'))
tfMm2.ms = as.data.frame(tfMm2.ms)

tfMm2.ms %>% arrange(OR)
tfMm2.ms %>% filter(BY < 0.1)  # non significant

#####################
#### number of TFs  for DiCo vs DiDi genes
bgN = attrM %>%
  right_join(ddgenes)

bgN %>%
  group_by(pattern) %>%
  summarise(n= length(unique(tf_sym)))
# pattern     n
# <fct>   <int>
# 1 DiCo      158
# 2 DiDi      158
bgsum = bgN %>%
  group_by(pattern, gene_id) %>%
  summarise(n= ifelse(is.na(tf_sym), 0, length(tf_sym)) ) %>% unique()

tfplot = bgsum %>%
  mutate(n = n+1) %>%
  ggplot(aes(x=pattern, y=n)) +
  scale_y_continuous(trans='log2') +
  geom_boxplot(outlier.color = 'gray30',outlier.size = 0.4) +
  #geom_violin() +
  xlab('') +
  ylab('Number of TFs (log2 scale)')
tfplot
ggsave('./results/tfbs/tf_N.pdf', tfplot, units='cm', height = 10, width = 8, useDingbats=F)
ggsave('./results/tfbs/tf_N.png', tfplot, units='cm', height = 10, width = 8)

# tfcorplot = bgsum %>% 
#   left_join(ddgenes) %>%
#   ggplot(aes(y=DiCoRank, x=n, color=pattern)) +
#   geom_point(size=0.6, alpha=0.4, show.legend = F) +
#   geom_smooth(aes(color=pattern), method='lm') +
#   geom_hline(yintercept=0, linetype=2, col='gray40') +
#   xlab('Number of miRNA targets') +
#   stat_cor(method='spearman', cor.coef.name = 'rho')
# tfcorplot
# ggsave('./results/tfbs/tf_corplot.pdf',tfcorplot, units='cm', height = 10, width = 10, useDingbats=F)
# ggsave('./results/tfbs/tf_corplot.png',tfcorplot, units='cm', height = 10, width = 10)

# permutation test:
xDiCo = bgsum$n[bgsum$pattern=='DiCo']
nDiCo = length(xDiCo)
xDiDi = bgsum$n[bgsum$pattern=='DiDi']
nDiDi = length(xDiDi)
pop = c(xDiCo, xDiDi)
obs = c('xDiCo' = mean(xDiCo), 'xDiDi' = mean(xDiDi) )
obstat = cohens_d(xDiCo, xDiDi)
#obstat = mean(xDiCo)-mean(xDiDi)
perms = sapply(1:1000, function(x){
  px = sample(pop)
  pDiCo = px[1:nDiCo]
  pDiDi = px[(nDiCo+1): length(pop)]
  pstat = cohens_d(pDiCo, pDiDi)
  #pstat = mean(pDiCo) - mean(pDiDi)
  return(pstat)
})
hist(perms)
mean(perms >= obstat)
# p = 0.426
mean(perms) / obstat
# fpr= -0.021

save(list=ls(),
     file='./results/tfbs/analysis.rdata')
####

# sinfo = readRDS('./data/processed/tidy/sample_info.rds') %>% select(tissue, age, sample_id)
# expr = readRDS('./data/processed/tidy/expression.rds')
# expr =  expr %>% left_join(sinfo)
# 
# mirplot = expr %>%
#   filter(gene_id%in%dicomir) %>% 
#   ggplot(aes(x=age, y=expression, color=tissue)) +
#   facet_wrap(gene_id~.,scale='free_y', ncol=4) +
#   geom_point(size=0.3, aes(group=tissue)) +
#   geom_smooth(method='loess', se=F, aes(group=tissue) ) +
#   scale_x_continuous(trans='log2') +
#   scale_color_manual(values=tissuecol) +
#   geom_vline(xintercept= 90, linetype='dashed', col='gray30')
