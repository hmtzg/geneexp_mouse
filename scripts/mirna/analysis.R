library(tidyverse)
library(ggpubr)
source('./scripts/functions.R')
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
attr = readRDS('./data/mirna/attr.rds')
ddgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')

ddgenes = reshape2::melt(ddgenes) %>% 
  rownames_to_column(var='gene_id') %>%
  rename('DiCoRank' = value) %>%
  mutate(pattern = ifelse(DiCoRank<0, 'DiCo', 'DiDi')) %>% tibble()

attrM = attr %>% 
  select(-GeneSym, -GeneID, -`NA`, -weight, -ENS_hs) %>%
  rename('gene_id' = ENS_mm)

###
bg = attrM %>%
  inner_join(ddgenes)

length(unique(bg$miRNA)) # 513 mirna to be tested.
length(unique(bg$gene_id)) # 5458 DiCo+DiDi genes to be tested.
table(sapply(unique(bg$miRNA), function(x){ sum(is.na(bg[bg$miRNA%in%x,'pattern'])) })) # background check.

# number of miRNA targets for each gene
bg %>%  
  group_by(gene_id) %>%
  summarise(n=unique(length(miRNA))) %>%
  select(n) %>% table()

# number of genes each miRNA targets
bg %>%  
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  select(n) %>% table() 
# 156 miRNAs target only one gene
# 69 miRNAs target only two gene
# 39 miRNAs target only three gene

# cutoff for at least 5 genes targeted by each miRNA
bgX = bg %>%  
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>4 & n<501) %>%
  left_join(bg) %>%
  mutate(pattern = factor(pattern))

# bgX = bg %>%  
#   group_by(gene_id) %>%
#   summarise(n=unique(length(miRNA))) %>%
#   filter(n>9) %>% select(gene_id) %>%
#   inner_join(bgX, by='gene_id')

# bg = bg %>% 
#   mutate(pattern = factor(pattern))

mirTest = function(mat, miRcol = 1, genecol=3, patcol=5){
  colnames(mat)[c(miRcol, genecol, patcol)] = c('miRNA', 'gene_id', 'pattern')
  # data frame with columns;
  # miRNA: miRNA names
  # gene_id: miRNA target gene names
  # pattern: factor with two levels to be tested, first level will be used in fisher test, OR>1
  
  sapply(unique(mat$miRNA), function(x){
    print(x)
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

#miRtestMm = mirTest(bg,miRcol = 1, genecol = 3,patcol = 5)
miRtestMm = mirTest(bgX, miRcol = 1, genecol = 4, patcol = 6)
miRtestMmresult = t(sapply(miRtestMm, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(miRtestMmresult) = c('OR', 'pvalue')
miRtestMmresult = cbind(miRtestMmresult, 'BY' = p.adjust(miRtestMmresult[,'pvalue'], method='BY') )
miRtestMmresult = as.data.frame(miRtestMmresult)

miRtestMmresult %>% arrange(OR)
miRtestMmresult %>% filter(BY < 0.1)  # non significant

#####################
#####################
gtexcov = readRDS('./results/GTEx/covcor.rds') %>%
  mutate(pattern = ifelse(rho<0,'Co', 'Di')) %>%
  select(-p, -adjusted_p)

attrG = attr %>% 
  select(-GeneSym, -GeneID, -`NA`, -weight, -ENS_mm) %>%
  rename('GeneID' = ENS_hs)
bg2 = attrG %>%
  inner_join(gtexcov)

# cutoff for at least 5 genes targeted by each miRNA
bgX2 = bg2 %>%  
  group_by(miRNA) %>%
  summarise(n=unique(length(GeneID))) %>%
  filter(n>4) %>%
  left_join(bg2) %>%
  mutate(pattern = factor(pattern))

miRtestGtex = mirTest(bgX2, miRcol = 1, genecol = 4, patcol = 6)
miRtestG2 = t(sapply(miRtestGtex, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(miRtestG2) = c('OR', 'pvalue')
miRtestG2 = cbind(miRtestG2, 'BY' = p.adjust(miRtestG2[,'pvalue'], method='BY'))
miRtestG2 = as.data.frame(miRtestG2)

miRtestG2 %>% filter(BY < 0.1)  # 0 significant
miRtestG2 %>% filter(pvalue < 0.1 & OR>1)

#####################
#####################
jonkercov = readRDS('./data/other_datasets/jonker/processed/covch.rds') %>% as.data.frame() %>%
  mutate(pattern = ifelse(rho<0, 'Co', 'Di')) %>% 
  select(-pval, -`BH`) 

bg3 = attrM %>%
  inner_join(jonkercov)

bgX3 = bg3 %>%  
  group_by(miRNA) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>4) %>%
  left_join(bg3) %>%
  mutate(pattern = factor(pattern))

miRtestJ = mirTest(bgX3, miRcol = 1, genecol = 4, patcol = 6)
miRtestJ2 = t(sapply(miRtestJ, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(miRtestJ2) = c('OR', 'pvalue')
miRtestJ2 = cbind(miRtestJ2, 'BY' = p.adjust(miRtestJ2[,'pvalue'], method='BY'))
miRtestJ2 = as.data.frame(miRtestJ2)

miRtestJ2 %>% filter(BY < 0.1)  # 1 significant
#                  OR       pvalue         BY
# hsa-miR-193b-3p 0.690087 1.401889e-05 0.02348303

miRtestJ2 %>% filter(pvalue < 0.1 & OR >1) # 

####
sinfo = readRDS('./data/processed/tidy/sample_info.rds') %>% select(tissue, age, sample_id)
expr = readRDS('./data/processed/tidy/expression.rds')
expr =  expr %>% left_join(sinfo)

dicomir = bg %>% filter(miRNA=='hsa-miR-193b-3p') %>% pull(gene_id) %>% unique()

mirplot = expr %>%
  filter(gene_id%in%dicomir) %>%
  ggplot(aes(x=age, y=expression, color=tissue)) +
  facet_wrap(gene_id~.,scale='free_y', ncol=4) +
  geom_point(size=0.3, aes(group=tissue)) +
  geom_smooth(method='loess', se=F, aes(group=tissue) ) +
  scale_x_continuous(trans='log2') +
  scale_color_manual(values=tissuecol) +
  geom_vline(xintercept= 90, linetype='dashed', col='gray30')
mirplot
ggsave('./results/mirna/mirplot.pdf',mirplot, units='cm', height = 25, width = 16, useDingbats=F)
ggsave('./results/mirna/mirplot.png',mirplot, units='cm', height = 25, width = 16)

#### number of miRNA targets for DiCo vs DiDi genes
bgN = attrM %>%
  right_join(ddgenes)

bgN %>%
  group_by(pattern) %>%
  summarise(n= length(unique(miRNA)))
# pattern     n
# <fct>   <int>
# 1 DiCo      406
# 2 DiDi      415

bgsum = bgN %>%
  group_by(pattern, gene_id) %>%
  summarise(n= ifelse(is.na(miRNA), 0, length(miRNA)) ) %>% unique()

mirNplot = bgsum %>%
  mutate(n=n+1) %>%
  ggplot(aes(x=pattern, y=n)) +
  #geom_boxplot(outlier.color = 'gray30',outlier.size = 0.4) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_y_continuous(trans='log2') +
  xlab('') +
  ylab('Number of miRNAs (log2 scale)')
mirNplot
ggsave('./results/mirna/miRNA_N.pdf',mirNplot, units='cm', height = 10, width = 8, useDingbats=F)
ggsave('./results/mirna/miRNA_N.png',mirNplot, units='cm', height = 10, width = 8)

# mircorplot = bgsum %>% 
#   left_join(ddgenes) %>%
#   ggplot(aes(y=DiCoRank, x=n, color=pattern)) +
#   geom_point(size=0.6, alpha=0.4, show.legend = T) +
#   geom_smooth(aes(color=pattern), method='lm', show.legend = T) +
#   scale_x_continuous(trans='log2') +
#   geom_hline(yintercept=0, linetype=2, col='gray40') +
#   xlab('Number of miRNA targets (log2 scale)') +
#   stat_cor(method='spearman', cor.coef.name = 'rho', show.legend = F)
# mircorplot
# ggsave('./results/mirna/miRNA_corplot.pdf',mircorplot, units='cm', height = 10, width = 10, useDingbats=F)
# ggsave('./results/mirna/miRNA_corplot.png',mircorplot, units='cm', height = 10, width = 10)


# permutation test:
xDiCo = bgsum$n[bgsum$pattern=='DiCo']
nDiCo = length(xDiCo)
xDiDi = bgsum$n[bgsum$pattern=='DiDi']
nDiDi = length(xDiDi)
pop = c(xDiCo, xDiDi)
obs = c('xDiCo' = mean(xDiCo), 'xDiDi' = mean(xDiDi) )
obstat = cohens_d(xDiCo, xDiDi)
#obstat = mean(xDiCo) -mean(xDiDi)
perms = sapply(1:1000, function(x){
  px = sample(pop)
  pDiCo = px[1:nDiCo]
  pDiDi = px[(nDiCo+1): length(pop)]
  pstat = cohens_d(pDiCo, pDiDi)
  #pstat = mean(pDiCo) - mean(pDiDi)
               
  return(pstat)
})
hist(perms)
mean(perms <= obstat)
# p = 0.023
# DiDi genes have more miRNA regulation
mean(perms) / obstat
# fpr= - 0.002

save(list=ls(),
     file='./results/mirna/analysis.rdata')
