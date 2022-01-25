library(tidyverse)
library(ggpubr)
source('./scripts/functions.R')
attr = readRDS('./data/tfbs/attr.rds')
colnames(attr)[c(3,4)] = c('tf_sym', 'tf_id')
ddgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')

ddgenes = reshape2::melt(ddgenes) %>% 
  rownames_to_column(var='gene_id') %>%
  rename('DiCoRank' = value) %>%
  mutate(pattern = ifelse(DiCoRank<0, 'DiCo', 'DiDi')) %>% tibble()

attrM = attr %>% 
  select(-GeneSym, -GeneID, -weight, -ENS_hs) %>%
  rename('gene_id' = ENS_mm)
#########
#########
#########
attr %>% 
  summarise(tf = unique(tf_sym))
library(biomaRt)
martx = biomaRt::useMart(biomart = 'ensembl')
marths_= biomaRt::useDataset('hsapiens_gene_ensembl', mart=martx)
listAttributes(marths_)[grep('entrez', listAttributes(marths_)[,1]),]
genemap = getBM(attributes = c('entrezgene_id', 'ensembl_gene_id'), filters = c('entrezgene_id'),
                values = unique(attrM$tf_id), mart = marths_)
sum(duplicated(genemap$entrezgene_id))
dups = genemap$entrezgene_id[duplicated(genemap$entrezgene_id)]
genemap = genemap[!genemap$entrezgene_id%in%dups,]

martmm_ = biomaRt::useDataset('mmusculus_gene_ensembl', mart=martx)
ensmap = biomaRt::getLDS(attributes = c('ensembl_gene_id'), filters = 'ensembl_gene_id', 
                         values = unique(genemap$ensembl_gene_id), mart = marths_, 
                         attributesL = c('ensembl_gene_id'), martL = martmm_)
colnames(ensmap) = c('tf_hsid', 'tf_mmid')
sum(duplicated(ensmap$tf_hsid)) # 0
sum(duplicated(ensmap$tf_mmid)) # 0 

colnames(genemap) = c('tf_id', 'tf_hsid')

saveRDS(ensmap, file='./data/tfbs/analysis.rdata')
ensmap = ensmap %>% 
  left_join(genemap, by='tf_hsid')
head(ensmap)
ensmap$tf_id = as.character(ensmap$tf_id)

ensmap %>%
  rename('gene_id'=tf_mmid) %>%
  inner_join(ddgenes) %>%
  group_by(pattern) %>%
  summarise(n = n())
  
tfdico = ensmap %>%
  rename('gene_id'=tf_mmid) %>%
  inner_join(ddgenes) %>% 
  arrange(DiCoRank)

head(tfdico)
head(attrM)

attrM %>% 
  filter(tf_id%in%'6667') %>% 
  inner_join(ddgenes) %>%
  group_by(pattern) %>%
  summarise(n=n())

#########
#########
#########
###
bg = attrM %>%
  inner_join(ddgenes)

length(unique(bg$tf_sym)) # 158 tf to be tested.
length(unique(bg$gene_id)) # 7247 DiCo+DiDi genes to be tested.
table(sapply(unique(bg$tf_sym), function(x){ sum(is.na(bg[bg$tf_sym%in%x,'pattern'])) })) # background check.

## number of TF targets for each gene
bg %>%  
  group_by(gene_id) %>%
  summarise(n=unique(length(tf_sym))) %>%
  dplyr::select(n) %>% table() # 

# number of genes each TF targets
bg %>%  
  group_by(tf_sym) %>%
  summarise(n=unique(length(gene_id))) %>%
  dplyr::select(n) %>% table() 
# 1 TF targets only 2 genes
# 1 TF targets only 4 gene
# 1 TF targets only 10 gene
# 2 TFs target only 20 genes

# cutoff for at least 5 genes targeted by each TF
bgX = bg %>%  
  group_by(tf_sym) %>%
  summarise(n=unique(length(gene_id))) %>%
  filter(n>=5) %>%
  left_join(bg) %>%
  mutate(pattern = factor(pattern))

# 
# bg = bg %>% 
#   mutate(pattern = factor(pattern))

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

tfMm = mirTest(bgX, miRcol = 1, genecol = 4, patcol = 6)
tfMm2 = t(sapply(tfMm, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(tfMm2) = c('OR', 'pvalue')
tfMm2 = cbind(tfMm2, 'BH' = p.adjust(tfMm2[,'pvalue'], method='BH'))
tfMm2 = as.data.frame(tfMm2)

tfMm2 %>% arrange(OR) # 158
tfMm2 %>% filter(BH < 0.1)  # non significant

saveRDS(tfMm2,'./results/source_data/f4/tf.rds')
#####################
#####################
gtexcov = readRDS('./results/GTEx/covcor.rds') %>%
  mutate(pattern = ifelse(rho<0,'Co', 'Di')) %>%
  select(-p, -adjusted_p)
attrG = attr %>% 
  select(-GeneSym, -GeneID, -weight, -ENS_mm) %>%
  rename('GeneID' = ENS_hs)
bg2 = attrG %>%
  inner_join(gtexcov)
tfG = mirTest(bg2, miRcol = 1, genecol = 3, patcol = 5)
tfG2 = t(sapply(tfG, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(tfG2) = c('OR', 'pvalue')
tfG2 = cbind(tfG2, 'BH' = p.adjust(tfG2[,'pvalue'], method='BH'))
tfG2 = as.data.frame(tfG2)

tfG2 %>% filter(BH < 0.1)  # 5 significant
# OR       pvalue         BH
# SNAI1  0.7881239 0.0003971847 0.01568879
# TCF3   0.7881239 0.0003971847 0.01568879
# SNAI2  0.7881239 0.0003971847 0.01568879
# TFAP2A 0.8436784 0.0003080968 0.01568879
# STAT5B 1.1660940 0.0025694414 0.08119435

tfG2 %>% filter(pvalue < 0.1 & OR>1) # 9 tfbs
# OR      pvalue         BH
# BCL6   1.218817 0.003864888 0.10177538
# TBP    1.126003 0.006360155 0.14355778
# HMGA1  1.106819 0.055635579 0.61887988
# STAT5B 1.166094 0.002569441 0.08119435
# NRF1   1.105333 0.095430027 0.61887988
# USF1   1.176772 0.029736375 0.47966558
# HIVEP1 1.545733 0.011146449 0.22014237
# RELB   1.117205 0.092380721 0.61887988
# POU3F1 1.587134 0.033394439 0.47966558

#####################
#####################
jonkercov = readRDS('./data/other_datasets/jonker/processed/covch.rds') %>% as.data.frame() %>%
  mutate(pattern = ifelse(rho<0, 'Co', 'Di')) %>% 
  select(-pval, -` BH`) %>%
  rownames_to_column(var='gene_id')
bg3 = attrM %>%
  inner_join(jonkercov)
tfJ = mirTest(bg3, miRcol = 1, genecol = 3, patcol = 5)
tfJ2 = t(sapply(tfJ, function(x) c(x$Fistest$est, x$Fistest$p.val ) ))
colnames(tfJ2) = c('OR', 'pvalue')
tfJ2 = cbind(tfJ2, 'BH' = p.adjust(tfJ2[,'pvalue'], method='BH'))
tfJ2 = as.data.frame(tfJ2)

tfJ2 %>% filter(BH < 0.1)  # 0 significant

tfJ2 %>% filter(pvalue < 0.1 & OR >1) # 11 tfbs
# OR       pvalue        BH
# TEAD2  1.105106 0.0350094557 0.8427443
# LTF    1.155023 0.0333980299 0.8427443
# NFIA   1.146694 0.0418427011 0.8427443
# TCF4   1.087183 0.0441925724 0.8427443
# ETV4   1.099051 0.0629716489 0.9044818
# POU2F1 1.135970 0.0159574481 0.8427443
# POU2F2 1.140249 0.0533382500 0.8427443
# ETS2   1.173458 0.0207586387 0.8427443
# SND1   1.134226 0.0481157628 0.8427443
# MAPK14 1.279514 0.0008205434 0.1296459
# GTF2I  1.327388 0.0839241031 0.9044818

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

# # permutation test:
# xDiCo = bgsum$n[bgsum$pattern=='DiCo']
# nDiCo = length(xDiCo)
# xDiDi = bgsum$n[bgsum$pattern=='DiDi']
# nDiDi = length(xDiDi)
# pop = c(xDiCo, xDiDi)
# obs = c('xDiCo' = mean(xDiCo), 'xDiDi' = mean(xDiDi) )
# obstat = cohens_d(xDiCo, xDiDi)
# #obstat = mean(xDiCo)-mean(xDiDi)
# perms = sapply(1:1000, function(x){
#   px = sample(pop)
#   pDiCo = px[1:nDiCo]
#   pDiDi = px[(nDiCo+1): length(pop)]
#   pstat = cohens_d(pDiCo, pDiDi)
#   #pstat = mean(pDiCo) - mean(pDiDi)
#   return(pstat)
# })
# hist(perms)
# mean(perms >= obstat)
# # p = 0.426
# mean(perms) / obstat
# # fpr= -0.021

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
