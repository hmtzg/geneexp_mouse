library(tidyverse)
library(ggpubr)
dat1 = readRDS('./data/processed/tidy/CoV_change.rds') %>% 
  tibble()  %>%
  rename('BH' = FDR) %>% relocate(gene_id)

# mouse our data, background is development divergent genes
dat1 = dat1 %>%
  filter(period=='development' & CoV_change > 0 ) %>%
  select(gene_id) %>%
  left_join(dat1, by ='gene_id') %>%
  filter(period=='aging')
bg = as.character(dat1$gene_id)

# jonker data:
dat2 = readRDS('./data/other_datasets/jonker/processed/covch.rds') %>% as.data.frame() %>%
  #rownames_to_column(var = 'gene_id') %>%
  set_names(c('gene_id','CoV_change','pval', 'BH')) %>% tibble() %>%
  filter(gene_id%in%bg)

# GTEx:
dat3 = readRDS('./results/GTEx/covcor.rds') %>%
  set_names(c('gene_id', 'CoV_change', 'pval', 'BH')) #  human, GTEx

## convert human ENS to mouse:
martx = biomaRt::useMart(biomart = 'ensembl')
martHs_ = biomaRt::useDataset('hsapiens_gene_ensembl', mart=martx)
martMm_ = biomaRt::useDataset('mmusculus_gene_ensembl', mart=martx)
ensmap = biomaRt::getLDS(attributes = c('ensembl_gene_id'), filters = 'ensembl_gene_id', 
                         values = dat3$gene_id, mart = martHs_, 
                         attributesL = c('ensembl_gene_id'), martL = martMm_)
colnames(ensmap) = c('ENS_hs', 'ENS_mm')
dup1 = unique(ensmap$ENS_hs[duplicated(ensmap$ENS_hs)]) #  508 duplicate genes
ensmap = ensmap[!ensmap$ENS_hs%in%dup1,]
dup2 = unique(ensmap$ENS_mm[duplicated(ensmap$ENS_mm)]) # 102 duplicate genes
ensmap = ensmap[!ensmap$ENS_mm%in%dup2,] 
colnames(ensmap)[1]='gene_id'

dat3 = dat3 %>% right_join(ensmap)
dat3backup = dat3

dat3 = dat3 %>%
  filter(ENS_mm%in%bg)


##########
########## Number of genes overlap:
##########

## our data and GTEx:
ovrlp13 = dat1 %>% select(-pval, -BH, -period) %>%
  inner_join(rename(select(dat3, -pval, -BH, -gene_id ), 'gene_id'=ENS_mm), by='gene_id'  ) %>%
  filter(gene_id%in%bg) %>%
  filter( CoV_change.x!=0 & CoV_change.y!=0 ) %>%
  summarise(tb = table(sign(CoV_change.x), sign(CoV_change.y)))
ovrlp13
# tb[,"-1"] [,"1"]
# <int>  <int>
# 1      2514   1727
# 2      2134   1537
fisher.test(ovrlp13)
# OR = 1.048
# p = 0.3029

## our data and Jonker:
ovrlp12 = dat1 %>% select(-pval, -BH, -period) %>%
  inner_join(select(dat2, -pval, -BH), by='gene_id'  ) %>%
  filter(gene_id%in%bg) %>%
  filter( CoV_change.x!=0 & CoV_change.y!=0 ) %>%
  summarise(tb = table(sign(CoV_change.x), sign(CoV_change.y)))
ovrlp12
# tb[,"-1"] [,"1"]
# <int>  <int>
# 1      2607   1856
# 2      2042   1804
fisher.test(ovrlp12)
# OR = 1.230 
# p = 1.21e-6

## Jonker and GTEx:
ovrlp23 = dat2 %>% select(-pval, -BH) %>%
  inner_join(rename(select(dat3, -pval, -BH, -gene_id), 'gene_id'=ENS_mm), by='gene_id'  ) %>%
  filter(gene_id%in%bg) %>%
  filter( CoV_change.x!=0 & CoV_change.y!=0 ) %>%
  summarise(tb = table(sign(CoV_change.x), sign(CoV_change.y)))
ovrlp23
# tb[,"-1"] [,"1"]
# <int>  <int>
# 1      2541   1707
# 2      1933   1434
fisher.test(ovrlp23)
# OR = 1.104
# p = 0.03495

##############################

##########
########## Correlation among datasets:
##########

dat1 %>% select(-pval, -BH, -period) %>%
  inner_join(rename(select(dat3, -pval, -BH, -gene_id ), 'gene_id'=ENS_mm), by='gene_id'  ) %>%
  summarise(rho = cor.test(CoV_change.x, CoV_change.y, m='s')$est,
            pval = cor.test(CoV_change.x, CoV_change.y, m='s')$p.val)
# rho     pval
# <dbl>    <dbl>
#   1 0.0408 0.000268
## our data and GTEx:

mm_gtex = dat1 %>% select(-pval, -BH, -period) %>%
  inner_join(rename(select(dat3, -pval, -BH, -gene_id ), 'gene_id'=ENS_mm), by='gene_id'  ) %>%
  ggplot(aes(x=CoV_change.x, y=CoV_change.y))  +
  geom_point(alpha= 0.4, color = 'gray30') +
  geom_smooth(method = 'lm') +
  xlab('our data') +
  ylab('GTEx') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size=2)

## our data and Jonker:
mm_jonker = dat1 %>% select(-pval, -BH, -period) %>%
  inner_join(select(dat2, -pval, -BH), by='gene_id'  ) %>%
  ggplot(aes(x=CoV_change.x, y=CoV_change.y))  +
  geom_point(alpha= 0.4, color = 'gray30') +
  geom_smooth(method = 'lm') +
  xlab('our data') +
  ylab('Jonker') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size=2)
# rho = 0.084, p = 1.9e-14  

## Jonker and GTEx:
jonker_gtex = dat2 %>% select(-pval, -BH) %>%
  inner_join(rename(select(dat3, -pval, -BH, -gene_id), 'gene_id'=ENS_mm), by='gene_id'  ) %>%
  ggplot(aes(x=CoV_change.x, y=CoV_change.y))  +
  geom_point(alpha= 0.4, color = 'gray30') +
  geom_smooth(method = 'lm') +
  ylab('GTEx') +
  xlab('Jonker') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size=2)
# rho = 0.02, p = 0.078

co_cors = ggarrange(mm_gtex, mm_jonker, jonker_gtex, ncol=3)
ggsave('./results/SI_figures/dataset_overlap/co_cors.pdf', co_cors,units='cm', height = 8, width = 16,
       useDingbats=F)
ggsave('./results/SI_figures/dataset_overlap/co_cors.png', co_cors, units='cm', height = 8, width = 16)
##########
########## with significant genes in either of datasets:
##########

dat1 %>% select(-pval, -period) %>%
  inner_join(rename(select(dat3, -pval, -gene_id ), 'gene_id'=ENS_mm), by='gene_id'  ) %>%
  filter(BH.x < 0.3 & BH.y < 0.3) %>%
  ggplot(aes(x=CoV_change.x, y=CoV_change.y))  +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('our data') +
  ylab('GTEx') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho')

dat1 %>% select(-pval, -period) %>%
  inner_join(select(dat2, -pval), by='gene_id'  ) %>%
  filter(BH.x < 0.1 & BH.y < 0.1) %>%
  ggplot(aes(x=CoV_change.x, y=CoV_change.y))  +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('our data') +
  ylab('Jonker') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho')

dat2 %>% select(-pval) %>%
  inner_join(rename(select(dat3, -pval, -gene_id), 'gene_id'=ENS_mm), by='gene_id'  ) %>%
  filter(BH.x < 0.1 | BH.y < 0.1) %>%
  ggplot(aes(x=CoV_change.x, y=CoV_change.y))  +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('Jonker') +
  ylab('GTEx') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho')

##########
########## Pathway overlap
##########
mm_gse = readRDS('./data/processed/raw/dc_gse.rds')
mmsig = mm_gse@result[mm_gse@result$qvalues<0.1, 1:10]

covch = readRDS('./data/processed/tidy/CoV_change.rds')
devdiv = covch %>% 
  filter(period=='development' & CoV_change > 0 ) %>% 
  select(gene_id) %>%
  left_join( filter(covch, period=='aging') ) %>% 
  select(gene_id, CoV_change)
divgenes = devdiv$CoV_change
names(divgenes) = devdiv$gene_id
library(clusterProfiler)
dc_gse = gseGO(geneList = sort(divgenes, decreasing = T), OrgDb = org.Mm.eg.db, ont = "BP", 
               pvalueCutoff = 1,
               keyType = "ENSEMBL", nPerm = 1000, minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BY',
               verbose = F)
dc_gse@result[1:5,1:10]
sum(dc_gse@result$qvalues<0.1)
saveRDS(dc_gse, file = './results/SI_figures/dataset_overlap/mm_co_gse.rds')

jonker_gse = readRDS('./data/other_datasets/jonker/dico_gse.rds')
jonker_gse@result[1:5,1:10]
jonkersig = jonker_gse@result[jonker_gse@result$qvalues<0.1, 1:10]
length(intersect(mmsig$ID, jonkersig$ID)) # 56 categories
gtex_gse = readRDS('./results/GTEx/dico_gse.rds')

#### sig categories overlap:
mmsig %>% select(ID, NES) %>%
  inner_join(select(jonkersig, ID, NES), by='ID') %>% # 56 categories common
  summarise(rho = cor.test(NES.x, NES.y, m='s')$est,
            p = cor.test(NES.x, NES.y, m='s')$p.val)
# rho = 0.15, p =0.27
mmsig %>% select(ID, NES, Description) %>%
  inner_join(select(jonkersig, ID, NES), by='ID') %>%
  filter(sign(NES.x)!= sign(NES.y))

mmsig %>% select(ID, NES) %>%
  inner_join(select(jonkersig, ID, NES), by='ID') %>%
  ggplot(aes(x=NES.x, y=NES.y)) +
  geom_point() +
  xlab('our data') +
  ylab('Jonker') +
  geom_smooth(method = 'lm') + 
  stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 'middle')

#### all GO correlation
select(mm_gse@result[,1:10], ID, NES, Description) %>%
  inner_join( select(jonker_gse@result[,1:10], ID, NES), by='ID' ) %>%
  summarise(rho = cor.test(NES.x, NES.y, m='s')$est,
            p = cor.test(NES.x, NES.y, m='s')$p.val)
# rho = 0.26, = 4.3e-70

select(mm_gse@result[,1:10], ID, NES, Description) %>%
  inner_join( select(gtex_gse@result[,1:10], ID, NES), by='ID' ) %>%
  summarise(rho = cor.test(NES.x, NES.y, m='s')$est,
            p = cor.test(NES.x, NES.y, m='s')$p.val)
# rho = 0.08, p=2.16e-7

select(jonker_gse@result[,1:10], ID, NES, Description) %>%
  inner_join( select(gtex_gse@result[,1:10], ID, NES), by='ID' ) %>%
  summarise(rho = cor.test(NES.x, NES.y, m='s')$est,
            p = cor.test(NES.x, NES.y, m='s')$p.val)
# rho = 0.17, p=8.1e-29

select(mm_gse@result[,1:10], ID, NES, Description) %>%
  inner_join( select(jonker_gse@result[,1:10], ID, NES), by='ID' ) %>%
  ggplot(aes(x=NES.x, y=NES.y)) +
  geom_point() +
  xlab('our data') +
  ylab('Jonker') +
  geom_smooth(method = 'lm') + 
  stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 'middle')

#### with our data analysed the same way as jonker and gtex:
dplyr::select(dc_gse@result[,1:10], ID, NES, Description) %>%
  inner_join( dplyr::select(jonker_gse@result[,1:10], ID, NES), by='ID' ) %>%
  summarise(rho = cor.test(NES.x, NES.y, m='s')$est,
            p = cor.test(NES.x, NES.y, m='s')$p.val)
# rho = 0.27, = 8.3e-83

dplyr::select(dc_gse@result[,1:10], ID, NES, Description) %>%
  inner_join( dplyr::select(gtex_gse@result[,1:10], ID, NES), by='ID' ) %>%
  summarise(rho = cor.test(NES.x, NES.y, m='s')$est,
            p = cor.test(NES.x, NES.y, m='s')$p.val)
# rho = 0.088, = 6.7e-9

save(list=ls(), file='./data/Co_overlap_among_datasets.rdata')

