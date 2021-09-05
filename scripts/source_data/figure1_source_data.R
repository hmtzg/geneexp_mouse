library(tidyverse)
library(openxlsx)
##########
##########
expr_ch = readRDS('data/processed/tidy/expression_change.rds')
########## Figure 1- source data 1
########## a,b
pca_data = readRDS('./data/processed/tidy/pca_data.rds') 
sample_info = readRDS('./data/processed/tidy/sample_info.rds')

pca_data = pca_data %>% # sdata1
  left_join(sample_info) %>%
  select(-log2age)

pc_age_cors  = pca_data %>% #  sdata2
  rename(`Scale period` = period) %>%
  left_join(select(sample_info,-ind_id,-log2age)) %>%
  mutate(`age period` = ifelse(age< 90,'Dev','Aging') ) %>%
  group_by(`Scale period`, type, PC, tissue, `age period`) %>%
  summarise(rho = cor.test(age,value,m='s')$est,
            p = cor.test(age,value,m='s')$p.val)

mdist = readRDS('./data/processed/tidy/mean_euclidean_dist.rds') # sdata3
cors  = readRDS('./data/processed/raw/pwise_expch_cors.rds') # sdata4

pwise_cor_perm_test = readRDS('./data/processed/tidy/pwise_expch_cor_perm_test.rds') # sdata5

n_exp_ch = expr_ch %>% # sdata6
  drop_na() %>%
  filter(`Expression Change` !=0 ) %>% # drop expression - age rho = 0 values
  mutate(direction = `Expression Change` > 0) %>%
  mutate(direction = ifelse( direction == TRUE, 'Up', 'Down')) %>%
  mutate(period = gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels = c('Development','Ageing'))) 

perm_overlaps = readRDS('./data/processed/raw/tissue_gene_overlaps_perm.rds') # sdata7
perm_overlaps_fdr0.1 = readRDS('./data/processed/raw/tissue_siggene_fdr01_overlaps_perm.rds') # sdata8

cors.06 = expr_ch %>%
  filter(abs(`Expression Change`)> 0.6 ) %>%
  select(-p, -FDR) %>%
  mutate(tissue = factor(tissue, levels= c('Cortex','Liver', 'Lung','Muscle'))) %>%
  pivot_wider(names_from = c(tissue, period), values_from=`Expression Change`) %>% 
  select(-gene_id) %>% 
  as.matrix() %>%
  cor(method='s', use ='pairwise.complete.obs')
# colnames(cors.06) = sub('\\_.*','',colnames(cors.06))
# rownames(cors.06) = sub('\\_.*','',rownames(cors.06))
cors.06 = cors.06[c(1,2,3,4,7,6,5,8),c(1,2,3,4,7,6,5,8)] # sdata9

revgenes = readRDS('data/processed/tidy/revgenes.tissues.rds') # sdata10

UD.tissue = readRDS(file = './data/processed/raw/updown_perm_each.rds')
DU.tissue = readRDS(file = './data/processed/raw/downup_perm_each.rds')
rev.tissue = list(UpDown = lapply(UD.tissue,function(x) x$dist),
                  DownUp = lapply(DU.tissue, function(x)x$dist))
revdat = reshape2::melt(rev.tissue) %>% 
  set_names(c('Proportion', 'tissue', 'Reversal')) %>% 
  mutate(Reversal = factor(Reversal, levels= c('UpDown', 'DownUp')) ) # sdata11

shared_rev_perm_test = readRDS('./data/processed/tidy/shared_rev_perm_test.rds') # sdata12

write.xlsx(list(sdata0 = NULL,
                sdata1 = pca_data, 
                sdata2 = pc_age_cors, 
                sdata3 = mdist, 
                sdata4 = cors,
                sdata5 = pwise_cor_perm_test,
                sdata6 = n_exp_ch,
                sdata7 = perm_overlaps,
                sdata8 = perm_overlaps_fdr0.1, 
                sdata9 = cors.06,
                sdata10 = revgenes, 
                sdata11 = revdat,
                sdata12 = shared_rev_perm_test), 
           file = './results/source_data/figure_1_source_data_1.xlsx', row.names=T)


# corperm.dev = readRDS('./data/processed/raw/pairwise_tissue_cors_perm_dev.rds') # sheet 5
# corperm.aging = readRDS('./data/processed/raw/pairwise_tissue_cors_perm_aging.rds') # sheet 6












