library(tidyverse)
library(openxlsx)
##########
##########
dc = readRDS('./data/other_datasets/scRNA-seq/processed/deconvolution_dc_genes.rds')
nondc = readRDS('./data/other_datasets/scRNA-seq/processed/deconvolution_nondc_genes.rds')
age = readRDS('./data/processed/raw/ages.rds')

coef = reshape2::melt(list(DiCo=dc,`Non-DiCo`=nondc)) %>%
  set_names(c('cell type','id','proportion','tissue','Geneset')) %>%
  left_join(data.frame(id=names(age), age= age))

top_celltype = coef %>%  # sdata1 
  group_by(tissue, `cell type`, Geneset) %>%
  summarise(mprop = mean(proportion)) %>%
  group_by(tissue, Geneset) %>%
  top_n(n=1, wt=mprop) %>%
  left_join(coef, by=c('tissue','cell type','Geneset')) %>%
  mutate(period= ifelse(age<90,'development','ageing')) 

prop_changes = readRDS('./data/other_datasets/scRNA-seq/processed/deconvolution_prop_change.rds')
prop_ch = reshape2::melt(prop_changes, id.vars=c('tissue','cell type','period')) %>% 
  rename('Proportion Change (rho)' = value, 'geneset' = L1) %>%
  spread(key = 'variable', value='Proportion Change (rho)') %>%
  rename('pval'= 'propch_p') #  sdata 2

coeffs = readRDS('./data/other_datasets/scRNA-seq/processed/deconvolutions_combined.rds') #  sdata 3

perm_result_tab = readRDS('./data/other_datasets/scRNA-seq/processed/dico_vs_nondico_prop_test.rds') # sdata 4

maxmin_density = readRDS('./data/other_datasets/scRNA-seq/processed/max_min_cors_density.rds') # sdata 5
maxmin_ch = readRDS('./data/other_datasets/scRNA-seq/processed/max_min_cors_ch.rds') # sdata 6

intra_ts_cov = readRDS('./data/other_datasets/scRNA-seq/processed/intra_tissue_mean_cov.rds')
intra_ts_cov = intra_ts_cov %>%
  mutate(`age group (month)` = factor(`age group`, levels = c('m3', 'm18', 'm24')) ) %>%
  mutate(age = as.numeric( gsub('[a-z]','',`age group`) ) )

write.xlsx(list(sdata0 = NULL,
                sdata1 = top_celltype,
                sdata2 = prop_ch,
                sdata3 = coeffs,
                sdata4 = perm_result_tab,
                sdata5 = maxmin_density,
                sdata6 = maxmin_ch,
                sdata7 = intra_ts_cov), 
           file = './results/source_data/figure_5_source_data_1.xlsx' )

#nondc_perms= readRDS('./data/other_datasets/scRNA-seq/processed/deconvolution_perm.rds')










