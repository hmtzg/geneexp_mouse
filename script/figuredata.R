library(tidyverse)

age = readRDS('data/processed/figures/raw/age.rds')
tis = readRDS('data/processed/figures/raw/ts.ord.rds')
exp = readRDS('data/processed/figures/raw/exp.rds')
tissue_names = setNames(c('Cortex','Lung','Liver','Muscle'), c('ctx','lng','lv','ms'))
samp_id = paste('s',1:length(age),sep='')


sample_info = data.frame(ind_id = colnames(exp), age = unname(age), tissue = tis) %>%
  mutate(tissue = tissue_names[as.character(tissue)]) %>%
  mutate(log2age = log2(age)) %>%
  mutate(sample_id = samp_id) 

colnames(exp) = samp_id

saveRDS(sample_info,'./data/processed/figures/tidy/sample_info.rds')

exp = reshape2::melt(exp) %>%
  set_names(c('gene_id','sample_id','expression')) 

saveRDS(exp,'./data/processed/figures/tidy/expression.rds')

pca_data_all = list(readRDS('data/processed/figures/raw/pc.aging.scaled.rds'),
                readRDS('data/processed/figures/raw/pc.all.rds'),
                readRDS('data/processed/figures/raw/pc.scaled.rds'),
                readRDS('data/processed/figures/raw/pc.aging.rds'),
                readRDS('data/processed/figures/raw/pc.dev.rds'),
                readRDS('data/processed/figures/raw/pc.dev.scaled.rds'))
names(pca_data_all) = c('aging_notissue','all_raw','all_notissue','aging_raw','development_raw','development_notissue')

pca_data = lapply(pca_data_all,function(x)x$x[,1:4])
rownames(pca_data[[1]]) = samp_id[age>=93]
rownames(pca_data[[2]]) = samp_id
rownames(pca_data[[3]]) = samp_id
rownames(pca_data[[4]]) = samp_id[age>=93]
rownames(pca_data[[5]]) = samp_id[age<93]
rownames(pca_data[[6]]) = samp_id[age<93]
pca_data = reshape2::melt(pca_data) %>%
  set_names(c('sample_id','PC','value','type')) %>%
  separate(type,into=c('period','type'))

pca_data = sapply(pca_data_all,function(x)summary(x)$imp[2,1:4]) %>%
  reshape2::melt() %>%
  set_names(c('PC','type','varExp')) %>%
  separate(type,into=c('period','type')) %>%
  right_join(pca_data) 

saveRDS(pca_data,'./data/processed/figures/tidy/pca_data.rds')

expch_dev = readRDS('./data/processed/figures/raw/dev.expage.cor.rds')
names(expch_dev) = tissue_names[names(expch_dev)]

expch_age = readRDS('./data/processed/figures/raw/aging.expage.cor.rds')
names(expch_age) = tissue_names[names(expch_age)]

expch_dev = reshape2::melt(expch_dev) %>%
  spread(key= `Var2`, value = `value`) %>%
  mutate(period = 'development')

expch_age = reshape2::melt(expch_age) %>%
  spread(key= `Var2`, value = `value`) %>%
  mutate(period = 'aging')

expch = rbind(expch_dev,expch_age) %>%
  set_names(c('gene_id','tissue','Expression Change','p','FDR','period'))

saveRDS(expch,'./data/processed/figures/tidy/expression_change.rds')

sdmat = readRDS('./data/processed/figures/raw/genesd.rds')
sdmat = reshape2::melt(sdmat) %>%
  set_names('gene_id','ind_id','CoV')
saveRDS(sdmat,'./data/processed/figures/tidy/CoV.rds')

covch_dev = readRDS('./data/processed/figures/raw/dev.sdage.cor.rds')
covch_aging = readRDS('./data/processed/figures/raw/aging.sdage.cor.rds')

covch_dev = data.frame(covch_dev, period = 'development') %>%
  rename(CoV_change = rho, 
         p = V2) %>%
  mutate(FDR = p.adjust(p, method='BH'))

covch_aging = data.frame(covch_aging, period = 'aging') %>%
  rename(CoV_change = rho, 
         p = V2) %>%
  mutate(FDR = p.adjust(p, method='BH'))

covch = rbind(covch_dev,covch_aging)
saveRDS(covch,'./data/processed/figures/tidy/CoV_change.rds')

aging_gse = readRDS('./data/processed/figures/raw/aging.gse.rds')@result[,1:6] %>%
  mutate(p.adj = p.adjust(pvalue,method='BH')) %>%
  mutate(period = 'aging')
dev_gse = readRDS('./data/processed/figures/raw/dev.gse.rds')@result[,1:6] %>%
  mutate(p.adj = p.adjust(pvalue,method='BH'))%>%
  mutate(period = 'development')
cov_gse = rbind(aging_gse,dev_gse)
saveRDS(cov_gse,'./data/processed/figures/tidy/CoV_GSEA.rds')

rev_gse = readRDS('./data/processed/figures/raw/rev.gse.rds')@result[,1:6] %>%
  mutate(p.adj = p.adjust(pvalue,method='BH'))
saveRDS(rev_gse,'./data/processed/figures/tidy/divcon_GSEA.rds')

