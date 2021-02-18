# Deconvolution:
library(stringr)
library(tidyverse)
library(ggpubr)
library(ggforce)
library(reshape2)
library(RColorBrewer)

bulk_expr = readRDS('./data/processed/raw/expression.rds')
age = readRDS('./data/processed/raw/ages.rds')
samp_tissue = readRDS('./data/processed/raw/tissue_ids.rds')
sc_expr = readRDS('./data/other_datasets/scRNA-seq/processed/celltype_expr.rds')[['m3']]
id_convert = readRDS('./data/uniprot_to_ENS.rds')

####################
#################### convert uniprot gene ids to ENS
####################
####################
for(i in 1:4){
  sgenes = intersect(id_convert[,1], rownames(sc_expr[[i]]))
  sc_expr[[i]] = sc_expr[[i]][sgenes,]
  rownames(sc_expr[[i]]) = id_convert[id_convert[,1]%in%sgenes,2]
}

#################### 
bulk_expr = sapply(unique(samp_tissue), function(x){ bulk_expr[,samp_tissue==x] })
bulk_expr_all = bulk_expr
age.each = sapply(unique(samp_tissue),function(x) age[samp_tissue==x])

names(sc_expr) = str_to_title(names(sc_expr))
sc_expr = sc_expr[c(4,1,2,3)]
names(sc_expr)[1] = names(bulk_expr)[1]


# take same genes in bulk and cell type objects, (same genes for each tissue):
for(i in 1:4){
  sgenes = intersect(rownames(bulk_expr[[i]]), rownames(sc_expr[[i]]))
  bulk_expr[[i]] = bulk_expr[[i]][sgenes,]
  sc_expr[[i]] = sc_expr[[i]][sgenes,]
}

####################
#################### deconvolution with all genes ####################
####################
####################

coeffs = sapply(names(bulk_expr), function(a){
  coeffs.tsX = sapply(colnames(bulk_expr[[a]]), function(x){
    coef(lm(bulk_expr[[a]][,x]~sc_expr[[a]]))
  })
  coeffs.tsX = coeffs.tsX[-1,]
  rownames(coeffs.tsX) = colnames(sc_expr[[a]])
  return(coeffs.tsX)
})

# reorder rows by the cell type proportions:
coeffs = lapply(coeffs, function(x) x[order(rowMeans(x),decreasing = T),] )
saveRDS(coeffs,'./data/other_datasets/scRNA-seq/processed/deconvolution_allgenes.rds')

####################
#################### deconvolution with DC genes  ####################
####################
####################

dgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
dcgenes = names(dgenes[dgenes<0])

bulk_expr_dc  = list()
sc_expr_dc  = list()
for(i in 1:4){
  sgenes = intersect(dcgenes, rownames(bulk_expr[[i]]))
  bulk_expr_dc[[i]] = bulk_expr[[i]][sgenes,]
  sc_expr_dc[[i]] = sc_expr[[i]][sgenes,]
}
names(bulk_expr_dc) = names(bulk_expr)
names(sc_expr_dc) = names(sc_expr)

coeffs_dc = sapply(names(bulk_expr_dc), function(a){
  coeffs.tsX = sapply(colnames(bulk_expr_dc[[a]]), function(x){
    coef(lm(bulk_expr_dc[[a]][,x]~sc_expr_dc[[a]]))
  })
  coeffs.tsX = coeffs.tsX[-1,]
  rownames(coeffs.tsX) = colnames(sc_expr_dc[[a]])
  return(coeffs.tsX)
})
# reorder rows by the cell type proportions:
coeffs_dc = lapply(coeffs_dc, function(x) x[order(rowMeans(x),decreasing = T),] )
saveRDS(coeffs_dc, './data/other_datasets/scRNA-seq/processed/deconvolution_dc_genes.rds')

####################
#################### test cell type proportion changes of dc genes with age ####################
####################
####################

obs_dc = reshape2::melt(coeffs_dc) %>% 
  set_names(c('cell type','id','proportion','tissue')) %>%
  left_join(data.frame(id = names(age), age=age) ) %>%
  mutate(period= ifelse(age < 90,'Development','Ageing'))  %>%
  group_by(tissue, `cell type`,period) %>% 
  summarise(`Proportion Change` = cor(age, proportion, m='s'),
            `propch_p` = cor.test(age, proportion, m='s')$p.val)
  

####################
#################### deconvolution with non-DC genes ####################
####################
####################

bulk_expr_nondc  = list()
sc_expr_nondc  = list()
for(i in 1:4){
  sgenes = intersect(dcgenes, rownames(bulk_expr[[i]]))
  bulk_expr_nondc[[i]] = bulk_expr[[i]][!rownames(bulk_expr[[i]])%in%sgenes,]
  sc_expr_nondc[[i]] = sc_expr[[i]][!rownames(sc_expr[[i]])%in%sgenes,]
}
names(bulk_expr_nondc) = names(bulk_expr)
names(sc_expr_nondc) = names(sc_expr)

coeffs_nondc = sapply(names(bulk_expr_nondc), function(a){
  coeffs.tsX = sapply(colnames(bulk_expr_nondc[[a]]), function(x){
    coef(lm(bulk_expr_nondc[[a]][,x]~sc_expr_nondc[[a]]))
  })
  coeffs.tsX = coeffs.tsX[-1,]
  rownames(coeffs.tsX) = colnames(sc_expr_nondc[[a]])
  return(coeffs.tsX)
})
# reorder rows by the cell type proportions:
coeffs_nondc = lapply(coeffs_nondc, function(x) x[order(rowMeans(x),decreasing = T),] )
saveRDS(coeffs_nondc, './data/other_datasets/scRNA-seq/processed/deconvolution_nondc_genes.rds')

####################
#################### test cell type proportion changes of non-dc genes with age
####################
####################

obs_nondc = reshape2::melt(coeffs_nondc) %>% 
  set_names(c('cell type','id','proportion','tissue')) %>%
  left_join(data.frame(id = names(age), age=age) ) %>%
  mutate(period= ifelse(age < 90,'Development','Ageing'))  %>%
  group_by(tissue, `cell type`,period) %>% 
  summarise(`Proportion Change` = cor(age, proportion, m='s'),
            `propch_p` = cor.test(age, proportion, m='s')$p.val)

obs_dc %>% filter(tissue=='Liver' & `cell type` =='hepatocyte')
obs_nondc %>% filter(tissue=='Liver' & `cell type` =='hepatocyte')
obs_dc %>% filter(tissue=='Lung' & `cell type` =='bronchial smooth muscle cell')
obs_nondc %>% filter(tissue=='Lung' & `cell type` =='bronchial smooth muscle cell')
obs_dc %>% filter(tissue=='Muscle')
obs_nondc %>% filter(tissue=='Muscle')
obs_dc %>% filter(tissue=='Cortex')
obs_nondc %>% filter(tissue=='Cortex')

####################
#################### 
####################  test dc vs non-dc difference with permutations ####################
####################

nondc_perms = list()
for(perm in 1:1000){
  print(perm/1000)
  bulk_expr_nondc_perm  = list()
  sc_expr_nondc_perm  = list()
  for(i in 1:4){
    sgenes = intersect(dcgenes, rownames(bulk_expr[[i]]))
    nondc = rownames(bulk_expr[[i]])[!rownames(bulk_expr[[i]])%in%sgenes]
    chooseg = sample(nondc, size = length(sgenes))
    bulk_expr_nondc_perm[[i]] = bulk_expr[[i]][chooseg,]
    sc_expr_nondc_perm[[i]] = sc_expr[[i]][chooseg,]
  }
  names(bulk_expr_nondc_perm) = names(bulk_expr)
  names(sc_expr_nondc_perm) = names(sc_expr)
  
  coeffs_nondc_perm = sapply(names(bulk_expr_nondc_perm), function(a){
    coeffs.tsX = sapply(colnames(bulk_expr_nondc_perm[[a]]), function(x){
      coef(lm(bulk_expr_nondc_perm[[a]][,x]~sc_expr_nondc_perm[[a]]))
    })
    coeffs.tsX = coeffs.tsX[-1,]
    rownames(coeffs.tsX) = colnames(sc_expr_nondc_perm[[a]])
    return(coeffs.tsX)
  })
  # nondc_perms[[perm]] = reshape2::melt(coeffs_nondc_perm) %>% 
  #   set_names(c('cell type','id','proportion','tissue')) %>%
  #   left_join(data.frame(id = names(age), age=age) ) %>%
  #   mutate(period= ifelse(age < 90,'Development','Ageing'))  %>%
  #   group_by(tissue, `cell type`, period) %>%  
  #   summarise(`Proportion Change` = cor(age, proportion, m='s'))
  nondc_perms[[perm]] = reshape2::melt(coeffs_nondc_perm) %>% 
    set_names(c('cell type','id','proportion','tissue')) %>%
    left_join(data.frame(id = names(age), age=age) ) %>%
    mutate(period= ifelse(age < 90,'Development','Ageing'))  %>%
    group_by(tissue, `cell type`, period) %>%  
    summarise(`Proportion Change` = coefficients(lm(proportion~age))[2] )
}
names(nondc_perms) = paste0('perm',1:1000)
#saveRDS(nondc_perms, './data/processed/raw/deconvolution_perm.rds')


nondc_perms = reshape2::melt(nondc_perms) %>% 
  rename('Proportion Change' = value, perm = L1) %>%
  select(-variable)

#######
#################### test cell type proportion changes of nondc genes with age observed values with lm
obs_nondc = reshape2::melt(coeffs_nondc) %>% 
  set_names(c('cell type','id','proportion','tissue')) %>%
  left_join(data.frame(id = names(age), age=age) ) %>%
  mutate(period= ifelse(age < 90,'Development','Ageing'))  %>%
  group_by(tissue, `cell type`,period) %>% 
  summarise(`Proportion Change` = coefficients(lm(proportion~age))[2],
            `propch_p` = summary(lm(proportion~age))$coeff[2,4] )

kpr = obs_nondc %>% filter(tissue=='Cortex' & `cell type`=='neuron', period=='Ageing') %>% pull(`Proportion Change`)
kk = nondc_perms %>%
  filter(tissue=='Cortex' & `cell type`=='neuron', period=='Ageing') %>% pull(`Proportion Change`)
mean(kk >= kpr)

perm_results = nondc_perms %>% 
  left_join(mutate(obs_nondc, obs='obs'), by=c('tissue','cell type','period') ) %>%
  rename('Perm Change'=`Proportion Change.x`,
         'Obs Change'=`Proportion Change.y`) %>%
  group_by(period, tissue, `cell type`) %>%
  summarise(pval = mean(`Perm Change` >= unique(`Obs Change`) ),
            FPR = median(`Perm Change`)/ unique(`Obs Change`) )

sum(perm_results$pval< 0.05)

obs_nondc %>% 
  mutate(period= factor(period, levels=c('Development','Ageing'))) %>%
  ggplot(aes(x=`cell type`, y=`Proportion Change`)) +
  facet_grid(tissue~period) +
  geom_bar(stat='identity', position = 'dodge')


####
####################
#################### 
####################
####################




