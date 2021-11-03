# Deconvolution with scRNA-seq data:
library(stringr)
library(tidyverse)
library(ggpubr)
library(ggforce)
library(reshape2)
library(RColorBrewer)
library(openxlsx)
library(ggridges)

theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
periodcol = c(Development = "#FE6100", Ageing ="#648FFF")

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
#################### range = [12492, 12849]
range(sapply(bulk_expr,nrow))

coeffs = sapply(names(bulk_expr), function(a){
  coeffs.tsX = sapply(colnames(bulk_expr[[a]]), function(x){
    coef(lm(bulk_expr[[a]][,x]~sc_expr[[a]]))
  })
  coeffs.tsX = coeffs.tsX[-1,]
  rownames(coeffs.tsX) = colnames(sc_expr[[a]])
  return(coeffs.tsX)
})

# reorder rows by the cell type proportions:
coeffs_allgenes = coeffs
coeffs = lapply(coeffs, function(x) x[order(rowMeans(x),decreasing = T),] )
saveRDS(coeffs, './data/other_datasets/scRNA-seq/processed/deconvolution_allgenes.rds')

obs_allgenes = reshape2::melt(coeffs) %>% 
  set_names(c('cell type','id','proportion','tissue')) %>%
  left_join(data.frame(id = names(age), age=age) ) %>%
  mutate(period= ifelse(age < 90,'Development','Ageing'))  %>%
  group_by(tissue, `cell type`, period) %>% 
  summarise(`Proportion Change` = cor(age, proportion, m='s'),
            `propch_p` = cor.test(age, proportion, m='s')$p.val)

####################
#################### deconvolution with DC genes  ####################
####################
#################### range = [4007, 4106]

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

range(sapply(bulk_expr_dc, nrow))
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
#################### range = [8485, 8743]

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

range(sapply(bulk_expr_nondc, nrow))
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

# obs_nondc %>% filter(tissue=='Liver' & `cell type` =='hepatocyte')
# obs_nondc %>% filter(tissue=='Lung' & `cell type` =='bronchial smooth muscle cell')
# obs_nondc %>% filter(tissue=='Cortex' & `cell type` =='neuron')
# obs_nondc %>% filter(tissue=='Muscle' & `cell type` =='skeletal muscle satellite cell')
# 
# obs_dc %>% filter(tissue=='Liver' & `cell type` =='hepatocyte')
# obs_dc %>% filter(tissue=='Lung' & `cell type` =='bronchial smooth muscle cell')
# obs_dc %>% filter(tissue=='Cortex' & `cell type` =='neuron')
# obs_dc %>% filter(tissue=='Muscle' & `cell type` =='skeletal muscle satellite cell')
# 
# obs_allgenes %>% filter(tissue=='Liver' & `cell type` =='hepatocyte')
# obs_allgenes %>% filter(tissue=='Lung' & `cell type` =='bronchial smooth muscle cell')
# obs_allgenes %>% filter(tissue=='Cortex' & `cell type` =='neuron')
# obs_allgenes %>% filter(tissue=='Muscle' & `cell type` =='skeletal muscle satellite cell')

prop_changes = list(allgenes = obs_allgenes, dc_genes = obs_dc, nondc_genes = obs_nondc)
saveRDS(prop_changes, './data/other_datasets/scRNA-seq/processed/deconvolution_prop_change.rds')

prop_ch = reshape2::melt(prop_changes, id.vars=c('tissue','cell type','period')) %>% 
  rename('Proportion Change (rho)' = value, 'geneset' = L1) %>%
  spread(key = 'variable', value='Proportion Change (rho)') %>%
  rename('pval'= 'propch_p')

write.xlsx(prop_ch, './results/SI_tables/TableS12.xlsx')

###########

coeffs = reshape2::melt(list(`All-genes` = coeffs_allgenes ,DiCo=coeffs_dc,`Non-DiCo` = coeffs_nondc)) %>%
  set_names(c('cell type', 'id', 'proportion', 'tissue','Geneset')) %>%
  left_join(data.frame(id=names(age), age= age))

saveRDS(coeffs,'./data/other_datasets/scRNA-seq/processed/deconvolutions_combined.rds')

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
    #nondc = rownames(bulk_expr[[i]])[!rownames(bulk_expr[[i]])%in%sgenes]
    # choose among all genes:
    allg = rownames(bulk_expr[[i]])
    chooseg = sample(allg, size = length(sgenes))
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
  nondc_perms[[perm]] = reshape2::melt(coeffs_nondc_perm) %>%
    set_names(c('cell type','id','proportion','tissue')) %>%
    left_join(data.frame(id = names(age), age=age) ) %>%
    mutate(period= ifelse(age < 90,'Development','Ageing'))  %>%
    group_by(tissue, `cell type`, period) %>%
    summarise(`Proportion Change` = cor(age, proportion, m='s'))
}
names(nondc_perms) = paste0('perm',1:1000)
#nondc_perms = nondc_perms[1:1000]
saveRDS(nondc_perms, './data/other_datasets/scRNA-seq/processed/deconvolution_perm.rds')


######################## perm results
nondc_perms2 = reshape2::melt(nondc_perms, id.vars=c('tissue','cell type','period')) %>% 
  rename('Proportion Change' = value, perm = L1) %>%
  select(-variable)

perm_result_dev = nondc_perms2 %>%
  filter(period == 'Development') %>% select(-period) %>%
  rename(Perms = `Proportion Change`) %>% 
  left_join( select(filter(obs_dc, period == 'Development'),-propch_p, -period),
             by=c('tissue','cell type') ) %>% 
  group_by(tissue, `cell type`) %>%
  summarise(pval = ifelse(unique(`Proportion Change`)>0,
                          mean(Perms >= unique(`Proportion Change`)),
                          mean(Perms <= unique(`Proportion Change`))),
            FPR = median(Perms) / unique(`Proportion Change`) ) %>%
  ungroup()

# perm_result_dev = nondc_perms2 %>%
#   filter(period == 'Development') %>% select(-period) %>%
#   rename(Perms = `Proportion Change`) %>% 
#   left_join( select(filter(obs_dc, period == 'Development'),-propch_p, -period),
#              by=c('tissue','cell type') ) %>% 
#   group_by(tissue, `cell type`) %>%
#   summarise(pval =  mean(Perms >= unique(`Proportion Change`)),
#             FPR = median(Perms) / unique(`Proportion Change`) ) %>%
#   ungroup()

perm_result_dev %>% filter(pval<0.05)
# tissue `cell type`                              pval    FPR
# <chr>  <fct>                                   <dbl>  <dbl>
# 1 Cortex interneuron                          0.0254    0.596
# 2 Liver  endothelial cell of hepatic sinusoid 0.000300 -5.17 
# 3 Lung   club cell of bronchiole              0.026     0.667
# 4 Muscle endothelial cell                     0        -1.37 
# 5 Muscle T cell                               0.0138    0.5  
# 6 Muscle B cell                               0.0132   -0.333

sum(perm_result_dev$pval<0.05) # 19 ->  6

perm_result_aging = nondc_perms2 %>%
  filter(period == 'Ageing') %>% select(-period) %>%
  rename(Perms = `Proportion Change`) %>% 
  left_join( select(filter(obs_dc, period == 'Ageing'),-propch_p, -period),
             by=c('tissue','cell type') ) %>% 
  group_by(tissue, `cell type`) %>%
  summarise(pval = ifelse(unique(`Proportion Change`)<0,
                          mean(`Perms` <= unique(`Proportion Change`) ),
                          mean(`Perms` >= unique(`Proportion Change`) )),
            FPR = median(`Perms`) / unique(`Proportion Change`) ) %>%
  ungroup()

# perm_result_aging = nondc_perms2 %>%
#   filter(period == 'Ageing') %>% select(-period) %>%
#   rename(Perms = `Proportion Change`) %>% 
#   left_join( select(filter(obs_dc, period == 'Ageing'),-propch_p, -period),
#              by=c('tissue','cell type') ) %>% 
#   group_by(tissue, `cell type`) %>%
#   summarise(pval = mean(`Perms` <= unique(`Proportion Change`) ),
#             FPR = median(`Perms`) / unique(`Proportion Change`) ) %>%
#   ungroup()

perm_result_aging %>% filter(pval<0.05)
# tissue `cell type`                     pval    FPR
# <chr>  <fct>                          <dbl>  <dbl>
#   1 Cortex interneuron                    0.014  0.260
# 2 Lung   non-classical monocyte         0.003  0.218
# 3 Lung   respiratory basal cell         0.003 -0.25 
# 4 Lung   regulatory T cell              0.025 -3    
# 5 Muscle skeletal muscle satellite cell 0.044  0.885

perm_results = list(Development = perm_result_dev,
                    Ageing = perm_result_aging)

perm_result2 = reshape2::melt(perm_results, id.vars = c('tissue','cell type')) %>%
  spread(key = 'variable', value = 'value') %>% 
  rename('period' = L1) %>%
  mutate(period = factor(period, levels = c('Development','Ageing') ) )

perm_result_tab = perm_result2 %>%
  left_join(select(obs_dc, -propch_p), by=c('tissue','cell type','period')) %>%
  rename('DiCo: Proportion Change' = `Proportion Change`) %>%
  left_join(select(obs_nondc, -propch_p), by=c('tissue','cell type','period')) %>%
  rename('non-DiCo: Proportion Change' = `Proportion Change`)

saveRDS(perm_result_tab, './data/other_datasets/scRNA-seq/processed/dico_vs_nondico_prop_test.rds')
write.xlsx(perm_result_tab, './results/SI_tables/TableS15.xlsx')

#########################
perm_dat  = nondc_perms2 %>% 
  rename(Perms = `Proportion Change`) %>% 
  mutate(period = factor(period, levels =c('Development','Ageing') ) )  #%>%  
  # left_join( select(obs_dc, -propch_p), by=c('tissue','cell type','period') ) %>%
  # rename('DiCo'= `Proportion Change`) %>% 
  # left_join( select(obs_nondc, -propch_p), by=c('tissue','cell type','period') ) %>% 
  # rename('Non-DiCo'= `Proportion Change`)

obs_dc = obs_dc %>% 
  mutate(period = factor(period, levels = c('Development','Ageing') ) )

obs_nondc = obs_nondc %>% 
  mutate(period = factor(period, levels = c('Development','Ageing') ) )

perm_result2 = perm_result2 %>%
  mutate(pval = ifelse(pval==0,'<0.001', pval))

#####
### cortex
ord1 = coeffs %>%
  filter(tissue=='Cortex') %>%
  group_by(tissue, `cell type`) %>%
  summarise(m=mean(proportion)) %>%
  arrange(desc(`m`)) %>% pull(`cell type`) %>% as.character()

# sig_ones = perm_result2 %>% 
#   filter(tissue=='Cortex' & pval < 0.15) %>% 
#   select(-pval, -FPR) %>%
#   left_join(nondc_perms2) %>%
#   rename(Perms = `Proportion Change`)
# head(sig_ones)

# label1 = perm_result2 %>% 
#   filter(tissue=='Cortex' & pval < 0.15) %>% 
#   select(-pval, -FPR) %>%
#   mutate(y_pos = 1)

fill1 = perm_result2 %>% 
  filter(tissue=='Cortex') %>% 
  mutate(`cell type` = factor(`cell type`, levels = ord1)) %>%
  mutate(cols = paste(`cell type`,period, sep='_') )  %>%
  mutate(col =  ifelse(pval < 0.05, 'goldenrod1', 'gray50') ) %>%
  mutate(texcol =  ifelse(pval < 0.05, 'darkred', 'gray30') ) %>%
  select(-FPR) 

p1dat = perm_dat %>%
  mutate(period = factor(period, levels = c('Development','Ageing') ) )  %>%
  filter(tissue=='Cortex') %>%
  mutate(`cell type` = factor(`cell type`, levels = ord1)) %>%
  mutate(cols = paste(`cell type`, period, sep='_') ) 
p1 = perm_dat %>%
  mutate(period = factor(period, levels = c('Development','Ageing') ) )  %>%
  filter(tissue=='Cortex') %>%
  mutate(`cell type` = factor(`cell type`, levels = ord1)) %>%
  mutate(cols = paste(`cell type`, period, sep='_') )  %>%
  ggplot(aes(y=`cell type`, x = Perms, fill=cols)) +
  facet_grid(period~tissue) +
  #geom_density_ridges(stat='density', trim=T,scale=5) 
  geom_density_ridges(scale = 1.1,  show.legend = F, rel_min_height=0.01) +
  xlim(-1, 1) +
  geom_vline(xintercept=0, linetype='dashed', color=adjustcolor('gray60', alpha.f = 0.5)) +
  geom_segment(aes(x= `Proportion Change`, xend= `Proportion Change`,y = as.numeric(`cell type`),
                   yend = as.numeric(`cell type`) + 0.9),
               mutate(filter(obs_dc, tissue=='Cortex'), `cell type` = factor(`cell type`, levels = ord1 )),
               inherit.aes = F, colour='darkred' , lineend = 'butt') +
  geom_segment(aes(x= `Proportion Change`, xend= `Proportion Change`,y = as.numeric(`cell type`),
                   yend = as.numeric(`cell type`) + 0.9),
               mutate(filter(obs_nondc, tissue=='Cortex'), `cell type` = factor(`cell type`, levels = ord1 )),
               inherit.aes = F, colour='dodgerblue' , lineend = 'butt', linetype=1) +
  geom_text(data = fill1, aes(x=-1, y = as.numeric(`cell type`) + 0.3, label=paste0('p=', pval)), 
            hjust = 0.5, inherit.aes = F, size=2,
            color= setNames(fill1$texcol, nm = fill1$cols)) +
  scale_fill_manual(values = setNames(fill1$col, nm = fill1$cols) ) +
  ylab('') +
  xlab('Proportion - Age Correlation') +
  theme_bw()
p1

saveRDS(p1dat,'results/source_data/f5/fs2.rds')
# ggsave('./results/SI_figures/Figure_S26.pdf', p1, units = 'cm', width = 16, height = 18, useDingbats = F)
# ggsave('./results/SI_figures/Figure_S26.png', p1, units = 'cm', width = 16, height = 18)

ggsave('./results/figure_supplements/fs5/FS2.pdf', p1, units = 'cm', width = 16, height = 18,
       useDingbats = F)
ggsave('./results/figure_supplements/fs5/FS2.png', p1, units = 'cm', width = 16, height = 18)

#####
### liver
ord2 = coeffs %>%
  filter(tissue=='Liver') %>%
  group_by(tissue, `cell type`) %>%
  summarise(m=mean(proportion)) %>%
  arrange(desc(`m`)) %>% pull(`cell type`) %>% as.character()

fill2 = perm_result2 %>% 
  filter(tissue=='Liver') %>% 
  mutate(`cell type` = factor(`cell type`, levels = ord2)) %>%
  mutate(cols = paste(`cell type`,period, sep='_') )  %>%
  mutate(col =  ifelse(pval < 0.05, 'goldenrod1', 'gray50') ) %>%
  mutate(texcol =  ifelse(pval < 0.05, 'darkred', 'gray30') ) %>%
  select(-FPR) 

p2dat = perm_dat %>%
  mutate(period = factor(period, levels = c('Development','Ageing') ) )  %>%
  filter(tissue=='Liver') %>%
  mutate(`cell type` = factor(`cell type`, levels = ord2)) %>%
  mutate(cols = paste(`cell type`, period, sep='_') ) 
p2 = perm_dat %>%
  mutate(period = factor(period, levels = c('Development','Ageing') ) )  %>%
  filter(tissue=='Liver') %>%
  mutate(`cell type` = factor(`cell type`, levels = ord2)) %>%
  mutate(cols = paste(`cell type`, period, sep='_') )  %>%
  ggplot(aes(y=`cell type`, x = Perms, fill=cols)) +
  facet_grid(period~tissue) +
  #geom_density_ridges(stat='density', trim=T,scale=5) 
  geom_density_ridges(scale = 1.1,  show.legend = F, rel_min_height=0.01) +
  xlim(-1, 1) +
  geom_vline(xintercept=0, linetype='dashed', color=adjustcolor('gray60', alpha.f = 0.5)) +
  geom_segment(aes(x= `Proportion Change`, xend= `Proportion Change`,y = as.numeric(`cell type`),
                   yend = as.numeric(`cell type`) + 0.9),
               mutate(filter(obs_dc, tissue=='Liver'), `cell type` = factor(`cell type`, levels = ord2 )),
               inherit.aes = F, colour='darkred' , lineend = 'butt') +
  geom_segment(aes(x= `Proportion Change`, xend= `Proportion Change`,y = as.numeric(`cell type`),
                   yend = as.numeric(`cell type`) + 0.9),
               mutate(filter(obs_nondc, tissue=='Liver'), `cell type` = factor(`cell type`, levels = ord2 )),
               inherit.aes = F, colour='dodgerblue' , lineend = 'butt', linetype=1) +
  geom_text(data = fill2, aes(x=-1, y = as.numeric(`cell type`) + 0.3, label=paste0('p=', pval)), 
            hjust = 0.5, inherit.aes = F, size=2,
            color= setNames(fill2$texcol, nm = fill2$cols)) +
  scale_fill_manual(values = setNames(fill2$col, nm = fill2$cols) ) +
  ylab('') +
  xlab('Proportion - Age Correlation') +
  theme_bw()
p2

saveRDS(p2dat,'results/source_data/f5/fs3.rds')
# ggsave('./results/SI_figures/Figure_S27.pdf', p2, units = 'cm', width = 16, height = 18, useDingbats = F)
# ggsave('./results/SI_figures/Figure_S27.png', p2, units = 'cm', width = 16, height = 18)

ggsave('./results/figure_supplements/fs5/FS3.pdf', p2, units = 'cm', width = 16, height = 18, 
       useDingbats = F)
ggsave('./results/figure_supplements/fs5/FS3.png', p2, units = 'cm', width = 16, height = 18)

#####
### lung
ord3 = coeffs %>%
  filter(tissue=='Lung') %>%
  group_by(tissue, `cell type`) %>%
  summarise(m=mean(proportion)) %>%
  arrange(desc(`m`)) %>% pull(`cell type`) %>% as.character()

fill3 = perm_result2 %>% 
  filter(tissue=='Lung') %>% 
  mutate(`cell type` = factor(`cell type`, levels = ord3)) %>%
  mutate(cols = paste(`cell type`,period, sep='_') )  %>%
  mutate(col =  ifelse(pval < 0.05, 'goldenrod1', 'gray50') ) %>%
  mutate(texcol =  ifelse(pval < 0.05, 'darkred', 'gray30') ) %>%
  select(-FPR) 

p3dat = perm_dat %>%
  mutate(period = factor(period, levels = c('Development','Ageing') ) )  %>%
  filter(tissue=='Lung') %>%
  mutate(`cell type` = factor(`cell type`, levels = ord3)) %>%
  mutate(cols = paste(`cell type`, period, sep='_') ) 
p3 = perm_dat %>%
  mutate(period = factor(period, levels = c('Development','Ageing') ) )  %>%
  filter(tissue=='Lung') %>%
  mutate(`cell type` = factor(`cell type`, levels = ord3)) %>%
  mutate(cols = paste(`cell type`, period, sep='_') )  %>%
  ggplot(aes(y=`cell type`, x = Perms, fill=cols)) +
  facet_grid(period~tissue) +
  #geom_density_ridges(stat='density', trim=T,scale=5) 
  geom_density_ridges(scale = 1.1,  show.legend = F, rel_min_height=0.01) +
  xlim(-1, 1) +
  geom_vline(xintercept=0, linetype='dashed', color=adjustcolor('gray60', alpha.f = 0.5)) +
  geom_segment(aes(x= `Proportion Change`, xend= `Proportion Change`,y = as.numeric(`cell type`),
                   yend = as.numeric(`cell type`) + 0.9),
               mutate(filter(obs_dc, tissue=='Lung'), `cell type` = factor(`cell type`, levels = ord3 )),
               inherit.aes = F, colour='darkred' , lineend = 'butt') +
  geom_segment(aes(x= `Proportion Change`, xend= `Proportion Change`,y = as.numeric(`cell type`),
                   yend = as.numeric(`cell type`) + 0.9),
               mutate(filter(obs_nondc, tissue=='Lung'), `cell type` = factor(`cell type`, levels = ord3 )),
               inherit.aes = F, colour='dodgerblue' , lineend = 'butt', linetype=1) +
  geom_text(data = fill3, aes(x=-1, y = as.numeric(`cell type`) + 0.3, label=paste0('p=', pval)), 
            hjust = 0.5, inherit.aes = F, size=2,
            color= setNames(fill3$texcol, nm = fill3$cols)) +
  scale_fill_manual(values = setNames(fill3$col, nm = fill3$cols) ) +
  ylab('') +
  xlab('Proportion - Age Correlation') +
  theme_bw()
p3

saveRDS(p3dat,'results/source_data/f5/fs4.rds')
# ggsave('./results/SI_figures/Figure_S28.pdf', p3, units = 'cm', width = 16, height = 18, useDingbats = F)
# ggsave('./results/SI_figures/Figure_S28.png', p3, units = 'cm', width = 16, height = 18)
ggsave('./results/figure_supplements/fs5/FS4.pdf', p3, units = 'cm', width = 16, height = 18, 
       useDingbats = F)
ggsave('./results/figure_supplements/fs5/FS4.png', p3, units = 'cm', width = 16, height = 18)

#####
### muscle
ord4 = coeffs %>%
  filter(tissue=='Muscle') %>%
  group_by(tissue, `cell type`) %>%
  summarise(m=mean(proportion)) %>%
  arrange(desc(`m`)) %>% pull(`cell type`) %>% as.character()

fill4 = perm_result2 %>% 
  filter(tissue=='Muscle') %>% 
  mutate(`cell type` = factor(`cell type`, levels = ord4)) %>%
  mutate(cols = paste(`cell type`,period, sep='_') )  %>%
  mutate(col =  ifelse(pval < 0.05, 'goldenrod1', 'gray50') ) %>%
  mutate(texcol =  ifelse(pval < 0.05, 'darkred', 'gray30') ) %>%
  select(-FPR) 

p4dat = perm_dat %>%
  mutate(period = factor(period, levels = c('Development','Ageing') ) )  %>%
  filter(tissue=='Muscle') %>%
  mutate(`cell type` = factor(`cell type`, levels = ord4)) %>%
  mutate(cols = paste(`cell type`, period, sep='_') )
p4 = perm_dat %>%
  mutate(period = factor(period, levels = c('Development','Ageing') ) )  %>%
  filter(tissue=='Muscle') %>%
  mutate(`cell type` = factor(`cell type`, levels = ord4)) %>%
  mutate(cols = paste(`cell type`, period, sep='_') )  %>%
  ggplot(aes(y=`cell type`, x = Perms, fill=cols)) +
  facet_grid(period~tissue) +
  #geom_density_ridges(stat='density', trim=T,scale=5) 
  geom_density_ridges(scale = 1.1,  show.legend = F, rel_min_height=0.01) +
  xlim(-1, 1) +
  geom_vline(xintercept=0, linetype='dashed', color=adjustcolor('gray60', alpha.f = 0.5)) +
  geom_segment(aes(x= `Proportion Change`, xend= `Proportion Change`,y = as.numeric(`cell type`),
                   yend = as.numeric(`cell type`) + 0.9),
               mutate(filter(obs_dc, tissue=='Muscle'), `cell type` = factor(`cell type`, levels = ord4 )),
               inherit.aes = F, colour='darkred' , lineend = 'butt') +
  geom_segment(aes(x= `Proportion Change`, xend= `Proportion Change`,y = as.numeric(`cell type`),
                   yend = as.numeric(`cell type`) + 0.9),
               mutate(filter(obs_nondc, tissue=='Muscle'), `cell type` = factor(`cell type`, levels = ord4 )),
               inherit.aes = F, colour='dodgerblue' , lineend = 'butt', linetype=1) +
  geom_text(data = fill4, aes(x=-1, y = as.numeric(`cell type`) + 0.3, label=paste0('p=', pval)), 
            hjust = 0.5, inherit.aes = F, size=2,
            color= setNames(fill4$texcol, nm = fill4$cols)) +
  scale_fill_manual(values = setNames(fill4$col, nm = fill4$cols) ) +
  ylab('') +
  xlab('Proportion - Age Correlation') +
  theme_bw()
p4

saveRDS(p4dat,'results/source_data/f5/fs5.rds')

# ggsave('./results/SI_figures/Figure_S29.pdf', p4, units = 'cm', width = 16, height = 18, useDingbats = F)
# ggsave('./results/SI_figures/Figure_S29.png', p4, units = 'cm', width = 16, height = 18)

ggsave('./results/figure_supplements/fs5/FS5.pdf', p4, units = 'cm', width = 16, height = 18, 
       useDingbats = F)
ggsave('./results/figure_supplements/fs5/FS5.png', p4, units = 'cm', width = 16, height = 18)

# perm_plot = ggarrange(p1,p2,p3,p4, ncol=2, nrow=2, common.legend = F, font.label = list(size=8), align = 'h')
# perm_plot
# 
 # ggsave('./results/SI_figures/Figure_S26.pdf', perm_plot, units = 'cm', width = 16, height = 18, 
 #        useDingbats = F)
# ggsave('./results/SI_figures/Figure_S26.png', perm_plot, units = 'cm', width = 16, height = 18)

##########

####################
#################### 
####################
####################

