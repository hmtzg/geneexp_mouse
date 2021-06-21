library(tidyverse)
library(reshape2)
library(openxlsx)
ct.exp  = readRDS('./data/other_datasets/scRNA-seq/processed/celltype_expr.rds')

####### get same genes for all age groups in each tissue, get same cell type for all age groups:
for(x in names(ct.exp$m3) ){
  sgenes = Reduce(intersect, list(rownames(ct.exp$m3[[x]]), 
                                  rownames(ct.exp$m18[[x]]), rownames(ct.exp$m24[[x]])))
  scelltype = Reduce(intersect,list(colnames(ct.exp$m3[[x]]),
                                    colnames(ct.exp$m18[[x]]), colnames(ct.exp$m24[[x]])))
  ct.exp$m3[[x]] = ct.exp$m3[[x]][sgenes, scelltype]
  ct.exp$m18[[x]] = ct.exp$m18[[x]][sgenes, scelltype]
  ct.exp$m24[[x]] = ct.exp$m24[[x]][sgenes, scelltype]
}

### function to calculate pairwise cell type correlations among tissues :
celltype_cors_f = function(ct.exp){
  ct_cors = sapply(names(ct.exp), function(x){
    m3 = sapply(names(ct.exp[[x]]), function(b){
      pairts = sapply(colnames(ct.exp[[x]][[b]]), function(y){
        othr = names(ct.exp[[x]])[!names(ct.exp[[x]])%in%b ]
        pairct = sapply(othr, function(c) {
          sgenes = intersect(rownames(ct.exp[[x]][[b]]), rownames(ct.exp[[x]][[c]]))
          paircor = sapply(colnames(ct.exp[[x]][[c]]), function(a){
            cor(ct.exp[[x]][[b]][sgenes,y], ct.exp[[x]][[c]][sgenes,a], m='s')
          }, simplify = F)
          paircor = data.frame(name = names(paircor), rho = unlist(paircor))
          return( paircor )
        }, simplify = F)
      }, simplify = F)
    })
  }, simplify = F)
}

### calculate pairwise cell type correlations among tissues with all genes:
ct_cors = celltype_cors_f(ct.exp)

saveRDS(ct_cors$m3, './data/other_datasets/scRNA-seq/processed/celltype_3m_pairwise_cors_allgenes.rds')
saveRDS(ct_cors, './data/other_datasets/scRNA-seq/processed/celltype_pairwise_cors_allgenes.rds')

##############
##############
##############
################# pairwise correlations of cell types among tissues within each age group (all genes):
#################

### all cell type correlations within each age group
cell_type_cors = reshape2::melt(ct_cors,c('name', 'rho')) %>%
  set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue','age gr') %>%
  mutate(`age gr`= factor(`age gr`, levels=c('m3', 'm18', 'm24'))) %>%
  mutate(age= as.numeric(gsub('[a-z]','', `age gr`)) )

##### maximally correlated cell types in 3m group
maxcors_3m = reshape2::melt(ct_cors$m3, id.vars = c('name', 'rho')) %>% 
  set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
  group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  top_n(n=1, wt=rho) %>% ungroup()

##### choose maximally correlated cell types from 3m group and add other age groups
allmaxcors = cell_type_cors %>%
  right_join(select(maxcors_3m, -rho),
             by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )

##### maximally correlated cell types correlation change with age
maxcorchange = allmaxcors %>%
  select(-`age gr`) %>%
  group_by(`1st tissue`,`1st cell type`) %>% 
  summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
  ungroup() %>%
  arrange(-`rho ch`)

maxcors_table_sx = maxcors_3m %>% 
  left_join(maxcorchange) %>%
  relocate(`1st tissue`, `1st cell type`,`rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

#################
#################
##### minimally correlated cell types in 3m group
mincors_3m = reshape2::melt(ct_cors$m3, id.vars = c('name', 'rho')) %>% 
  set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
  group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  top_n(n=1, wt = -rho) %>% ungroup()

##### choose minimally correlated cell types from 3m group and add other age groups
allmincors = cell_type_cors %>%
  right_join(select(mincors_3m, -rho),
             by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )

##### minimally correlated cell types correlation change with age
mincorchange = allmincors %>%
  select(-`age gr`) %>%
  group_by(`1st tissue`,`1st cell type`) %>%
  summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
  ungroup() %>%
  arrange(-`rho ch`)

#####
mincors_table_sx = mincors_3m %>% 
  left_join(mincorchange) %>%
  relocate(`1st tissue`, `1st cell type`, `rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

##############
############## 
##############
#################
#################
#################
################# repeat with DC genes
#################

# get same genes for all age groups in each tissue, get same cell type for all age groups:
dgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
dcgenes = names(which(dgenes<0))

## convert ENS DC genes to uniprot:
convert = readRDS('./data/uniprot_to_ENS.rds')
choosedc = convert[convert[,2]%in%dcgenes, 1]

ct.exp.dc = list()
for(x in names(ct.exp$m3) ){
  sgenes = intersect(choosedc, rownames(ct.exp$m3[[x]]))
  ct.exp.dc$m3[[x]] = ct.exp$m3[[x]][sgenes, ]
  ct.exp.dc$m18[[x]] = ct.exp$m18[[x]][sgenes,]
  ct.exp.dc$m24[[x]] = ct.exp$m24[[x]][sgenes,]
}

ct_cors_dc = celltype_cors_f(ct.exp.dc)

saveRDS(ct_cors_dc, './data/other_datasets/scRNA-seq/processed/celltype_pairwise_cors_dcgenes.rds')

### all cell type correlations within each age group
cell_type_cors_dc = reshape2::melt(ct_cors_dc, c('name', 'rho')) %>%
  set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue','age gr') %>%
  mutate(`age gr`= factor(`age gr`, levels=c('m3', 'm18', 'm24'))) %>%
  mutate(age= as.numeric(gsub('[a-z]','', `age gr`)) )

##### maximally correlated cell types in 3m group
maxcors_3m_dc = reshape2::melt(ct_cors_dc$m3, id.vars = c('name', 'rho')) %>% 
  set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
  group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  top_n(n=1, wt=rho) %>% ungroup()

##### choose maximally correlated cell types from 3m group and add other age groups
allmaxcors_dc = cell_type_cors_dc %>%
  right_join(select(maxcors_3m_dc, -rho),
             by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )

##### maximally correlated cell types correlation change with age
maxcorchange_dc = allmaxcors_dc %>%
  select(-`age gr`) %>%
  group_by(`1st tissue`,`1st cell type`) %>% 
  summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
  ungroup() %>%
  arrange(-`rho ch`)

#####
maxcors_table_sx_dc = maxcors_3m_dc %>% 
  left_join(maxcorchange_dc) %>%
  relocate(`1st tissue`, `1st cell type`,`rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

#################
#################
#################
##### minimally correlated cell types in 3m group among dc genes

##### minimally correlated cell types in 3m group
mincors_3m_dc = reshape2::melt(ct_cors_dc$m3, id.vars = c('name', 'rho')) %>% 
  set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
  group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  top_n(n=1, wt = -rho) %>% ungroup()

##### choose minimally correlated cell types from 3m group and add other age groups
allmincors_dc = cell_type_cors_dc %>%
  right_join(select(mincors_3m_dc, -rho),
             by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )

##### minimally correlated cell types correlation change with age
mincorchange_dc = allmincors_dc %>%
  select(-`age gr`) %>%
  group_by(`1st tissue`,`1st cell type`) %>%
  summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
  ungroup() %>%
  arrange(-`rho ch`)

#####
mincors_table_sx_dc = mincors_3m_dc %>% 
  left_join(mincorchange_dc) %>%
  relocate(`1st tissue`, `1st cell type`, `rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

############## 
##############
##############
##############
############## 
##############
##############
##############
#################
################# repeat with non-DC genes
#################

# get same genes for all age groups in each tissue, get same cell type for all age groups:
ct.exp.nondc = list()
for(x in names(ct.exp$m3) ){
  dgenes = setdiff(rownames(ct.exp$m3[[x]]), choosedc)
  ct.exp.nondc$m3[[x]] = ct.exp$m3[[x]][dgenes, ]
  ct.exp.nondc$m18[[x]] = ct.exp$m18[[x]][dgenes,]
  ct.exp.nondc$m24[[x]] = ct.exp$m24[[x]][dgenes,]
}

ct_cors_nondc = celltype_cors_f(ct.exp.nondc)

saveRDS(ct_cors_nondc, './data/other_datasets/scRNA-seq/processed/celltype_pairwise_cors_nondcgenes.rds')

### all cell type correlations within each age group
cell_type_cors_nondc = reshape2::melt(ct_cors_nondc, c('name', 'rho')) %>%
  set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue','age gr') %>%
  mutate(`age gr`= factor(`age gr`, levels=c('m3', 'm18', 'm24'))) %>%
  mutate(age= as.numeric(gsub('[a-z]','', `age gr`)) )

##### maximally correlated cell types in 3m group
maxcors_3m_nondc = reshape2::melt(ct_cors_nondc$m3, id.vars = c('name', 'rho')) %>% 
  set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
  group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  top_n(n=1, wt=rho) %>% ungroup()

##### choose maximally correlated cell types from 3m group and add other age groups
allmaxcors_nondc = cell_type_cors_nondc %>%
  right_join(select(maxcors_3m_nondc, -rho),
             by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )

##### maximally correlated cell types correlation change with age
maxcorchange_nondc = allmaxcors_nondc %>%
  select(-`age gr`) %>%
  group_by(`1st tissue`,`1st cell type`) %>% 
  summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
  ungroup() %>%
  arrange(-`rho ch`)

#####
maxcors_table_sx_nondc = maxcors_3m_nondc %>% 
  left_join(maxcorchange_nondc) %>%
  relocate(`1st tissue`, `1st cell type`,`rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

#################
#################
#################
##### minimally correlated cell types in 3m group among dc genes

##### minimally correlated cell types in 3m group
mincors_3m_nondc = reshape2::melt(ct_cors_nondc$m3, id.vars = c('name', 'rho')) %>% 
  set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
  group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  top_n(n=1, wt = -rho) %>% ungroup()

##### choose minimally correlated cell types from 3m group and add other age groups
allmincors_nondc = cell_type_cors_nondc %>%
  right_join(select(mincors_3m_nondc, -rho),
             by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )

##### minimally correlated cell types correlation change with age
mincorchange_nondc = allmincors_nondc %>%
  select(-`age gr`) %>%
  group_by(`1st tissue`,`1st cell type`) %>%
  summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
  ungroup() %>%
  arrange(-`rho ch`)

#####
mincors_table_sx_nondc = mincors_3m_nondc %>% 
  left_join(mincorchange_nondc) %>%
  relocate(`1st tissue`, `1st cell type`, `rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

saveRDS(list(allgenes = cell_type_cors, dcgenes=cell_type_cors_dc, nondcgenes = cell_type_cors_nondc),
        './data/other_datasets/scRNA-seq/processed/celltype_pairwise_cors.rds')

##########################

######## write max cors to table Sx (all, dc, non-dc)
maxcors_table = list(allgenes = maxcors_table_sx, DiCo_genes = maxcors_table_sx_dc,
                     nonDiCo_genes  = maxcors_table_sx_nondc)
write.xlsx(maxcors_table, './results/SI_tables/TableS13.xlsx')

######## write min cors to table S8 (all, dc, non-dc)
mincors_table = list(allgenes = mincors_table_sx, DiCo_genes = mincors_table_sx_dc,
                     nonDiCo_genes = mincors_table_sx_nondc)
write.xlsx(mincors_table, './data/TableS14.xlsx')

##############
##############
############## test max/min cor changes with permutation
##############

table(sign(maxcorchange$`rho ch`))
table(sign(maxcorchange_dc$`rho ch`))
table(sign(maxcorchange_nondc$`rho ch`))

table(sign(mincorchange$`rho ch`))
table(sign(mincorchange_dc$`rho ch`))
table(sign(mincorchange_nondc$`rho ch`))

##############
############## Permutation test to test DC and non-DC cor change differences
##############
##############

nondc_perms = list()

perm=1
x='lung'

nondc_perms_max = list()
nondc_perms_min = list()
for(perm in 1:1000){
  print(perm/1000)
  ct.exp.nondc.perm = list()
  for(x in names(ct.exp$m3) ){
    nondcgenes = setdiff(rownames(ct.exp$m3[[x]]), choosedc)
    chooseg = sample(nondcgenes, size = length(choosedc) ) 
    ct.exp.nondc.perm$m3[[x]] = ct.exp$m3[[x]][chooseg, ]
    ct.exp.nondc.perm$m18[[x]] = ct.exp$m18[[x]][chooseg,]
    ct.exp.nondc.perm$m24[[x]] = ct.exp$m24[[x]][chooseg,]
  }
  names(ct.exp.nondc.perm) = names(ct.exp)

  ct_cors_nondc_perm = celltype_cors_f(ct.exp.nondc.perm)
  
  ### all cell type correlations within each age group
  cell_type_cors_nondc_perm = reshape2::melt(ct_cors_nondc_perm, c('name', 'rho')) %>%
    set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue','age gr') %>%
    mutate(`age gr`= factor(`age gr`, levels=c('m3', 'm18', 'm24'))) %>%
    mutate(age= as.numeric(gsub('[a-z]','', `age gr`)) )
  
  ##### maximally correlated cell types in 3m group
  maxcors_3m_nondc_perm = reshape2::melt(ct_cors_nondc_perm$m3, id.vars = c('name', 'rho')) %>% 
    set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
    group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
    top_n(n=1, wt=rho) %>% ungroup()
  
  ##### choose maximally correlated cell types from 3m group and add other age groups
  allmaxcors_nondc_perm = cell_type_cors_nondc_perm %>%
    right_join(select(maxcors_3m_nondc_perm, -rho),
               by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )
  
  ##### maximally correlated cell types correlation change with age
  maxcorchange_nondc_perm = allmaxcors_nondc_perm %>%
    select(-`age gr`) %>%
    group_by(`1st tissue`,`1st cell type`) %>% 
    summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
    ungroup() %>%
    arrange(-`rho ch`)
  
  ##### minimally correlated cell types in 3m group
  mincors_3m_nondc_perm = reshape2::melt(ct_cors_nondc_perm$m3, id.vars = c('name', 'rho')) %>% 
    set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
    group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
    top_n(n=1, wt = -rho) %>% ungroup()
  
  ##### choose minimally correlated cell types from 3m group and add other age groups
  allmincors_nondc_perm = cell_type_cors_nondc_perm %>%
    right_join(select(mincors_3m_nondc_perm, -rho),
               by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )
  
  ##### minimally correlated cell types correlation change with age
  mincorchange_nondc_perm = allmincors_nondc_perm %>%
    select(-`age gr`) %>%
    group_by(`1st tissue`,`1st cell type`) %>%
    summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
    ungroup() %>%
    arrange(-`rho ch`)
  
  nondc_perms_max[[perm]] = maxcorchange_nondc_perm
  nondc_perms_min[[perm]] = mincorchange_nondc_perm
}
names(nondc_perms_max) = paste0('perm',1:1000)
names(nondc_perms_min) = paste0('perm',1:1000)

saveRDS(nondc_perms_max, './data/other_datasets/scRNA-seq/processed/celltype_pairwise_maxcors_nondc_permds.rds')
saveRDS(nondc_perms_min, './data/other_datasets/scRNA-seq/processed/celltype_pairwise_mincors_nondc_permds.rds')

nondc_perms_max$perm1$`rho ch`
maxcorchange_nondc$`rho ch`

table(sign(maxcorchange_nondc$`rho ch`))['-1']
table(sign(maxcorchange_dc$`rho ch`))['-1']

mean(sapply(nondc_perms_max, function(x) table(sign(x$`rho ch`))['-1'] ) >= 35)
# 0.895
hist(sapply(nondc_perms_max, function(x) table(sign(x$`rho ch`))['-1'] ))
abline(v=35)

table(sign(mincorchange_nondc$`rho ch`))['1']
table(sign(mincorchange_dc$`rho ch`))['1']
mean(sapply(nondc_perms_min, function(x) table(sign(x$`rho ch`))['1'] )>= 44)
# 0.794
hist(sapply(nondc_perms_min, function(x) table(sign(x$`rho ch`))['1'] ))
abline(v=44)

maxcorchange_nondc








