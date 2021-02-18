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

maxcors_table_s7 = maxcors_3m %>% 
  left_join(maxcorchange) %>%
  relocate(`1st tissue`, `1st cell type`,`rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

write.xlsx(maxcors_table_s7, file='./data/Table_S7_allgenes.xlsx', row.names=F)

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
mincors_table_s8 = mincors_3m %>% 
  left_join(mincorchange) %>%
  relocate(`1st tissue`, `1st cell type`, `rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

write.xlsx(mincors_table_s8, file='./data/Table_S8_allgenes.xlsx', row.names=F)

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
maxcors_table_s7_dc = maxcors_3m_dc %>% 
  left_join(maxcorchange_dc) %>%
  relocate(`1st tissue`, `1st cell type`,`rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

write.xlsx(maxcors_table_s7_dc, file='./data/Table_S7_dc.xlsx', row.names=F)

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
mincors_table_s8_dc = mincors_3m_dc %>% 
  left_join(mincorchange_dc) %>%
  relocate(`1st tissue`, `1st cell type`, `rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

write.xlsx(mincors_table_s8_dc, file='./data/Table_S8_dc.xlsx', row.names=F)

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
maxcors_table_s7_nondc = maxcors_3m_nondc %>% 
  left_join(maxcorchange_nondc) %>%
  relocate(`1st tissue`, `1st cell type`,`rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

write.xlsx(maxcors_table_s7_nondc, file='./data/Table_S7_nondc.xlsx', row.names=F)

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
mincors_table_s8_nondc = mincors_3m_nondc %>% 
  left_join(mincorchange_nondc) %>%
  relocate(`1st tissue`, `1st cell type`, `rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B'))

write.xlsx(mincors_table_s8_nondc, file='./data/Table_S8_nondc.xlsx', row.names=F)

saveRDS(list(allgenes = cell_type_cors, dcgenes=cell_type_cors_dc, nondcgenes = cell_type_cors_nondc),
        './data/other_datasets/scRNA-seq/processed/celltype_pairwise_cors.rds')

##########################

######## write max cors to table S7 (all, dc, non-dc)
maxcors_table = list(allgenes = maxcors_table_s7, DiCo_genes = maxcors_table_s7_dc,
                     nonDiCo_genes  = maxcors_table_s7_nondc)
write.xlsx(maxcors_table, './data/Table_S7.xlsx')

######## write min cors to table S8 (all, dc, non-dc)
mincors_table = list(allgenes = mincors_table_s8, DiCo_genes = mincors_table_s8_dc,
                     nonDiCo_genes = mincors_table_s8_nondc)
write.xlsx(mincors_table, './data/Table_S8.xlsx')

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






















