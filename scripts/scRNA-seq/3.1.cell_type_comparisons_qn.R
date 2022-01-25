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

xg = Reduce(intersect, lapply(ct.exp$m3, rownames ))
for(i in 1:3){
  for(j in 1:4){
    ct.exp[[i]][[j]] = ct.exp[[i]][[j]][rownames(ct.exp[[i]][[j]])%in%xg,] 
  }
}
sapply(ct.exp, function(x)sapply(x,ncol))
# order: names(ct.exp$m3)
expr = cbind(do.call(cbind, ct.exp$m3 ),
      do.call(cbind, ct.exp$m18 ),
      do.call(cbind, ct.exp$m24 ))
library(preprocessCore)
exprqn = normalize.quantiles(expr)
dimnames(exprqn) = dimnames(expr)

ct.expqn = list('m3' = list('lung' = exprqn[,1:23],
                            'liver'=exprqn[,24:33],
                            'muscle'=exprqn[,34:39],
                            'brain'=exprqn[,40:54]),
                'm18' = list('lung' = exprqn[,1:23+54],
                             'liver'=exprqn[,24:33+54],
                             'muscle'=exprqn[,34:39+54],
                             'brain'=exprqn[,40:54+54]),
                'm24' = list('lung' = exprqn[,1:23+54*2],
                             'liver'=exprqn[,24:33+54*2],
                             'muscle'=exprqn[,34:39+54*2],
                             'brain'=exprqn[,40:54]+54*2))
sapply(ct.expqn, function(x)sapply(x,ncol))
ct.exp = ct.expqn

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

# saveRDS(ct_cors$m3, './data/other_datasets/scRNA-seq/processed/celltype_3m_pairwise_cors_allgenes.rds')
# saveRDS(ct_cors, './data/other_datasets/scRNA-seq/processed/celltype_pairwise_cors_allgenes.rds')

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
  summarise(`rho ch` = round(cor(rho,age, m='s'),2),
            `rho p` = cor.test(rho,age, m='s')$p.val) %>%
  mutate(fdr = p.adjust(`rho p`,method = 'BH')) %>%
  ungroup() %>%
  arrange(-`rho ch`)

maxcorchange %>% filter(fdr<0.1)

table(sign(maxcorchange$`rho ch`))
# -1  0  1 
# 37  3 14 
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
  summarise(`rho ch` = round(cor(rho,age, m='s'),2),
            `rho p` = cor.test(rho,age, m='s')$p.val) %>%
  mutate(fdr = p.adjust(`rho p`,method = 'BH')) %>%
  ungroup() %>%
  arrange(-`rho ch`)

mincorchange %>% filter(fdr<0.5)
table(sign(mincorchange$`rho ch`))
# -1  1 
# 3 51 
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

# saveRDS(ct_cors_dc, './data/other_datasets/scRNA-seq/processed/celltype_pairwise_cors_dcgenes.rds')

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
  summarise(`rho ch` = round(cor(rho,age, m='s'),2),
            `rho p` = cor.test(rho,age, m='s')$p.val) %>%
  mutate(fdr = p.adjust(`rho p`,method = 'BH')) %>%
  ungroup() %>%
  arrange(-`rho ch`)

maxcorchange_dc %>% filter(fdr<0.1)
table(sign(maxcorchange_dc$`rho ch`))
# -1  0  1 
# 36  2 16 
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
  summarise(`rho ch` = round(cor(rho,age, m='s'),2),
            `rho p` = cor.test(rho,age, m='s')$p.val) %>%
  mutate(fdr = p.adjust(`rho p`,method = 'BH')) %>%
  ungroup() %>%
  arrange(-`rho ch`)

mincorchange_dc %>% filter(fdr<0.1)
table(sign(mincorchange_dc$`rho ch`))
# -1  0  1 
# 2  1 51 
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

# saveRDS(ct_cors_nondc, './data/other_datasets/scRNA-seq/processed/celltype_pairwise_cors_nondcgenes.rds')

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
  summarise(`rho ch` = round(cor(rho,age, m='s'),2),
            `rho p` = cor.test(rho,age, m='s')$p.val) %>%
  mutate(fdr = p.adjust(`rho p`,method = 'BH')) %>%
  ungroup() %>%
  arrange(-`rho ch`)

maxcorchange_nondc %>% filter(fdr<0.1)
table(sign(maxcorchange_nondc$`rho ch`))
# -1  0  1 
# 36  3 15 
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
  summarise(`rho ch` = round(cor(rho,age, m='s'),2),
            `rho p` = cor.test(rho,age, m='s')$p.val) %>%
  mutate(fdr = p.adjust(`rho p`,method = 'BH')) %>%
  ungroup() %>%
  arrange(-`rho ch`)

mincorchange_nondc %>% filter(fdr<0.1)
table(sign(mincorchange_nondc$`rho ch`))
# -1  0  1 
# 2  1 51 
#####
mincors_table_sx_nondc = mincors_3m_nondc %>% 
  left_join(mincorchange_nondc) %>%
  relocate(`1st tissue`, `1st cell type`, `rho ch`,
           `2nd tissue`, `2nd cell type`, rho) %>%
  set_names(c('TissueA','Cell-typeA', 'Similarity Change with Age',
              'TissueB','Cell-typeB', 'Similarity to B','Change in Similarity','fdr'))

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

table(sign(maxcorchange$`rho ch`)) # 36 dec
table(sign(maxcorchange_dc$`rho ch`)) # 30 dec
table(sign(maxcorchange_nondc$`rho ch`)) # 35 dec

table(sign(mincorchange$`rho ch`)) # 45 inc
table(sign(mincorchange_dc$`rho ch`)) # 47 inc
table(sign(mincorchange_nondc$`rho ch`)) # 44 inc

##############
############## Permutation test to test DC and non-DC cor change differences
##############
##############

#nondc_perms = list()

# perm=1
# x='lung'

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

maxcorchange_nondc$`rho ch`

table(sign(maxcorchange_nondc$`rho ch`))['-1'] #  35
table(sign(maxcorchange$`rho ch`))['-1'] #  36
table(sign(maxcorchange_dc$`rho ch`))['-1'] # 30 obs

mean(sapply(nondc_perms_max, function(x) table(sign(x$`rho ch`))['-1'] ) >= 30)
# 0.895 -> 1
hist(sapply(nondc_perms_max, function(x) table(sign(x$`rho ch`))['-1'] ))
abline(v=36)

table(sign(mincorchange_nondc$`rho ch`))['1'] # 44
table(sign(mincorchange_dc$`rho ch`))['1'] # 47 obs
table(sign(mincorchange$`rho ch`))['1'] # 45 

mean(sapply(nondc_perms_min, function(x) table(sign(x$`rho ch`))['1'] )>= 47)
# 0.794 -> 0.661
hist(sapply(nondc_perms_min, function(x) table(sign(x$`rho ch`))['1'] ))
abline(v=45)

########
######## repeat permutation with all genes:
########

allg_perms_max = list()
allg_perms_min = list()
for(perm in 1:1000){
  print(perm/1000)
  #ct.exp.nondc.perm = list()
  ct.exp.allg.perm = list()
  for(x in names(ct.exp$m3) ){
    #nondcgenes = setdiff(rownames(ct.exp$m3[[x]]), choosedc)
    allgenes = rownames(ct.exp$m3[[x]])
    chooseg = sample(allgenes, size = length(choosedc) ) 
    ct.exp.allg.perm$m3[[x]] = ct.exp$m3[[x]][chooseg, ]
    ct.exp.allg.perm$m18[[x]] = ct.exp$m18[[x]][chooseg,]
    ct.exp.allg.perm$m24[[x]] = ct.exp$m24[[x]][chooseg,]
  }
  names(ct.exp.allg.perm) = names(ct.exp)
  
  #ct_cors_nondc_perm = celltype_cors_f(ct.exp.allg.perm)
  ct_cors_allg_perm = celltype_cors_f(ct.exp.allg.perm)
  
  ### all cell type correlations within each age group
  # cell_type_cors_nondc_perm = reshape2::melt(ct_cors_allg_perm, c('name', 'rho')) %>%
  #   set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue','age gr') %>%
  #   mutate(`age gr`= factor(`age gr`, levels=c('m3', 'm18', 'm24'))) %>%
  #   mutate(age= as.numeric(gsub('[a-z]','', `age gr`)) )
  cell_type_cors_allg_perm = reshape2::melt(ct_cors_allg_perm, c('name', 'rho')) %>%
    set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue','age gr') %>%
    mutate(`age gr`= factor(`age gr`, levels=c('m3', 'm18', 'm24'))) %>%
    mutate(age= as.numeric(gsub('[a-z]','', `age gr`)) )
  
  
  ##### maximally correlated cell types in 3m group
  # maxcors_3m_nondc_perm = reshape2::melt(ct_cors_allg_perm$m3, id.vars = c('name', 'rho')) %>% 
  #   set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
  #   group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  #   top_n(n=1, wt=rho) %>% ungroup()
  maxcors_3m_allg_perm = reshape2::melt(ct_cors_allg_perm$m3, id.vars = c('name', 'rho')) %>% 
    set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
    group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
    top_n(n=1, wt=rho) %>% ungroup()
  
  ##### choose maximally correlated cell types from 3m group and add other age groups
  # allmaxcors_nondc_perm = cell_type_cors_allg_perm %>%
  #   right_join(select(maxcors_3m_allg_perm, -rho),
  #              by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )
  allmaxcors_allg_perm = cell_type_cors_allg_perm %>%
    right_join(select(maxcors_3m_allg_perm, -rho),
               by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )
  
  ##### maximally correlated cell types correlation change with age
  # maxcorchange_nondc_perm = allmaxcors_allg_perm %>%
  #   select(-`age gr`) %>%
  #   group_by(`1st tissue`,`1st cell type`) %>% 
  #   summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
  #   ungroup() %>%
  #   arrange(-`rho ch`)
  maxcorchange_allg_perm = allmaxcors_allg_perm %>%
    select(-`age gr`) %>%
    group_by(`1st tissue`,`1st cell type`) %>% 
    summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
    ungroup() %>%
    arrange(-`rho ch`)
  
  ##### minimally correlated cell types in 3m group
  # mincors_3m_nondc_perm = reshape2::melt(ct_cors_allg_perm$m3, id.vars = c('name', 'rho')) %>% 
  #   set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
  #   group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  #   top_n(n=1, wt = -rho) %>% ungroup()
  mincors_3m_allg_perm = reshape2::melt(ct_cors_allg_perm$m3, id.vars = c('name', 'rho')) %>% 
    set_names('2nd cell type', 'rho', '2nd tissue','1st cell type', '1st tissue') %>%
    group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
    top_n(n=1, wt = -rho) %>% ungroup()
  
  ##### choose minimally correlated cell types from 3m group and add other age groups
  # allmincors_nondc_perm = cell_type_cors_allg_perm %>%
  #   right_join(select(mincors_3m_allg_perm, -rho),
  #              by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )
  allmincors_allg_perm = cell_type_cors_allg_perm %>%
    right_join(select(mincors_3m_allg_perm, -rho),
               by= c("2nd cell type", "2nd tissue", "1st cell type", "1st tissue") )
  
  ##### minimally correlated cell types correlation change with age
  # mincorchange_nondc_perm = allmincors_allg_perm %>%
  #   select(-`age gr`) %>%
  #   group_by(`1st tissue`,`1st cell type`) %>%
  #   summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
  #   ungroup() %>%
  #   arrange(-`rho ch`)
  mincorchange_allg_perm = allmincors_allg_perm %>%
    select(-`age gr`) %>%
    group_by(`1st tissue`,`1st cell type`) %>%
    summarise(`rho ch` = round(cor(rho,age, m='s'),2)) %>%
    ungroup() %>%
    arrange(-`rho ch`)
  
  allg_perms_max[[perm]] = maxcorchange_allg_perm
  allg_perms_min[[perm]] = mincorchange_allg_perm
}
names(allg_perms_max) = paste0('perm',1:1000)
names(allg_perms_min) = paste0('perm',1:1000)

saveRDS(allg_perms_max, './data/other_datasets/scRNA-seq/processed/celltype_pairwise_maxcors_allg_permds.rds')
saveRDS(allg_perms_min, './data/other_datasets/scRNA-seq/processed/celltype_pairwise_mincors_allg_permds.rds')

table(sign(maxcorchange_nondc$`rho ch`))['-1'] #  35 dec
table(sign(maxcorchange$`rho ch`))['-1'] # 36 dec
table(sign(maxcorchange_dc$`rho ch`))['-1'] # 30 obs

mean(sapply(allg_perms_max, function(x) table(sign(x$`rho ch`))['-1'] ) >= 30)
# 0.895 -> 0.99
hist(sapply(allg_perms_max, function(x) table(sign(x$`rho ch`))['-1'] ), xlim = c(29,43))
abline(v=30)

table(sign(mincorchange_nondc$`rho ch`))['1'] # 44
table(sign(mincorchange$`rho ch`))['1'] # 45 
table(sign(mincorchange_dc$`rho ch`))['1'] # 47 obs

mean(sapply(allg_perms_min, function(x) table(sign(x$`rho ch`))['1'] )>= 47)
# 0.456
hist(sapply(allg_perms_min, function(x) table(sign(x$`rho ch`))['1'] ))
abline(v=45)

