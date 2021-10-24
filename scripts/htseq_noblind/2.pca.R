library(tidyverse)
library(openxlsx)

exp = readRDS('./data/htseq/expr_mat_no_blind.rds')
age = readRDS('./data/processed/raw/ages.rds')
ts.ord = readRDS('./data/processed/raw/tissue_ids.rds')
sample_info = readRDS('./data/processed/tidy/sample_info.rds')
samp_id = names(ts.ord)
## PCA all regions with expression data:
pc = prcomp(t(exp),scale = T)
# cor(pc$x[,1], age, m="s") # 0.23
# cor(pc$x[,2], age, m="s") # 0.27
# cor(pc$x[,3], age, m="s") # 0.019
# cor.test(pc$x[,4], age, m="s") # 0.60

agedev = age[age<93]
expdev = exp[, age<93]
pcdev = prcomp(t(expdev), scale = T)

ageag = age[age>=93]
expag = exp[, age>=93]
pcag = prcomp(t(expag), scale = T)

exp.sc = cbind(t(scale(t(exp[,ts.ord=='Cortex']))), t(scale(t(exp[,ts.ord=="Lung"]))),
               t(scale(t(exp[,ts.ord=="Liver"]))), t(scale(t(exp[,ts.ord=="Muscle"]))))
# gene: ENSMUSG00000070323 has same value in individuals of cortex, so remove it
exp.sc = exp.sc[complete.cases(exp.sc),]
pc.sc = prcomp(t(exp.sc), scale = T)

exp.sc.dev = cbind(t(scale(t(exp[,ts.ord=="Cortex" & age<93]))),
                   t(scale(t(exp[,ts.ord=="Lung" & age<93]))),
                   t(scale(t(exp[,ts.ord=="Liver" & age<93]))),
                   t(scale(t(exp[,ts.ord=="Muscle" & age<93]))))
exp.sc.dev = exp.sc.dev[complete.cases(exp.sc.dev),]
pc.sc.dev = prcomp(t(exp.sc.dev),scale = T)

exp.sc.ag = cbind(t(scale(t(exp[,ts.ord=="Cortex" & age>=93]))),
                  t(scale(t(exp[,ts.ord=="Lung" & age>=93]))),
                  t(scale(t(exp[,ts.ord=="Liver" & age>=93]))),
                  t(scale(t(exp[,ts.ord=="Muscle" & age>=93]))))
exp.sc.ag = exp.sc.ag[complete.cases(exp.sc.ag),]
pc.sc.ag = prcomp(t(exp.sc.ag),scale = T)

pca_data_all = list(pc.sc.ag,
                    pc,
                    pc.sc,
                    pcag,
                    pcdev,
                    pc.sc.dev)
names(pca_data_all) = c('aging_notissue','all_raw','all_notissue','aging_raw','development_raw',
                        'development_notissue')

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

pca_data = sapply(pca_data_all, function(x) summary(x)$imp[2,1:4]) %>%
  reshape2::melt() %>%
  set_names(c('PC','type','varExp')) %>%
  separate(type,into=c('period','type')) %>%
  right_join(pca_data) 

saveRDS(pca_data,'./data/htseq/no_blind/pca_data.rds')

# pc_age_cors = pca_data %>%
#   left_join(select(sample_info,-ind_id,-log2age)) %>%
#   group_by(period,type,PC, tissue) %>%
#   summarise(rho = cor.test(age,value,m='s')$est,
#             p = cor.test(age,value,m='s')$p.val)
# write.xlsx(pc_age_cors, './data/Table_S0.xlsx')

pc_age_cors  = pca_data %>%
  rename(`Scale period` = period) %>%
  left_join(select(sample_info,-ind_id,-log2age)) %>%
  mutate(`age period` = ifelse(age< 90,'Dev','Aging') ) %>%
  group_by(`Scale period`, type, PC, tissue, `age period`) %>%
  summarise(rho = cor.test(age,value,m='s')$est,
            p = cor.test(age,value,m='s')$p.val)

pc_age_cors %>%
  filter(PC=='PC4' & `Scale period` == 'all' & type =='raw') %>%
  group_by(`age period`) %>%
  summarise(range(rho))
# `age period` `range(rho)`
# <chr>               <dbl>
# 1 Aging              0.0240
# 2 Aging              0.882 
# 3 Dev               -0.991 
# 4 Dev               -0.631 

#################################
################################# does PC1,2,3 separate samples according to tissues? check with Anova
################################# and eclidean distance
#################################

aov_dat = pca_data %>%
  filter(type == 'raw', period == 'all', PC%in%c('PC1','PC2','PC3')) %>% 
  select(-period, -type, -varExp) %>%
  left_join(select(sample_info, -ind_id,  -log2age), by = c('sample_id'))

# summary(aov(value~age*tissue*PC, data=aov_dat)) # drop age
# summary(aov(value~tissue*PC, data=aov_dat)) # drop PC
# summary(aov(value~tissue*age, data=aov_dat)) # drop age
# summary(aov(value~PC*age, data=aov_dat)) # drop age and PC
# summary(aov(value~tissue, data=aov_dat)) # drop both
# ### for each PC between tissue and age
# summary(aov(value~tissue*age, data=filter(aov_dat,PC=='PC1' ) )) # mostly tissue
# summary(aov(value~tissue*age, data=filter(aov_dat,PC=='PC2' ) )) # tissue
# summary(aov(value~tissue*age, data=filter(aov_dat,PC=='PC3' ) )) # tissue
# summary(aov(value~tissue*age, data=filter(aov_dat,PC=='PC4' ) ))# age

for(i in unique(aov_dat$PC)){
  print(c(i,summary(aov(value~tissue, data = aov_dat[aov_dat$PC==i,])),
          TukeyHSD(aov(value~tissue, data = aov_dat[aov_dat$PC==i,]))))
}

###### function to calculate mean pairwise distance among rows:
# pwise_distm_f = function(df){ 
#   rownames(df) = NULL
#   stopifnot(ncol(df) == 3)
#   colnames(df) = c('x','y','id')
#   pwise_distm = df %>%
#     column_to_rownames('id') %>%
#     dist() %>%
#     as.matrix() %>%
#     as.data.frame() %>%
#     rownames_to_column(var = 'id.x') %>%
#     gather(key = 'id.y', value=dist, -id.x) %>% 
#     filter(id.x < id.y) %>%
#     summarise(mean(dist))
#   pwise_distm = unname(pwise_distm)
#   return(pwise_distm)
# }
## mean pairwise distance of tissues in PC3-PC4 space:
# dist_dat = pca_data %>%
#   filter(PC %in% c('PC3', 'PC4'), period == 'all', type =='raw') %>%
#   select(-period, -type, -varExp) %>% 
#   left_join(select(sample_info, -age,-log2age), by='sample_id') %>% 
#   mutate(ind_id = as.character(ind_id)) %>%
#   spread(key=PC, value=value) %>%
#   mutate(id = paste(ind_id, tissue, sep='-')) %>%
#   select(-sample_id,-tissue)
## calculate pairwise distance among tissues per individual:
# mdist = unlist(sapply(unique(dist_dat$ind_id), function(x){
#   pwise_distm_f(dist_dat[dist_dat$ind_id == x,-1])
# }))
# mdist = data.frame(ind_id= names(mdist), mdist = mdist, row.names = NULL) %>%
#   left_join(unique(select(sample_info, ind_id, age)), by = 'ind_id')
saveRDS(mdist, './data/htseq/no_blind/mean_euclidean_dist.rds')

pwise_distMean = function(mat, id_col=1){
  # mat: matrix or data frame with columns as coordinates to calculate distance, and one id column
  # distance is calculated pairwise among rows of same ids
  # id_col: index of the id column
  
  sapply(unique(mat[,id_col] ), function(x){
    ind = mat[mat[,id_col]==x,]
    mean(as.vector(dist(ind[, -id_col])))
  })
}

dist_data = pca_data %>%
  filter(PC %in% c('PC1', 'PC2','PC3', 'PC4'), period == 'all', type =='raw') %>%
  select(-period, -type, -varExp) %>% 
  left_join(select(sample_info, -age,-log2age), by='sample_id') %>% 
  mutate(ind_id = as.character(ind_id)) %>%
  spread(key=PC, value=value) %>%
  select(-sample_id,-tissue)

## calculate pairwise distance among tissues per individual:
mdist = pwise_distMean(dist_data)

mdist = data.frame(ind_id= names(mdist), mdist = mdist, row.names = NULL) %>%
  left_join(unique(select(sample_info, ind_id, age)), by = 'ind_id')

saveRDS(mdist, './data/htseq/no_blind/mean_euclidean_dist.rds')

# spearman test between age and mean pairwise distance in dev. and ageing:
mdist %>%
  mutate(period = ifelse(age <90, 'dev','ageing')) %>%
  group_by(period) %>%
  summarise(rho = cor(mdist,age, m='s'),
            p = cor.test(mdist,age, m='s')$p.val)
# period    rho         p
# <chr>   <dbl>     <dbl>
# 1 ageing -0.597 0.0899   
# 2 dev     0.991 0.0000146

#### dev only pairwise dist:
dist_data_dev = pca_data %>%
  filter(PC %in% c('PC1','PC2','PC3', 'PC4'), period == 'development', type =='raw') %>%
  select(-period, -type, -varExp) %>% 
  left_join(select(sample_info, -age,-log2age), by='sample_id') %>% 
  mutate(ind_id = as.character(ind_id)) %>%
  spread(key=PC, value=value) %>%
  select(-sample_id,-tissue)

## calculate pairwise distance among tissues per individual:
mdist_dev = pwise_distMean(dist_data_dev)

mdist_dev = data.frame(ind_id= names(mdist_dev), mdist = mdist_dev, row.names = NULL) %>%
  left_join(unique(select(sample_info, ind_id, age)), by = 'ind_id')

mdist_dev %>%
  summarise(rho = cor(mdist,age, m='s'),
            p = cor.test(mdist,age, m='s')$p.val)
# rho            p
# 1 0.9910312 1.456125e-05

#### ageing only pairwise dist:
dist_dat_aging = pca_data %>%
  filter(PC %in% c('PC1','PC2','PC3', 'PC4'), period == 'aging', type =='raw') %>%
  select(-period, -type, -varExp) %>% 
  left_join(select(sample_info, -age,-log2age), by='sample_id') %>% 
  mutate(ind_id = as.character(ind_id)) %>%
  spread(key=PC, value=value) %>%
  select(-sample_id,-tissue)

mdist_aging = pwise_distMean(dist_dat_aging)

mdist_aging = data.frame(ind_id= names(mdist_aging), mdist = mdist_aging, row.names = NULL) %>%
  left_join(unique(select(sample_info, ind_id, age)), by = 'ind_id')

mdist_aging %>%
  summarise(rho = cor(mdist,age, m='s'),
            p = cor.test(mdist,age, m='s')$p.val)
# rho         p
# 1 -0.5042195 0.1663127

