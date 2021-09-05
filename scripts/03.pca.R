library(tidyverse)
library(openxlsx)

exp = readRDS('./data/processed/raw/expression.rds')
age = readRDS('./data/processed/raw/ages.rds')
ts.ord = readRDS('./data/processed/raw/tissue_ids.rds')
sample_info = readRDS('./data/processed/tidy/sample_info.rds')
samp_id = names(ts.ord)
## PCA all regions with expression data:
pc = prcomp(t(exp),scale = T)
# cor(pc$x[,1], age, m="s") # 0.16
# cor(pc$x[,2], age, m="s") # 0.26
# cor(pc$x[,3], age, m="s") # 0.06
# cor.test(pc$x[,4], age, m="s") # 0.62

agedev = age[age<93]
expdev = exp[, age<93]
pcdev = prcomp(t(expdev), scale = T)

ageag = age[age>=93]
expag = exp[, age>=93]
pcag = prcomp(t(expag), scale = T)

exp.sc = cbind(t(scale(t(exp[,ts.ord=='Cortex']))), t(scale(t(exp[,ts.ord=="Lung"]))),
               t(scale(t(exp[,ts.ord=="Liver"]))), t(scale(t(exp[,ts.ord=="Muscle"]))))
pc.sc = prcomp(t(exp.sc), scale = T)

exp.sc.dev = cbind(t(scale(t(exp[,ts.ord=="Cortex" & age<93]))),
                   t(scale(t(exp[,ts.ord=="Lung" & age<93]))),
                   t(scale(t(exp[,ts.ord=="Liver" & age<93]))),
                   t(scale(t(exp[,ts.ord=="Muscle" & age<93]))))
pc.sc.dev = prcomp(t(exp.sc.dev),scale = T)

exp.sc.ag = cbind(t(scale(t(exp[,ts.ord=="Cortex" & age>=93]))),
                  t(scale(t(exp[,ts.ord=="Lung" & age>=93]))),
                  t(scale(t(exp[,ts.ord=="Liver" & age>=93]))),
                  t(scale(t(exp[,ts.ord=="Muscle" & age>=93]))))
which(apply(is.na(exp.sc.ag),1,function(x)sum(x)>0) )
exp.sc.ag = exp.sc.ag[!apply(is.na(exp.sc.ag),1,function(x)sum(x)>0),]
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

saveRDS(pca_data,'./data/processed/tidy/pca_data.rds')

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

write.xlsx(pc_age_cors, './results/SI_tables/TableS1.xlsx')

pc_age_cors %>%
  filter(type=='notissue' & `Scale period` == 'all' & `age period` == 'Dev') %>%
  filter(PC=='PC1')
# `Scale period` type     PC    tissue `age period`   rho         p
# <chr>          <chr>    <fct> <fct>  <chr>        <dbl>     <dbl>
#   1 all            notissue PC1   Cortex Dev          0.955 0.000806 
# 2 all            notissue PC1   Liver  Dev          0.991 0.0000146
# 3 all            notissue PC1   Lung   Dev          0.883 0.00845  
# 4 all            notissue PC1   Muscle Dev          0.955 0.000806 
pc_age_cors %>%
  filter(type=='notissue' & `Scale period` == 'all' & `age period` == 'Dev') %>%
  filter(PC=='PC2')
# `Scale period` type     PC    tissue `age period`    rho         p
# <chr>          <chr>    <fct> <fct>  <chr>         <dbl>     <dbl>
#   1 all            notissue PC2   Cortex Dev          -0.991 0.0000146
# 2 all            notissue PC2   Liver  Dev           0.883 0.00845  
# 3 all            notissue PC2   Lung   Dev          -0.901 0.00562  
# 4 all            notissue PC2   Muscle Dev           0.306 0.504  
pc_age_cors %>%
  filter(type=='raw' & `Scale period` == 'development') %>%
  filter(PC=='PC2')
# `Scale period` type  PC    tissue `age period`   rho       p
# <chr>          <chr> <fct> <fct>  <chr>        <dbl>   <dbl>
#   1 development    raw   PC2   Cortex Dev          0.793 0.0334 
# 2 development    raw   PC2   Liver  Dev          0.937 0.00185
# 3 development    raw   PC2   Lung   Dev          0.721 0.0676 
# 4 development    raw   PC2   Muscle Dev          0.775 0.0408 
pc_age_cors %>%
  filter(type=='raw' & `Scale period` == 'development') %>%
  filter(PC=='PC4')
# `Scale period` type  PC    tissue `age period`    rho         p
# <chr>          <chr> <fct> <fct>  <chr>         <dbl>     <dbl>
#   1 development    raw   PC4   Cortex Dev          -0.955 0.000806 
# 2 development    raw   PC4   Liver  Dev          -0.991 0.0000146
# 3 development    raw   PC4   Lung   Dev          -0.883 0.00845  
# 4 development    raw   PC4   Muscle Dev          -0.955 0.000806 

#################################
################################# does PC1,2,3 separate samples according to tissues? check with Anova
################################# and eclidean distance
#################################

aov_dat = pca_data %>%
  filter(type == 'raw', period == 'all', PC%in%c('PC1','PC2','PC3','PC4')) %>% 
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
pwise_distMean = function(mat, id_col=1){
  # mat: matrix or data frame with columns as coordinates to calculate distance, and one id column
  # distance is calculated pairwise among rows of same ids
  # id_col: index of the id column
  
  sapply(unique(mat[,id_col] ), function(x){
    ind = mat[mat[,id_col]==x,]
    mean(as.vector(dist(ind[, -id_col])))
  })
}

## mean pairwise distance of tissues in PC1,2,3,4 space:
dist_dat = pca_data %>%
  filter(PC %in% c('PC1','PC2','PC3', 'PC4'), period == 'all', type =='raw') %>%
  select(-period, -type, -varExp) %>% 
  left_join(select(sample_info, -age,-log2age), by='sample_id') %>% 
  mutate(ind_id = as.character(ind_id)) %>%
  spread(key=PC, value=value) %>%
  select(-sample_id,-tissue)

mdist = pwise_distMean(dist_dat)

mdist = data.frame(ind_id= names(mdist), mdist = mdist, row.names = NULL) %>%
  left_join(unique(select(sample_info, ind_id, age)), by = 'ind_id')

saveRDS(mdist, './data/processed/tidy/mean_euclidean_dist.rds')

# spearman test between age and mean pairwise distance in dev. and ageing:
mdist %>%
  mutate(period = ifelse(age <90, 'dev','ageing')) %>%
  group_by(period) %>%
  summarise(rho = cor(mdist,age, m='s'),
            p = cor.test(mdist,age, m='s')$p.val)

# dev: (PC3,4) 0.90, 0.00562 (old) -> (PC1,2,3,4) 0.99, 0.0000146
# ageing: (PC3,4) -0.303, 0.429 (old) -> (PC1,2,3,4) -0.87, 0.00256

######
#### dev only pairwise dist:
dist_dat_dev = pca_data %>%
  filter(PC %in% c('PC1','PC2','PC3', 'PC4'), period == 'development', type =='raw') %>%
  select(-period, -type, -varExp) %>% 
  left_join(select(sample_info, -age,-log2age), by='sample_id') %>% 
  mutate(ind_id = as.character(ind_id)) %>%
  spread(key=PC, value=value) %>%
  select(-sample_id,-tissue)

mdist_dev = pwise_distMean(dist_dat_dev)

mdist_dev = data.frame(ind_id= names(mdist_dev), mdist = mdist_dev, row.names = NULL) %>%
  left_join(unique(select(sample_info, ind_id, age)), by = 'ind_id')

saveRDS(mdist_dev, './data/processed/tidy/mean_euclidean_dist_dev.rds')

mdist_dev %>%
  summarise(rho = cor(mdist,age, m='s'),
            p = cor.test(mdist,age, m='s')$p.val)
#(old (PC3-4 only: rho=0.63, p=0.13)) --> rho = 0.95, p = 0.0008

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

#saveRDS(mdist_aging, './data/processed/tidy/mean_euclidean_dist_aging.rds')

mdist_aging %>%
  summarise(rho = cor(mdist,age, m='s'),
            p = cor.test(mdist,age, m='s')$p.val)
# rho= -0.64, p=0.05959

pc_age_cors %>%
  filter(type=='raw' & `Scale period` == 'aging' & `age period` == 'Aging') %>%
  filter(PC=='PC1')
# `Scale period` type  PC    tissue `age period`     rho        p
# <chr>          <chr> <fct> <fct>  <chr>          <dbl>    <dbl>
#   1 aging          raw   PC1   Cortex Aging         0.683  0.0621  
# 2 aging          raw   PC1   Liver  Aging         0.924  0.000363
# 3 aging          raw   PC1   Lung   Aging         0.832  0.00541 
# 4 aging          raw   PC1   Muscle Aging        -0.0420 0.915  

pc_age_cors %>%
  filter(type=='raw' & `Scale period` == 'aging' & `age period` == 'Aging') %>%
  filter(PC=='PC4')
# `Scale period` type  PC    tissue `age period`    rho      p
# <chr>          <chr> <fct> <fct>  <chr>         <dbl>  <dbl>
#   1 aging          raw   PC4   Cortex Aging        -0.299 0.471 
# 2 aging          raw   PC4   Liver  Aging        -0.723 0.0278
# 3 aging          raw   PC4   Lung   Aging        -0.773 0.0145
# 4 aging          raw   PC4   Muscle Aging         0.109 0.780 





