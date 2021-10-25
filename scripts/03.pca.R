library(tidyverse)
library(openxlsx)

exp = readRDS('./data/processed/raw/expression.rds')
age = readRDS('./data/processed/raw/ages.rds')
ts.ord = readRDS('./data/processed/raw/tissue_ids.rds')
sample_info = readRDS('./data/processed/tidy/sample_info.rds')
samp_id = names(ts.ord)
## PCA all regions with expression data:
pc = prcomp(t(exp),scale = T)
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

pca_data_all = list(pc.sc.ag, #  aging scaled
                    pc, # all raw
                    pc.sc, # all scaled
                    pcag, # aging raw
                    pcdev, # dev raw
                    pc.sc.dev) # dev scaled
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

pc_age_cors = pca_data %>%
  rename(`PCA period` = period) %>%
  left_join(select(sample_info,-ind_id,-log2age)) %>%
  mutate(`age period` = ifelse(age< 90,'Dev','Aging') ) %>%
  group_by(`PCA period`, type, PC, tissue, `age period`) %>%
  summarise(rho = cor.test(age,value,m='s')$est,
            p = cor.test(age,value,m='s')$p.val)

saveRDS(pc_age_cors, file='./results/source_data/f1/pc_age_cors.rds')
write.xlsx(pc_age_cors, './results/SI_tables/TableS1.xlsx')

## PC4 - age cors (alldata - raw - dev): (fig 1b)
pc_age_cors %>%
  filter(type=='raw' & `PCA period` == 'all' & `age period` == 'Dev') %>%
  filter(PC=='PC4')
# `PCA period` type  PC    tissue `age period`    rho         p
# <chr>        <chr> <fct> <fct>  <chr>         <dbl>     <dbl>
# 1 all          raw   PC4   Cortex Dev          -0.955 0.000806 
# 2 all          raw   PC4   Liver  Dev          -0.991 0.0000146
# 3 all          raw   PC4   Lung   Dev          -0.883 0.00845  
# 4 all          raw   PC4   Muscle Dev          -0.955 0.000806 

## PC1 - age cors (alldata - scaled - dev): (fig 1 - fs2)
pc_age_cors %>%
  filter(type=='notissue' & `PCA period` == 'all' & `age period` == 'Dev') %>%
  filter(PC=='PC1')
# `PCA period` type     PC    tissue `age period`   rho         p
# <chr>        <chr>    <fct> <fct>  <chr>        <dbl>     <dbl>
# 1 all          notissue PC1   Cortex Dev          0.955 0.000806 
# 2 all          notissue PC1   Liver  Dev          0.991 0.0000146
# 3 all          notissue PC1   Lung   Dev          0.883 0.00845  
# 4 all          notissue PC1   Muscle Dev          0.955 0.000806 

## PC2 - age cors (alldata - scaled - dev): (fig 1 - fs2)
pc_age_cors %>%
  filter(type=='notissue' & `PCA period` == 'all' & `age period` == 'Dev') %>%
  filter(PC=='PC2')
# `PCA period` type     PC    tissue `age period`    rho         p
# <chr>        <chr>    <fct> <fct>  <chr>         <dbl>     <dbl>
# 1 all          notissue PC2   Cortex Dev          -0.991 0.0000146
# 2 all          notissue PC2   Liver  Dev           0.883 0.00845  
# 3 all          notissue PC2   Lung   Dev          -0.901 0.00562  
# 4 all          notissue PC2   Muscle Dev           0.306 0.504   

## PC2 - age cors (devdata - raw): (fig 1 - fs3ab)
pc_age_cors %>%
  filter(type=='raw' & `PCA period` == 'development') %>%
  filter(PC=='PC2')
# `PCA period` type  PC    tissue `age period`   rho       p
# <chr>        <chr> <fct> <fct>  <chr>        <dbl>   <dbl>
# 1 development  raw   PC2   Cortex Dev          0.793 0.0334 
# 2 development  raw   PC2   Liver  Dev          0.937 0.00185
# 3 development  raw   PC2   Lung   Dev          0.721 0.0676 
# 4 development  raw   PC2   Muscle Dev          0.775 0.0408 

## PC4 - age cors (devdata - raw): (fig 1 - fs3ab)
pc_age_cors %>%
  filter(type=='raw' & `PCA period` == 'development') %>%
  filter(PC=='PC4')
# `PCA period` type  PC    tissue `age period`    rho         p
# <chr>        <chr> <fct> <fct>  <chr>         <dbl>     <dbl>
# 1 development  raw   PC4   Cortex Dev          -0.955 0.000806 
# 2 development  raw   PC4   Liver  Dev          -0.991 0.0000146
# 3 development  raw   PC4   Lung   Dev          -0.883 0.00845  
# 4 development  raw   PC4   Muscle Dev          -0.955 0.000806 

## PC4 - age cors (agdata - raw): (fig 1 - fs3e)
pc_age_cors %>%
  filter(type=='raw' & `PCA period` == 'aging') %>%
  filter(PC=='PC4')
# `PCA period` type  PC    tissue `age period`    rho      p
# <chr>        <chr> <fct> <fct>  <chr>         <dbl>  <dbl>
# 1 aging        raw   PC4   Cortex Aging        -0.299 0.471 
# 2 aging        raw   PC4   Liver  Aging        -0.723 0.0278
# 3 aging        raw   PC4   Lung   Aging        -0.773 0.0145
# 4 aging        raw   PC4   Muscle Aging         0.109 0.780 

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
saveRDS(mdist, './results/source_data/f1/mdist.rds')

# spearman test between age and mean pairwise distance in dev. and ageing:
mdist %>%
  mutate(period = ifelse(age <90, 'dev','ageing')) %>%
  group_by(period) %>%
  summarise(rho = cor(mdist,age, m='s'),
            p = cor.test(mdist,age, m='s')$p.val)
# period    rho         p
# <chr>   <dbl>     <dbl>
# 1 ageing -0.866 0.00256  
# 2 dev     0.991 0.0000146

# mdistplot = mdist %>%
#   mutate(period=ifelse(age<90,'development', 'aging')) %>%
#   mutate(period = factor(period, levels=c('development', 'aging'))) %>%
#   ggplot(aes(x=age, y=mdist)) +
#   geom_point() +
#   scale_x_continuous(trans='log2') +
#   geom_smooth(method = 'loess', se=F, color='midnightblue') +
#   geom_vline(xintercept = 90, linetype='dashed', color='gray20') +
#   ggpubr::stat_cor(aes(group=period),method='spearman', cor.coef.name = 'rho',
#                    label.x.npc = c(0.1, 0.65), label.y.npc = c(0.9, 0.4), size=3) +
#   ylab('Mean Euclidean Distance') +
#   xlab('Age in days (in log2 scale)') + theme_bw()
# 
# ggsave('./results/figure1/mdis.pdf', mdistplot, units = 'cm', width = 12, height = 8,
#        useDingbats = F)
# ggsave('./results/figure1/mdist.png', mdistplot, units = 'cm', width = 12, height = 8)

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
saveRDS(mdist_dev, './results/source_data/f1/mdist_dev.rds')

mdist_dev %>%
  summarise(rho = cor(mdist,age, m='s'),
            p = cor.test(mdist,age, m='s')$p.val)
# rho            p
# 1 0.9549937 0.0008055353

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

saveRDS(mdist_aging, './data/processed/tidy/mean_euclidean_dist_aging.rds')
saveRDS(mdist_aging, './results/source_data/f1/mdist_aging.rds')

mdist_aging %>%
  summarise(rho = cor(mdist,age, m='s'),
            p = cor.test(mdist,age, m='s')$p.val)
# rho= -0.64, p=0.05959

pc_age_cors %>%
  filter(type=='raw' & `PCA period` == 'aging') %>%
  filter(PC=='PC1')
# `Scale period` type  PC    tissue `age period`     rho        p
# <chr>          <chr> <fct> <fct>  <chr>          <dbl>    <dbl>
# 1 aging          raw   PC1   Cortex Aging         0.683  0.0621  
# 2 aging          raw   PC1   Liver  Aging         0.924  0.000363
# 3 aging          raw   PC1   Lung   Aging         0.832  0.00541 
# 4 aging          raw   PC1   Muscle Aging        -0.0420 0.915  

pc_age_cors %>%
  filter(type=='raw' & `PCA period` == 'aging') %>%
  filter(PC=='PC4')
# `Scale period` type  PC    tissue `age period`    rho      p
# <chr>          <chr> <fct> <fct>  <chr>         <dbl>  <dbl>
# 1 aging          raw   PC4   Cortex Aging        -0.299 0.471 
# 2 aging          raw   PC4   Liver  Aging        -0.723 0.0278
# 3 aging          raw   PC4   Lung   Aging        -0.773 0.0145
# 4 aging          raw   PC4   Muscle Aging         0.109 0.780 

