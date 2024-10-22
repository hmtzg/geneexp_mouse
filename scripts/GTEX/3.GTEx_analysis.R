#MD
library(tidyverse)
library(ggpubr)
library(ggforce)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)
sexcolors <- setNames(c('#ffc9b5', '#8fb8de'), c('Female', 'Male'))
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('Up','Down'))
# revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))
revcol = setNames(c('gray25', 'gray25', 'gray25', 'gray25'), c('UpDown','DownUp','UpUp','DownDown'))
sexshape =  c('female' = 17, 'male' = 19)

#### Functions ####

plotsave <- function(ggobj, prefix, width, height, ...){
  path = strsplit(prefix,'/')[[1]]
  path = paste(path[-length(path)], collapse = '/')
  if(!file.exists(path)){
    system(paste('mkdir -p',path))
  }
  saveRDS(object = ggobj, file = paste(prefix,'.rds',sep=''))
  ggsave(file = paste(prefix, '.pdf', sep = ''), plot = ggobj, units = 'cm', 
         width = width, height = height, useDingbats = F, limitsize = F)
  ggsave(file = paste(prefix, '.png', sep = ''), plot = ggobj, units = 'cm', bg='white',
         width = width, height = height, limitsize = F)
}

tablesave <- function(tib, prefix, ...){
  path = strsplit(prefix,'/')[[1]]
  path = paste(path[-length(path)], collapse = '/')
  if(!file.exists(path)){
    system(paste('mkdir -p',path))
  }
  readr::write_csv(tib, path = paste(prefix,'.csv', sep = ''), append = F)
  readr::write_tsv(tib, path = paste(prefix,'.tsv', sep = ''), append = F)
  saveRDS(object = tib, file = paste(prefix,'.rds',sep=''))
}

#### PCA ####
attr = readRDS('./data/processed/GTEx/attr.rds') 

  # filter(!exc) #exclude individuals with more than one sample in a minor tissue

expvals = list(Cortex = readRDS('./data/processed/GTEx/expression/tpm/Brain-Cortex.rds'),
               Lung = readRDS('./data/processed/GTEx/expression/tpm/Lung.rds'),
               Liver = readRDS('./data/processed/GTEx/expression/tpm/Liver.rds'),
               Muscle = readRDS('./data/processed/GTEx/expression/tpm/Muscle-Skeletal.rds'))


samplelist = unique(unlist(lapply(expvals, function(x)colnames(x))))

# there are samples from same individual but not having expression:
attr = attr %>%
  filter(minor_tissue %in% c('Lung','Liver','Brain - Cortex','Muscle - Skeletal')) %>%
  filter(sample_id %in% samplelist) %>%
  group_by(id,minor_tissue) %>%
  summarise(n = length(unique(sample_id))) %>%
  ungroup() %>%
  mutate(exc = n>1) %>%
  right_join(attr) 
table(attr$exc) # no duplication
summary(attr$exc)
genelist = unique(unlist(lapply(expvals,function(expx){
  sapply(strsplit(as.character(rownames(expx)),'[.]'),function(x)x[1])
})))
martx = biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
#martx = biomaRt::useEnsembl('ensembl', 'hsapiens_gene_ensembl', mirror='asia')
biotype = biomaRt::getBM(attributes = c('ensembl_gene_id', 'gene_biotype'), filters = 'ensembl_gene_id', 
                         values = genelist, mart = martx)
codinggenes = unique(filter(biotype,gene_biotype == 'protein_coding')$ensembl_gene_id)

incids = setdiff(unique(( sapply(expvals, colnames) %>% 
  reshape2::melt() %>%
  set_names(c('sample_id','Tissue')) %>%
  left_join(attr) %>%
  group_by(id) %>%
  summarise( n = length(unique(Tissue))) %>%
  mutate(inc = n==4) %>%
  filter(inc))$id),NA) # compile the list of individuals with samples in all 4 tissues

idmap = setNames(attr$id, attr$sample_id)[samplelist]
inc_sampleIDs = names(idmap[idmap%in%incids]) # get the sample ids for the individuals to be included

# subset attr data for only samples common in four tissues:
attr = attr %>% 
  filter(sample_id %in% inc_sampleIDs ) 

expvals = lapply(expvals, function(expx){
  rownames(expx) = sapply(strsplit(as.character(rownames(expx)),'[.]'),function(x)x[1]) # remove gene version
  expx = expx[rownames(expx) %in% codinggenes,] #get only the coding genes
  dupgenes = unique(rownames(expx)[duplicated(rownames(expx))]) # to remove duplicated genes
  expx = expx[!(rownames(expx)%in%dupgenes), colnames(expx) %in% inc_sampleIDs] 
  # exc = apply(expx,1,function(x) (max(x)<1)|(mean(x==0)==1))
  # expx = expx[!exc,] 
  expx
})

sample2id = setNames(attr$id, attr$sample_id)
expx = expvals[[4]]
sum(duplicated(unname(sample2id[colnames(expx)])))

genelist = Reduce('intersect',lapply(expvals,rownames))
expvals = cbind(expvals[[1]][genelist,],expvals[[2]][genelist,],
                expvals[[3]][genelist,],expvals[[4]][genelist,]) #to make sure the gene order is the same
expvals = expvals [ ! (rowMeans(expvals==0)>0.25), ] # exclude the genes if expression value is 0 in more than
# 25% of the samples
dim(expvals)
#16211 - 188 (47 samples)  -> 16197

#normalization
exp_l2 = log2(expvals+1)
exp_l2_qn = preprocessCore::normalize.quantiles(exp_l2)
dimnames(exp_l2_qn) = dimnames(exp_l2)
# 16197 188

pcx = prcomp(t(exp_l2_qn),scale=T)
pca_data = data.frame(pcx$x[,1:4], sample_id = rownames(pcx$x)) %>%
  left_join(attr) %>%
  mutate(Tissue = ifelse(major_tissue=='Brain','Cortex',major_tissue)) %>%
  mutate(Age = c(25,35,45,55,65,75)[age])
pcimp = round(summary(pcx)$imp[2,1:4]*100,2)
pc1_2 = pca_data %>%
  #ggplot(aes( x = PC1, y= PC2 , color = Tissue, size = age, shape=sex)) +
  ggplot(aes( x = PC1, y= PC2 , color = Tissue, size = age)) +
  geom_point() +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range = c(0.2,1.5)) +
  scale_shape_manual(values = sexshape) +
  xlab(paste('PC 1 (',pcimp[1],'%)',sep=''))+
  ylab(paste('PC 2 (',pcimp[2],'%)',sep=''))

pc3_4 = pca_data %>%
  #ggplot(aes( x = PC3, y= PC4 , color = Tissue, size = age, shape=sex)) +
  ggplot(aes( x = PC3, y= PC4 , color = Tissue, size = age)) +
  geom_point() +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range = c(0.2,1.5)) +
  scale_shape_manual(values = sexshape) +
  xlab(paste('PC 3 (',pcimp[3],'%)',sep=''))+
  ylab(paste('PC 4 (',pcimp[4],'%)',sep=''))

pc1_age = pca_data %>%
  #ggplot(aes(x = Age, y = PC1, color = Tissue, linetype=sex)) +
  ggplot(aes(x = Age, y = PC1, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = T, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y') +
  ggtitle('PC1') +ylab(NULL)

pc2_age = pca_data %>%
  ggplot(aes(x = Age, y = PC2, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = F, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y')+
  ggtitle('PC2') +ylab(NULL)

pc3_age = pca_data %>%
  ggplot(aes(x = Age, y = PC3, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = F, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y')+
  ggtitle('PC3') +ylab(NULL)

pc4_age = pca_data %>%
  ggplot(aes(x = Age, y = PC4, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = F, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y')+
  ggtitle('PC4') +ylab(NULL)

##########
##########
########## all PC-age cors plot:
##########
##########
# pcdatall = data.frame(pcx$x, sample_id = rownames(pcx$x)) %>%
#   left_join(select(attr,sample_id, age, major_tissue)) %>%
#   mutate(Tissue = ifelse(major_tissue=='Brain','Cortex', major_tissue)) %>%
#   mutate(Age = c(25,35,45,55,65,75)[age]) %>%
#   select(-age, -major_tissue)
# 
# pccors = pcdatall %>% 
#   gather(key = 'PC', value = 'value', -Age, -sample_id, -Tissue) %>%
#   group_by(PC, Tissue) %>%
#   summarise(rho = abs(cor.test(Age, value, method='spearman')$est),
#             pval = cor.test(Age, value, method='spearman')$p.val ) %>% ungroup() %>%
#   mutate(PC = factor(PC,  levels = paste0('PC', 1:188) ) ) %>%
#   ggplot(aes(x = PC, y = rho, color=Tissue)) +
#   geom_bar(stat='identity') +
#   facet_wrap(~Tissue) +
#   ylab( bquote('Abs('*rho[age-PC]~')') ) +
#   xlab('PC1-PC188') +
#   scale_color_manual(values = tissuecol) +
#   theme(axis.text.x  = element_blank(),
#         axis.ticks.x = element_blank())
# pccors
# 
# plotsave(ggobj = pccors, prefix = './results/GTEx/pc_cors', width = 16, height = 12)

#calculate the Eucledian distance using all PCs because there is no clear relationship between PC and age

# has bug, repeating values:
# meanEuc = reshape2::melt(as.matrix(dist(pcx$x))) %>%
#   mutate(id1 = idmap[as.character(Var1)],
#          id2 = idmap[as.character(Var2)]) %>%
#   filter(id1 == id2) %>%
#   filter(Var1!=Var2) %>%
#   group_by(id1) %>%
#   summarise(dist = mean(unique(value))) %>%
#   ungroup() %>% rename(id = id1) %>%
#   left_join(attr)

####
pwise_distMean = function(mat, id_col=1){
  # mat: matrix or data frame with columns as coordinates to calculate distance, and one id column
  # distance is calculated pairwise among rows of same ids
  # id_col: index of the id column
  
  sapply(unique(mat[,id_col] ), function(x){
    ind = mat[mat[,id_col]==x,]
    mean(as.vector(dist(ind[, -id_col])))
  })
}

meanEuc = reshape2::melt( as.data.frame(pcx$x) %>%
  mutate(id = idmap[rownames(pcx$x)]) %>%
  relocate(id) %>% 
  pwise_distMean(., id_col = 1) ) %>%
  mutate(id = rownames(.)) %>% 
  left_join(attr) %>%
  distinct(id, value, age, sex) %>%
  mutate(Age = c(25,35,45,55,65,75)[age]) %>%
  mutate(sex = str_to_title(sex))

eucdist = meanEuc %>%
  #ggplot(aes(x = Age, y= value,  shape=sex)) +
  ggplot(aes(x = Age, y= value)) +
  geom_point(aes(size = age)) +
  #geom_point() +
  scale_size_discrete(range = c(0.2,1.5)) +
  geom_smooth(method = 'lm', se = F, color='darkred') +
  #scale_color_manual(values=c('Female'= 'darkslategrey' ,'Male'= 'olivedrab3')) +
  stat_cor(method = 'spearman', size = 6/pntnorm, cor.coef.name = 'rho') +
  ylab('Mean Pairwise Euclidean Distance')

meanEuc %>% count(sex)
## 11 female
## 36 male

####
# eucdist = meanEuc %>%
#   mutate(Age = c(25,35,45,55,65,75)[age]) %>%
#   ggplot(aes(x = Age, y= dist)) +
#   geom_point(aes(size = age)) +
#   scale_size_discrete(range = c(0.2,1.5)) +
#   geom_smooth(method = 'lm', se = F, color = 'darkred') +
#   stat_cor(method = 'spearman', size = 6/pntnorm, cor.coef.name = 'rho') +
#   ylab('Mean Pairwise Euclidean Distance')

pcaplots1 = ggarrange(pc1_2, pc3_4, eucdist, nrow = 1, ncol =3 , labels = c('a.','b.','c.'), 
                      font.label = list(size = 8), widths = c(2,2,1.5), common.legend = T, legend = 'right')
pcaplots2 = ggarrange(pc1_age, pc2_age, pc3_age,pc4_age, nrow = 1, ncol =4 , 
                      labels = c('d.','e.','f.','g.'), font.label = list(size = 8), legend = 'none')

pca_plots = ggarrange(pcaplots1, pcaplots2, ncol=1,nrow=2, heights = c(2,2.5))

plotsave(ggobj = pca_plots, prefix = './results/GTEx/pca', width = 16, height = 12)
plotsave(ggobj = pca_plots, prefix = './results/figure_supplements/fs2/FS8', width = 16, height = 12)

saveRDS(pca_data,'results/source_data/f2/fs8_pca.rds')
saveRDS(meanEuc, 'results/source_data/f2/fs8_euc.rds')

#### Gene Expression Analysis ####

expvals = reshape2::melt(exp_l2_qn) %>%
  set_names(c('GeneID','sample_id','Expression')) %>%
  left_join(attr) %>%
  mutate(Tissue = ifelse(major_tissue=='Brain','Cortex',major_tissue)) %>%
  mutate(Age = c(25,35,45,55,65,75)[age]) 

#Expression - Age correlation
agecor = expvals %>%
  group_by(Tissue, GeneID) %>%
  summarise(rho = cor.test(Expression, Age, method = 's')$est,
            p = cor.test(Expression, Age, method = 's')$p.val)

agecor = ungroup(agecor) %>%
  mutate(adjusted_p = p.adjust(p, method = 'fdr'))
saveRDS(agecor, 'results/GTEx/agecor.rds')

#### CoV analysis #### 

#CoV calculation for each gene of each individual (taking sd and mean of 4 tissues)
expvals = expvals %>%
  group_by(GeneID, id, sex) %>%
  summarise(CoV = sd(Expression)/mean(Expression)) %>%
  ungroup() %>%
  right_join(expvals)

# summarise CoV across genes by taking the mean or median of CoVs per gene
sumCov = expvals %>%
  select(id,CoV,age,Age,GeneID) %>%
  unique() %>%
  group_by(id, age, Age) %>%
  summarise(meancov = mean(CoV),
            medcov = median(CoV)) %>%
  ungroup()

meancovplot = ggplot(sumCov,aes(x = age, y = meancov)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
  geom_jitter(width = 0.1, size = 0.3) +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) +
  xlab(NULL) + ylab('Mean CoV') +
  theme(axis.text.x = element_text(angle=90))

medcovplot = ggplot(sumCov,aes(x = age, y = medcov)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
  geom_jitter(width = 0.1, size = 0.3) +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', label.y = 0.31, size = 6/pntnorm) +
  xlab(NULL) + ylab('Median CoV') +
  theme(axis.text.x = element_text(angle=90))

# Calculate pairwise correlation coefficients between tissues
pairwisedat = expvals %>%
  select(Tissue, Expression, GeneID, id, age, Age) %>%
  unique()

pairwisedat = pairwisedat %>%
  spread(key = Tissue, value = Expression) %>%
  group_by(id, age, Age) %>%
  summarise(`Cortex-Liver` = cor(Cortex, Liver, method = 'spearman'),
            `Cortex-Lung` = cor(Cortex, Lung, method = 'spearman'),
            `Cortex-Muscle` = cor(Cortex, Muscle, method = 'spearman'),
            `Liver-Lung` = cor(Liver, Lung, method = 'spearman'),
            `Liver-Muscle` = cor(Liver, Muscle, method = 'spearman'),
            `Lung-Muscle` = cor(Lung, Muscle, method = 'spearman'))

pairwiseplot = pairwisedat %>%
  gather(key = 'type', value = 'Correlation', -id, -age, -Age) %>%
  ggplot(aes(x = age, y = Correlation)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
  geom_jitter(width = 0.1, size = 0.3) +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) +
  facet_wrap(~type, scales = 'free_y') +
  ylab('Spearman\'s Correlation Coefficient') + xlab(NULL) +
  theme(axis.text.x = element_text(angle=90))
  

p1 = ggarrange(meancovplot, medcovplot, ncol = 1, nrow = 2, labels = c('a.','b.'), font.label = list(size = 8))
covresplot = ggarrange(p1, pairwiseplot, ncol =2, nrow = 1, labels = c(NA,'c.'), font.label = list(size = 8), 
                       widths = c(1,3))

plotsave(ggobj = covresplot, prefix = './results/GTEx/CoV',width = 16, height = 8)
plotsave(ggobj = covresplot, prefix = './results/figure_supplements/fs2/FS9',width = 16, height = 8)

saveRDS(sumCov,'results/source_data/f2/fs9_CoV_mean_med.rds')
saveRDS(pairwisedat,'results/source_data/f2/fs9_pairwisecors.rds')

## mean pairwise exp cors:
mpwise= pairwisedat %>% 
  gather(key = 'pairs', value='rho', -id, -age,-Age) %>%
  group_by(id, Age) %>%
  summarise(mean = mean(rho),
            median = median(rho)) %>%
  gather(key='method', value = 'rho', mean, median) %>%
  ggplot(aes(x=Age,y=rho)) +
  facet_wrap(~method, strip.position = 'left',
             labeller = as_labeller(c(mean = 'Mean P.wise Expr. Corr.',
                                      median = 'Median P.wise Expr. Corr.'))) +
  geom_point(size=1.5, color="steelblue", alpha=0.9) +
  geom_smooth(se=T, method = 'lm', color = 'midnightblue', fill='lightblue') +
  scale_x_continuous(trans = 'log2') +
  xlab('Age in years (in log2 scale)') + ylab(NULL) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size=4) +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 8)) 

ggsave('./results/GTEx/meanpwisecor.pdf', mpwise, units = 'cm', width = 12, height = 8,
       useDingbats =F)
ggsave('./results/GTEx/meanpwisecor.png', mpwise, units = 'cm', width = 12, height = 8)


# Calculate CoV change with age per each gene
covcor = expvals %>%
  select(GeneID, Age, CoV, id) %>%
  unique() %>%
  group_by(GeneID) %>%
  summarise(rho = cor.test(CoV, Age, method = 's')$est,
            p = cor.test(CoV, Age, method = 's')$p.val)

covcor = ungroup(covcor) %>%
  mutate(adjusted_p = p.adjust(p, method = 'fdr'))

covcor %>%
  filter(adjusted_p<0.1)
#nothing significant
table(covcor$rho<0)
# FALSE  TRUE 
# 6977  9220 
# 9220 convergent
saveRDS(covcor, file='./results/GTEx/covcor.rds')

############
############ repeat with different sexes
############

# per sex (47 individuals):
agecorS = expvals %>%
  group_by(Tissue, GeneID, sex) %>%
  summarise(rho = cor.test(Expression, Age, method = 's')$est,
            p = cor.test(Expression, Age, method = 's')$p.val)

agecorS = ungroup(agecorS) %>%
  mutate(adjusted_p = p.adjust(p, method = 'fdr'))

# Expression Change correlations between female and male (47 samples total):
agecorS %>%
  select(-p, -adjusted_p) %>% 
  spread(key= sex, value = rho) %>%
  group_by(Tissue) %>%
  summarise(rho = cor.test(female, male, method='s')$est,
            p = cor.test(female, male, method='s')$p.val)
# 1 Cortex  0.295  8.  e-323
# 2 Liver   0.0984 3.88e- 36
# 3 Lung   -0.0377 1.61e-  6
# 4 Muscle  0.196  1.34e-139

sexcorplot = agecorS %>%
  select(-p, -adjusted_p) %>% 
  spread(key= sex, value = rho) %>%
  ggplot(aes(x=female, y=male)) +
  facet_grid(~Tissue) +
  geom_point(size=0.4, alpha=0.3, color='gray30') +
  geom_smooth(method = 'lm', se=F, color='darkred') +
  stat_cor(method = 'spearman', size = 6/pntnorm, cor.coef.name = 'rho') 

sexcorplot
plotsave(ggobj = sexcorplot, prefix = './results/GTEx/sex_agecor',width = 16, height = 12)

## using all samples for each tissue:
##....

# pca_data %>%
#   ggplot(aes( x = PC1, y= PC2 , color = Tissue, size = age, shape=sex)) +
#   geom_point() +
#   scale_color_manual(values = tissuecol) +
#   scale_size_discrete(range = c(0.2,1.5)) +
#   scale_shape_manual(values = c('female' = 17, 'male' = 19)) +
#   xlab(paste('PC 1 (',pcimp[1],'%)',sep='')) +
#   ylab(paste('PC 2 (',pcimp[2],'%)',sep=''))


#### CoV analysis #### 

############# CoV by excluding each tissue:
ts = unique(expvals$Tissue)
#ts=ts[1]

### CoV by sex:
# summarise CoV across genes by taking the mean or median of CoVs per gene
sumCovS = expvals %>%
  select(id,CoV,age,Age,GeneID, sex) %>%
  unique() %>%
  group_by(id, age, Age, sex) %>%
  summarise(meancov = mean(CoV),
            medcov = median(CoV)) %>%
  ungroup()

# mean CoV change per sex:
sumCovS %>%
  group_by(sex) %>%
  summarise(rho = cor.test(Age, meancov, m='s', na.rm=T)$est,
            p = cor.test(Age, meancov, m='s', na.rm=T)$p.val)
# sex        rho      p
# <fct>    <dbl>  <dbl>
# 1 female -0.584  0.0593
# 2 male    0.0516 0.765 

# median CoV change per sex:
sumCovS %>%
  group_by(sex) %>%
  summarise(rho = cor.test(Age, medcov, m='s', na.rm=T)$est,
            p = cor.test(Age, medcov, m='s', na.rm=T)$p.val)
# sex        rho      p
# <fct>    <dbl>  <dbl>
# 1 female -0.673  0.0231
# 2 male   -0.0540 0.754

# meancovplotsex = ggplot(sumCovS, aes(x = age, y = meancov)) +
#   geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
#   facet_grid(sex~.) +
#   geom_jitter(width = 0.1, size = 0.3) +
#   #geom_smooth(aes(x=as.numeric(age)), se = F, size = 0.4, show.legend = T, method = 'lm') +
#   stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm,
#            label.x = c(0.5,0.5), label.y = c(0.51,0.51) ) +
#   xlab(NULL) + ylab('Mean CoV') +
#   theme(axis.text.x = element_text(angle=90))

meancovplotsex = ggplot(sumCovS, aes(x = age, y = meancov)) +
  #geom_boxplot(width = 0.4, outlier.shape = NA, fill = adjustcolor('gray70', alpha.f = 0.3)) +
  facet_grid(sex~.) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_smooth(aes(x=as.numeric(age)), se = F, size = 0.8, show.legend = T, method = 'lm') +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm,
           label.x = c(0.5,0.5), label.y = c(0.51,0.51) ) +
  xlab(NULL) + ylab('Mean CoV') +
  theme(axis.text.x = element_text(angle=90))

medcovplotsex = ggplot(sumCovS, aes(x = age, y = medcov)) +
  facet_grid(sex~.) +
  #geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_smooth(aes(x=as.numeric(age)), se = F, size = 0.8, show.legend = T, method = 'lm') +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho',  size = 6/pntnorm, 
           label.x = c(0.5,0.5) ) +
  xlab(NULL) + ylab('Median CoV') +
  theme(axis.text.x = element_text(angle=90))

covsex = ggarrange(meancovplotsex, medcovplotsex, nrow=1, ncol=2, labels=c('a.','b.'), 
                   font.label = list(size=8) )

plotsave(ggobj = covsex, prefix = './results/GTEx/CoVsex',width = 16, height = 14)

# Calculate CoV change with age per each gene
covcorS = expvals %>%
  select(GeneID, Age, CoV, id, sex) %>%
  unique() %>%
  group_by(GeneID, sex) %>%
  summarise(rho = cor.test(CoV, Age, method = 's')$est,
            p = cor.test(CoV, Age, method = 's')$p.val)

covcorS = ungroup(covcorS) %>%
  mutate(adjusted_p = p.adjust(p, method = 'fdr'))

## correlation between female and male CoV change values:
covcorS %>%
  select(-p, -adjusted_p) %>%
  spread(key=sex, value=rho) %>%
  summarise(rho = cor.test(female, male, method='s')$est,
            p = cor.test(female, male, method='s')$p.val)
# rho        p
# <dbl>    <dbl>
#  1 0.135 5.07e-67

# Calculate pairwise correlation coefficients between tissues
pairwisedatsex = expvals %>%
  select(Tissue, Expression, GeneID, id, age, Age, sex) %>%
  unique()

pairwisedatsex = pairwisedatsex %>%
  spread(key = Tissue, value = Expression) %>%
  group_by(id, age, Age, sex) %>%
  summarise(`Cortex-Liver` = cor(Cortex, Liver, method = 'spearman'),
            `Cortex-Lung` = cor(Cortex, Lung, method = 'spearman'),
            `Cortex-Muscle` = cor(Cortex, Muscle, method = 'spearman'),
            `Liver-Lung` = cor(Liver, Lung, method = 'spearman'),
            `Liver-Muscle` = cor(Liver, Muscle, method = 'spearman'),
            `Lung-Muscle` = cor(Lung, Muscle, method = 'spearman'))

# pairwiseplotsex1 = pairwisedatsex %>%
#   gather(key = 'type', value = 'Correlation', -id, -age, -Age, -sex) %>%
#   filter(sex=='female') %>%
#   ggplot(aes(x = age, y = Correlation)) +
#   facet_wrap(~type, scales='free_y', ncol=6)+
#   #geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
#   geom_jitter(width = 0.1, size = 0.3) +
#   #geom_smooth(aes(as.numeric(x=age)), method='lm', se=F ) +
#   stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) +
#   ylab('Spearman\'s Correlation Coefficient') + xlab(NULL) +
#   theme(axis.text.x = element_text(angle=90))

pairwiseplotsex1 = pairwisedatsex %>%
  gather(key = 'type', value = 'Correlation', -id, -age, -Age, -sex) %>%
  filter(sex=='female') %>%
  ggplot(aes(x = age, y = Correlation)) +
  facet_wrap(~type, scales='free_y', ncol=6) +
  #geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_smooth(aes(x=as.numeric(age) ), method='lm', se=F, size=0.8) +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm,
           label.x = 0.8, label.sep = '\n') +
  ylab('Sp. Correlation Coefficient') + xlab(NULL) +
  theme(axis.text.x = element_text(angle=90))
pairwiseplotsex1

xpos = pairwisedatsex %>%
  gather(key = 'type', value = 'Correlation', -id, -age, -Age, -sex) %>%
  filter(sex=='male')
pairwiseplotsex2 = pairwisedatsex %>%
  gather(key = 'type', value = 'Correlation', -id, -age, -Age, -sex) %>%
  filter(sex=='male') %>%
  ggplot(aes(x = age, y = Correlation)) +
  facet_wrap(~type, scales='free_y', ncol=6) +
  #geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_smooth(aes(x = c(1,2,3,4,5)[xpos$age]-1 ), method='lm', se=F, size=0.8 ) +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm, 
           label.x=0.8, label.sep = '\n') +
  ylab('Sp. Correlation Coefficient') + xlab(NULL) +
  theme(axis.text.x = element_text(angle=90))
pairwiseplotsex2

pairwiseplotsex = ggarrange(pairwiseplotsex1, pairwiseplotsex2, nrow=2, labels = c('c.','d.'), 
                            font.label = list(size=8))
pairwiseplotsex 
plotsave(ggobj = pairwiseplotsex, prefix = './results/GTEx/pwisesex', width = 16, height = 12)

sexplots = ggarrange(covsex, pairwiseplotsex, nrow=2, heights = c(0.8,1) )
sexplots

plotsave(ggobj = sexplots, prefix = './results/GTEx/sexplots', width = 16, height = 14)
plotsave(ggobj = sexplots, prefix = './results/figure_supplements/fs2/FS16', width = 16, height = 14)

saveRDS(sumCovS,'results/source_data/f2/fs16_cov_mean_med.rds')
saveRDS(pairwisedatsex,'results/source_data/f2/fs16_pairwisecors.rds')

# Calculate CoV change with age per each gene
covcorsex = expvals %>%
  select(GeneID, Age, CoV, id,sex) %>%
  unique() %>%
  group_by(GeneID,sex) %>%
  summarise(rho = cor.test(CoV, Age, method = 's')$est,
            p = cor.test(CoV, Age, method = 's')$p.val)

covcorsex = ungroup(covcorsex) %>%
  mutate(adjusted_p = p.adjust(p, method = 'fdr'))

covcorsex %>%
  filter(adjusted_p<0.1)
#nothing significant
covcorsex %>%
  group_by(sex) %>% 
  summarise(table(rho<0))
# 1 female 6259            
# 2 female 9938            
# 3 male   7934            
# 4 male   8263 
# 9938 female 8263 male convergent

# meanEuc %>%
#   ggplot(aes(x = Age, y= value)) +
#   facet_grid(sex~.) +
#   geom_point(aes(size = age)) +
#   scale_size_discrete(range = c(0.2,1.5)) +
#   geom_smooth(method = 'lm', se = F, color = 'darkred') +
#   stat_cor(method = 'spearman', size = 6/pntnorm, cor.coef.name = 'rho') +
#   ylab('Mean Pairwise Euclidean Distance')

#####
#####


attr %>% group_by(minor_tissue, sex) %>% summarise(n=length(sample_id))
# 1 Brain - Cortex    female    11
# 2 Brain - Cortex    male      36
# 3 Liver             female    11
# 4 Liver             male      36
# 5 Lung              female    11
# 6 Lung              male      36
# 7 Muscle - Skeletal female    11
# 8 Muscle - Skeletal male      36
attr %>% group_by(sex, age) %>% summarise(n=length(sample_id))

##########
########## dico enrichment:
##########

ddc_genes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
## convert human ENS to mouse:
martx = biomaRt::useMart(biomart = 'ensembl')
martHs_ = biomaRt::useDataset('hsapiens_gene_ensembl', mart=martx)
martMm_ = biomaRt::useDataset('mmusculus_gene_ensembl', mart=martx)
ensmap = biomaRt::getLDS(attributes = c('ensembl_gene_id'), filters = 'ensembl_gene_id', 
                         values = names(ddc_genes), mart = martMm_, 
                         attributesL = c('ensembl_gene_id'), martL = martHs_)
colnames(ensmap) = c('ENS_mm', 'ENS_hs')
dup1 = unique(ensmap$ENS_hs[duplicated(ensmap$ENS_hs)]) #  108 duplicate genes
ensmap = ensmap[!ensmap$ENS_hs%in%dup1,]
dup2 = unique(ensmap$ENS_mm[duplicated(ensmap$ENS_mm)]) # 258 duplicate genes
ensmap = ensmap[!ensmap$ENS_mm%in%dup2,] 
colnames(ensmap)[2]='GeneID'

covdiv = covcor %>%
  inner_join(ensmap)
divgenes = covdiv$rho
names(divgenes) = covdiv$GeneID # 7929 genes

library(clusterProfiler)
#BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
dc_gse = gseGO(geneList = sort(divgenes, decreasing = T), OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "BP", 
               pvalueCutoff = 1,
               keyType = "ENSEMBL", nPerm = 1000, minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BH',
               verbose = F)
dc_gse@result[1:5,1:10]
sum(dc_gse@result$p.adjust<0.1) # 0
saveRDS(dc_gse, './results/GTEx/dico_gse.rds')

##### gse without di background:
glist = covcor$rho
names(glist) = covcor$GeneID
dc_gse2 = gseGO(geneList = sort(glist, decreasing = T), OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "BP", 
               pvalueCutoff = 1,
               keyType = "ENSEMBL", nPerm = 1000, minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BH',
               verbose = F)
dc_gse@result[1:5,1:10]
sum(dc_gse@result$p.adjust<0.1) # 0
saveRDS(dc_gse, './results/GTEx/co_gse.rds')

###
#### Effect size:
source("scripts/functions.R")

unique(expvals$Age)
expvals %>% filter(Age==25) %>%
  group_by(major_tissue) %>%
  summarise(n=length(unique(id)) ) #  1 sample from each tissue
expvals %>% filter(Age%in%c(25, 35)) %>%
  group_by(major_tissue) %>%
  summarise(n=length(unique(id)) ) #  4 sample from each tissue


expp = expvals %>% 
  filter(Age%in%c(25, 35)) %>%
  dplyr::select(GeneID, Expression, id, major_tissue) %>%
  mutate(Grp = paste(id, major_tissue, sep='-')) %>%
  dplyr::select(-id, -major_tissue) %>%
  spread(key='Grp', value='Expression') %>%
  column_to_rownames(var = 'GeneID') %>%
  as.matrix

ts.ord = sapply(strsplit(colnames(expp),'-'),`[[`,3)

ES = sapply(unique(ts.ord), function(y){
  sapply(rownames(expp), function(x){
    cohens_d(expp[x, ts.ord == y], expp[x, ts.ord != y] )
  })
})
saveRDS(ES, file='data/other_datasets/GTEx/effectsize.rds')
# get tissue with highest ES for each gene:
ts.spec = colnames(ES)[apply(ES, 1, which.max)]
names(ts.spec) = rownames(ES)

# get tissue specific genes using >3Q of genes assigned to a tissue:
ts.specQ3 = sapply(unique(ts.spec), function(x){
  ts.genes = names(which(ts.spec==x))
  cutoff = summary(ES[ts.genes,x])['3rd Qu.']
  q3.genes = names(which(ES[ts.genes,x] > cutoff))
  return(q3.genes)
})
saveRDS(ts.specQ3,'./data/other_datasets/GTEx/ts.specQ3.genes.rds')

# get ES values for tissue specific genes:
ts.spec.ES = sapply(names(ts.specQ3),function(x) {
  ES[ts.specQ3[[x]],]
}, simplify = F)

ts.spec.ES2 = reshape2::melt(ts.spec.ES) %>% 
  set_names(c('gene id', 'tissue', 'ES','spec'))

saveRDS(ts.spec.ES2,'./data/other_datasets/GTEx/ts.spec.ES.rds')

ts.specQ3.genes = unlist(ts.specQ3)
names(ts.specQ3.genes) = gsub('[0-9]','', names(ts.specQ3.genes))
sapply(ts.specQ3, length)
##### in which tissue the highest expression change occurs for each gene:

### use beta from linear regression:

cortex = expvals %>%
  filter(major_tissue=='Brain') %>%
  select(GeneID, Expression, id) %>%
  spread(key='id', value = 'Expression') %>%
  column_to_rownames(var = 'GeneID') %>% as.matrix()
lung = expvals %>%
  filter(major_tissue=='Lung') %>%
  select(GeneID, Expression, id) %>%
  spread(key='id', value = 'Expression') %>%
  column_to_rownames(var = 'GeneID') %>% as.matrix()
liver = expvals %>%
  filter(major_tissue=='Liver') %>%
  select(GeneID, Expression, id) %>%
  spread(key='id', value = 'Expression') %>%
  column_to_rownames(var = 'GeneID') %>% as.matrix()
muscle = expvals %>%
  filter(major_tissue=='Muscle') %>%
  select(GeneID, Expression, id) %>%
  spread(key='id', value = 'Expression') %>%
  column_to_rownames(var = 'GeneID') %>% as.matrix()

cortage = attr %>% 
  filter(id%in%colnames(cortex) & major_tissue=='Brain') %>% 
  pull(age)
cortage = log2(c(25,35,45,55,65,75)[cortage])
lungage = attr %>% 
  filter(id%in%colnames(lung) & major_tissue=='Lung') %>% 
  pull(age)
lungage = log2(c(25,35,45,55,65,75)[lungage])
liverage = attr %>% 
  filter(id%in%colnames(liver) & major_tissue=='Liver') %>% 
  pull(age)
liverage = log2(c(25,35,45,55,65,75)[liverage])
muscleage = attr %>% 
  filter(id%in%colnames(muscle) & major_tissue=='Muscle') %>% 
  pull(age)
muscleage = log2(c(25,35,45,55,65,75)[muscleage])

cortbeta = t(apply(cortex,1,function(x){
  summary(lm(x~cortage))$coef[2,c(1,4)]
}))
lungbeta = t(apply(lung,1,function(x){
  summary(lm(x~lungage))$coef[2,c(1,4)]
}))
liverbeta = t(apply(liver,1,function(x){
  summary(lm(x~liverage))$coef[2,c(1,4)]
}))
musclebeta = t(apply(muscle,1,function(x){
  summary(lm(x~muscleage))$coef[2,c(1,4)]
}))

expbeta = cbind(cortbeta[,1], lungbeta[,1], liverbeta[,1], musclebeta[,1])
colnames(expbeta) = c('Cortex', 'Lung', 'Liver', 'Muscle')

ts.expr.ch = colnames(expbeta)[apply(expbeta, 1, function(x) which.max(abs(x)))]
names(ts.expr.ch) = rownames(expbeta)

# direction of expression change for those genes :
ts.expr.ch.dir = sapply(1:length(ts.expr.ch), function(x){ sign(expbeta[x, ts.expr.ch[x]]) } )
names(ts.expr.ch.dir) = names(ts.expr.ch)

########
######## Expr change in native tissue :
########
# for tissue-specific genes (>Q3) :
mat = data.frame(sameness = names(ts.specQ3.genes) == ts.expr.ch[ts.specQ3.genes],
                 exp_dir = ts.expr.ch.dir[ts.specQ3.genes] )
table(mat)[,c(2,1)]
fisher.test(table(mat)[,c(2,1)])
# OR: 0.8552, p = 0.02 (eski)
# OR: 1.628626, p =6.743e-12

saveRDS(list(tbl = table(mat)[,c(2,1)],
             fisher = fisher.test(table(mat)[,c(2,1)])),
        file = 'data/other_datasets/GTEx/specloss_fisher.rds')

### using only Co genes:
## losing expr in native, gain in other tissues:
cogenes = covcor %>% filter(rho<0) %>% pull(GeneID)

specsub = ts.specQ3.genes[ts.specQ3.genes%in%cogenes]
expchsub = ts.expr.ch[cogenes]
expchdirsub = ts.expr.ch.dir[cogenes]
matsub = data.frame(sameness = names(specsub) == expchsub[specsub],
                    expdir = expchdirsub[specsub])
table(matsub)[,c(2,1)]
sum(table(matsub)[,c(2,1)]) # 2407 genes
fisher.test(table(matsub)[,c(2,1)])
fisher.test(table(matsub)[,c(2,1)])$p.val
# OR = 7.20875, p = 7.04019e-87

saveRDS(list(tbl = table(matsub)[,c(2,1)],
             fisher = fisher.test(table(matsub)[,c(2,1)])),
        file = 'data/other_datasets/GTEx/specloss_fisher_co.rds')

covcor %>% 
  filter(adjusted_p<0.1)
# no sig genes

## DiCo vs DiDi and tis spec vs non-tis spec:
devdiv = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
divg = names(devdiv)
ensmap = readRDS('./results/GTEx/ensmap.rds')
divg = ensmap[ensmap$ENS_mm%in%divg,1] # 7976

convg = as.character(covcor %>% filter(rho < 0) %>% pull(GeneID))
dicog = intersect(divg, convg) # 4681

spec.pat.mat = data.frame(pat = divg%in%dicog,
                          ts.spec = divg%in%ts.specQ3.genes)
table(spec.pat.mat)
fisher.test(table(spec.pat.mat ))
fisher.test(table(spec.pat.mat ))$p.val
# OR = 0.8956, p = 1.072e-8 (old)
# OR: 1.107809, p = 0.05218

saveRDS(list(table(spec.dc.mat),fisher.test(table(spec.dc.mat))),
        file = 'data/other_datasets/GTEx/tisspec_dico_fisher.rds')

save(list=ls(),file = './results/GTEx/data.RData')

