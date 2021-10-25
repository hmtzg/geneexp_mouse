library(Rvislib)
library(tidyverse)
library(ggpubr)
library(ggforce)
theme_set(theme_rvis(base_size = 6, legend.pos = 'right'))
pntnorm <- (1/0.352777778)
sexcolors <- setNames(c('#ffc9b5', '#8fb8de'), c('Female', 'Male'))
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('Up','Down'))
# revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))
revcol = setNames(c('gray25', 'gray25', 'gray25', 'gray25'), c('UpDown','DownUp','UpUp','DownDown'))

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
  ggsave(file = paste(prefix, '.png', sep = ''), plot = ggobj, units = 'cm', 
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
attr = attr %>%
  filter(minor_tissue %in% c('Lung','Liver','Brain - Cortex','Muscle - Skeletal')) %>%
  filter(sample_id %in% samplelist) %>%
  group_by(id,minor_tissue) %>%
  summarise(n = length(unique(sample_id))) %>%
  ungroup() %>%
  mutate(exc = n>1) %>%
  right_join(attr) 
table(attr$exc) # no duplication
genelist = unique(unlist(lapply(expvals,function(expx){sapply(strsplit(as.character(rownames(expx)),'[.]'),function(x)x[1])})))
martx = biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
biotype = biomaRt::getBM(attributes = c('ensembl_gene_id', 'gene_biotype'), filters = 'ensembl_gene_id', values = genelist, mart = martx)
codinggenes = unique(filter(biotype,gene_biotype == 'protein_coding')$ensembl_gene_id)

incids = setdiff(unique((sapply(expvals, colnames) %>% 
  reshape2::melt() %>%
  set_names(c('sample_id','Tissue')) %>%
  left_join(attr) %>%
  group_by(id) %>%
  summarise( n = length(unique(Tissue))) %>%
  mutate(inc = n==4) %>%
  filter(inc))$id),NA) # compile the list of individuals with samples in all 4 tissues

inc_sampleIDs = unique((attr %>% filter(id %in% incids))$sample_id) # get the sample ids for the individuals to be included

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
duplicated(unname(sample2id[colnames(expx)]))


genelist = Reduce('intersect',lapply(expvals,rownames))
expvals = cbind(expvals[[1]][genelist,],expvals[[2]][genelist,],expvals[[3]][genelist,],expvals[[4]][genelist,]) #to make sure the gene order is the same
expvals = expvals [ ! (rowMeans(expvals==0)>0.25), ] # exclude the genes if expression value is 0 in more than  25% of the samples
dim(expvals)
#16211 - 188 (47 samples)

#normalization
exp_l2 = log2(expvals+1)
exp_l2_qn = preprocessCore::normalize.quantiles(exp_l2)
dimnames(exp_l2_qn) = dimnames(exp_l2)


pcx = prcomp(t(exp_l2_qn),scale=T)
pca_data = data.frame(pcx$x[,1:4], sample_id = rownames(pcx$x)) %>%
  left_join(attr) %>%
  mutate(Tissue = ifelse(major_tissue=='Brain','Cortex',major_tissue)) %>%
  mutate(Age = c(25,35,45,55,65,75)[age]) 
pcimp = round(summary(pcx)$imp[2,1:4]*100,2)
pc1_2 = pca_data %>%
  ggplot(aes( x = PC1, y= PC2 , color = Tissue, size = age)) +
  geom_point() +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range = c(0.2,1.5)) +
  xlab(paste('PC 1 (',pcimp[1],'%)',sep=''))+
  ylab(paste('PC 2 (',pcimp[2],'%)',sep=''))

pc3_4 = pca_data %>%
  ggplot(aes( x = PC3, y= PC4 , color = Tissue, size = age)) +
  geom_point() +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range = c(0.2,1.5)) +
  xlab(paste('PC 3 (',pcimp[3],'%)',sep=''))+
  ylab(paste('PC 4 (',pcimp[4],'%)',sep=''))

pc1_age = pca_data %>%
  ggplot(aes(x = Age, y = PC1, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = F, method = 'lm') +
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

idmap = setNames(attr$id, attr$sample_id)
#calculate the Eucledian distance using all PCs because there is no clear relationship between PC and age
meanEuc = reshape2::melt(as.matrix(dist(pcx$x))) %>%
  mutate(id1 = idmap[as.character(Var1)],
         id2 = idmap[as.character(Var2)]) %>%
  filter(id1 == id2) %>%
  filter(Var1!=Var2) %>%
  group_by(id1) %>%
  summarise(dist = mean(value)) %>%
  ungroup() %>% rename(id = id1) %>%
  left_join(attr)

eucdist = meanEuc %>%
  mutate(Age = c(25,35,45,55,65,75)[age]) %>%
  select(id, dist, Age, age) %>% unique( ) %>% 
  ggplot(aes(x = Age, y= dist)) +
  geom_point(aes(size = age)) +
  scale_size_discrete(range = c(0.2,1.5)) +
  geom_smooth(method = 'lm', se = F, color = 'darkred') +
  stat_cor(method = 'spearman', size = 6/pntnorm, cor.coef.name = 'rho') +
  ylab('Mean Pairwise Euclidean Distance')

fig2_supp9abdefg = pca_data %>% select(PC1, PC2, PC3, PC4, Tissue, age) %>% unique()
fig2_supp9c = eucdist$data
pcaplots1 = ggarrange(pc1_2, pc3_4, eucdist, nrow = 1, ncol =3 , labels = c('a','b','c'), font.label = list(size = 8), widths = c(2,2,1.5), common.legend = T, legend = 'right')
pcaplots2 = ggarrange(pc1_age, pc2_age, pc3_age,pc4_age, nrow = 1, ncol =4 , labels = c('d','e','f','g'), font.label = list(size = 8), legend = 'none')
pca_plots = ggarrange(pcaplots1, pcaplots2, ncol=1,nrow=2, heights = c(2,2.5))
plotsave(ggobj = pca_plots, prefix = './results/GTEx/fig2supp9',width = 16, height = 12)
tablesave(fig2_supp9abdefg, prefix = './results/GTEx/fig2supp9abdefg')
tablesave(fig2_supp9c, prefix = './results/GTEx/fig2supp9c')

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

#### CoV analysis #### 

#CoV calculation for each gene of each individual (takin sd and mean of 4 tissues)
expvals = expvals %>%
  group_by(GeneID, id) %>%
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
  theme_rvis(base_size = 6, x.text.angle = 90)
fig2_supp10ab = meancovplot$data
medcovplot = ggplot(sumCov,aes(x = age, y = medcov)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
  geom_jitter(width = 0.1, size = 0.3) +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', label.y = 0.31, size = 6/pntnorm) +
  xlab(NULL) + ylab('Median CoV')+
  theme_rvis(base_size = 6, x.text.angle = 90)


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
  ylab('Spearman Correlation Coefficient') + xlab(NULL) +
  theme_rvis(base_size = 6, x.text.angle = 90)
  
fig2_supp10c = pairwiseplot$data%>% select(-`.group`)

p1 = ggarrange(meancovplot, medcovplot, ncol = 1, nrow = 2, labels = c('a','b'), font.label = list(size = 8))
covresplot = ggarrange(p1, pairwiseplot, ncol =2, nrow = 1, labels = c(NA,'c'), font.label = list(size = 8), widths = c(1,3))

plotsave(ggobj = covresplot, prefix = './results/GTEx/CoV',width = 16, height = 8)
tablesave(fig2_supp10ab, prefix = './results/GTEx/fig2supp10ab')
tablesave(fig2_supp10c, prefix = './results/GTEx/fig2supp10c')

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
# 6985  9226 
#9226 convergent

save(list=ls(),file = './results/GTEx/data.RData')
