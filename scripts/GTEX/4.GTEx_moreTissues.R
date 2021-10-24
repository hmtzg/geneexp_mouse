library(tidyverse)
library(ggpubr)
library(ggforce)
pntnorm <- (1/0.352777778)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
sexcolors <- setNames(c('#ffc9b5', '#8fb8de'), c('Female', 'Male'))
tissuecol = c(`Adipose Tissue` = '#B2BFC2',
              `Blood` = '#A3392E',
              `Blood Vessel` = '#F75E4D',
              `Brain` = '#233789',
              `Lung` = '#f49e92',
              `Liver` = '#801008',
              `Muscle` = '#dbb32e',
              `Nerve` = '#074A28',
              Pituitary = '#556110',
              Skin = '#C28A46',
              Thyroid = '#7B1791')
varcol = setNames(c('dodgerblue','firebrick3'),c('Divergent','Convergent'))
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

attr = attr %>%
  group_by(major_tissue, minor_tissue) %>%
  summarise(n = length(unique(id))) %>%
  ungroup() %>%
  group_by(major_tissue) %>%
  top_n(1,n) %>%
  ungroup() %>%
  left_join(attr)

indids = lapply(unique(attr$minor_tissue),function(x)unique(filter(attr,minor_tissue==x)$id))
names(indids) = unique(attr$minor_tissue)
indids = reshape2::melt(indids) %>%
  set_names(c('id','minor_tissue')) %>%
  mutate(val = 1) %>%
  spread(key = minor_tissue,value = val, fill = 0) %>%
  as.data.frame()
rownames(indids) = indids$id
indids$id = NULL

pheatmap::pheatmap(t(indids), cellwidth = 2,cellheight = 10, cutree_rows = 3, show_colnames = F,
                   color = c('gray99','lightblue'),legend = F,
                   filename = './results/GTEx/alltissues/sample_tissue_heatmap.png')

pheatmap::pheatmap(t(indids), cellwidth = 2,cellheight = 10, cutree_rows = 3, show_colnames = F,
                   color = c('gray99','lightblue'),legend = F,
                   filename = './results/GTEx/alltissues/sample_tissue_heatmap.pdf')

pheatmap::pheatmap(t(indids), cellwidth = 2,cellheight = 10, cutree_rows = 3, show_colnames = F,
                   color = c('gray99','lightblue'),legend = F,
                   filename = './results/figure_supplements/fs2/FS13.png')

pheatmap::pheatmap(t(indids), cellwidth = 2,cellheight = 10, cutree_rows = 3, show_colnames = F,
                   color = c('gray99','lightblue'),legend = F,
                   filename = './results/figure_supplements/fs2/FS13.pdf')

tislist = names(which(cutree(hclust(dist(t(indids))), 3)==1))
fn = paste('./data/processed/GTEx/expression/tpm/',
           sapply(strsplit(gsub(' ','',tislist),'[(]'),function(x)x[[1]]),'.rds',sep='')

expvals = lapply(fn, function(fnx){readRDS(fnx)})
names(expvals) = tislist
samplelist = unique(unlist(lapply(expvals, function(x)colnames(x))))
attr = attr %>%
  # filter(minor_tissue %in% c('Lung','Liver','Brain - Cortex','Muscle - Skeletal')) %>%
  filter(sample_id %in% samplelist) %>%
  group_by(id,minor_tissue) %>%
  summarise(nsamp = length(unique(sample_id))) %>%
  ungroup() %>%
  mutate(exc = nsamp>1) %>%
  right_join(attr)  %>%
  filter(sample_id %in% samplelist)
table(attr$exc) # no duplication

genelist = unique(unlist(lapply(expvals,function(expx){
  sapply(strsplit(as.character(rownames(expx)),'[.]'),function(x)x[1])})))
martx = biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
biotype = biomaRt::getBM(attributes = c('ensembl_gene_id', 'gene_biotype'), filters = 'ensembl_gene_id',
                         values = genelist, mart = martx)
codinggenes = unique(filter(biotype,gene_biotype == 'protein_coding')$ensembl_gene_id)
saveRDS(codinggenes, file='./results/GTEx/alltissues/codinggenes.rds')
incids = setdiff(unique((sapply(expvals, colnames) %>% 
                           reshape2::melt() %>%
                           set_names(c('sample_id','Tissue')) %>%
                           left_join(attr) %>%
                           group_by(id) %>%
                           summarise( n = length(unique(Tissue))) %>%
                           mutate(inc = n==length(expvals)) %>%
                           filter(inc))$id),NA) # compile the list of individuals with samples in all 4 tissues

inc_sampleIDs = unique((attr %>% filter(id %in% incids))$sample_id) 
# get the sample ids for the individuals to be included

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
expvals = Reduce('cbind',lapply(expvals,function(x)x[genelist,]))
expvals = expvals [ ! (rowMeans(expvals==0)>0.25), ] 
# exclude the genes if expression value is 0 in more than  25% of the samples
dim(expvals)
#16305 - 350 (35 samples) -> 16290 350

#normalization
exp_l2 = log2(expvals+1)
exp_l2_qn = preprocessCore::normalize.quantiles(exp_l2)
dimnames(exp_l2_qn) = dimnames(exp_l2)

pcx = prcomp(t(exp_l2_qn),scale=T)
pca_data = data.frame(pcx$x[,1:4], sample_id = rownames(pcx$x)) %>%
  left_join(attr) %>%
  mutate(Tissue = major_tissue) %>%
  mutate(Age = c(25,35,45,55,65,75)[age]) 
pcimp = round(summary(pcx)$imp[2,1:4]*100,2)
pc1_2 = pca_data %>%
  ggplot(aes( x = PC1, y= PC2 , color = Tissue, size = age)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range = c(0.2,1.5)) +
  xlab(paste('PC 1 (',pcimp[1],'%)',sep='')) +
  ylab(paste('PC 2 (',pcimp[2],'%)',sep='')) +
  theme(legend.spacing.x = unit(0.1,'cm'),
        legend.key.size = unit(0.1,'cm'),
        legend.background = element_rect(color='gray70'))

pc3_4 = pca_data %>%
  ggplot(aes( x = PC3, y= PC4 , color = Tissue, size = age)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range = c(0.2,1.5)) +
  xlab(paste('PC 3 (',pcimp[3],'%)',sep='')) +
  ylab(paste('PC 4 (',pcimp[4],'%)',sep='')) +
  theme(legend.spacing.x = unit(0.1,'cm'),
        legend.key.size = unit(0.1,'cm'),
        legend.background = element_rect(color='gray70'))

pc1_age = pca_data %>%
  ggplot(aes(x = Age, y = PC1, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = F, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y') +
  ggtitle('PC1') + ylab(NULL) +
  theme(strip.text = element_text(size=6))

pc2_age = pca_data %>%
  ggplot(aes(x = Age, y = PC2, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = F, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y')+
  ggtitle('PC2') + ylab(NULL)+
  theme(strip.text = element_text(size=6))

pc3_age = pca_data %>%
  ggplot(aes(x = Age, y = PC3, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = F, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y')+
  ggtitle('PC3') +ylab(NULL) +
  theme(strip.text = element_text(size=6))

pc4_age = pca_data %>%
  ggplot(aes(x = Age, y = PC4, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = F, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y') +
  ggtitle('PC4') + ylab(NULL) +
  theme(strip.text = element_text(size=6))

idmap = setNames(attr$id, attr$sample_id)
#calculate the Eucledian distance using all PCs because there is no clear relationship between PC and age
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
  ggplot(aes(x = Age, y= value)) +
  geom_point(aes(size = age)) +
  scale_size_discrete(range = c(0.2,1.5)) +
  geom_smooth(method = 'lm', se = F, color = 'darkred') +
  stat_cor(method = 'spearman', size = 6/pntnorm, cor.coef.name = 'rho') +
  ylab('Mean Pairwise Euclidean Distance') +
  theme(axis.title.y = element_text(size=5))

pcaplots1 = ggarrange(pc1_2, pc3_4, eucdist, nrow = 1, ncol =3 , labels = c('a.','b.','c.'),
                      font.label = list(size = 8), widths = c(2,2,1.5), common.legend = T, legend = 'right')
pcaplots2 = ggarrange(pc1_age, pc2_age, pc3_age,pc4_age, nrow = 1, ncol =4 , labels = c('d.','e.','f.','g.'),
                      font.label = list(size = 8), legend = 'none')

pca_plots = ggarrange(pcaplots1, pcaplots2, ncol=1,nrow=2, heights = c(1.5,6))

plotsave(ggobj = pca_plots, prefix = './results/GTEx/alltissues/pca', width = 16, height = 26)
plotsave(ggobj = pca_plots, prefix = './results/figure_supplements/fs2/FS10', width = 16, height = 20)

#### Gene Expression Analysis ####

expvals = reshape2::melt(exp_l2_qn) %>%
  set_names(c('GeneID','sample_id','Expression')) %>%
  left_join(attr) %>%
  mutate(Tissue = major_tissue) %>%
  mutate(Age = c(25,35,45,55,65,75)[age]) 

#Expression - Age correlation
agecor = expvals %>%
  group_by(Tissue, GeneID) %>%
  summarise(rho = cor.test(Expression, Age, method = 's')$est,
            p = cor.test(Expression, Age, method = 's')$p.val)

agecor = ungroup(agecor) %>%
  mutate(adjusted_p = p.adjust(p, method = 'fdr'))

#### CoV analysis #### 

#CoV calculation for each gene of each individual (taking sd and mean of 10 tissues)
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
  theme(axis.text.x = element_text(angle = 90))

medcovplot = ggplot(sumCov,aes(x = age, y = medcov)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
  geom_jitter(width = 0.1, size = 0.3) +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) +
  xlab(NULL) + ylab('Median CoV') +
  theme(axis.text.x = element_text(angle = 90))


# Calculate pairwise correlation coefficients between tissues
pairwisedat = expvals %>%
  select(Tissue, Expression, GeneID, id, age, Age) %>%
  unique()
xx = t(combn(unique(pairwisedat$Tissue),2))
indvids = unique(pairwisedat$id)

pairwisedat = pairwisedat %>%
  spread(key = Tissue, value = Expression) 

cors = lapply(indvids, function(ind){
  indexp = filter(pairwisedat, id == ind) %>% select(-c(1:4)) %>% as.data.frame()
  cor(as.matrix(indexp), method = 'spearman')
})

names(cors) = indvids

pairwisedat = reshape2::melt(cors) %>%
  set_names(c('tis1','tis2','Spearman Correlation Coefficient','id')) %>%
  filter(tis1!=tis2) %>%
  left_join(unique( select(attr,id,age) ))

pairwisedat2 = pairwisedat %>%
  mutate(Age = c(25,35,45,55,65,75)[age]) %>%
  group_by(tis1,tis2) %>%
  # head()
  summarise(corx = list(cor.test(`Spearman Correlation Coefficient`, Age, method = 'spearman'))) %>%
  ungroup() 

meancordat = pairwisedat %>%
  mutate(Age = c(25,35,45,55,65,75)[age]) %>%
  group_by(tis1, tis2) %>% 
  summarise(meancor = mean(`Spearman Correlation Coefficient`)) %>%
  ungroup() 

pairwisedat2$rho = sapply(pairwisedat2$corx, function(x) x$est)
pairwisedat2$p = sapply(pairwisedat2$corx,function(x)x$p.val)
pairwisedat2$adjusted_p = p.adjust(pairwisedat2$p, method = 'BH')
table(pairwisedat2$adjusted_p<0.1)
# none significant

pairwiseplot = pairwisedat2 %>%
  unique() %>%
  full_join(select(meancordat, meancor, tis1, tis2)) %>%
  ggplot(aes(x = tis1, y = tis2, color = rho, size=meancor)) +
  geom_point() +
  scale_color_gradient2(low = varcol['Divergent'], mid = 'gray90', midpoint = 0, high = varcol['Convergent'])+
  xlab(NULL) + ylab(NULL) +
  scale_size_continuous(range= c(0.5,5), guide = guide_legend('Mean\nSimilarity')) +
  theme(legend.position = 'right', axis.text.x = element_text(angle=90),
        legend.key.size = unit(0.2,'cm'))

# pairwiseplot = pairwisedat  %>%
#   mutate(Age = c(25,35,45,55,65,75)[age]) %>%
#   ggplot(aes(x = age, y = `Spearman Correlation Coefficient`)) +
#   geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70', size = 0.3) +
#   geom_jitter(width = 0.1, size = 0.05) +
#   stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) +
#   facet_wrap(~type, scales = 'free_y', strip.position = 'top', ncol = 5) +
#   ylab('Spearman Correlation Coefficient') + xlab(NULL) +
#   theme_rvis(base_size = 6, x.text.angle = 90) 

p1 = ggarrange(meancovplot, medcovplot, ncol = 1, nrow = 2, labels = c('a.','b.'), font.label = list(size = 8))
covresplot = ggarrange(p1, pairwiseplot, ncol =2, nrow = 1, labels = c(NA,'c.'), font.label = list(size = 8), 
                       widths = c(1,3))
  
plotsave(ggobj = covresplot, prefix = './results/GTEx/alltissues/CoV',width = 16, height = 10)

plotsave(ggobj = covresplot, prefix = './results/figure_supplements/fs2/FS11',width = 16, height = 10)
  
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
# GeneID             rho          p adjusted_p
# <fct>            <dbl>      <dbl>      <dbl>
# 1 ENSG00000214102 -0.685 0.00000554     0.0451
# 2 ENSG00000183621 -0.657 0.0000180      0.0976
# 3 ENSG00000156968 -0.697 0.00000337     0.0451
table(covcor$rho<0)
# FALSE  TRUE 
# 6737  9568 
#9568 convergent

mostsignifplot = expvals %>%
  filter(GeneID %in% 'ENSG00000214102') %>%
  ggplot(aes(x = age, y = Expression, color = Tissue)) +
  geom_jitter(width = 0.1, size = 0.2) +
  geom_smooth(method = 'lm' , se= F, aes(x = as.numeric(age)-2), size = 0.5) +
  scale_color_manual(values = tissuecol)+
  ggtitle('WEE2') + xlab('Age') +
  #theme_rvis(base_size = 6, legend.pos = c(1,1), just.par = c(1,1)) +
  theme(plot.title = element_text(face = 'italic')) 
  
plotsave(ggobj = mostsignifplot, prefix = './results/GTEx/alltissues/signif',width = 8, height = 8)

save(list=ls(),file = './results/GTEx/alltissues/data.RData')
