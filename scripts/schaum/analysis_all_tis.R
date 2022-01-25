library(tidyverse)
library(ggpubr)
library(RColorBrewer)

expr = readRDS('data/other_datasets/schaum/rawexp.rds')
sinfo = readRDS('data/other_datasets/schaum/sinfo.rds')
pntnorm <- (1/0.352777778)
sexcol = c('f' = '#EE6352', 'm' = '#2b2b45')
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))

sinfo = sinfo %>%
  rename('Tissue'=Organ) %>%
  mutate(Tissue = ifelse(Tissue=='Limb_Muscle', 'Muscle', Tissue)) %>%
  mutate(age = factor(Age, levels=sort(unique(sinfo$Age))))

unique(sinfo$Tissue) # 17 different tissue
sinfo = sinfo %>% 
  filter(Age!=1) 

sinfo %>% group_by(Mouse_ID) %>% summarise(n=n()) %>% unique() %>% nrow() #  50 individuals

## remove outliers 21m5, 3m7:
sinfo = sinfo %>% filter(!Mouse_ID%in%c('3m7', '21m5'))

sinfo
expr = expr[,colnames(expr)%in%sinfo$Sample_Name]
dim(expr) # 799

libsize = reshape2::melt(colSums(expr)) %>%
  rownames_to_column() %>%
  set_names('Sample_Name', 'Lib_Size')
sinfo %>%
  left_join(libsize) %>%
  group_by(Tissue) %>%
  summarise(cor = cor.test(Lib_Size, Age, m='s')$est,
            pval = cor.test(Lib_Size, Age, m='s')$p.val)

## get only protein coding genes:
pcoding=readRDS('data/other_datasets/schaum/pcodinggenes.rds')
expr = expr[rownames(expr)%in%pcoding,] # 20562

### count filtration
expr = expr[, colSums(expr) > 4e6] # 657

sinfo = sinfo %>%
  filter(Sample_Name%in%colnames(expr))
  
sinfo %>%
  filter(Sample_Name%in%colnames(expr)) %>%
  group_by(Mouse_ID) %>%
  summarise(n=length(Tissue )) %>% 
  filter(n==17) # only 1 sample has 17 tissue samples

sinfo %>%
  group_by(Mouse_ID) %>%
  summarise(n=n())

samp_ids = lapply(unique(sinfo$Tissue),function(x) unique(filter(sinfo, Tissue==x)$Mouse_ID ) )
names(samp_ids) = unique(sinfo$Tissue)

samp_ids = samp_ids[!names(samp_ids)%in%c('Brown_Fat', 'Gonadal_Fat', 'Mesenteric_Fat')]

samp_ids = reshape2::melt(samp_ids) %>%
  set_names('Mouse_ID', 'Tissue') %>%
  mutate(val = 1) %>%
  spread(key = Tissue, value = val, fill = 0) %>%
  as.data.frame()
rownames(samp_ids) = samp_ids$Mouse_ID
samp_ids$Mouse_ID = NULL
table(rowSums(samp_ids))

pheatmap::pheatmap(t(samp_ids), cellwidth = 2, cellheight = 10, show_colnames = F,
                   color = c('gray99','lightblue'), legend = F, cutree_rows = 2, 
                   filename = 'results/schaum/8tis/sample_tissue_heatmap.png')
pheatmap::pheatmap(t(samp_ids), cellwidth = 2, cellheight = 10, show_colnames = F,
                   color = c('gray99','lightblue'), legend = F, cutree_rows = 2, 
                   filename = 'results/schaum/8tis/sample_tissue_heatmap.pdf')

sampage = as.numeric(sapply(strsplit(rownames(samp_ids),'[a-z]'),`[`,1))
order(sampage)
pheatmap::pheatmap(t(samp_ids[order(sampage),]), cellwidth = 4, cellheight = 10, show_colnames = T, 
                   fontsize_col = 5, gaps_col = c(5,11,17,23,29,35,41,45),
                   color = c('gray99','lightblue'), legend = F, cutree_rows = 7, cluster_cols = F,
                   filename = 'results/schaum/8tis/sample_tissue_heatmap2.pdf')

pheatmap::pheatmap(t(samp_ids[order(sampage),]), cellwidth = 4, cellheight = 10, show_colnames = T, 
                   fontsize_col = 5, gaps_col = c(5,11,17,23,29,35,41,45),
                   color = c('gray99','lightblue'), legend = F, cutree_rows = 7, cluster_cols = F,
                   filename = 'results/schaum/8tis/sample_tissue_heatmap2.png')

# plot(hclust(dist(t(samp_ids))))
# cutree(hclust(dist(t(samp_ids))),7)
#tislist = names(which(cutree(hclust(dist(t(samp_ids))),8)==2)) # 10 tissue
#tislist = tislist[!tislist%in%c('Brown_Fat', 'Spleen','Muscle')]
tislist = c('Liver', 'Heart', 'Kidney', 'Brain', 'Lung', 'Muscle', 'Subcutaneous_Fat', 'Spleen')

sinfo = sinfo %>%
  filter(Tissue%in%tislist)
expr = expr[,colnames(expr)%in%sinfo$Sample_Name]

incids = unique((sinfo %>%
  group_by(Mouse_ID) %>%
  summarise(n = length(unique(Tissue))) %>%
  mutate(inc = n==length(tislist)) %>%
  filter(inc))$Mouse_ID) # compile the list of individuals with samples in included tissues

# 26 individuals
incSamps = (sinfo %>% filter(Mouse_ID%in% incids))$Sample_Name
# get the sample ids for the individuals to be included

sinfo = sinfo %>% filter(Sample_Name %in% incSamps) # 208 sample
expr = expr[,colnames(expr)%in%incSamps]

age_dist = sinfo %>%
  group_by(age,Sex) %>%
  summarise(n=length(unique(Mouse_ID))) %>%
  ggplot(aes(x=age, y=n, fill=Sex)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=sexcol) +
  ylab('') + xlab('Age (in months)')

ggsave('results/schaum/8tis/age_dist.pdf', age_dist, width = 10, height = 10, useDingbats=F)
ggsave('results/schaum/8tis/age_dist.png', age_dist, width = 10, height = 10)

sinfo %>%
  group_by(Mouse_ID) %>%
  summarise(n=n())

sinfo %>%
  group_by(Mouse_ID) %>%
  summarise(n=unique(Sex)) %>% count(n)

# filter 0 genes in more than 25% of samples:
sum(rowMeans(expr==0)>0.25) # 2943
exprf = expr[!rowMeans(expr==0)>0.25,] # 17619, 208

## quantile normalisation:
expqn = preprocessCore::normalize.quantiles(as.matrix(log2(exprf+1)))
dimnames(expqn) = dimnames(exprf)
boxplot(expqn[,sample(56,10)])
rm(exprf)

## PCA
pcx = prcomp(t(expqn), scale=T)
pca_dat = data.frame(pcx$x[,1:4], Sample_Name=rownames(pcx$x)) %>%
  left_join(sinfo)
pcimp = round(summary(pcx)$imp[2,1:4]*100,0)

tissuecol = colorRampPalette(brewer.pal(8,'Dark2'))( length(tislist) )

pc12 = pca_dat %>%
  ggplot(aes(x=PC1, y=PC2, color=Tissue, size=age)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range=c(1, 2.5)) +
  #coord_fixed(ratio = pcimp[2]/pcimp[1], clip='off') +
  xlab(paste('PC 1 (',pcimp[1],'%)',sep='')) +
  ylab(paste('PC 2 (',pcimp[2],'%)',sep='')) +
  guides(color = guide_legend('Tissue'),
         size = guide_legend('Age')) +
  theme(legend.position = 'right')

#pca_dat %>% filter(PC2 > 50 & PC2 < 70) #  3m7 lung

pc34 = pca_dat %>%
  ggplot(aes(x=PC3, y=PC4, color=Tissue, size=age)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range=c(1, 2.5)) +
  #coord_fixed(ratio = pcimp[4]/pcimp[3], clip='off') +
  xlab(paste('PC 3 (',pcimp[3],'%)',sep='')) +
  ylab(paste('PC 4 (',pcimp[4],'%)',sep='')) +
  theme(legend.position = 'right')

pca_dat %>% filter( PC3 < 10 & PC3 > -10 ) #  3m7 lung

pc1_age = pca_dat %>%
  ggplot(aes(x = Age, y = PC1, color = Tissue, size=age)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range=c(1, 2.5)) +
  geom_smooth(se = F, size = 0.3, show.legend = T, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y') +
  ggtitle('PC1') + ylab(NULL) +
  theme(strip.text = element_text(size=8/pntnorm),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.1))

pc2_age = pca_dat %>%
  ggplot(aes(x = Age, y = PC2, color = Tissue, size=age)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range=c(1, 2.5)) +
  geom_smooth(se = F, size = 0.3, show.legend = T, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~.,  scales = 'free_y') +
  ggtitle('PC2') + ylab(NULL) +
  theme(strip.text = element_text(size=8/pntnorm),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.1))

pc3_age = pca_dat %>%
  ggplot(aes(x = Age, y = PC3, color = Tissue, size=age)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range=c(1, 2.5)) +
  geom_smooth(se = F, size = 0.3, show.legend = T, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~.,  scales = 'free_y') +
  ggtitle('PC3') +ylab(NULL) +
  theme(strip.text = element_text(size=8/pntnorm),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.1))

pc4_age = pca_dat %>%
  ggplot(aes(x = Age, y = PC4, color = Tissue, size=age)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  scale_size_discrete(range=c(1, 2.5)) +
  geom_smooth(se = F, size = 0.3, show.legend = T, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~.,  scales = 'free_y') +
  ggtitle('PC4') +ylab(NULL) +
  theme(strip.text = element_text(size=8/pntnorm),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.1))

pwise_distMean = function(mat, id_col=1){
  # mat: matrix or data frame with columns as coordinates to calculate distance, and one id column
  # distance is calculated pairwise among rows of same ids
  # id_col: index of the id column
  
  sapply(unique(mat[,id_col] ), function(x){
    ind = mat[mat[,id_col]==x,]
    mean(as.vector(dist(ind[, -id_col])))
  })
}

dist_dat = pca_dat %>% 
  dplyr::select(PC1,PC2,PC3,PC4, Mouse_ID)
mdist = pwise_distMean(dist_dat, id_col = 5)
mdist = reshape2::melt(mdist) %>%
  rownames_to_column('Mouse_ID') %>%
  left_join(unique(dplyr::select(sinfo, Mouse_ID, Age, Sex )))

cor.test(mdist$value, mdist$Age, m='s')

eucdist = mdist %>%
  ggplot(aes(x = Age, y= value)) +
  geom_point() +
  scale_size_continuous(range = c(0.2,1.5)) +
  geom_smooth(method = 'lm', se = F, color='darkred') +
  stat_cor(method = 'spearman', size = 6/pntnorm, cor.coef.name = 'rho') +
  ylab('Mean Pairwise Euclidean Distance')
eucdist

#mdist %>% filter(value>150) # 21m5

pcaplots1 = ggarrange(pc12, pc34, eucdist, nrow = 1, ncol = 3 , labels = c('a.','b.','c.'), 
                      font.label = list(size = 8), widths = c(2,2,2), common.legend = T, legend = 'right')

pcaplots2 = ggarrange(pc1_age, pc2_age, pc3_age,pc4_age, nrow = 1, ncol =4, hjust = c(-0.5),
                      labels = c('d.','e.','f.','g.'), font.label = list(size = 8), legend = 'none')

pca_plots = ggarrange(pcaplots1, pcaplots2, ncol=1,nrow=2, heights = c(1.5, 2.5))

ggsave('./results/schaum/8tis/pca.pdf', pca_plots, units='cm', width = 16, height = 16, 
       useDingbats=F)
ggsave('./results/schaum/8tis/pca.png', pca_plots, units='cm', width = 16, height = 16,
       bg='white')

saveRDS(pca_dat, 'results/source_data/f2/fs19_pca.rds')
saveRDS(mdist, 'results/source_data/f2/fs19_eucdist.rds')

### Age-related expression change ###

#expmat (expqn) final processed data to be used
expmat = reshape2::melt(expqn) %>%
  set_names('gene_id','Sample_Name', 'expression') %>%
  left_join(sinfo)

agecor = expmat %>%
  group_by(Tissue, gene_id) %>%
  summarise(rho = cor.test(expression, Age, m='s')$est,
            p = cor.test(expression, Age, m='s')$p.val)

saveRDS(agecor, 'data/other_datasets/schaum/8tis/age_exp_cors.rds')
### CoV Analysis ###

# calculate CoV:

covdat = expmat %>%
  group_by(gene_id, Mouse_ID) %>%
  summarise(CoV = sd(expression)/mean(expression))

covdat = covdat %>%  ungroup()

covdat = covdat %>% left_join(sinfo)

# mean CoV:
sumcovdat = covdat %>%
  dplyr::select(Mouse_ID, CoV, Age, gene_id)

sumcovdat = sumcovdat %>%  unique()

# CoV change:
covch = sumcovdat %>%
  group_by(gene_id) %>%
  summarise(rho = cor.test(CoV, Age, m='s')$est,
            pval = cor.test(CoV, Age, m='s')$p.val )
covch = covch %>%
  mutate(fdr = p.adjust(pval, method = 'BH'))

covch %>%
  filter(fdr<0.1) %>%
  mutate(pattern = ifelse(rho<0 , 'conv', 'div')) %>%
  group_by(pattern) %>%
  summarise(n = n())
# pattern     n
# <chr>   <int>
# 1 conv      132
# 2 div       112

saveRDS(covch, 'data/other_datasets/schaum/8tis/covch.rds')

covch %>%
  mutate(dir = ifelse(rho<0 ,'conv','div') ) %>%
  group_by(dir) %>%
  summarise(n = n())
# dir       n
# <chr> <int>
# 1 conv   8656
# 2 div    8963

sumCoV = sumcovdat %>%
  group_by(Mouse_ID, Age) %>%
  summarise(meancov = mean(CoV),
            medcov = median(CoV)) %>% ungroup()

cor.test(sumCoV$Age, sumCoV$meancov, m='s')
cor.test(sumCoV$Age, sumCoV$medcov, m='s')

meancovplot = ggplot(sumCoV, aes(x = Age, y = meancov)) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_smooth(method='lm') +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) +
  xlab(NULL) + ylab('Mean CoV') +
  theme(axis.text.x = element_text(angle=90))
meancovplot

medcovplot = ggplot(sumCoV,aes(x = Age, y = medcov)) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_smooth(method='lm') +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) +
  xlab(NULL) + ylab('Median CoV') +
  theme(axis.text.x = element_text(angle=90))
medcovplot

# pairwise correlation coefficients between tissues:
pairwisedat = expmat %>%
  dplyr::select(Tissue, expression, gene_id, Mouse_ID, Age) %>%
  unique()
indvids = unique(pairwisedat$Mouse_ID)

pexpdat = pairwisedat %>% 
  spread(key=Tissue, value=expression)

cors = lapply(indvids, function(ind){
  indexp = filter(pexpdat, Mouse_ID == ind) %>% select(-c(1:3)) %>% as.data.frame()
  cor(as.matrix(indexp), method = 'spearman')
})
names(cors) = indvids
cors

pairwisedat = reshape2::melt(cors) %>%
  set_names(c('tis1','tis2','Spearman Correlation Coefficient','Mouse_ID')) %>%
  filter(tis1!=tis2) %>%
  left_join(unique( dplyr::select(sinfo,Mouse_ID, Age) ))

####
pcors = pairwisedat %>%
  ggplot(aes(y=`Spearman Correlation Coefficient`, x=Age)) +
  geom_point(size=0.4) +
  geom_smooth(method = 'lm') +
  #facet_grid(tis1~tis2, scales = 'free_y') +
  facet_wrap(tis1~tis2, scales = 'free_y') +
  stat_cor(method='spearman', cor.coef.name = 'rho', size=6/pntnorm)
pcors
ggsave('results/schaum/8tis/pcors.pdf', pcors, width = 20, height = 20, units = 'cm')
ggsave('results/schaum/8tis/pcors.png', pcors, width = 20, height = 20, units = 'cm')

pairwisecor_ch = pairwisedat %>%
  #mutate(Age = c(25,35,45,55,65,75)[age]) %>%
  group_by(tis1,tis2) %>%
  summarise(corx = list(cor.test(`Spearman Correlation Coefficient`, Age, method = 'spearman'))) %>%
  ungroup() 

meancordat = pairwisedat %>%
  #mutate(Age = c(25,35,45,55,65,75)[age]) %>%
  group_by(tis1, tis2) %>% 
  summarise(meancor = mean(`Spearman Correlation Coefficient`, na.rm=T)) %>%
  ungroup() 

pairwisecor_ch$rho = sapply(pairwisecor_ch$corx, function(x) x$est)
pairwisecor_ch$p = sapply(pairwisecor_ch$corx,function(x)x$p.val)
pairwisecor_ch$adjusted_p = p.adjust(pairwisecor_ch$p, method = 'BH')
table(pairwisecor_ch$adjusted_p<0.1)
# XX significant
pairwisecor_ch %>% filter(adjusted_p<0.1) %>%
  summarise(sign(rho)) %>% table() # XX significant increasing correlation

pairwiseplotdat = pairwisecor_ch %>%
  unique() %>%
  full_join(select(meancordat, meancor, tis1, tis2))

pairwisecor_ch %>%
  mutate(tis1= as.character(tis1),
         tis2=as.character(tis2)) %>%
  filter( !duplicated(paste0(pmax(tis1,tis2),pmin(tis1,tis2) ) ) ) %>%
  pull(rho) %>% sign() %>% table() # 16 positive

pairwisecor_ch %>%
  mutate(tis1= as.character(tis1),
         tis2=as.character(tis2)) %>%
  filter( !duplicated(paste0(pmax(tis1,tis2),pmin(tis1,tis2) ) ) & adjusted_p<0.1 ) %>%
  pull(rho) %>% sign() %>% table() # 5 positive, 2 negative sig change

pairwisecor_ch %>%
  mutate(tis1= as.character(tis1),
         tis2=as.character(tis2)) %>%
  filter( !duplicated(paste0(pmax(tis1,tis2),pmin(tis1,tis2) ) ) & rho<0)
## 9/12 involves muscle or fat, 75%


varcol = setNames(c('dodgerblue','firebrick3'),c('Divergent','Convergent'))

pairwiseplot = pairwisecor_ch %>%
  unique() %>%
  full_join(select(meancordat, meancor, tis1, tis2)) %>%
  ggplot(aes(x = tis1, y = tis2, color = rho, size=meancor)) +
  geom_point() +
  geom_point(data=pairwisecor_ch %>%
               filter(adjusted_p<0.1),
             pch=21, size=5, color = 'purple') +
  scale_color_gradient2(low = varcol['Divergent'], mid = 'gray90', midpoint = 0, high = varcol['Convergent'])+
  xlab(NULL) + ylab(NULL) +
  scale_size_continuous(range= c(0.5,5), guide = guide_legend('Mean\nSimilarity')) +
  theme(legend.position = 'right', axis.text.x = element_text(angle=90),
        legend.key.size = unit(0.2,'cm'))
pairwiseplot

p1 = ggarrange(meancovplot, medcovplot, ncol = 1, nrow = 2, labels = c('a.','b.'), font.label = list(size = 8))
covresplot = ggarrange(p1, pairwiseplot, ncol =2, nrow = 1, labels = c(NA,'c.'), font.label = list(size = 8), 
                       widths = c(1,3))

ggsave('results/schaum/8tis/CoV.pdf', covresplot, units='cm', width = 16, height = 12, useDingbats=F)
ggsave('results/schaum/8tis/CoV.png', covresplot, units='cm', width = 16, height = 12)

saveRDS(sumCoV, 'results/source_data/f2/fs20_sumcov.rds')
saveRDS(pairwisecor_ch, 'results/source_data/f2/fs20_pwisecors.rds')

## ES:
source("scripts/functions.R")

expmat %>%
  filter(Age==3) %>%
  group_by(Tissue) %>%
  summarise(length(unique(Mouse_ID)))

expp = expmat %>%
  filter(Age==3) %>%
  dplyr::select(gene_id, expression, Mouse_ID, Tissue ) %>%
  mutate(Grp = paste(Mouse_ID, Tissue, sep='-' )) %>%
  dplyr::select(-Mouse_ID, -Tissue) %>%
  spread(key='Grp', value='expression' ) %>%
  column_to_rownames(var = 'gene_id') %>%
  as.matrix()

ts.ord = sapply(strsplit(colnames(expp),'-'),`[[`,2)

ES = sapply(unique(ts.ord), function(y){
  sapply(rownames(expp), function(x){
    cohens_d(expp[x, ts.ord == y], expp[x, ts.ord != y] )
  })
})
saveRDS(ES, file='data/other_datasets/schaum/8tis/effectsize.rds')

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

saveRDS(ts.specQ3,'./data/other_datasets/schaum/8tis/ts.specQ3.genes.rds')

# get ES values for tissue specific genes:
ts.spec.ES = sapply(names(ts.specQ3),function(x) {
  ES[ts.specQ3[[x]],]
}, simplify = F)

ts.spec.ES2 = reshape2::melt(ts.spec.ES) %>% 
  set_names(c('gene id', 'tissue', 'ES','spec'))

saveRDS(ts.spec.ES2,'./data/other_datasets/schaum/8tis/ts.spec.ES.rds')

ts.specQ3.genes = unlist(ts.specQ3)
names(ts.specQ3.genes) = gsub('[0-9]','', names(ts.specQ3.genes))

##### in which tissue the highest expression change occurs for each gene:

### use beta from linear regression:
unique(ts.ord)
cortex = expmat %>%
  filter(Tissue=='Brain') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()
heart = expmat %>%
  filter(Tissue=='Heart') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()
kidney = expmat %>%
  filter(Tissue=='Kidney') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()
liver = expmat %>%
  filter(Tissue=='Liver') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()
lung = expmat %>%
  filter(Tissue=='Lung') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()
muscle = expmat %>%
  filter(Tissue=='Muscle') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()
spleen = expmat %>%
  filter(Tissue=='Spleen') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()
fat = expmat %>%
  filter(Tissue=='Subcutaneous_Fat') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()

age = setNames(sinfo$Age, sinfo$Sample_Name)
cortage = log2(age[colnames(cortex)]*30)
heartage = log2(age[colnames(heart)]*30)
kidneyage = log2(age[colnames(kidney)]*30)
liverage = log2(age[colnames(liver)]*30)
lungage = log2(age[colnames(lung)]*30)
muscleage = log2(age[colnames(muscle)]*30)
spleenage = log2(age[colnames(spleen)]*30)
fatage = log2(age[colnames(fat)]*30)

cortbeta = t(apply(cortex,1,function(x){
  summary(lm(x~cortage))$coef[2,c(1,4)]
}))
heartbeta = t(apply(heart,1,function(x){
  summary(lm(x~heartage))$coef[2,c(1,4)]
}))
kidneybeta = t(apply(kidney,1,function(x){
  summary(lm(x~kidneyage))$coef[2,c(1,4)]
}))
liverbeta = t(apply(liver,1,function(x){
  summary(lm(x~liverage))$coef[2,c(1,4)]
}))
lungbeta = t(apply(lung,1,function(x){
  summary(lm(x~lungage))$coef[2,c(1,4)]
}))
musclebeta = t(apply(muscle,1,function(x){
  summary(lm(x~muscleage))$coef[2,c(1,4)]
}))
spleenbeta = t(apply(spleen,1,function(x){
  summary(lm(x~spleenage))$coef[2,c(1,4)]
}))
fatbeta = t(apply(fat,1,function(x){
  summary(lm(x~fatage))$coef[2,c(1,4)]
}))

expbeta = cbind(cortbeta[,1], heartbeta[,1],  kidneybeta[,1], liverbeta[,1],
                lungbeta[,1], musclebeta[,1], spleenbeta[,1], fatbeta[,1])
colnames(expbeta) = c('Cortex', 'Heart','Kidney', 'Liver', 'Lung', 'Muscle', 'Spleen', 'Fat')

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
# OR: 3.718371, p = 2.2e-16
saveRDS(list(tbl = table(mat)[,c(2,1)],
             fisher = fisher.test(table(mat)[,c(2,1)])),
        file = 'data/other_datasets/schaum/8tis/specloss_fisher.rds')

### using only Co genes:
# losing expr. in native tissue, gain expr. in other tissues:
cogenes = covch %>% filter(rho<0) %>% pull(gene_id)

specsub = ts.specQ3.genes[ts.specQ3.genes%in%cogenes]
expchsub = ts.expr.ch[cogenes]
expchdirsub = ts.expr.ch.dir[cogenes]
matsub = data.frame(sameness = names(specsub) == expchsub[specsub],
                    expdir = expchdirsub[specsub])
table(matsub)[,c(2,1)]
sum(table(matsub)[,c(2,1)]) # 2380 genes
fisher.test(table(matsub)[,c(2,1)])
fisher.test(table(matsub)[,c(2,1)])$p.val
# OR = 84.20069, p = 9.733053e-96

saveRDS(list(tbl = table(matsub)[,c(2,1)],
             fisher = fisher.test(table(matsub)[,c(2,1)])),
        file = 'data/other_datasets/schaum/8tis/specloss_fisher_co.rds')

covch %>% 
  filter(fdr<0.1) %>%
  mutate(pattern = ifelse(rho<0, 'conv', 'div')) %>%
  group_by(pattern) %>%
  summarise(n=n())
# pattern     n
# <chr>   <int>
# 1 conv      132
# 2 div       112
132/(132+112) # 0.5409

save(list=ls(), file = 'data/other_datasets/schaum/8tis/analysis.rdata')

