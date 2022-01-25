library(tidyverse)
library(ggpubr)
library(RColorBrewer)

## 1 ens genler (denendi ayni sonuclar)
#  2 protein kodlayan genler
#  3 count filtration

expr = readRDS('data/other_datasets/schaum/rawexp.rds')
sinfo = readRDS('data/other_datasets/schaum/sinfo.rds')
pntnorm <- (1/0.352777778)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))

sexcol = c('f' = '#EE6352', 'm' = '#2b2b45')
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'), c('Brain','Lung','Liver','Muscle'))

sinfo = sinfo %>%
  rename('Tissue'=Organ) %>%
  mutate(Tissue = ifelse(Tissue=='Limb_Muscle', 'Muscle', Tissue)) %>%
  mutate(age = factor(Age, levels=sort(unique(sinfo$Age))))

sinfo4 = sinfo %>%
  filter(Tissue%in%c('Muscle','Brain','Liver','Lung')) %>%
  filter(Age!=1)

sinfo4 = sinfo4 %>%
  group_by(Mouse_ID) %>%
  summarise(n=n()) %>%
  filter(n==4) %>%
  left_join(sinfo4) %>%
  select(-n)

# 3m7 was outlier, remove it
sinfo4 = sinfo4 %>% filter(!Mouse_ID%in%c('3m7'))

sinfo4 %>%
  group_by(Mouse_ID) %>%
  summarise(n=n()) %>%
  filter(n!=4) # checked.

age_dist = sinfo4 %>% 
  group_by(age, Sex) %>%
  summarise(n= length(unique(Mouse_ID)) ) %>%
  ggplot(aes(x=age, y=n, fill=Sex)) +
  geom_bar(stat='identity', position = position_stack()) +
  scale_fill_manual(values=sexcol) +
  ylab('') + xlab('Age in months')

ggsave('results/schaum/4tissue/age_dist.pdf', age_dist, units='cm', width = 10, height = 10, useDingbats=F)
ggsave('results/schaum/4tissue/age_dist.png', age_dist, units='cm', width = 10, height = 10)

rm(sinfo)
sinfo4
expr = expr[,colnames(expr)%in%sinfo4$Sample_Name]
dim(expr) # 180

libsize = reshape2::melt(colSums(expr)) %>%
  rownames_to_column() %>%
  set_names('Sample_Name', 'Lib_Size')
sinfo4 %>%
  left_join(libsize) %>%
  group_by(Tissue) %>%
  summarise(cor = cor.test(Lib_Size, Age, m='s')$est,
            pval = cor.test(Lib_Size, Age, m='s')$p.val)
# Tissue     cor  pval
# <chr>    <dbl> <dbl>
# 1 Brain   0.0152 0.921
# 2 Liver  -0.0152 0.921
# 3 Lung    0.154  0.312
# 4 Muscle  0.0652 0.670

## get ens genes:
# library(biomaRt)
# mart = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
# convert = getBM(filters = "uniprot_gn_symbol", attributes = c("uniprot_gn_symbol","ensembl_gene_id"),
#                 mart = mart, values = rownames(expr) )
# dup1 = unique(convert[duplicated(convert[,1]),1])
# dup2 = unique(convert[duplicated(convert[,2]),2])
# uconvert = convert[!convert[,1]%in%dup1,]
# uuconvert = uconvert[!uconvert[,2]%in%dup2,]
# head(uuconvert)
# idmap = uuconvert$uniprot_gn_symbol
# names(idmap) = uuconvert$ensembl_gene_id
# 
# sgenes = intersect(uuconvert$uniprot_gn_symbol, rownames(expr))
# expr = expr[sgenes,]
# identical(unname(idmap[idmap%in%sgenes]), rownames(expr)) # true
# rownames(expr) = names(idmap[idmap%in%sgenes])
# dim(expr) # 20549
 
# get only protein coding genes:
library(biomaRt)
listMarts()
mus.mart = useMart(biomart = 'ensembl')
mus.ensembl = useDataset('mmusculus_gene_ensembl', mart = mus.mart )
grep('uniprot',listAttributes(mart = mus.ensembl)[,1])
listAttributes(mart = mus.ensembl)[grep('uniprot',listAttributes(mart = mus.ensembl)[,1]),]
genetype = getBM(attributes = c('uniprot_gn_symbol', 'gene_biotype'), filters = 'uniprot_gn_symbol',
                  values = rownames(expr), mart = mus.ensembl)
# genetype = getBM(attributes = c('ensembl_gene_id', 'gene_biotype'), filters = 'ensembl_gene_id',
#                  values = '', mart = mus.ensembl)
pcoding = unique(genetype$uniprot_gn_symbol[genetype$gene_biotype=='protein_coding'])
#pcoding = unique(genetype$ensembl_gene_id[genetype$gene_biotype=='protein_coding'])
# saveRDS(pcoding, 'data/other_datasets/schaum/pcodinggenes_ens.rds')
# pcoding = readRDS('data/other_datasets/schaum/pcodinggenes_ens.rds')
saveRDS(pcoding, 'data/other_datasets/schaum/pcodinggenes.rds')
pcoding = readRDS('data/other_datasets/schaum/pcodinggenes.rds')
#pcoding=readRDS('data/other_datasets/schaum/pcodinggenes.rds')
exprp = expr[rownames(expr)%in%pcoding,] # 20562

### count filtration
exprp = exprp[, colSums(exprp) > 4e6]

sinfo4 = sinfo4 %>%
  filter(Sample_Name%in%colnames(exprp)) %>%
  group_by(Mouse_ID) %>%
  summarise(n=n()) %>%
  filter(n==4) %>%
  left_join(sinfo4) %>%
  dplyr::select(-n)

sinfo4 %>%
  group_by(Mouse_ID) %>%
  summarise(n=n()) %>%
  filter(n!=4) # checked.

exprp = exprp[, colnames(exprp)%in%sinfo4$Sample_Name] #  148 samples

libsize = reshape2::melt(colSums(exprp)) %>%
  rownames_to_column() %>%
  set_names('Sample_Name', 'Lib_Size')
sinfo4 %>%
  left_join(libsize) %>%
  group_by(Tissue) %>%
  summarise(cor = cor.test(Lib_Size, Age, m='s')$est,
            pval = cor.test(Lib_Size, Age, m='s')$p.val)
# Tissue      cor  pval
# <chr>     <dbl> <dbl>
# 1 Brain   0.00967 0.955
# 2 Liver  -0.0940  0.580
# 3 Lung    0.0699  0.681
# 4 Muscle -0.146   0.389

# filter 0 genes in more than 25% of samples:
sum(rowMeans(exprp==0)>0.25) # 3756
exprf = exprp[!rowMeans(exprp==0)>0.25,] # 16806

## quantile normalisation:
expqn = preprocessCore::normalize.quantiles(as.matrix(log2(exprf+1)))
dimnames(expqn) = dimnames(exprf)
boxplot(expqn[,sample(144,10)])

## PCA
pcx = prcomp(t(expqn), scale=T)
pca_dat = data.frame(pcx$x[,1:4], Sample_Name=rownames(pcx$x)) %>%
  left_join(sinfo4)
pcimp = round(summary(pcx)$imp[2,1:4]*100,0)

pc12 = pca_dat %>%
  ggplot(aes(x=PC1, y=PC2, color=Tissue, size=age)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = tissuecol) +
  #scale_size_continuous(range=c(1, 2.5), trans='log2') +
  scale_size_discrete(range=c(0.5, 2.5)) +
  #coord_fixed(ratio = pcimp[2]/pcimp[1], clip='off') +
  xlab(paste('PC 1 (',pcimp[1],'%)',sep='')) +
  ylab(paste('PC 2 (',pcimp[2],'%)',sep='')) +
  guides(color = guide_legend('Tissue'),
         size = guide_legend('Age'))

pca_dat %>% filter(Tissue=='Lung' & PC1>25) # 3m7 seems to be an outlier

pc34 = pca_dat %>%
  ggplot(aes(x=PC3, y=PC4, color=Tissue, size=age)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = tissuecol) +
  #scale_size_continuous(range=c(0.2,1.5), trans='log2') +
  scale_size_discrete(range=c(0.5, 2.5)) +
  #coord_fixed(ratio = pcimp[4]/pcimp[3], clip='off') +
  xlab(paste('PC 3 (',pcimp[3],'%)',sep='')) +
  ylab(paste('PC 4 (',pcimp[4],'%)',sep=''))

pca_dat %>% filter(Tissue=='Lung' & PC3>0) # 3m7 seems to be an outlier, same as above

pc1_age = pca_dat %>%
  ggplot(aes(x = Age, y = PC1, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = T, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y') +
  ggtitle('PC1') +ylab(NULL) 

pc2_age = pca_dat %>%
  ggplot(aes(x = Age, y = PC2, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = T, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~.,  scales = 'free_y') +
  ggtitle('PC2') + ylab(NULL) 

pc3_age = pca_dat %>%
  ggplot(aes(x = Age, y = PC3, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = T, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y') +
  ggtitle('PC3') +ylab(NULL)

pc4_age = pca_dat %>%
  ggplot(aes(x = Age, y = PC4, color = Tissue)) +
  geom_jitter(size = 0.1, width = 0.3) +
  scale_color_manual(values = tissuecol) +
  geom_smooth(se = F, size = 0.3, show.legend = T, method = 'lm') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
           size = 6/pntnorm, show.legend = F) +
  facet_grid(Tissue~., scales = 'free_y') +
  ggtitle('PC4') +ylab(NULL) 

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
  left_join(unique(dplyr::select(sinfo4, Mouse_ID, Age, Sex )))

cor.test(mdist$value, mdist$Age, m='s')
# rho = 0.1259, p=0.4575

eucdist = mdist %>%
  ggplot(aes(x = Age, y= value)) +
  geom_point() +
  scale_size_continuous(range = c(0.2,1.5)) +
  geom_smooth(method = 'lm', se = F, color='darkred') +
  stat_cor(method = 'spearman', size = 6/pntnorm, cor.coef.name = 'rho') +
  ylab('Mean Pairwise Euclidean Distance')
eucdist

mdist %>% filter(value<165) # 3m7 seem to be an outlier, same as above

pcaplots1 = ggarrange(pc12, pc34, eucdist, nrow = 1, ncol =3 , labels = c('a.','b.','c.'), 
                      font.label = list(size = 8), widths = c(2,2,1.5), common.legend = T, legend = 'right')
pcaplots2 = ggarrange(pc1_age, pc2_age, pc3_age,pc4_age, nrow = 1, ncol =4 , 
                      labels = c('d.','e.','f.','g.'), font.label = list(size = 8), legend = 'none')

pca_plots = ggarrange(pcaplots1, pcaplots2, ncol=1,nrow=2, heights = c(2,2.5))

ggsave('./results/schaum/4tissue/pca.pdf', pca_plots , units = 'cm', width = 14, height = 10, useDingbats=F)
ggsave('./results/schaum/4tissue/pca.png', pca_plots , units = 'cm', width = 14, height = 10, bg='white')

saveRDS(pca_dat, 'results/source_data/f2/fs17_pca.rds')
saveRDS(mdist, 'results/source_data/f2/fs17_eucdist.rds')

# rm Mouse_ID: 3m7

# ^
## go back to start

### Age-related expression change ###

#expmat (expqn) final processed data to be used
expmat = reshape2::melt(expqn) %>%
  set_names('gene_id','Sample_Name', 'expression') %>%
  left_join(sinfo4)

expmat %>%
  group_by(gene_id, Mouse_ID) %>%
  summarise(n=length(Tissue)) %>%
  filter(n<4)

agecor = expmat %>%
  group_by(Tissue, gene_id) %>%
  summarise(rho = cor.test(expression, Age, m='s')$est,
            p = cor.test(expression, Age, m='s')$p.val)

saveRDS(agecor, 'data/other_datasets/schaum/4tissue/age_exp_cors.rds')
saveRDS(expmat, 'data/other_datasets/schaum/4tissue/expr.rds')

agecor_pearson = expmat %>%
  group_by(Tissue, gene_id) %>%
  summarise(cor = cor.test(expression, Age, m='pearson')$est,
            p = cor.test(expression, Age, m='pearson')$p.val)

saveRDS(agecor_pearson, 'data/other_datasets/schaum/4tissue/age_exp_cors_pearson.rds')

### CoV Analysis ###

# calculate CoV:
covdat = expmat %>%
  group_by(gene_id, Mouse_ID) %>%
  summarise(CoV = sd(expression)/mean(expression)) %>%
  ungroup()

covdat = covdat %>%
  left_join(sinfo4)

sumCoV = covdat %>%
  dplyr::select(Mouse_ID, CoV, Age, gene_id) %>%
  unique()

saveRDS(sumCoV, 'data/other_datasets/schaum/4tissue/cov.rds')

# CoV change:
covch = sumCoV %>%
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
# 1 conv      170
# 2 div       149
# total = 319
# 170/(170+149)

saveRDS(covch, 'data/other_datasets/schaum/4tissue/covch.rds')

covch %>%
  mutate(dir = ifelse(rho<0 ,'conv','div') ) %>%
  group_by(dir) %>%
  summarise(n = n())
# dir       n
# <chr> <int>
# 1 conv   7761
# 2 div    9045

# mean CoV:
summaryCoV = sumCoV %>%
  group_by(Mouse_ID, Age) %>%
  summarise(meancov = mean(CoV),
            medcov = median(CoV)) %>% ungroup()

cor.test(summaryCoV$Age, summaryCoV$meancov, m='s')
cor.test(summaryCoV$Age, summaryCoV$medcov, m='s')

meancovplot = ggplot(summaryCoV, aes(x = Age, y = meancov)) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_smooth(method='lm') +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) +
  xlab(NULL) + ylab('Mean CoV') +
  theme(axis.text.x = element_text(angle=90))
meancovplot

medcovplot = ggplot(summaryCoV,aes(x = Age, y = medcov)) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_smooth(method='lm') +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho',  size = 6/pntnorm) +
  xlab(NULL) + ylab('Median CoV') +
  theme(axis.text.x = element_text(angle=90))
medcovplot

# pairwise correlation coefficients between tissues:
pairwisedat = expmat %>%
  dplyr::select(Tissue, expression, gene_id, Mouse_ID, Age) %>%
  unique()

pairwisedat = pairwisedat %>%
  spread(key = Tissue, value = expression) %>%
  group_by(Mouse_ID, Age) %>%
  summarise(`Brain-Liver` = cor(Brain, Liver, method = 'spearman'),
            `Brain-Lung` = cor(Brain, Lung, method = 'spearman'),
            `Brain-Muscle` = cor(Brain, Muscle, method = 'spearman'),
            `Liver-Lung` = cor(Liver, Lung, method = 'spearman'),
            `Liver-Muscle` = cor(Liver, Muscle, method = 'spearman'),
            `Lung-Muscle` = cor(Lung, Muscle, method = 'spearman'))

cortests = list(`Brain-Liver` = cor.test(pairwisedat$Age, pairwisedat$`Brain-Liver`, m='s'),
     `Brain-Lung` = cor.test(pairwisedat$Age, pairwisedat$`Brain-Lung`, m='s'),
     `Brain-Muscle` = cor.test(pairwisedat$Age, pairwisedat$`Brain-Muscle`, m='s'),
     `Liver-Lung` = cor.test(pairwisedat$Age, pairwisedat$`Liver-Lung`, m='s'),
     `Liver-Muscle` = cor.test(pairwisedat$Age, pairwisedat$`Liver-Muscle`, m='s'),
     `Lung-Muscle` = cor.test(pairwisedat$Age, pairwisedat$`Lung-Muscle`, m='s'))
cortests = t(sapply(cortests,function(x) c(x$est, x$p.val)))
colnames(cortests)[2] = 'p'
cortests = cbind(cortests, 'BH' = p.adjust(cortests[,2], method = 'BH') )
# rho           p          BH
# Brain-Liver   0.3655462 0.026091742 0.078275225
# Brain-Lung    0.3275705 0.047801439 0.095602879
# Brain-Muscle  0.2224804 0.185667258 0.227485387
# Liver-Lung    0.1247944 0.461776239 0.461776239
# Liver-Muscle -0.2205697 0.189571156 0.227485387
# Lung-Muscle  -0.5069401 0.001364814 0.008188883

pairwiseplot = pairwisedat %>%
  gather(key = 'type', value = 'Correlation', -Mouse_ID, -Age) %>%
  ggplot(aes(x = Age, y = Correlation)) +
  #geom_boxplot(width = 0.4, outlier.shape = NA, fill = 'gray70') +
  geom_smooth(method='lm') +
  geom_jitter(width = 0.1, size = 0.3) +
  stat_cor(aes(x = Age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) +
  facet_wrap(~type, scales = 'free_y') +
  ylab('Spearman\'s Correlation Coefficient') + xlab(NULL) +
  theme(axis.text.x = element_text(angle=90))
pairwiseplot

p1 = ggarrange(meancovplot, medcovplot, ncol = 1, nrow = 2, labels = c('a.','b.'), font.label = list(size = 8))
covresplot = ggarrange(p1, pairwiseplot, ncol =2, nrow = 1, labels = c(NA,'c.'), font.label = list(size = 8), 
                       widths = c(1,3))

ggsave('./results/schaum/4tissue/covchange.pdf', covresplot, units='cm', width = 16, height = 12, useDingbats=F)
ggsave('./results/schaum/4tissue/covchnage.png', covresplot, units='cm', width = 16, height = 12, bg='white')

saveRDS(summaryCoV, 'results/source_data/f2/fs18_sumcov.rds')
saveRDS(pairwisedat, 'results/source_data/f2/fs18_pwisecors.rds')

## mean pairwise exp cors:
mpwise = pairwisedat %>% 
  gather(key = 'pairs', value='rho', -Mouse_ID,-Age) %>%
  group_by(Mouse_ID, Age) %>%
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
mpwise

ggsave('./results/schaum/4tissue/meanpwisecor.pdf', mpwise, units='cm', width = 16, height = 12, useDingbats=F)
ggsave('./results/schaum/4tissue/meanpwisecor.png', mpwise, units='cm', width = 16, height = 12, bg='white')

#### Effect size:
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
ts = ts.ord[1]

ES = sapply(unique(ts.ord), function(y){
  sapply(rownames(expp), function(x){
    cohens_d(expp[x, ts.ord == y], expp[x, ts.ord != y] )
  })
})
saveRDS(ES, file='data/other_datasets/schaum/4tissue/effectsize.rds')
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

saveRDS(ts.specQ3,'./data/other_datasets/schaum/4tissue/ts.specQ3.genes.rds')

# get ES values for tissue specific genes:
ts.spec.ES = sapply(names(ts.specQ3),function(x) {
  ES[ts.specQ3[[x]],]
}, simplify = F)

ts.spec.ES2 = reshape2::melt(ts.spec.ES) %>% 
  set_names(c('gene id', 'tissue', 'ES','spec'))

saveRDS(ts.spec.ES2,'./data/other_datasets/schaum/4tissue/ts.spec.ES.rds')

ts.specQ3.genes = unlist(ts.specQ3)
names(ts.specQ3.genes) = gsub('[0-9]','', names(ts.specQ3.genes))

##### in which tissue the highest expression change occurs for each gene:

### use beta from linear regression:

cortex = expmat %>%
  filter(Tissue=='Brain') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()
lung = expmat %>%
  filter(Tissue=='Lung') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()
liver = expmat %>%
  filter(Tissue=='Liver') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()
muscle = expmat %>%
  filter(Tissue=='Muscle') %>%
  select(gene_id, expression, Sample_Name) %>%
  spread(key='Sample_Name', value = 'expression') %>%
  column_to_rownames(var = 'gene_id') %>% as.matrix()

age = setNames(sinfo4$Age, sinfo4$Sample_Name)
cortage = log2(age[colnames(cortex)]*30)
lungage = log2(age[colnames(lung)]*30)
liverage = log2(age[colnames(liver)]*30)
muscleage = log2(age[colnames(muscle)]*30)

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
# OR: 1.49678, p = 2.307e-7
saveRDS(list(tbl = table(mat)[,c(2,1)],
        fisher = fisher.test(table(mat)[,c(2,1)])),
        file = 'data/other_datasets/schaum/4tissue/specloss_fisher.rds')

### using only Co genes:
## loss of expr in native, gain expression in other tissues:

cogenes = covch %>% filter(rho<0) %>% pull(gene_id)

specsub = ts.specQ3.genes[ts.specQ3.genes%in%cogenes]
expchsub = ts.expr.ch[cogenes]
expchdirsub = ts.expr.ch.dir[cogenes]
matsub = data.frame(sameness = names(specsub) == expchsub[specsub],
                    expdir = expchdirsub[specsub])
table(matsub)[,c(2,1)]
sum(table(matsub)[,c(2,1)]) # among 2124 genes
fisher.test(table(matsub)[,c(2,1)])
fisher.test(table(matsub)[,c(2,1)])$p.val
# OR = 58.025, p = 1.505944e-197

saveRDS(list(tbl = table(matsub)[,c(2,1)],
             fisher = fisher.test(table(matsub)[,c(2,1)])),
        file = 'data/other_datasets/schaum/4tissue/specloss_fisher_co.rds')

covch %>% 
  filter(fdr<0.1) %>%
  mutate(pattern = ifelse(rho<0, 'conv', 'div')) %>%
  group_by(pattern) %>%
  summarise(n=n())
# pattern     n
# <chr>   <int>
# 1 conv      170
# 2 div       149
170/(170+149) # %53


## DiCo vs DiDi and tis spec vs non-tis spec:
devdiv = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
divg = names(devdiv)
idmap = readRDS('data/uniprot_to_ENS.rds')
divg = idmap[idmap$ensembl_gene_id%in%divg,1] # 7765

convg = covch %>% filter(rho < 0) %>% pull(gene_id)
dicog = intersect(divg, convg) # 3439

spec.pat.mat = data.frame(pat = divg%in%dicog,
                          ts.spec = divg%in%ts.specQ3.genes)
table(spec.pat.mat )
fisher.test(table(spec.pat.mat ))
fisher.test(table(spec.pat.mat ))$p.val
# OR = 1.33, p = 1.072e-8

saveRDS(list(table(spec.dc.mat),fisher.test(table(spec.dc.mat))),
        file = 'data/other_datasets/schaum/4tissue/tisspec_dico_fisher.rds')

######## DiCo enrichment for Schaum with dev. divergent genes from our dataset:
########

ddc_genes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
idconvert = readRDS('data/other_datasets/schaum/4tissue/idconvert.rds')

covdiv = covch %>%
  right_join(idconvert) %>%
  filter(ensembl_gene_id%in%(names(ddc_genes)))

divgenes = covdiv$rho
names(divgenes) = covdiv$ensembl_gene_id

library(clusterProfiler)
require(org.Mm.eg.db)
dc_gse = gseGO(geneList = sort(divgenes, decreasing = T), OrgDb = org.Mm.eg.db, ont = "BP", 
               pvalueCutoff = 1,
               keyType = "ENSEMBL", nPerm = 1000, minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BH',
               verbose = F)
dc_gse@result[1:5,1:10]
sum(dc_gse@result$p.adjust<0.1) # 0
saveRDS(dc_gse, './data/other_datasets/schaum/4tissue/dico_gse.rds')

save(list=ls(), file = 'data/other_datasets/schaum/4tissue/analysis.rdata')
