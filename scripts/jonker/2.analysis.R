library(tidyverse)
library(ggpubr)
library(RColorBrewer)
exp = readRDS('./data/other_datasets/jonker/raw/exp_qn.rds')

expr = reshape2::melt(exp) %>%
  set_names('gene_id','sample_id','expression')
saveRDS(expr,'./data/other_datasets/jonker/processed/expression.rds')

age = readRDS('./data/other_datasets/jonker/raw/age.rds')
tissue.ord = readRDS('./data/other_datasets/jonker/raw/tissue_id.rds')
sample_info = readRDS('./data/other_datasets/jonker/raw/sample_info.rds')
sample_info = sample_info %>% 
  mutate(Tissue= str_to_title(tissue)) %>% 
  select(-tissue)

chstwo = brewer.pal(9,'Paired')[c(8,4)]
Tissuecol = setNames(c('#233789', '#f49e92', '#801008',chstwo),c('Brain','Lung','Liver','Kidney','Spleen'))
pntnorm <- (1/0.352777778)
theme_set(theme_pubr(base_size = 6, legend = 'right') +
            theme(legend.key.size = unit(2,'pt'),
                  axis.text = element_text(size=6),
                  axis.title = element_text(size=6)))

############### PCA
pc = prcomp(t(exp),scale = T)
sample_info = sample_info %>% 
  mutate(agew= age) %>%
  mutate(age= age*7) %>%
  rename('Age'=age, 'Tissue'=Tissue) %>%
  select(-log2age)

varpc = reshape2::melt(summary(pc)$imp[2,1:4]) %>%
  rownames_to_column() %>%
  rename('varExp'=value, 'PC'=rowname)

pca_data = reshape2::melt(pc$x[,1:4]) %>%
  set_names(c('sample_id','PC','value')) %>%
  left_join(varpc)

pcdat = pca_data %>%
  select(-varExp) %>%
  spread(key='PC', value='value') %>%
  left_join(sample_info, by='sample_id')

#####
pc12 = pcdat %>%
  ggplot(aes(x=PC1,y=PC2, size=Age, color=Tissue)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = Tissuecol) +
  scale_size_continuous(range = c(0.5,2), trans= 'log2') +
  xlab(paste('PC1 (', round(varpc[1,2]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(varpc[2,2]*100),'%)',sep=''))
 
pc34 = pcdat %>%
  ggplot(aes(x=PC3,y=PC4, size=Age, color=Tissue)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = Tissuecol) +
  scale_size_continuous(range = c(0.5,2), trans= 'log2') +
  #coord_fixed(ratio = varpc[4,2]/varpc[3,2], clip = 'off') +
  xlab(paste('PC3 (', round(varpc[3,2]*100),'%)',sep='')) +
  ylab(paste('PC4 (', round(varpc[4,2]*100),'%)',sep=''))

pc1age = pcdat %>%
  ggplot(aes(x = Age, y = PC1, color = Tissue)) +
  geom_smooth(alpha = 0.1, se=F, show.legend = F, method='lm') +
  geom_point(size = 0.5) +
  facet_grid(Tissue~., scales='free_y') +
  scale_color_manual(values = Tissuecol) +
  guides(color = F) +
  xlab('Age') +
  ylab(NULL) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F) +
  ggtitle('PC1')

pc2age = pcdat %>%
  ggplot(aes(x = Age, y = PC2, color = Tissue)) +
  geom_smooth(alpha = 0.1, se=F, show.legend = F, method='lm') +
  geom_point(size = 0.5) +
  facet_grid(Tissue~., scales='free_y') +
  scale_color_manual(values = Tissuecol) +
  guides(color = F) +
  xlab('Age') +
  ylab(NULL) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F) +
  ggtitle('PC2')

pc3age = pcdat %>%
  ggplot(aes(x = Age, y = PC3, color = Tissue)) +
  geom_smooth(alpha = 0.1, se=F, show.legend = F, method='lm') +
  geom_point(size = 0.5) +
  facet_grid(Tissue~., scales='free_y') +
  scale_color_manual(values = Tissuecol) +
  guides(color = F) +
  xlab('Age') +
  ylab(NULL) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F) +
  ggtitle('PC3')

pc4age = pcdat %>%
  ggplot(aes(x = Age, y = PC4, color = Tissue)) +
  geom_smooth(alpha = 0.1, se=F, show.legend = F, method='lm') +
  geom_point(size = 0.5) +
  facet_grid(Tissue~., scales='free_y') +
  scale_color_manual(values = Tissuecol) +
  guides(color = F) +
  xlab('Age') +
  ylab(NULL) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F) +
  ggtitle('PC4')

pc1234 = ggarrange(pc12, pc34, common.legend = T, legend = 'right', labels=c('a.','b.'),
                   font.label = list(size=8), widths = c(1,1))

pcages = ggarrange(pc1age,pc2age,pc3age,pc4age, nrow=1, ncol=4, labels=c('c.','d.','e.','f.'),
                   font.label = list(size=8), legend='none')

pcaplots = ggarrange(pc1234, pcages, ncol=1, nrow=2, heights = c(1,1))

ggsave('./results/Jonker/pca.pdf', pcaplots, units='cm', width = 16, height = 12, useDingbats=F)
ggsave('./results/Jonker/pca.png', pcaplots, units='cm', width = 16, height = 12)

pwise_distMean = function(mat, id_col=1){
  # mat: matrix or data frame with columns as coordinates to calculate distance, and one id column
  # distance is calculated pairwise among rows of same ids
  # id_col: index of the id column
  
  sapply(unique(mat[,id_col] ), function(x){
    ind = mat[mat[,id_col]==x,]
    mean(as.vector(dist(ind[, -id_col])))
  })
}

dist_dat = pcdat %>% 
  dplyr::select(PC1,PC2,PC3,PC4, ind_id)
mdist = pwise_distMean(dist_dat, id_col = 5)
mdist = reshape2::melt(mdist) %>%
  rownames_to_column('ind_id') %>%
  left_join(unique(dplyr::select(sample_info, ind_id, Age)))

cor.test(mdist$value, mdist$Age, m='s')

eucdist = mdist %>%
  ggplot(aes(x = Age, y= value)) +
  geom_point() +
  scale_size_continuous(range = c(0.2,1.5)) +
  geom_smooth(method = 'lm', se = F, color='darkred') +
  stat_cor(method = 'spearman', size = 6/pntnorm, cor.coef.name = 'rho') +
  ylab('Mean Pairwise Euclidean Distance')
eucdist

pc1234 = ggarrange(pc12, pc34, eucdist, ncol=3, common.legend = T, legend = 'right', labels=c('a.','b.','c.'),
                   font.label = list(size=8), widths = c(1,1))

############# pairwise distance.....

## age cors:

agecor = expr %>%
  left_join(sample_info) %>%
  group_by(Tissue, gene_id) %>%
  summarise(rho = cor.test(Age, expression, m='s')$est,
            pval = cor.test(Age, expression, m='s')$p.val)

saveRDS(agecor, 'data/other_datasets/jonker/processed/agecor.rds')

#################
############ CoV analysis
ids = as.character(sample_info %>% pull(ind_id))
ageu = sample_info %>% pull(Age)
ageu = ageu[1:18] 

genecov = sapply(unique(ids), function(x){
  sapply(rownames(exp),function(y){
    sd(exp[y, ids == x]) / mean(exp[y, ids == x])
  })
})

saveRDS(genecov,'./data/other_datasets/jonker/processed/CoV.rds')

genecovtidy = reshape2::melt(genecov) %>%
  set_names('gene_id','ind_id','CoV')

# CoV change with age:
covch = t(sapply(rownames(genecov),function(x){
  x = cor.test(genecov[x,], ageu, m ="s")
  c(x$est, x$p.val)
}))
covch = cbind(covch, p.adjust(covch[,2], method = "BH"))
colnames(covch)[c(2,3)] = c("pval", "BH")
covch = as.data.frame(covch) %>% rownames_to_column('gene_id') 
saveRDS(covch,'./data/other_datasets/jonker/processed/covch.rds')

covch %>% filter(BH<0.1) %>% mutate(pat = ifelse(rho<0,'cov', 'div')) %>% group_by(pat) %>% summarise(n=n())
##### 1735 significant cov changes
##### 1144 convergent
##### 591 divergent

##### mean CoV:
cov_dat_sum = genecovtidy %>%
  left_join(unique(select(sample_info,-Tissue,-sample_id)), by='ind_id') %>%
  group_by(ind_id, Age) %>%
  summarise(meanCoV = mean(CoV),
            medianCoV = median(CoV)) %>% ungroup()

saveRDS(cov_dat_sum,file = './data/other_datasets/jonker/processed/mean_cov.rds')

# mean CoV change with age:
cov_cordat = cov_dat_sum %>%
  summarise(cor = cor.test(meanCoV, Age, method = 's')$est,
            cor.p = cor.test(meanCoV, Age, method = 's')$p.val)
cov_cordat
# cor  cor.p
# <dbl>  <dbl>
# 1 -0.480 0.0440

# median CoV change with age:
cov_cordat_median = cov_dat_sum %>%
  summarise(cor = cor.test(medianCoV, Age, method = 's')$est,
            cor.p = cor.test(medianCoV, Age, method = 's')$p.val)
cov_cordat_median
# cor cor.p
# <dbl> <dbl>
# 1 -0.0282 0.912

# mean CoV change plot:
meancov = cov_dat_sum %>%
  ggplot(aes(x = Age, y= meanCoV)) +
  geom_smooth(se=T, method='lm', color = 'midnightblue', fill = 'lightblue', show.legend = F) +
  geom_point(size = 1, color='steelblue') +
  scale_color_manual(values = Tissuecol) +
  xlab('Age') +
  ylab('Mean CoV') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F)

# median CoV change plot:
mediancov = cov_dat_sum %>%
  ggplot(aes(x = Age, y= medianCoV)) +
  geom_smooth(se=T, method='lm', color = 'midnightblue', fill = 'lightblue', show.legend = F) +
  geom_point(size = 1, color='steelblue') +
  scale_color_manual(values = Tissuecol) +
  xlab('Age') +
  ylab('Median CoV') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F)

# top DiCo gene: ENSMUSG00000041560, ENSMUSG00000029552
covch %>% 
  top_n(n = 1, wt = -`rho`) %>%
  slice(2) %>% 
  left_join(expr) %>%
  left_join(sample_info) %>% 
  ggplot(aes(x=Age, y=expression, col=Tissue)) +
  geom_point() +
  geom_smooth(method='lm')

#################### pairwise Tissue corelations
tissue.ord =  str_to_title(tissue.ord)
exp2 = exp
colnames(exp2) = ids

pcors = list()
chs = combn(5,2)
for(i in 1:10){
  tsx = chs[,i]
  chsts = unique(tissue.ord)
  pname = paste(chsts[tsx[1]], chsts[tsx[2]], sep = '-')
  e1 = exp2[, tissue.ord == chsts[tsx[1]] ]
  e2 = exp2[, tissue.ord == chsts[tsx[2]] ]
  sameind = intersect(colnames(e1), colnames(e2))
  pair = sapply(sameind, function(x){ cor(e1[,x], e2[,x], m='s')})
  pair = data.frame(rho = pair, pair=names(pair), row.names=NULL)
  pcors[[pname]] = pair
}

# pairwise expression correlations change with age:
pcorch = sapply(names(pcors), function(x){
  c(round(cor.test(pcors[[x]][,1], ageu, m='s')$est,2),
    round(cor.test(pcors[[x]][,1], ageu, m='s')$p.val,3))
})
rownames(pcorch)[2] = 'p.val'

names(ageu) = ids[1:18]
pexpcors = reshape2::melt(pcors) %>%
  set_names('ind_id', 'variable','rho','pair') %>%
  select(-variable) %>%
  left_join(data.frame(Age=ageu,ind_id=names(ageu)))

saveRDS(pexpcors, './data/other_datasets/jonker/processed/pwise_exp_cors.rds')

pcortests = pexpcors %>%
  group_by(pair) %>%
  summarise(p = cor.test(rho, Age, m='s')$p.val,
            rho = cor.test(rho, Age, m='s')$est)

pcortests = cbind(pcortests, 'BH' = p.adjust(pcortests$p, method = 'BH'))
# pair          p         rho        BH
# 1   Kidney-Brain 0.50993343  0.16615726 0.7284763
# 2   Kidney-Liver 0.93112332 -0.02194530 0.9311233
# 3    Kidney-Lung 0.14162955  0.36052991 0.2832591
# 4  Kidney-Spleen 0.08011082  0.42323077 0.2832591
# 5    Liver-Brain 0.06205187  0.44831111 0.2832591
# 6     Liver-Lung 0.23000491  0.29782906 0.3833415
# 7   Liver-Spleen 0.09598649  0.40442051 0.2832591
# 8     Lung-Brain 0.59381084 -0.13480684 0.7422635
# 9   Spleen-Brain 0.89194359 -0.03448547 0.9311233
# 10   Spleen-Lung 0.12052850  0.37934017 0.2832591
#####
p2 = ggplot(pexpcors) +
  aes(x = Age, y = rho) +
  facet_wrap(~pair, scales = 'free', ncol=2) +
  geom_point(size=1, color="steelblue", alpha=0.9) +
  geom_smooth(se=T, method='lm', color = 'midnightblue', fill='lightblue') +
  scale_color_manual(values = Tissuecol) +
  ylab('Sperman correlation coefficient') +
  xlab('Age') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F)

p1  = ggarrange(meancov, mediancov, nrow=2 , labels=c('a.','b.'), font.label = list(size=8),vjust = c(1,-1) )

pwisecorplot = ggarrange(p1, p2, ncol=2, widths = c(2,3), labels= c(NA,'c.'),font.label = list(size=8))

ggsave('./results/figure_supplements/fs2/FS7.pdf', pwisecorplot, units='cm', width = 16, height = 11, 
       useDingbats=F)
ggsave('./results/figure_supplements/fs2/FS7.png', pwisecorplot, units='cm', width = 16, height = 11)

p2.2 = ggplot(pexpcors) +
  aes(x = Age, y = rho) +
  facet_wrap(~pair, scales = 'free', ncol=4) +
  geom_point(size=1, color="steelblue", alpha=0.9) +
  geom_smooth(se=T, method='lm', color = 'midnightblue', fill='lightblue') +
  scale_color_manual(values = Tissuecol) +
  ylab('Sperman correlation coefficient') +
  xlab('Age') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F)

p1  = ggarrange(meancov, mediancov, nrow=2 , labels=c('d.','e.'), font.label = list(size=8),vjust = c(1,-1) )
pwisecorplot2 = ggarrange(p1, p2.2, ncol=2, widths = c(1,3), labels= c(NA,'f.'),font.label = list(size=8))
fs7 = ggarrange(pc1234, pwisecorplot2, nrow=2, heights = c(1,1.3))

ggsave('./results/figure_supplements/fs2/FS7.pdf', fs7, units='cm', width = 16, height = 13, 
       useDingbats=F)
ggsave('./results/figure_supplements/fs2/FS7.png', fs7, units='cm', width = 16, height = 13, bg='white')

saveRDS(cov_dat_sum,'results/source_data/f2/fs7_mean_median.rds')
saveRDS(pexpcors,'results/source_data/f2/fs7_pexpcors.rds')

###
saveRDS(pcdat, 'results/source_data/f2/fs7pca.rds')
saveRDS(mdist, 'results/source_data/f2/fs7_eucdist.rds')


## mean pairwise corrs:
meancors = pexpcors %>% 
  group_by(ind_id, Age) %>%
  summarise(mean = mean(rho),
            median = median(rho)) %>%
  gather(key='method', value = 'rho', mean, median)

szx = 4
mcorsplot = meancors %>%
  ggplot(aes(x=Age,y=rho)) +
  facet_wrap(~method, strip.position = 'left', scales = 'free_y',
             labeller = as_labeller(c(mean = 'Mean P.wise Expr. Correlation',
                                      median = 'Median P.wise Expr. Correlation'))) +
  geom_point(size=1.5, color="steelblue", alpha=0.9) +
  geom_smooth(se=T, method = 'lm',color = 'midnightblue', fill='lightblue') +
  scale_x_continuous(trans = 'log2') +
  xlab('Age in days (in log2 scale)') + ylab(NULL) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 8)) + theme_bw()
ggsave('./results/Jonker/meanexpcors.pdf', mcorsplot, units = 'cm', width = 16, height = 10,
       useDingbats =F)
ggsave('./results/Jonker/meanexpcors.png', mcorsplot, units = 'cm', width = 16, height = 10)

########
######## DiCo enrichment for Jonker with dev. divergent genes from our dataset:
########

ddc_genes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')

covdiv =  covch[covch$gene_id%in%names(ddc_genes),] # 8377 genes
divgenes = covdiv$rho
names(divgenes) = covdiv$gene_id

library(clusterProfiler)
require(org.Mm.eg.db)
dc_gse = gseGO(geneList = sort(divgenes, decreasing = T), OrgDb = org.Mm.eg.db, ont = "BP", 
               pvalueCutoff = 1,
               keyType = "ENSEMBL", nPerm = 1000, minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BH',
               verbose = F)
dc_gse@result[1:5,1:10]
sum(dc_gse@result$p.adjust<0.1) # 705
saveRDS(dc_gse, './data/other_datasets/jonker/dico_gse.rds')

# dc_gse_genelist =  strsplit(dc_gse@result[,11], split = '/')
# names(dc_gse_genelist) = dc_gse@result[,'ID']
# dc_gse_genelist =reshape2::melt(dc_gse_genelist) %>%
#   set_names(c('gene_id','GO_ID'))
# dc_gse_table = list(enrichment = dc_gse@result[,1:10],
#                     genelist = dc_gse_genelist)
#write.xlsx(dc_gse_table, './results/Jonker/dico_gse_table.xlsx')

####
#####
#### Effect size:
source("scripts/functions.R")

age
expp = exp[,age==13]
ts = setNames(sample_info$Tissue, sample_info$sample_id)
ts.ord = ts[colnames(expp)]

ES = sapply(unique(tissue.ord), function(y){
  sapply(rownames(expp), function(x){
    cohens_d(expp[x, ts.ord == y], expp[x, ts.ord != y] )
  })
})
saveRDS(ES, file='data/other_datasets/jonker/effectsize.rds')
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

saveRDS(ts.specQ3,'./data/other_datasets/jonker/ts.specQ3.genes.rds')

# get ES values for tissue specific genes:
ts.spec.ES = sapply(names(ts.specQ3),function(x) {
  ES[ts.specQ3[[x]],]
}, simplify = F)

ts.spec.ES2 = reshape2::melt(ts.spec.ES) %>% 
  set_names(c('gene id', 'tissue', 'ES','spec'))

saveRDS(ts.spec.ES2,'./data/other_datasets/jonker/ts.spec.ES.rds')

ts.specQ3.genes = unlist(ts.specQ3)
names(ts.specQ3.genes) = gsub('[0-9]','', names(ts.specQ3.genes))

##### in which tissue the highest expression change occurs for each gene:
### use beta from linear regression:

cortex = exp[, tissue.ord=='Brain']
lung = exp[, tissue.ord=='Brain']
liver = exp[, tissue.ord=='Liver']
spleen = exp[, tissue.ord=='Spleen']
kidney = exp[, tissue.ord=='Kidney']

ages = setNames(sample_info$Age, sample_info$sample_id)

cortage = log2(ages[tissue.ord=='Brain'])
lungage = log2(ages[tissue.ord=='Lung'])
liverage = log2(ages[tissue.ord=='Liver'])
spleenage = log2(ages[tissue.ord=='Spleen'])
kidneyage = log2(ages[tissue.ord=='Kidney'])
  
cortbeta = t(apply(cortex,1,function(x){
  summary(lm(x~cortage))$coef[2,c(1,4)]
}))
lungbeta = t(apply(lung,1,function(x){
  summary(lm(x~lungage))$coef[2,c(1,4)]
}))
liverbeta = t(apply(liver,1,function(x){
  summary(lm(x~liverage))$coef[2,c(1,4)]
}))
spleenbeta = t(apply(spleen,1,function(x){
  summary(lm(x~spleenage))$coef[2,c(1,4)]
}))
kidneybeta = t(apply(kidney,1,function(x){
  summary(lm(x~kidneyage))$coef[2,c(1,4)]
}))

expbeta = cbind(cortbeta[,1], lungbeta[,1], liverbeta[,1], spleenbeta[,1], kidneybeta[,1])
colnames(expbeta) = c('Cortex', 'Lung', 'Liver', 'Spleen', 'Kidney')


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
# OR = 1.08, p = 0.1957 (old)
# OR = 1.673539, p = 2.2e-16

saveRDS(list(tbl = table(mat)[,c(2,1)],
             fisher = fisher.test(table(mat)[,c(2,1)])),
        file = 'data/other_datasets/jonker/processed/specloss_fisher.rds')

### using only Co genes:
# losing expr in native, gain in other tissues:
cogenes = covch %>% filter(rho<0) %>% pull(gene_id)

specsub = ts.specQ3.genes[ts.specQ3.genes%in%cogenes]
expchsub = ts.expr.ch[cogenes]
expchdirsub = ts.expr.ch.dir[cogenes]
matsub = data.frame(sameness = names(specsub) == expchsub[specsub],
                    expdir = expchdirsub[specsub])
table(matsub)[,c(2,1)]
sum(table(matsub)[,c(2,1)]) # 2967
fisher.test(table(matsub)[,c(2,1)])
fisher.test(table(matsub)[,c(2,1)])$p.val
# OR = 6.413, p = 2.2e-16 (old)
# OR = 7.518605, p = 6.464742e-109

saveRDS(list(tbl = table(matsub)[,c(2,1)],
             fisher = fisher.test(table(matsub)[,c(2,1)])),
        file = 'data/other_datasets/jonker/processed/specloss_fisher_co.rds')

covch %>% 
  filter(BH<0.1) %>%
  mutate(pattern = ifelse(rho<0, 'conv', 'div')) %>%
  group_by(pattern) %>%
  summarise(n=n())
# pattern     n
# <chr>   <int>
# 1 conv     1144
# 2 div       591
1144/(1144+591) # %66

#### Do DiCo genes enriched among tissue spec genes:
# devdiv = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
# divg = names(which(devdiv>0))
# idmap = readRDS('data/uniprot_to_ENS.rds')
# divg = idmap[idmap$ensembl_gene_id%in%divg,1]
# convg = covch %>% filter(rho < 0) %>% pull(gene_id)
# dicog = intersect(divg, convg)
# 
# spec.dc.mat = data.frame(spec = rownames(exp)%in%ts.specQ3.genes ,
#                          divconv = rownames(exp)%in%dicog )
# # spec 1: not ts specific, 2: ts specific
# # divconv 1: not dc gene, 2: dc gene
# table(spec.dc.mat)
# fisher.test(table(spec.dc.mat)) # 
# ## OR = 1.559732, p = 2.2e-16
# 
# dc_vs_tisspec = table(spec.dc.mat)[c(2,1),c(2,1)]
# colnames(dc_vs_tisspec) = c('DC', ' Non-DC')
# rownames(dc_vs_tisspec) = c('Specific to a tissue', 'Not tissue-specific')
# dc_vs_tisspec
# 
# saveRDS(list(table(spec.dc.mat),fisher.test(table(spec.dc.mat))),
#         file = 'data/other_datasets/jonker/processed/tisspec_dico_fisher.rds')


save(list=ls(), file='./data/other_datasets/jonker/analysis.rdata')

