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
# pc12 = pcdat %>%
#   ggplot(aes(x=PC1,y=PC2, size=Age, color=Tissue)) +
#   geom_point(alpha = 0.5) +
#   scale_color_manual(values = Tissuecol) +
#   scale_size_continuous(range = c(0.5,2), trans= 'log2') +
#   xlab(paste('PC1 (', round(varpc[1,2]*100),'%)',sep='')) +
#   ylab(paste('PC2 (', round(varpc[2,2]*100),'%)',sep='')) 
# 
# pc34 = pcdat %>% 
#   ggplot(aes(x=PC3,y=PC4, size=Age, color=Tissue)) +
#   geom_point(alpha = 0.5) +
#   scale_color_manual(values = Tissuecol) +
#   scale_size_continuous(range = c(0.5,2), trans= 'log2') +
#   #coord_fixed(ratio = varpc[4,2]/varpc[3,2], clip = 'off') +
#   xlab(paste('PC3 (', round(varpc[3,2]*100),'%)',sep='')) +
#   ylab(paste('PC4 (', round(varpc[4,2]*100),'%)',sep='')) 
# 
# pc1age = pcdat %>%
#   ggplot(aes(x = Age, y = PC1, color = Tissue)) +
#   geom_smooth(alpha = 0.1, se=F, show.legend = F, method='lm') +
#   geom_point(size = 0.5) +
#   facet_grid(Tissue~., scales='free_y') +
#   scale_color_manual(values = Tissuecol) +
#   guides(color = F) +
#   xlab('Age') +
#   ylab(NULL) +
#   stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
#            size = 6/pntnorm, show.legend = F) +
#   ggtitle('PC1') 
# 
# pc2age = pcdat %>%
#   ggplot(aes(x = Age, y = PC2, color = Tissue)) +
#   geom_smooth(alpha = 0.1, se=F, show.legend = F, method='lm') +
#   geom_point(size = 0.5) +
#   facet_grid(Tissue~., scales='free_y') +
#   scale_color_manual(values = Tissuecol) +
#   guides(color = F) +
#   xlab('Age') +
#   ylab(NULL) +
#   stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
#            size = 6/pntnorm, show.legend = F) +
#   ggtitle('PC2') 
# 
# pc3age = pcdat %>%
#   ggplot(aes(x = Age, y = PC3, color = Tissue)) +
#   geom_smooth(alpha = 0.1, se=F, show.legend = F, method='lm') +
#   geom_point(size = 0.5) +
#   facet_grid(Tissue~., scales='free_y') +
#   scale_color_manual(values = Tissuecol) +
#   guides(color = F) +
#   xlab('Age') +
#   ylab(NULL) +
#   stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
#            size = 6/pntnorm, show.legend = F) +
#   ggtitle('PC3') 
# 
# pc4age = pcdat %>%
#   ggplot(aes(x = Age, y = PC4, color = Tissue)) +
#   geom_smooth(alpha = 0.1, se=F, show.legend = F, method='lm') +
#   geom_point(size = 0.5) +
#   facet_grid(Tissue~., scales='free_y') +
#   scale_color_manual(values = Tissuecol) +
#   guides(color = F) +
#   xlab('Age') +
#   ylab(NULL) +
#   stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black', 
#            size = 6/pntnorm, show.legend = F) +
#   ggtitle('PC4') 
# 
# pc1234 = ggarrange(pc12, pc34, common.legend = T, legend = 'right', labels=c('a.','b.'), 
#                    font.label = list(size=8), widths = c(1,1))
# 
# pcages = ggarrange(pc1age,pc2age,pc3age,pc4age, nrow=1, ncol=4, labels=c('c.','d.','e.','f.'),
#                    font.label = list(size=8), legend='none')
# 
# pcaplots = ggarrange(pc1234, pcages, ncol=1, nrow=2, heights = c(1,1))

#ggsave('./results/other_datasets/jonker/pca.pdf', pcaplots, units='cm', width = 16, height = 12, useDingbats=F)
#ggsave('./results/other_datasets/jonker/pca.png', pcaplots, units='cm', width = 16, height = 12)

############# pairwise distance.....


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

pexpcors %>% 
  group_by(pair) %>%
  summarise(corrho = cor.test(rho,Age, m='s')$est,
            corp = cor.test(rho,Age, m='s')$p)
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

saveRDS(cov_dat_sum,'results/source_data/f2/fs7_mean_median.rds')
saveRDS(pexpcors,'results/source_data/f2/fs7_pexpcors.rds')

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

save(list=ls(), file='./data/other_datasets/jonker/analysis.rdata')

