library(tidyverse)
library(ggpubr)
library(RColorBrewer)
exp = readRDS('./data/other_datasets/jonker/raw/expression.rds')
age = readRDS('./data/other_datasets/jonker/raw/age.rds')
tissue.ord = readRDS('./data/other_datasets/jonker/raw/tissue_id.rds')
sample_info = readRDS('./data/other_datasets/jonker/processed/sample_info.rds')
sample_info = sample_info %>% 
  mutate(Tissue= str_to_title(tissue)) %>% select(-tissue)

chstwo= brewer.pal(9,'Paired')[c(8,4)]
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

pc1234 = ggarrange(pc12, pc34, common.legend = T, legend = 'right', labels=c('a.','b.'), font.label = list(size=8),
          widths = c(1,1))

pcages = ggarrange(pc1age,pc2age,pc3age,pc4age, nrow=1, ncol=4, labels=c('c.','d.','e.','f.'),
                   font.label = list(size=8), legend='none')

pcaplots = ggarrange(pc1234, pcages, ncol=1, nrow=2, heights = c(1,1))
 
ggsave('./results/other_datasets/jonker/pca.pdf', pcaplots, units='cm', width = 16, height = 12, useDingbats=F)
ggsave('./results/other_datasets/jonker/pca.png', pcaplots, units='cm', width = 16, height = 12)

############# pairwise distance.....


###################################
###########################
#######################
####################
##################
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

covch = t(sapply(rownames(genecov),function(x){
  x = cor.test(genecov[x,], ageu, m ="s")
  c(x$est, x$p.val)
}))
covch = cbind(covch, p.adjust(covch[,2], method = "BH"))
colnames(covch)[c(2,3)] = c("pval"," BH")
saveRDS(covch,'./data/other_datasets/jonker/processed/covch.rds')
##### 1774 significant cov changes
##### 1215 convergent
##### 559 divergent
dcranks = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
length(intersect(names(dcranks), names(which(covch[,3]<0.1) )) ) # 937
length(intersect(names(which(dcranks<0)), names(which(covch[,3]<0.1) )) ) #492
length(intersect(names(which(dcranks<0)), names(which(covch[covch[,1]<0,3]<0.1) )) ) # 370
length(intersect(names(which(dcranks<0)), names(which(covch[covch[,1]>0,3]<0.1) )) ) # 122

#####
covch2 = data.frame(covch, gene_id = rownames(covch), row.names = NULL) %>%
  rename(CoV_change = rho, 
         FDR = X.BH)

cov_dat_sum = genecovtidy %>%
  left_join(unique(select(sample_info,-Tissue,-sample_id)), by='ind_id') %>%
  group_by(ind_id, Age) %>%
  summarise(meanCoV = mean(CoV),
            medianCoV = median(CoV)) %>% ungroup()

saveRDS(cov_dat_sum,file = './data/other_datasets/jonker/processed/mean_cov.rds')

cov_cordat = cov_dat_sum %>%
  summarise(cor = cor.test(meanCoV, Age, method = 's')$est,
            cor.p = cor.test(meanCoV, Age, method = 's')$p.val)

cov_cordat_median = cov_dat_sum %>%
  summarise(cor = cor.test(medianCoV, Age, method = 's')$est,
            cor.p = cor.test(medianCoV, Age, method = 's')$p.val)

meancov = cov_dat_sum %>%
  ggplot(aes(x = Age, y= meanCoV)) +
  geom_smooth(se=T, method='lm', color = 'midnightblue', fill = 'lightblue', show.legend = F) +
  geom_point(size = 1, color='steelblue') +
  scale_color_manual(values = Tissuecol) +
  xlab('Age') +
  ylab('Mean CoV') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F)

mediancov = cov_dat_sum %>%
  ggplot(aes(x = Age, y= medianCoV)) +
  geom_smooth(se=T, method='lm', color = 'midnightblue', fill = 'lightblue', show.legend = F) +
  geom_point(size = 1, color='steelblue') +
  scale_color_manual(values = Tissuecol) +
  xlab('Age') +
  ylab('Median CoV') +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', color = 'black',
           size = 6/pntnorm, show.legend = F)

covch2 %>% 
  top_n(n = 1, wt = -`CoV_change`) %>%
  left_join(reshape2::melt(exp) %>% set_names(c('gene_id','sample_id','expression'))) %>%
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
p2
p1  = ggarrange(meancov, mediancov, nrow=2 , labels=c('a.','b.'), font.label = list(size=8),vjust = c(1,-1) )

pwisecorplot = ggarrange(p1, p2, ncol=2, widths = c(2,3), labels= c(NA,'c.'),font.label = list(size=8))
ggsave('./results/SI_figures/Figure_S15.pdf', units='cm', width = 16, height = 11, useDingbats=F)
ggsave('./results/SI_figures/Figure_S15.png', units='cm', width = 16, height = 11)
#
save(list=ls(), file = './data/other_datasets/jonker/processed/analysis.rdata')
### or 
# a = ggarrange(pcaplots, pwisecorplot, nrow=2)
# ggsave('./results/SI_figures/Figure_S16.pdf',a, units='cm', width = 16, height = 16, useDingbats=F)

############### age related expression change

######## calculate CoV with excluding each tissue:
idmap = ids
names(idmap) = colnames(exp)
expr = exp
colnames(expr) = ids
ts = unique(tissue.ord)
covex = lapply(ts, function(i){
  print(paste('Calculating CoV without:',i))
  cov3tsX = sapply(unique(colnames(expr)), function(x){
    expX = expr[, !tissue.ord%in%i]
    sapply(rownames(expX), function(y){
      sd(expX[y, colnames(expX)==x]) / mean( expX[y, colnames(expX) == x] )
    })
  })
  cov3tsX
})

names(covex) = paste0('wo_', ts)

names(ageu) = unique(ids)
covextidy  = reshape2::melt(covex) %>%
  set_names(c('gene_id', 'ind_id', 'CoV', 'Excluded')) %>%
  mutate(ind_id = as.character(ind_id)) %>%
  left_join(data.frame(age=ageu, ind_id=names(ageu)))

saveRDS(covextidy, file='./data/other_datasets/jonker/processed/CoV_wo_eachtissue.rds')

covexch = covextidy %>%
  group_by(Excluded, gene_id) %>%
  summarise(rho  = cor.test(CoV, age, m='s')$est,
            pval = cor.test(CoV, age, m='s')$p.val)

covexch = ungroup(covexch) %>%
  mutate(BH = p.adjust(pval, method = 'BH'))

saveRDS(covexch, './data/other_datasets/jonker/processed/CoV_change_wo_eachtissue.rds')

cov_sum3ts = covextidy %>%
  mutate(Excluded = gsub('wo_','', Excluded)) %>%
  mutate(Excluded = paste0('Exclude: ', Excluded) ) %>%
  mutate(ind_id = factor(ind_id)) %>%
  #left_join(unique(select(sample_info,-tissue,-sample_id, -log2age))) %>%
  group_by(ind_id, age, Excluded) %>%
  summarise(meanCoV = mean(CoV),
            medianCoV = median(CoV))

cov_sum3ts

cov_sumch3ts = cov_sum3ts %>%
  group_by(Excluded) %>%
  summarise(mean_rho = cor.test(meanCoV, age, method = 's')$est,
            mean_p = cor.test(meanCoV, age, method = 's')$p.val,
            median_rho = cor.test(medianCoV, age, method = 's')$est,
            median_p = cor.test(medianCoV, age, method = 's')$p.val)


cov_exc_mean = cov_sum3ts %>%
  ggplot( aes(x = age, y = meanCoV)) +
  facet_wrap(~Excluded, ncol = 5) +
  geom_point(size=1.5, color="steelblue", alpha=0.9) +
  geom_smooth(method = 'lm', se=T,color = 'midnightblue', fill='lightblue') +
  scale_x_continuous(trans = 'log2') +
  xlab('Age in days (in log2 scale)') + 
  ylab('Mean CoV') +
  stat_cor(aes(x = age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) 

cov_exc_median = cov_sum3ts %>%
  ggplot( aes(x = age, y = medianCoV)) +
  facet_wrap(~Excluded, ncol = 5) +
  geom_point(size=1.5, color="steelblue", alpha=0.9) +
  geom_smooth(method = 'lm', se=T,color = 'midnightblue', fill='lightblue') +
  scale_x_continuous(trans = 'log2') +
  xlab('Age in days (in log2 scale)') + 
  ylab('Median CoV') +
  stat_cor(aes(x = age), method = 'spearman', cor.coef.name = 'rho', size = 6/pntnorm) 

jex = ggarrange(cov_exc_mean, cov_exc_median, nrow=2, labels=c('a.','b.'), align='hv',
                  vjust=c(1.1, -0.2))
jex
ggsave('./results/Jonker/covch_exc.pdf', jex, units = 'cm', width = 16, height = 12,
       useDingbats = F)
ggsave('./results/Jonker/covch_exc.png', jex, units = 'cm', width = 16, height = 12)


save(list=ls(), file='results/Jonker/data.rdata')
