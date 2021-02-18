library(tidyverse)
library(reshape2)
library(openxlsx)
library(ggpubr)
library(gridExtra)
library(grid)
library(cluster)
source('./scripts/functions.R')

tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))

cov = readRDS('./data/processed/tidy/CoV.rds')
covch = readRDS('./data/processed/tidy/CoV_change.rds')
expr = readRDS('./data/processed/tidy/expression.rds')
sample_info = readRDS('./data/processed/tidy/sample_info.rds')
dgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
dcgenes = names(which(dgenes<0))

### gap statistics to find optimum K:
dccov = cov %>% 
  right_join(data.frame(gene_id=dcgenes)) %>%
  group_by(gene_id) %>%
  mutate(scCoV = scale(CoV)[,1] )

inds = sample_info %>% 
  select(ind_id, age) %>% distinct() %>%
  filter(ind_id !='465') %>%
  mutate(ind_id = factor(ind_id))

###############
scmat = dccov %>%
  select(-CoV) %>%
  acast(gene_id~ind_id, value.var = 'scCoV')

gap = clusGap(scmat, FUNcluster = kmeans , K.max = 16, B = 100, d.power = 2 )

#print(scmat, method='firstSEmax', SE.factor=1)

######### kmeans 9 clusters
km9 = dccov %>%
  select(-CoV) %>%
  acast(gene_id~ind_id, value.var = 'scCoV') %>%
  kmeans(., centers = 9)

cl9 = table(km9$cluster)
cl9 = paste('Cl:',names(cl9), ' (n=',cl9,')', sep='')
names(cl9)= names(table(km9$cluster))

dccovkm9 = dccov %>%
  select(-CoV) %>%
  mutate(ind_id= factor(ind_id)) %>%
  left_join(inds) %>%
  left_join( data.frame( gene_id = names(km9$cluster), cluster = km9$cluster, row.names=NULL ) )

####### plot 9 clusters
dc_cluster9 = dccovkm9 %>%
  mutate(cluster = factor(cluster, levels=1:9)) %>%
  ggplot(aes(x= age, y= scCoV )) +
  facet_wrap(~cluster, ncol = 3,  scales = 'free_y', labeller=labeller(cluster=cl9)) +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept=90, linetype='dashed', alpha=0.3) +
  geom_line(stat='smooth', method='loess', aes(group=gene_id), alpha=0.05, col='midnightblue' ) +
  xlab('Age in days (in log2 scale)') +
  ylab('CoV scaled across genes') +
  theme(legend.position = 'none',
        axis.text = element_text(size=6),
        axis.title = element_text(size=8),
        strip.text = element_text(size=8, margin=margin(l=1,b=1,t=1)))

ggsave('./results/SI_figures/Figure_S24.pdf', dc_cluster9, units='cm', width = 10, height = 10, useDingbat=F)
ggsave('./results/SI_figures/Figure_S24.png', dc_cluster9, units='cm', width = 10, height = 10)

######### kmeans 12 clusters
km12 = dccov %>%
  select(-CoV) %>%
  acast(gene_id~ind_id, value.var = 'scCoV') %>%
  kmeans(., centers = 12)

cl12 = table(km12$cluster)
cl12 = paste('Cl:',names(cl12), ' (n=',cl12,')', sep='')
names(cl12)= names(table(km12$cluster))

dccovkm12 = dccov %>%
  select(-CoV) %>%
  mutate(ind_id= factor(ind_id)) %>%
  left_join(inds) %>%
  left_join( data.frame( gene_id = names(km12$cluster), cluster = km12$cluster, row.names=NULL ) )

####### plot 12 clusters
dc_cluster12 = dccovkm12 %>%
  mutate(cluster = factor(cluster, levels=1:12)) %>%
  ggplot(aes(x= age, y= scCoV )) +
  facet_wrap(~cluster, ncol = 4,  scales = 'free_y', labeller=labeller(cluster=cl12)) +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept=90, linetype='dashed', alpha=0.3) +
  geom_line(stat='smooth', method='loess', aes(group=gene_id), alpha=0.05, col='midnightblue' ) +
  xlab('Age in days (in log2 scale)') +
  ylab('CoV scaled across genes') +
  theme(legend.position = 'none',
        axis.text = element_text(size=6),
        axis.title = element_text(size=8),
        strip.text = element_text(size=8, margin=margin(l=1,b=1,t=1)))

ggsave('./results/SI_figures/Figure_S25.pdf', dc_cluster12, units='cm', width = 10, height = 10, useDingbat=F)
ggsave('./results/SI_figures/Figure_S25.png', dc_cluster12, units='cm', width = 10, height = 10)


######### kmeans 6 clusters
km6 = dccov %>%
  select(-CoV) %>%
  acast(gene_id~ind_id, value.var = 'scCoV') %>%
  kmeans(., centers = 6)

cl6 = table(km6$cluster)
cl6 = paste('Cl:',names(cl6), ' (n=',cl6,')', sep='')
names(cl6)= names(table(km6$cluster))

dccovkm6 = dccov %>%
  select(-CoV) %>%
  mutate(ind_id= factor(ind_id)) %>%
  left_join(inds) %>%
  left_join( data.frame( gene_id = names(km6$cluster), cluster = km6$cluster, row.names=NULL ) )

####### plot 6 clusters
dc_cluster6 = dccovkm6 %>%
  mutate(cluster = factor(cluster, levels=1:6)) %>%
  ggplot(aes(x= age, y= scCoV )) +
  facet_wrap(~cluster, ncol = 3,  scales = 'free_y', labeller=labeller(cluster=cl6)) +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept=90, linetype='dashed', alpha=0.3) +
  geom_line(stat='smooth', method='loess', aes(group=gene_id), alpha=0.05, col='midnightblue' ) +
  xlab('Age in days (in log2 scale)') +
  ylab('CoV scaled across genes') +
  theme(legend.position = 'none',
        axis.text = element_text(size=6),
        axis.title = element_text(size=8),
        strip.text = element_text(size=8, margin=margin(l=1,b=1,t=1)))

ggsave('./results/SI_figures/Figure_S26.pdf', dc_cluster6, units='cm', width = 10, height = 10, useDingbat=F)
ggsave('./results/SI_figures/Figure_S26.png', dc_cluster6, units='cm', width = 10, height = 10)

########## save clusters

km_clusters = list(km6 = km6, km9 = km9, km12 = km12)
saveRDS(km_clusters, './data/processed/raw/dc_km_clusters.rds')


#####################
##################### gora enrichment of dc cluster genes:
##################### 
#####################

############ gora of km9 :
bg  = setNames(rep(0,length(km9$cluster)), nm = names(km9$cluster))
dc_cluster9_gora = lapply(as.numeric(names(table(km9$cluster))), function(i){
  geneset = bg
  geneset[km9$cluster==i] = 1
  gora = go_bp_enrich.test.Mm(genelist = geneset, selection = 1)
  return(gora)
})
names(dc_cluster9_gora) = paste0('cl', names(cl9))

############ gora of km6 :
bg  = setNames(rep(0,length(km6$cluster)), nm = names(km6$cluster))
dc_cluster6_gora = lapply(as.numeric(names(table(km6$cluster))), function(i){
  geneset = bg
  geneset[km6$cluster==i] = 1
  gora = go_bp_enrich.test.Mm(genelist = geneset, selection = 1)
  return(gora)
})
names(dc_cluster6_gora) = paste0('cl', names(cl6))

############ gora of km12 :
bg  = setNames(rep(0,length(km12$cluster)), nm = names(km12$cluster))
dc_cluster12_gora = lapply(as.numeric(names(table(km12$cluster))), function(i){
  geneset = bg
  geneset[km12$cluster==i] = 1
  gora = go_bp_enrich.test.Mm(genelist = geneset, selection = 1)
  return(gora)
})
names(dc_cluster12_gora) = paste0('cl', names(cl12))

dc_clusters_gora = list(km6 = dc_cluster6_gora, km9 = dc_cluster9_gora, km12 = dc_cluster12_gora)

saveRDS(dc_clusters_gora, './data/processed/raw/dc_km_clusters_gora.rds')

dc_clusters_table = unlist(dc_clusters_gora, recursive = F)
write.xlsx(dc_clusters_table, './data/Table_S12.xlsx')

################
################  
################
################

##################
################## cluster expression levels of dc clusters (km9) gene:

sc_exp = expr %>%
  #filter(gene_id%in%names(which(km9$cluster==i)))  %>%
  left_join(sample_info, by='sample_id') %>%
  right_join(inds, by=c('ind_id','age')) %>% # remove one individual from other tissues that lacks expression in cortex
  group_by(gene_id) %>%
  mutate(scexp = scale(expression)[,1]) %>%
  dplyr::select(-log2age, -ind_id, -age, -tissue, -expression)

dccl9 = table(km9$cluster)
dccl9 = paste('DC_cl:',names(dccl9), ' (n=',dccl9,')', sep='')
names(dccl9)= names(table(km9$cluster))

exp_km_clusters_dc9 = list()

for(i in 1:length(cl9)){
  kms = sc_exp %>%
    filter(gene_id%in%names(which(km9$cluster==i)))  %>%
    acast(gene_id~sample_id, value.var = 'scexp') %>%
    kmeans(., center=12)
  exp_km_clusters_dc9[[i]] = kms
  names(exp_km_clusters_dc9)[i] = paste0('expcl',i)
  
  expcl = sc_exp %>%
    left_join(sample_info) %>%
    dplyr::select(-log2age) %>%
    right_join( data.frame( gene_id = names(kms$cluster), cluster = kms$cluster, row.names=NULL ) )
  
  clx = table(kms$cluster)
  clx = paste('Cl:',names(clx), ' (n=',clx,')', sep='')
  names(clx)= names(table(kms$cluster))
  
  dc_exp_x = expcl %>%
    ggplot(aes(x= age, y= scexp, col = tissue, group = interaction(tissue, gene_id) )) +
    facet_wrap(~cluster, ncol = 4,  scales = 'free_y', labeller = labeller(cluster=clx)) +
    #scale_y_continuous(limits = c(-4, 4)) +
    scale_x_continuous(trans = 'log2') +
    geom_vline(xintercept=90, linetype='dashed', alpha=0.4) +
    geom_line(stat='smooth', method = 'loess', alpha = 0.3, show.legend = T, size = 0.1 ) +
    scale_color_manual(values = tissuecol) +
    xlab('Age in days (in log2 scale)') +
    ylab('Scaled Expression') +
    guides(colour = guide_legend(override.aes = list(size=2, alpha=1 )) ) +
    labs(colour = 'Tissue') + 
    theme(legend.position = 'right')
  ###
  dc_cluster_x = dccovkm9 %>%
    filter(cluster==i) %>%
    mutate(cluster = factor(cluster)) %>%
    ggplot(aes(x= age, y= scCoV )) +
    facet_wrap(~cluster, labeller=labeller(cluster=dccl9[i])) +
    scale_x_continuous(trans = 'log2') +
    geom_vline(xintercept=90, linetype='dashed', alpha=0.4) +
    geom_line(stat='smooth', method='loess', aes(group=gene_id), alpha=0.1, col='midnightblue' ) +
    xlab('Age in days (in log2 scale)') +
    ylab('CoV scaled') +
    theme(axis.title.y = element_text(size=6))
  
  leg = get_legend(dc_exp_x)
  blank = grid.rect(gp=gpar(col='white'))
  px = ggarrange(ggarrange(blank, dc_cluster_x, leg, ncol=3, widths = c(0.4, 1, 0.6)),
            dc_exp_x, nrow = 2, heights = c(0.5,1), legend = 'none')
  px
  ###
  ix = i + 26
  ggsave(paste0('./results/SI_figures/Figure_S',ix,'.pdf'), px, units='cm', width = 10, height = 10, useDingbats=F)
  ggsave(paste0('./results/SI_figures/Figure_S',ix,'.png'), px, units='cm', width = 10, height = 10)  
}

##### save expression clusters of dc9 cluster genes:
saveRDS(exp_km_clusters_dc9, './data/processed/raw/dc_km9_exp_clusters.rds')

##############################
names(exp_km_clusters_dc9) = paste0('CoV_cl', 1:9)

dcexpcl_genes  = lapply(exp_km_clusters_dc9, function(x) {
  data.frame(cl = x$cluster, gene_id = names(x$cluster), row.names = NULL)}) %>%
  reshape2::melt() %>% 
  set_names(c('gene_id','var','Exp_cl','CoV_cl')) %>%
  dplyr::select(-var) %>%
  mutate(CoV_cl = as.numeric(factor(CoV_cl)) )

dcexpcl_genes = dcexpcl_genes %>%
  left_join(covch, by='gene_id')

write.xlsx(dcexpcl_genes, './data/Table_S11.xlsx')

################
################  gora enrichment of expression clusters among dc cluster 7 (km9) genes:
################
################


i=1
bg = setNames(rep(0,length(exp_km_clusters_dc9$CoV_cl1$cluster)),
              nm = names(exp_km_clusters_dc9$CoV_cl1$cluster))
exp_cl12_dcx_gora = lapply(as.numeric(names(table(exp_km_clusters_dc9$CoV_cl1$cluster))), function(i){
  geneset = bg
  geneset[exp_km_clusters_dc9$CoV_cl1$cluster==i] = 1
  gora = go_bp_enrich.test.Mm(genelist = geneset, selection = 1)
  return(gora)
})
names(exp_cl12_dcx_gora) = paste0('exp_cl', 1:12)

expcl_gora = lapply(exp_km_clusters_dc9, function(x){
  bg = setNames(rep(0,length(x$cluster)),
                nm = names(x$cluster))
  expx_cl12 = lapply(as.numeric(names(table(x$cluster))), function(i){
    geneset = bg
    geneset[x$cluster==i] = 1
    gora = go_bp_enrich.test.Mm(genelist = geneset, selection = 1)
    return(gora)
  })
  names(expx_cl12) = paste0('exp_cl', 1:12)
  return(expx_cl12)
})

saveRDS(expcl_gora, './data/processed/raw/dc_km9_exp_clusters_gora.rds')

################
################
################
################
################
################
