library(tidyverse)
library(reshape2)
library(openxlsx)
library(ggpubr)
library(gridExtra)
library(grid)
library(cluster)
library(factoextra)
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

############### get dc genes and scale
dccov = cov %>% 
  right_join(data.frame(gene_id=dcgenes)) %>%
  group_by(gene_id) %>%
  mutate(scCoV = scale(CoV)[,1] )

inds = sample_info %>% 
  select(ind_id, age) %>% distinct() %>%
  filter(ind_id !='465') %>%
  mutate(ind_id = factor(ind_id))

############### get dc gene matrix
scmat = dccov %>%
  select(-CoV) %>%
  acast(gene_id~ind_id, value.var = 'scCoV')


### gap statistics to find optimum K for kmeans:
# gap1 = clusGap(scmat, FUNcluster = kmeans , K.max = 15, B = 500, d.power = 2, spaceH0 = 'original' )
# gap2 = clusGap(scmat, FUNcluster = kmeans , K.max = 15, B = 500, d.power = 2, spaceH0 = 'original' )
gap3 = clusGap(scmat, FUNcluster = kmeans , K.max = 15, B = 500, d.power = 2, spaceH0 = 'original' )
# gap4 = clusGap(scmat, FUNcluster = kmeans , K.max = 15, B = 500, d.power = 2, spaceH0 = 'original' )
# gap5 = clusGap(scmat, FUNcluster = kmeans , K.max = 15, B = 500, d.power = 2, spaceH0 = 'original' )

# maxSE(gap1$Tab[,'gap'], gap1$Tab[,'SE.sim'], method = 'Tibs2001SEmax')
# maxSE(gap2$Tab[,'gap'], gap2$Tab[,'SE.sim'], method = 'Tibs2001SEmax')
maxSE(gap3$Tab[,'gap'], gap3$Tab[,'SE.sim'], method = 'Tibs2001SEmax')
# maxSE(gap4$Tab[,'gap'], gap4$Tab[,'SE.sim'], method = 'Tibs2001SEmax')
# maxSE(gap5$Tab[,'gap'], gap5$Tab[,'SE.sim'], method = 'Tibs2001SEmax')
factoextra::fviz_gap_stat(gap3)

######### kmeans 7 clusters
kmX = dccov %>%
  select(-CoV) %>%
  acast(gene_id~ind_id, value.var = 'scCoV') %>%
  kmeans(., centers = 7, iter.max = 20, nstart = 50)

clX = table(kmX$cluster)
clX = paste('Cl:',names(clX), ' (n=',clX,')', sep='')
names(clX)= names(table(kmX$cluster))

dccovkmX = dccov %>%
  select(-CoV) %>%
  mutate(ind_id= factor(ind_id)) %>%
  left_join(inds) %>%
  left_join( data.frame( gene_id = names(kmX$cluster), cluster = kmX$cluster, row.names=NULL ) )

####### plot 7 clusters
dc_clusterX = dccovkmX %>%
  mutate(cluster = factor(cluster, levels=1:length(clX))) %>%
  ggplot(aes(x= age, y= scCoV )) +
  facet_wrap(~cluster, ncol = 3,  scales = 'free_y', labeller=labeller(cluster=clX)) +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept=90, linetype='dashed', alpha=0.3) +
  geom_line(stat='smooth', method='loess', aes(group=gene_id), alpha=0.05, col='midnightblue' ) +
  xlab('Age in days (in log2 scale)') +
  ylab('CoV scaled across genes') +
  theme(legend.position = 'none',
        axis.text = element_text(size=6),
        axis.title = element_text(size=8),
        strip.text = element_text(size=8, margin=margin(l=1,b=1,t=1)))

ggsave('./results/SI_figures/Figure_S24.pdf', dc_clusterX, units='cm', width = 10, height = 10, useDingbat=F)
ggsave('./results/SI_figures/Figure_S24.png', dc_clusterX, units='cm', width = 10, height = 10)

########## save clusters
saveRDS(kmX, './data/processed/raw/dc_km_clusters.rds')

#####################
##################### gora enrichment of dc cluster genes:
##################### 
#####################

############ gora of kmX :
bg  = setNames(rep(0,length(kmX$cluster)), nm = names(kmX$cluster))
dc_clusterX_gora = lapply(as.numeric(names(table(kmX$cluster))), function(i){
  geneset = bg
  geneset[kmX$cluster==i] = 1
  gora = go_bp_enrich.test.Mm(genelist = geneset, selection = 1)
  return(gora)
})
names(dc_clusterX_gora) = paste0('cl', names(clX))

saveRDS(dc_clusterX_gora, './data/processed/raw/dc_km_clusters_gora.rds')

write.xlsx(dc_clusterX_gora, './results/SI_tables/TableS7.xlsx')

################
################  
################
################

##################
################## cluster expression levels of dc genes:

sc_exp = expr %>%
  filter(gene_id%in%dcgenes )  %>%
  left_join(sample_info, by='sample_id') %>%
  right_join(inds, by=c('ind_id','age')) %>% # remove one individual from other tissues that lacks expression in cortex
  group_by(gene_id) %>%
  mutate(scexp = scale(expression)[,1]) %>%
  dplyr::select(-log2age, -ind_id, -age, -tissue, -expression)

############### get expression gene matrix
expmat = sc_exp %>%
  acast(gene_id~sample_id, value.var = 'scexp')

gapexp = clusGap(expmat, FUNcluster = kmeans, K.max = 15, B = 500, d.power = 2, spaceH0 = 'original' )
gapexp2 = clusGap(expmat, FUNcluster = kmeans, K.max = 20, B = 500, d.power = 2, spaceH0 = 'original' )
gapexp3 = clusGap(expmat, FUNcluster = kmeans, K.max = 30, B = 100, d.power = 2, spaceH0 = 'original' )

maxSE(gapexp$Tab[,'gap'], gapexp$Tab[,'SE.sim'], method = 'Tibs2001SEmax')
maxSE(gapexp2$Tab[,'gap'], gapexp2$Tab[,'SE.sim'], method = 'Tibs2001SEmax')
maxSE(gapexp3$Tab[,'gap'], gapexp3$Tab[,'SE.sim'], method = 'Tibs2001SEmax')
factoextra::fviz_gap_stat(gapexp)
factoextra::fviz_gap_stat(gapexp2)
factoextra::fviz_gap_stat(gapexp3)

kmexp = kmeans(expmat, centers = 25, iter.max = 20, nstart = 50)

saveRDS(kmexp, './data/processed/raw/dc_km_exp_clusters.rds')

cle = table(kmexp$cluster)
cle = paste('Cl:',names(cle), ' (n=',cle,')', sep='')
names(cle)= names(table(kmexp$cluster))

expcl = sc_exp %>%
  left_join(sample_info) %>%
  dplyr::select(-log2age) %>%
  right_join( data.frame( gene_id = names(kmexp$cluster), cluster = kmexp$cluster, row.names=NULL ) )

exp_clusters = expcl %>%
    ggplot(aes(x= age, y= scexp, col = tissue, group = interaction(tissue, gene_id) )) +
    facet_wrap(~cluster, ncol = 4,  scales = 'free_y', labeller = labeller(cluster=cle)) +
    scale_x_continuous(trans = 'log2') +
    geom_vline(xintercept=90, linetype='dashed', alpha=0.4) +
    geom_line(stat='smooth', method = 'loess', alpha = 0.1, show.legend = T, size = 0.1 ) +
    scale_color_manual(values = tissuecol) +
    xlab('Age in days (in log2 scale)') +
    ylab('Scaled Expression') +
    guides(colour = guide_legend(override.aes = list(size=2, alpha=1 )) ) +
    labs(colour = 'Tissue') +
    theme(legend.position = 'bottom')

ggsave('./results/SI_figures/Figure_S25.pdf', exp_clusters, units='cm', width = 15, height = 18, useDingbats=F)
ggsave('./results/SI_figures/Figure_S25.png', exp_clusters, units='cm', width = 15, height = 18)


##### save expression clusters of dc genes:

# dcexpcl_genes  = reshape2::melt(kmexp$cluster) %>% 
#   rownames_to_column(var = 'gene_id') %>%
#   set_names(c('gene_id','exp_Cl')) %>% 
#   left_join(data.frame(gene_id = names(kmX$cluster), DC_Cl = kmX$cluster, row.names=NULL ) ) %>% 
#   left_join(covch)
# 
# write.xlsx(dcexpcl_genes, './data/Table_S11.xlsx')

########### gora of expression clusters :
bg  = setNames(rep(0,length(kmexp$cluster)), nm = names(kmexp$cluster))
dc_exp_clX_gora = lapply(as.numeric(names(table(kmexp$cluster))), function(i){
  geneset = bg
  geneset[kmexp$cluster==i] = 1
  gora = go_bp_enrich.test.Mm(genelist = geneset, selection = 1)
  return(gora)
})
names(dc_exp_clX_gora) = paste0('cl', names(cle))

saveRDS(dc_exp_clX_gora, './data/processed/raw/dc_exp_km_clusters_gora.rds')

write.xlsx(dc_exp_clX_gora, './results/SI_tables/TableS8.xlsx')

######### save DC and exp cluster numbers of DC genes
# kmX_dat = reshape2::melt(kmX$cluster) %>% 
#   rownames_to_column(var = 'gene_id') %>%
#   set_names(c('gene_id','CoV_Cl'))
# 
# DiCo_cl_genes = covch %>% 
#   right_join(kmX_dat)
# 
# write.xlsx(DiCo_cl_genes, './data/Table_S11.xlsx')

DiCo_cl_genes = covch %>%
  right_join(
    reshape2::melt(kmX$cluster) %>% 
      rownames_to_column() %>% 
      set_names(c('gene_id','CoV_Cl')) %>%
      left_join(
        reshape2::melt(kmexp$cluster) %>%
          rownames_to_column() %>% 
          set_names(c('gene_id','Exp_Cl'))
      )
  ) %>%
  relocate(gene_id, period, CoV_change, pval, FDR, CoV_Cl, Exp_Cl)

saveRDS(DiCo_cl_genes, './data/processed/tidy/DiCo_cl_genes.rds')
write.xlsx(DiCo_cl_genes, './results/SI_tables/TableS6.xlsx')

##############################

save(list=ls(), file = './data/processed/raw/dc_cluster_alldata.rdata')
################
################
################
################
################
################
