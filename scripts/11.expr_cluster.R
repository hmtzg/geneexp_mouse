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

#cov = readRDS('./data/processed/tidy/CoV.rds')
#covch = readRDS('./data/processed/tidy/CoV_change.rds')
expr = readRDS('./data/processed/tidy/expression.rds')
sample_info = readRDS('./data/processed/tidy/sample_info.rds')
#dgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
#dcgenes = names(which(dgenes<0))
load('./results/expr_cluster_gap.rdata')
gapstat
sapply(gapstat, function(x){ maxSE(x$Tab[,'gap'], x$Tab[,'SE.sim'], method = 'Tibs2001SEmax') })
# Cortex  Liver   Lung Muscle 
# 15     14     17     17 
############### scale expression values:
scexpr = expr %>%
  left_join(select(sample_info,-log2age)) %>%
  group_by(gene_id, tissue) %>%
  mutate(scExpr = scale(expression)[,1]) %>% ungroup()

scexpT = split(scexpr, scexpr$tissue)

######### kmeans clusters:
cluT = sapply(gapstat, function(x){ maxSE(x$Tab[,'gap'], x$Tab[,'SE.sim'], method = 'Tibs2001SEmax') })
k = 0
kmT = lapply(scexpT, function(x) {
  k <<- k+1
  print(paste(k, names(cluT[k]),':',cluT[k]))
  x %>%
    acast(gene_id~sample_id, value.var = 'scExpr') %>%
    kmeans(., centers = cluT[k], iter.max = 20, nstart = 50)
})
clT = lapply(kmT, function(x) table(x$cluster))
clT = lapply(clT, function(x) paste('Cl:', names(x), ' (n=', x, ')', sep=''))
for(i in 1:4) names(clT[[i]]) = 1:cluT[i]

scexpT = sapply(names(scexpT), function(x){
  scexpT[[x]] %>% 
    mutate(sample_id = factor(sample_id)) %>%
    left_join( data.frame(gene_id = names(kmT[[x]]$cluster), cluster = kmT[[x]]$cluster, row.names=NULL) )
}, simplify = F, USE.NAMES=T)

########## save clusters
saveRDS(kmT, './data/processed/raw/expression_clusters_each_tissue.rds')

dgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
dico = names(which(dgenes<0))

scexpT = lapply(scexpT, function(x)
  x %>% mutate(pattern = ifelse(gene_id%in%dico,'DiCo', 'non-DiCo')) )

enrichmentf = function(dd, clcol = 8, genecol=1, patcol=9){
  colnames(dd)[c(clcol, genecol, patcol)] = c('cluster', 'gene_id', 'pattern')
  # clcol =  column number indicating which cluster the gene belongs to
  # genecol = column number of the gene ids
  # patcol = pattern column to compare categories in the column (2 groups, fisher)
  
  sapply(unique(dd$cluster), function(x){
    gr1 = dd %>%
      filter(cluster==x) %>%
      group_by(pattern) %>%
      summarise(n=length(unique(gene_id)))
    gr2 = dd %>%
      filter(!cluster==x) %>%
      group_by(pattern) %>%
      summarise(n=length(unique(gene_id)))
    
    fisherX = gr1 %>% full_join(gr2, by='pattern') %>%
      rename( 'clX'= n.x, 'clOthers' = n.y) %>%
      arrange(pattern) %>% 
      column_to_rownames('pattern') %>% 
      t() 
    if(sum(is.na(fisherX))>0) fisherX = fisherX %>% replace_na(replace = 0)
    fisherXtest = fisher.test(fisherX)
    list(table= fisherX, Fistest = fisherXtest)
  }, USE.NAMES=T, simplify=F)
}

dicoenrich = lapply(scexpT, function(x) enrichmentf(x, clcol = 8, genecol = 1, patcol=9) )
enrichresult = lapply(dicoenrich, function(x) {
  sapply(x, function(y) c('OR' = unname(y$Fistest$est), 'p.val' = y$Fistest$p.val))
})
enrichresult = lapply(enrichresult, function(x){
  rbind(x, 'BH' = p.adjust(x[2,], method='BH'))
})
enrichresult = lapply(enrichresult, function(x) {colnames(x) = 1:ncol(x); return(x) })

# significant results:
lapply(enrichresult, function(x) x[,x['BH',] < 0.1]  )
# DiCo depleted groups:
lapply(enrichresult, function(x) x[,x['BH',] < 0.1  & x['OR',] < 1 ]  )
# DiCo enriched groups:
lapply(enrichresult, function(x) x[,x['BH',] < 0.1  & x['OR',] > 1 ]  )


##################### 
#####################

##### Plot clusters:
clT = lapply(kmT, function(x) table(x$cluster))
clT = lapply(clT, function(x) paste('Cl:', names(x), ' (n=', x, ')', sep=''))
for(i in 1:4) names(clT[[i]]) = 1:cluT[i]

# for(i in names(enrichresult)){
#   dep = which( enrichresult[[i]]['BH',]<0.1 & enrichresult[[i]]['OR',] < 1 )
#   enr = which( enrichresult[[i]]['BH',]<0.1 & enrichresult[[i]]['OR',] > 1 )
#   clT[[i]][names(clT[[i]])%in%dep] = paste(clT[[i]][names(clT[[i]])%in%dep], '-')
#   clT[[i]][names(clT[[i]])%in%enr] = paste(clT[[i]][names(clT[[i]])%in%enr], '+')
# }

enrcol = c('+1' = 'firebrick', '0' = 'gray50' , '-1' = 'midnightblue')

xxct = scexpT$Cortex
enrct = which(enrichresult$Cortex['OR',]>1 & enrichresult$Cortex['BH',] < 0.1)
depct = which(enrichresult$Cortex['OR',]<1 & enrichresult$Cortex['BH',] < 0.1)

xxct = xxct %>%
  mutate(enr = ifelse(cluster%in%enrct, '+1', '0'))
xxct$enr[xxct$cluster%in%depct] = '-1'

# cortex plot : 15 clusters
expcl_Cortex = xxct %>%
  mutate(cluster = factor(cluster, levels=1:length(unique(cluster)))) %>%
  ggplot(aes(x= age, y= scExpr, color=enr )) +
  facet_wrap(~cluster, ncol = 4,  scales = 'free_y', labeller=labeller(cluster=clT$Cortex)) +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept=90, linetype='dashed', alpha=0.3) +
  geom_line(stat='smooth', method='loess', aes(group=gene_id), alpha=0.05) +
  scale_color_manual(values = enrcol) +
  xlab('Age in days (in log2 scale)') +
  ylab('Expression scaled across genes') +
  theme(legend.position = 'none',
        axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        strip.text = element_text(size=6, margin=margin(l=1,b=1,t=1)))

# expcl_Cortex = scexpT$Cortex %>%
#   mutate(cluster = factor(cluster, levels=1:length(unique(cluster)))) %>%
#   ggplot(aes(x= age, y= scExpr )) +
#   facet_wrap(~cluster, ncol = 4,  scales = 'free_y', labeller=labeller(cluster=clT$Cortex)) +
#   scale_x_continuous(trans = 'log2') +
#   geom_vline(xintercept=90, linetype='dashed', alpha=0.3) +
#   geom_line(stat='smooth', method='loess', aes(group=gene_id), alpha=0.05, col='midnightblue' ) +
#   xlab('Age in days (in log2 scale)') +
#   ylab('Expression scaled across genes') +
#   theme(legend.position = 'none',
#         axis.text = element_text(size=6),
#         axis.title = element_text(size=6),
#         strip.text = element_text(size=6, margin=margin(l=1,b=1,t=1)))
#expcl_Cortex
ggsave('./results/figure_supplements/f1s/FS12.pdf', expcl_Cortex, units='cm', width = 10, height = 10,
       useDingbat=F)
ggsave('./results/figure_supplements/f1s/FS12.png', expcl_Cortex, units='cm', width = 10, height = 10)

xxln = scexpT$Lung
enrln = which(enrichresult$Lung['OR',]>1 & enrichresult$Lung['BH',] < 0.1)
depln = which(enrichresult$Lung['OR',]<1 & enrichresult$Lung['BH',] < 0.1)

xxln = xxln %>%
  mutate(enr = ifelse(cluster%in%enrln, '+1', '0'))
xxln$enr[xxln$cluster%in%depln] = '-1'
# lung plot : 17 clusters
expcl_Lung = xxln %>%
  mutate(cluster = factor(cluster, levels=1:length(unique(cluster)))) %>%
  ggplot(aes(x= age, y= scExpr, color=enr )) +
  facet_wrap(~cluster, ncol = 4,  scales = 'free_y', labeller=labeller(cluster=clT$Lung)) +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept=90, linetype='dashed', alpha=0.3) +
  geom_line(stat='smooth', method='loess', aes(group=gene_id), alpha=0.05) +
  scale_color_manual(values = enrcol) +
  xlab('Age in days (in log2 scale)') +
  ylab('Expression scaled across genes') +
  theme(legend.position = 'none',
        axis.text = element_text(size=6),
        axis.title = element_text(size=8),
        strip.text = element_text(size=8, margin=margin(l=1,b=1,t=1)))
#expcl_Lung
ggsave('./results/figure_supplements/f1s/FS13.pdf', expcl_Lung, units='cm', width = 10, height = 10,
       useDingbat=F)
ggsave('./results/figure_supplements/f1s//FS13.png', expcl_Lung, units='cm', width = 10, height = 10)

xxlv = scexpT$Liver
enrlv = which(enrichresult$Liver['OR',]>1 & enrichresult$Liver['BH',] < 0.1)
deplv = which(enrichresult$Liver['OR',]<1 & enrichresult$Liver['BH',] < 0.1)

xxlv = xxlv %>%
  mutate(enr = ifelse(cluster%in%enrlv, '+1', '0'))
xxlv$enr[xxlv$cluster%in%deplv] = '-1'
# liver plot : 14 clusters
expcl_Liver = xxlv %>%
  mutate(cluster = factor(cluster, levels=1:length(unique(cluster)))) %>%
  ggplot(aes(x= age, y= scExpr, color=enr )) +
  facet_wrap(~cluster, ncol = 4,  scales = 'free_y', labeller=labeller(cluster=clT$Liver)) +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept=90, linetype='dashed', alpha=0.3) +
  geom_line(stat='smooth', method='loess', aes(group=gene_id), alpha=0.05) +
  scale_color_manual(values = enrcol) +
  xlab('Age in days (in log2 scale)') +
  ylab('Expression scaled across genes') +
  theme(legend.position = 'none',
        axis.text = element_text(size=6),
        axis.title = element_text(size=8),
        strip.text = element_text(size=8, margin=margin(l=1,b=1,t=1)))
#expcl_Liver
ggsave('./results/figure_supplements/f1s/FS14.pdf', expcl_Liver, units='cm', width = 10, height = 10,
       useDingbat=F)
ggsave('./results/figure_supplements/f1s/FS14.png', expcl_Liver, units='cm', width = 10, height = 10)


xxms = scexpT$Muscle
enrms = which(enrichresult$Muscle['OR',]>1 & enrichresult$Muscle['BH',] < 0.1)
depms = which(enrichresult$Muscle['OR',]<1 & enrichresult$Muscle['BH',] < 0.1)

xxms = xxms %>%
  mutate(enr = ifelse(cluster%in%enrms, '+1', '0'))
xxms$enr[xxms$cluster%in%depms] = '-1'
# muscle plot : 17 clusters
expcl_Muscle = xxms %>%
  mutate(cluster = factor(cluster, levels=1:length(unique(cluster)))) %>%
  ggplot(aes(x= age, y= scExpr, color=enr )) +
  facet_wrap(~cluster, ncol = 4,  scales = 'free_y', labeller=labeller(cluster=clT$Muscle)) +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept=90, linetype='dashed', alpha=0.3) +
  geom_line(stat='smooth', method='loess', aes(group=gene_id), alpha=0.05) +
  scale_color_manual(values = enrcol) +
  xlab('Age in days (in log2 scale)') +
  ylab('Expression scaled across genes') +
  theme(legend.position = 'none',
        axis.text = element_text(size=6),
        axis.title = element_text(size=8),
        strip.text = element_text(size=8, margin=margin(l=1,b=1,t=1)))
#expcl_Muscle
ggsave('./results/figure_supplements/f1s/FS15.pdf', expcl_Muscle, units='cm', width = 10, height = 10,
       useDingbat=F)
ggsave('./results/figure_supplements/f1s/FS15.png', expcl_Muscle, units='cm', width = 10, height = 10)

saveRDS(xxct,'results/source_data/f1/fs12.rds')
saveRDS(xxln,'results/source_data/f1/fs13.rds')
saveRDS(xxlv,'results/source_data/f1/fs14.rds')
saveRDS(xxms,'results/source_data/f1/fs15.rds')

########### gora of expression clusters :
bg  = setNames(rep(0,length(kmT$Cortex$cluster)), nm = names(kmT$Cortex$cluster))
exp_cl_gora = lapply(as.numeric(names(table(kmT$Cortex$cluster))), function(i){
  print(i)
  geneset = bg
  geneset[kmT$Cortex$cluster==i] = 1
  gora = go_bp_enrich.test.Mm(genelist = geneset, selection = 1)
  return(gora)
})

dico_cl_gora = lapply(as.numeric(names(table(kmT$Cortex$cluster))), function(i){
    print(i)
    clg = names(which(kmT$Cortex$cluster==i))
    geneset = setNames(rep(0, length(clg)), nm = clg)
    geneset[names(geneset)%in%dico] = 1
    #geneset[kmT$Cortex$cluster==i] = 1
    gora = go_bp_enrich.test.Mm(genelist = geneset, selection = 1)
    return(gora) 
})
dico_cl_gora
lapply(dico_cl_gora, head)
save(list=ls(), file = './results/SI_figures/expr_cluster/analysis.rdata')
################
################
################
################
################
################

# example plot for tik
# expcplot = expr %>% left_join(sample_info) %>%
#   filter(tissue%in%c('Cortex','Lung') & ind_id%in%'841') %>%
#   dplyr::select(-log2age, -age, -ind_id, -sample_id) %>% 
#   spread(key = 'tissue', value = expression) %>%
#   ggplot(aes(x=Cortex, y=Lung)) +
#   geom_point(alpha = 0.4, color='gray30') +
#   geom_smooth(method='lm') +
#   stat_cor(method='spearman', cor.coef.name = 'rho')
# ggsave('../../PhD/Tik_1/expr_cor.pdf', expcplot, units='cm', width = 10, height = 10, useDingbats=F)
# ggsave('../../PhD/Tik_1/expr_cor.png', expcplot, units='cm', width = 10, height = 10)
# 
# expcplot2 = expr %>% left_join(sample_info) %>%
#   filter(tissue%in%c('Cortex','Lung') & ind_id%in%'2067') %>%
#   dplyr::select(-log2age, -age, -ind_id, -sample_id) %>% 
#   spread(key = 'tissue', value = expression) %>%
#   ggplot(aes(x=Cortex, y=Lung)) +
#   geom_point(alpha = 0.4, color='gray30') +
#   geom_smooth(method='lm') +
#   stat_cor(method='spearman', cor.coef.name = 'rho')
# ggsave('../../PhD/Tik_1/expr_cor2.pdf', expcplot2, units='cm', width = 10, height = 10, useDingbats=F)
# ggsave('../../PhD/Tik_1/expr_cor2.png', expcplot2, units='cm', width = 10, height = 10)
