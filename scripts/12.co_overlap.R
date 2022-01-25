library(tidyverse)
library(ggpubr)
library(ggforce)
library(RColorBrewer)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)

Our_dat = readRDS('./data/processed/tidy/CoV_change.rds') %>% 
  tibble() %>% 
  filter(period=='aging') %>%
  relocate(gene_id) %>%
  mutate(pattern = ifelse(CoV_change<0, 'convergent', 'divergent')) %>%
  dplyr::select(-period, -CoV_change, -pval, -FDR)

Jonker = readRDS('./data/other_datasets/jonker/processed/covch.rds') %>% tibble() %>%
  dplyr::rename(FDR='BH', CoV_change='rho') %>%
  mutate(pattern = ifelse(CoV_change<0, 'convergent', 'divergent')) %>%
  dplyr::select(-CoV_change, -pval, -FDR)

ensmap = readRDS('./results/GTEx/ensmap.rds')
GTEx4 = readRDS('./results/GTEx/covcor.rds') %>%
  set_names(c('gene_id','CoV_change', 'pval', 'FDR')) %>%
  inner_join(ensmap) %>%
  dplyr::select(-gene_id) %>%
  rename(gene_id='ENS_mm') %>% relocate(gene_id) %>%
  mutate(pattern = ifelse(CoV_change<0, 'convergent', 'divergent')) %>%
  dplyr::select(-CoV_change, -pval, -FDR)

GTEx10 = readRDS('./results/GTEx/alltissues/covch.rds') %>%
  set_names(c('gene_id','CoV_change', 'pval', 'FDR')) %>%
  inner_join(ensmap) %>%
  dplyr::select(-gene_id) %>%
  rename(gene_id='ENS_mm') %>% relocate(gene_id) %>%
  mutate(pattern = ifelse(CoV_change<0, 'convergent', 'divergent')) %>%
  dplyr::select(-CoV_change, -pval, -FDR)

idconv = readRDS('data/other_datasets/schaum/4tissue/idconvert.rds')
Schaum4 = readRDS('data/other_datasets/schaum/4tissue/covch.rds') %>%
  set_names(c('gene_id','CoV_change', 'pval', 'FDR')) %>%
  inner_join(idconv) %>%
  dplyr::select(-gene_id) %>%
  rename(gene_id='ensembl_gene_id') %>% relocate(gene_id) %>%
  mutate(pattern = ifelse(CoV_change<0, 'convergent', 'divergent')) %>%
  dplyr::select(-CoV_change, -pval, -FDR)

Schaum8 = readRDS('data/other_datasets/schaum/8tis/covch.rds') %>%
  set_names(c('gene_id','CoV_change', 'pval', 'FDR')) %>%
  inner_join(idconv) %>%
  dplyr::select(-gene_id) %>%
  rename(gene_id='ensembl_gene_id') %>% relocate(gene_id) %>%
  mutate(pattern = ifelse(CoV_change<0, 'convergent', 'divergent')) %>%
  dplyr::select(-CoV_change, -pval, -FDR)

length(Reduce(intersect, list(Our_dat$gene_id, Jonker$gene_id, GTEx4$gene_id, Schaum4$gene_id)))

alldat = reshape2::melt(list(Our_dat = Our_dat, Jonker = Jonker, GTEx4 = GTEx4, GTEx10 = GTEx10,
     Schaum4 = Schaum4, Schaum8 = Schaum8)) %>%
  rename(dataset='L1') %>% tibble()
alldat

datasets = unique(alldat$dataset)
pairs = t(combn(datasets,2))

pwisetests = lapply(1:nrow(pairs), function(x){
  alldat %>% filter(dataset%in%pairs[x,]) %>%
    spread(key='dataset', value='pattern') %>%
    na.omit() %>%
    group_by_at(pairs[x,] ) %>%
    summarise(n=length(unique(gene_id))) %>% ungroup() %>%
    reshape2::acast(.[[1]]~.[[2]], value.var = 'n') %>%
    fisher.test()
})
##### check:
xx = Our_dat %>%
  inner_join(Jonker, by='gene_id')
coco = nrow(xx %>% filter(pattern.x=='convergent' & pattern.y=='convergent'))
codi = nrow(xx %>% filter(pattern.x=='convergent' & pattern.y=='divergent'))
dico = nrow(xx %>% filter(pattern.x=='divergent' & pattern.y=='convergent'))
didi = nrow(xx %>% filter(pattern.x=='divergent' & pattern.y=='divergent'))
fisher.test(matrix(c(coco, codi, dico, didi), ncol=2, byrow = T))

#####
names(pwisetests) = paste(pairs[,1], pairs[,2],sep = '-')
pOR = reshape2::melt(lapply(pwisetests, function(x) data.frame(x$est, x$p.val)) ) %>%
  spread(key='variable', value='value') %>%
  set_names(c('datasets', 'OR', 'pval'))

pORtest = pOR %>%
  mutate(dataset1 = unlist(lapply(strsplit(datasets,split = '-'), `[[`, 1)),
         dataset2 = unlist(lapply(strsplit(datasets,split = '-'), `[[`, 2))) %>%
  dplyr::select(-datasets)

pORtest
##### convert to matrix similar to correlation matrix:
pOR = pORtest %>%
  reshape2::acast(dataset1~dataset2,value.var = 'OR')
pOR = cbind(pOR[,1:3], Our_dat = c(pOR[4,1:3],NA, pOR[4,4]), pOR[,4:5])
pOR = rbind(pOR, Schaum8 = c(pOR[,6], NA))
pOR[1:3,1:3][upper.tri(pOR[1:3,1:3])] = pOR[1:3,1:3][lower.tri(pOR[1:3,1:3])]
pOR[5,1:3] =  pOR[1:3,5]
pOR
#           GTEx10    GTEx4    Jonker   Our_dat  Schaum4   Schaum8
# GTEx10         NA 1.637976 1.0847708 0.9703041 1.072130 1.0015795
# GTEx4   1.6379756       NA 1.1036949 1.0589178 1.106714 1.0673532
# Jonker  1.0847708 1.103695        NA 1.2210456 1.243193 0.9572973
# Our_dat 0.9703041 1.058918 1.2210456        NA 1.883678 1.3245659
# Schaum4 1.0721303 1.106714 1.2431934 1.8836783       NA 5.7638882
# Schaum8 1.0015795 1.067353 0.9572973 1.3245659 5.763888        NA

pORtest %>%
  mutate(fdr = p.adjust(pval, method = 'BH')) %>%
  filter(fdr<0.1)
  
pORpval = pORtest %>%
  mutate(fdr = p.adjust(pval, method = 'BH')) %>%
  reshape2::acast(dataset1~dataset2,value.var = 'fdr')
pORpval = cbind(pORpval[,1:3], Our_dat = c(pORpval[4,1:3],NA, pORpval[4,4]), pORpval[,4:5])
pORpval = rbind(pORpval, Schaum8 = c(pORpval[,6], NA))
pORpval[1:3,1:3][upper.tri(pORpval[1:3,1:3])] = pORpval[1:3,1:3][lower.tri(pORpval[1:3,1:3])]
pORpval[5,1:3] =  pORpval[1:3,5]
pORpval

simcol = setNames(c('dodgerblue','#E4572E'), c('depleted','enriched'))

pORpval[lower.tri(pORpval)] = NA
pstar = reshape2::melt(pORpval) %>%
  mutate(star = ifelse(value<0.1 & value > 0.01, '*', '')) %>%
  mutate(star = ifelse(value<0.01 & value > 0.001, '**', star)) %>%
  mutate(star = ifelse(value<0.001, '***', star))

cooverlap = reshape2::melt(pOR) %>% 
  ggplot(aes(x=Var1, y=Var2, fill=log2(value))) +
  geom_tile() +
  scale_fill_gradient2(low = simcol['depleted'], mid = 'white', midpoint = 0, high = simcol['enriched']) +
  #scale_fill_continuous(mycol, breaks=mybreak) +
  geom_text(data = pstar,mapping = aes(x=Var1, y=Var2), label= pstar$star, inherit.aes = F) +
  xlab(NULL) + ylab(NULL) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 6, title='Log2(OR)',reverse = T)) + 
  theme(legend.direction = 'vertical',
        legend.position = 'right')
cooverlap

ggsave('results/co_overlaps.pdf', cooverlap, units = 'cm', height = 7, width = 10, useDingbats=F)
ggsave('results/co_overlaps.png', cooverlap, units = 'cm', height = 7, width = 10, bg='white')

ggsave('results/figure_supplements/fs4/co_overlaps.pdf', cooverlap, units = 'cm', height = 7, width = 10, 
       useDingbats=F)
ggsave('results/figure_supplements/fs4/co_overlaps.png', cooverlap, units = 'cm', height = 7, width = 10, 
       bg='white')

#
range(log2(pOR), na.rm=T)
# -0.06296104  2.52704236
pallen = 25
mybreak = c(seq(-2.52704236, 0, length.out=pallen/2 + 1),
             seq(1/pallen, 2.52704236, length.out = pallen/2 ) )
mycol = colorRampPalette(c('dodgerblue', 'white','#E4572E'))(pallen)
#mycol[c(1,50)] = 'gray99'
#hdat = log2(pOR)
#diag(hdat) = 0
pheatmap::pheatmap(log2(pOR), cluster_rows = F, cluster_cols = F, cellwidth = 25, cellheight = 25,
                   color = mycol,  breaks = mybreak, show_rownames = T, show_colnames = T,
                   legend_breaks = c(-2, -1.25, 0, 1.25, 2), fontsize = 8, 
                   filename = 'results/co_overlaps2.pdf')

pheatmap::pheatmap(log2(pOR), cluster_rows = F, cluster_cols = F, cellwidth = 25, cellheight = 25,
                   color = mycol,  breaks = mybreak, show_rownames = T, show_colnames = T,
                   legend_breaks = c(-2, -1.25, 0, 1.25, 2), fontsize = 8, 
                   filename = 'results/figure_supplements/fs4/co_overlaps2.pdf')

hmap = pheatmap::pheatmap(log2(pOR), cluster_rows = F, cluster_cols = F, cellwidth = 25, cellheight = 25,
                   color = mycol,  breaks = mybreak, show_rownames = T, show_colnames = T,
                   legend_breaks = c(-2, -1.25, 0, 1.25, 2), fontsize = 8, silent = T)

saveRDS(log2(pOR), 'results/source_data/f4/fs2b.rds')
# as_ggplot(hmap$gtable) +
#   geom_text(data = pstar,mapping = aes(x=Var1, y=Var2), label= pstar$star, inherit.aes = F) 


#### correlation between datasets in terms of age-related expression change:

ourdat = readRDS('data/processed/tidy/expression_change.rds') %>% 
  tibble() %>% 
  filter(period=='aging') %>%
  relocate(gene_id) %>%
  dplyr::select(-period, -p, -FDR) %>%
  rename(rho='Expression Change', Tissue='tissue')
jonker = readRDS('data/other_datasets/jonker/processed/agecor.rds') %>%
  relocate(gene_id) %>% 
  dplyr::select(-pval) %>%
  mutate(Tissue = ifelse(Tissue=='Brain', 'Cortex', Tissue))
# gtex4 = readRDS('results/GTEx/agecor.rds') %>%
#   rename(gene_id='GeneID') %>%
#   inner_join(ensmap) %>%
#   dplyr::select(-gene_id, -p, -adjusted_p) %>%
#   rename(gene_id='ENS_mm') %>% relocate(gene_id)
gtex10 = readRDS('results/GTEx/alltissues/agecor.rds') %>%
  rename(gene_id='GeneID') %>%
  inner_join(ensmap) %>%
  dplyr::select(-gene_id, -p, -adjusted_p) %>%
  rename(gene_id='ENS_mm') %>% relocate(gene_id) %>%
  mutate(Tissue = ifelse(Tissue=='Brain', 'Cortex', Tissue) )
# schaum4 = readRDS('data/other_datasets/schaum/4tissue/age_exp_cors.rds') %>%
#   inner_join(idconv) %>%
#   dplyr::select(-gene_id, -p) %>%
#   rename(gene_id='ensembl_gene_id') %>% relocate(gene_id)
schaum8 = readRDS('data/other_datasets/schaum/8tis/age_exp_cors.rds') %>%
  inner_join(idconv) %>%
  dplyr::select(-gene_id, -p) %>%
  rename(gene_id='ensembl_gene_id') %>% relocate(gene_id) %>%
  mutate(Tissue = ifelse(Tissue=='Brain', 'Cortex', Tissue) )

agecors = list(our_dat=ourdat, jonker=jonker, GTEx=gtex10, Schaum=schaum8)
agecors = reshape2::melt(agecors) %>%
  dplyr::select(-variable) %>%
  set_names(c('gene_id','Tissue', 'rho', 'dataset'))

cors = agecors %>%
  mutate(val = paste(dataset, Tissue,sep='-')) %>%
  reshape2::acast(gene_id~val, value.var = 'rho') %>%
  cor(method = 'spearman', use='pairwise.complete.obs')

# reshape2::melt(cors) %>%
#   ggplot(aes(x=Var1, y=Var2, fill=value)) +
#   geom_tile() +
#   scale_fill_gradient2(low = simcol['depleted'], mid = 'white', midpoint = 0, high = simcol['enriched']) +
#   guides(fill = guide_colorbar(barwidth = 0.5, barheight = 6, title='Rho',reverse = T)) +
#   theme(legend.direction = 'vertical',
#         legend.position = 'right')

max(abs(cors))
palettelen = 50
# mycol = colorRampPalette(c('dodgerblue','white','#E4572E'))(palettelen)
# mybreak = c(seq(min(cors), 0, length.out=palettelen/2 + 1),
#             seq(max(cors)/palettelen, max(cors), length.out = palettelen/2 ) )

dset = sapply(strsplit(rownames(cors), split = '-'),`[[`,1)
dsettis = sapply(strsplit(rownames(cors), split = '-'),`[[`,2)

annots = data.frame(Dataset = factor(dset), row.names = rownames(cors),
                       Tissue = factor(dsettis))
dsetcol = c('#B8002B', '#325D59', '#34836C', '#112C2B')

tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle') )
othertiscol = setNames(brewer.pal(12,'Set3')[-6],
                       levels(annots$Tissue)[!levels(annots$Tissue)%in%names(tissuecol)])

dsettiscol = c(tissuecol, othertiscol)

annotcols = list(Dataset = setNames(dsetcol, levels(annots$Dataset)),
                 Tissue = dsettiscol)

cors2=cors
diag(cors2)=0
range(cors2)
palettelen = 50
mybreak2 = c(seq(-0.52, 0, length.out=palettelen/2 + 1),
            seq(1/palettelen, 0.52, length.out = palettelen/2 ) )
mycol2 = colorRampPalette(c('dodgerblue','white','#E4572E'))(palettelen)
mycol2[c(1,50)] = 'gray99'
library(pheatmap)
pheatmap::pheatmap(cors, cluster_rows = T, cluster_cols = T, cellwidth = 5, cellheight = 5,
                   color = mycol2, gaps_row = c(10,15,19), gaps_col = c(10,15,19), breaks = mybreak2,
                   annotation_row = annots, annotation_col = annots, annotation_colors = annotcols,
                   show_rownames = F, show_colnames = F, annotation_names_row = F,
                   cutree_rows = 3, cutree_cols = 3, legend_breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
                   treeheight_col = 0, treeheight_row = 20,
                   fontsize = 6, width = 12/2.54, height = 8/2.54,
                   filename = 'results/figure_supplements/fs4/expch_cors.pdf')
# filename = '
pheatmap::pheatmap(cors, cluster_rows = T, cluster_cols = T, cellwidth = 5, cellheight = 5,
                   color = mycol2, gaps_row = c(10,15,19), gaps_col = c(10,15,19), breaks = mybreak2,
                   annotation_row = annots, annotation_col = annots, annotation_colors = annotcols,
                   show_rownames = F, show_colnames = F, annotation_names_row = F,
                   cutree_rows = 3, cutree_cols = 3, legend_breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
                   treeheight_col = 0, treeheight_row = 20,
                   fontsize = 6, width = 12/2.54, height = 8/2.54,
                   filename = 'results/figure_supplements/fs4/expch_cors.png')

saveRDS(cors, 'results/source_data/f4/fs2a.rds')

corsnew = reshape2::melt(cors) %>%
  mutate(dset1 = sapply(strsplit(as.character(Var1), split = '-'), `[[`,1)) %>%
  mutate(dset1_tis = sapply(strsplit(as.character(Var1), split = '-'), `[[`,2)) %>%
  mutate(dset2 = sapply(strsplit(as.character(Var2), split = '-'), `[[`,1)) %>%
  mutate(dset2_tis = sapply(strsplit(as.character(Var2), split = '-'), `[[`,2)) %>%
  dplyr::select(-Var1, -Var2) %>% tibble()
corsnew = corsnew %>%
  mutate(dset1_tis = ifelse(dset1_tis=='Adipose Tissue', 'Subcutaneous_Fat', dset1_tis)) %>%
  mutate(dset2_tis = ifelse(dset1_tis=='Adipose Tissue', 'Fat', dset2_tis))
corsnew %>%
  filter(dset1 == 'GTEx' & dset2 !='GTEx') %>%
  summarise(mean = mean(value))
# mean rho = 0.0153  

corsnew %>% 
  group_by(dset1, dset2) %>%
  summarise(m = mean(value)) %>% ungroup() %>% slice(2:4)
# dset1 dset2         m
# <chr> <chr>     <dbl>
# 1 GTEx  jonker  0.00748
# 2 GTEx  our_dat 0.0260 
# 3 GTEx  Schaum  0.0148 

corsnew %>% 
  group_by(dset1, dset2) %>%
  filter(dset1_tis == dset2_tis ) %>%
  summarise(m = mean(value))
# dset1   dset2         m
# <chr>   <chr>     <dbl>
#   1 GTEx    GTEx     1     
# 2 GTEx    jonker   0.116 
# 3 GTEx    our_dat  0.133 
# 4 GTEx    Schaum  0.0104

###
corsnew %>%
  filter(dset1=='GTEx' & dset2=='jonker' & dset1_tis==dset2_tis) %>%
  summarise(mean = mean(value)) # 0.116
corsnew %>% 
  filter(dset1=='GTEx' & dset2=='our_dat' & dset1_tis==dset2_tis) %>%
  summarise(mean = mean(value)) # 0.133
corsnew %>%
  filter(dset1=='GTEx' & dset2=='Schaum' & dset1_tis==dset2_tis) %>%
  summarise(mean = mean(value)) # 0.0104
mean(c(0.116, 0.133, 0.0104))
#  0.086

# same:
corsnew %>% 
  group_by(dset1, dset2) %>%
  filter(dset1_tis == dset2_tis ) %>%
  summarise(m = mean(value)) %>% ungroup() %>%
  slice(2:4) %>%
  summarise(mean=mean(m))
# 0.086

corsnew %>%
  filter(dset1=='GTEx') %>% summarise(unique(dset1_tis))
corsnew %>%
  filter(dset1=='our_dat') %>% summarise(unique(dset1_tis))
corsnew %>%
  filter(dset1=='jonker') %>% summarise(unique(dset1_tis))
corsnew %>%
  filter(dset1=='Schaum') %>% summarise(unique(dset1_tis))

