library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(goseq)

theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)

dc = readRDS('./data/processed/raw/dc_gse.rds')
signifGO = dc@result[,1:10] %>%
  filter( p.adjust< 0.1 & NES < 0) %>%
  dplyr::select(ID, Description, NES) #  184 DiCo enriched

# write.table(signifGO, file = 'results/figure4/signifGO_DiCo.csv', row.names = F, quote = F,
#             sep = '\t')

# get dico genes
genes = unique(names(which(readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')<0)))
allgo = getgo(genes,"mm9","ensGene")
allgo = allgo[!sapply(allgo,is.null)] #  4741 DiCo genes present in gos
allgo = reshape2::melt(allgo) %>%
  set_names(c('ID','gene'))
#

signifGO_genes = left_join(signifGO,allgo)
gogenemat = signifGO_genes %>%
  dplyr::select(ID, gene) %>%
  unique() %>%
  mutate(value = 1) %>%
  spread(ID, value, fill = 0)

gogenemat = as.data.frame(gogenemat)
rownames(gogenemat) = gogenemat$gene
gogenemat$gene = NULL
gogenemat = as.matrix(gogenemat)
# gogenemat: list of all genes present in GO categories 

jaccardsim = apply(gogenemat, 2, function(x){
  apply(gogenemat,2,function(y){
    sum(x==1 & y==1) / sum(x==1 | y==1)
  })
})

gocatx = signifGO
simmat = jaccardsim[gocatx$ID,gocatx$ID] # change column/rows orders back
k = 25
treex = hclust(dist( t(gogenemat) ))
treecl = cutree(treex, k)

# choose representative GO groups based on max mean similarity to other groups in the same cluster:
reps = sapply(1:k, function(i){
  xx=names(which(treecl==i))
  if(length(xx)>1){
    xx = simmat[xx,xx]
    names(which.max(rowMeans(xx)))
  } else{
    xx
  }
})

repclus = setNames(lapply(1:k,function(i) names(which(treecl==i)) ),reps)
newdf = reshape2::melt(repclus) %>% 
  set_names(c('ID','rep')) %>%
  arrange(rep) %>% 
  left_join(signifGO) %>%
  unique() 

# representative categories:
newdf %>%
  mutate(rep = ID == rep) %>%
  filter(rep) %>%
  dplyr::select(ID, Description)

# check median jaccard similarities of categories in the same cluster:
streps = sort(sapply(names(repclus), function(i){ median(jaccardsim[repclus[[i]], i]) }))
sort(sapply(names(repclus), function(i){ mean(jaccardsim[repclus[[i]], i]) }))

newdf %>%
  filter(rep%in%names(streps)[1]) %>%
  dplyr::select(Description) #  unrelated groups

newdf %>%
  filter(rep%in%names(streps)[2]) # related groups

newdf %>%
  filter(rep%in%names(streps)[3]) # related groups

newdf %>%
  filter(rep%in%names(streps)[23]) # related groups

reprgenes = signifGO_genes %>%
  filter(ID %in%names(repclus)) %>%
  dplyr::select(ID, gene, Description) %>%
  set_names('ID','gene_id', 'Description')

reprg = reprgenes %>%
  mutate(repNames = ifelse(ID%in%'GO:0030193','Other GO', Description ) ) 

reprg = reprg %>%
  group_by(ID) %>%
  summarise(n=n()) %>%
  arrange(n) %>%
  right_join(reprg)

saveRDS(reprg, file='./results/figure4/gorepresentatives.rds')

#####
#####
#####
## re-cluster the out-group: GO:0030193

sort(sapply(names(repclus), function(i){ median(jaccardsim[repclus[[i]], i]) }))
# out-group: 'GO:0030193'
outg = repclus[['GO:0030193']]

gogenemat2 = gogenemat[,outg]
jaccardsim2 = apply(gogenemat2,2,function(x){
  apply(gogenemat,2,function(y){
    sum(x==1 & y==1) / sum(x==1 | y==1)
  })
})

treex2 = hclust(dist( t(gogenemat2)))
##
k2 = 20
treecl2 = cutree(treex2,k2)
reps2 = sapply(1:k2, function(i){
  xx=names(which(treecl2==i))
  if(length(xx)>1){
    xx = simmat[xx,xx]
    names(which.max(rowMeans(xx)))
  } else{
    xx
  }
})

simmat['GO:0032543','GO:0072655']

repclus2 = setNames(lapply(1:k2,function(i)names(which(treecl2==i))),reps2)
newdf2 = reshape2::melt(repclus2) %>% 
  set_names(c('ID','rep')) %>%
  arrange(rep) %>% 
  left_join(signifGO) %>%
  unique() 

newdf2 %>%
  mutate(rep = ID == rep) %>%
  filter(rep) %>%
  dplyr::select(ID, Description)

sort(sapply(names(repclus2), function(i){ median(jaccardsim2[repclus2[[i]], i]) }))

# outgroup: median jaccardsim: 0.00725
newdf2 %>%
  filter(rep=='GO:0072577') 

reprgenes2 = signifGO_genes %>%
  filter(ID %in%names(repclus2)) %>%
  dplyr::select(ID, gene, Description) %>%
  set_names('ID','gene_id', 'Description')

reprg2 = reprgenes2 %>%
  mutate(repNames = ifelse(ID%in%'GO:0072577','Other GO', Description ) ) 

reprg2 = reprg2 %>%
  group_by(ID) %>%
  summarise(n=n()) %>%
  arrange(n) %>%
  right_join(reprg2)

reprg2

saveRDS(reprg2, file='./results/figure4/gorepresentatives2.rds')

sort(sapply(names(repclus2), function(i){ median(jaccardsim[repclus2[[i]], i]) }))
# out-group: 'GO:0072577'
outg = repclus2[['GO:0072577']]

################
################
################
################
##
expch = readRDS('./data/processed/tidy/expression_change.rds') %>%
  mutate(period = gsub('aging', 'Ageing', period)) %>%
  mutate(period = str_to_title(period) ) %>%
  mutate(period = factor(period, levels = c('Development', 'Ageing'))) %>%
  dplyr::rename(Period = period)

expch %>%
  inner_join(reprg) %>% head


enricplot_dat = expch %>%
  inner_join(reprg) %>%
  group_by(Period, ID, tissue, Description, repNames, n) %>%
  summarise(mrho = mean(`Expression Change`),
            medrho = median(`Expression Change`))

range(enricplot_dat$n)
#GO:0009611 : repsonse to wounding
# gx = reprg %>% filter(ID%in%'GO:0009611') %>% pull(gene_id)
# 
# expch %>%
#   filter(gene_id%in%gx) %>%
#   group_by(Period, tissue) %>%
#   mutate(mrho = median(`Expression Change`)) %>%
#   ggplot(aes(fill=Period, y=mrho, x=tissue)) +
#   geom_bar(stat='identity', position= position_dodge()) +
#   geom_point( aes(x=tissue, y=`Expression Change`), inherit.aes = F )

enricplot_ogr_dat = expch %>%
  inner_join(reprg2) %>%
  group_by(Period, ID, tissue, Description, repNames, n) %>%
  summarise(mrho = mean(`Expression Change`),
            medrho = median(`Expression Change`))

range(unique(enricplot_ogr_dat$n))

enricplot_ogr = expch %>%
  inner_join(reprg2) %>%
  group_by(Period, ID, tissue, Description, repNames, n) %>%
  summarise(mrho = mean(`Expression Change`),
            medrho = median(`Expression Change`)) %>%
  ggplot(aes(fill=Period, y=mrho, x=reorder(repNames, n) ) ) +
  geom_bar(stat='identity', position=position_dodge()) +
  facet_wrap(~tissue, ncol=4) +
  scale_fill_manual(values=brewer.pal(3,"Set1")[c(2,1)]) +
  coord_flip() +
  geom_hline(yintercept = 0, size=0.3, linetype='solid', color='gray30') +
  geom_vline(xintercept = seq(1.5,25, by=1), linetype = 'dashed', size = 0.2, color = 'gray') +
  theme(legend.position = 'top',
        axis.text.x = element_text(size=4, vjust=2), 
        axis.ticks.length.x = unit(0,'pt'),
        axis.text.y = element_text(size=5),
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(vjust = -0.5),
        legend.background = element_rect(color='black', size=0.1),
        legend.key.size = unit(3, 'pt'),
        axis.title.x = element_text(size=6)) +
  xlab('') +
  ylab(bquote('Mean Expression Change ('*rho*')'))
enricplot_ogr

saveRDS(enricplot_ogr_dat, 'results/source_data/f4/fs1.rds')

ggsave('./results/figure_supplements/fs4/FS1.pdf', enricplot_ogr, units='cm', width = 16, height = 12,
       useDingbats=F)
ggsave('./results/figure_supplements/fs4/FS1.png', enricplot_ogr, units='cm', width = 16, height = 12)  
# 
