library(tidyverse)
#library(clusterProfiler)
sgo = readRDS('./data/processed/raw/dc_gse.rds')
#ddc = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')

#sum(sapply(sgo@geneSets, function(x) sum(is.na(x)))) # no empty go
## get all genes for all go groups:
allgo = reshape2::melt(sgo@geneSets) %>% 
  set_names('gene','ID')

## get only significant go:
signifGO = sgo@result %>% 
  filter(qvalues < 0.1 & NES<0) 

# filter only signif go :
signifGO_genes = signifGO %>% 
  select(ID) %>%
  left_join(allgo)

## gene category matrix:
gogenemat = signifGO_genes %>%
  mutate(value=1) %>%
  spread(ID, value, fill=0)

#gogenemat = as.data.frame(gogenemat)
rownames(gogenemat) = gogenemat$gene
gogenemat$gene = NULL
gogenemat = as.matrix(gogenemat)

# jaccard similarity:
jaccardsim = apply(gogenemat, 2, function(x){
  apply(gogenemat, 2, function(y){
    sum( x == 1 & y == 1) / sum( x == 1 | y == 1)
  })
})

## define go clusters with highest similarity with threshold:
#gocat = unique(signifGO$ID)

k = min(which(sapply(1:nrow(jaccardsim), function(i){
  cl = cutree(hclust(dist(jaccardsim)),i)
  #cl = cutree(hclust(as.dist(1-jaccardsim)),i)
  all(sapply(1:i, function(k){
    xx = names(which(cl==k))
    if(length(xx)==1){
      return(1)
    } else {
      xx = jaccardsim[xx,xx]
      median(xx[upper.tri(xx)])
    }
  }) >= 0.2)
})))
k

#k = ifelse(k==1, nrow(jaccardsim), k)
treex = hclust(dist(jaccardsim))
treecl = cutree(treex, k)

# get representative go group for go clusters:
reps = sapply(1:k, function(i){
  xx = names(which(treecl == i))
  if(length(xx)>1){
    xx = jaccardsim[xx,xx]
    names(which.max(rowMeans(xx)))
  } else{
    xx
  }
})
reps
repclus = setNames(lapply(1:k, function(i) names(which(treecl==i))), reps)

newdf = reshape2::melt(repclus) %>%
  set_names('ID', 'rep') %>% 
  arrange(rep) %>% 
  left_join(signifGO) %>%
  select(-core_enrichment)

head(newdf)

dico_go = newdf %>% 
  filter(ID==rep) %>% 
  mutate(termname = ifelse(nchar(Description)>40,
                           paste(substr(Description,1,37),'...',sep=''), Description)) %>%
  ggplot(aes(x=NES, y=reorder(termname, -NES))) +
  geom_bar(stat='identity') +
  theme(axis.text = element_text(size=6)) +
  ylab('')

ggsave('results/figure4/fig4g.pdf', dico_go, units = 'cm', width = 16, height = 16, useDingbats=F)

newdf[newdf$Description=='leukocyte chemotaxis',]

newdf[newdf$Description=='chemotaxis',] # GO:0006935
newdf[newdf$Description=='leukocyte chemotaxis',] # GO:0030595
newdf[newdf$Description=='granulocyte chemotaxis',] # GO:0071621

# repclus[names(repclus)%in%'GO:0006935']
# repclus[names(repclus)%in%'GO:0030595']
# repclus[names(repclus)%in%'GO:0071621']
# cor(jaccardsim[,'GO:0042330'],jaccardsim[,'GO:0006935'])
# cor(jaccardsim[,'GO:0030595'],jaccardsim[,'GO:0097529'])
# cor(jaccardsim[,'GO:0071621'],jaccardsim[,'GO:0030593'])
# cor(jaccardsim[,'GO:0071621'],jaccardsim[,'GO:0097530'])
# cor(jaccardsim[,'GO:0071621'],jaccardsim[,'GO:0042330'])
# head(repclus2)
# int = intersect(names(repclus), names(repclus2) )
# head(repclus[int])
# head(repclus2[int])
# a = 'GO:0006935'
# b = 'GO:0030595'
# sum(gogenemat[,a]==1 & gogenemat[, b] == 1) /
#   sum(gogenemat[,a]==1 | gogenemat[, b] == 1)
# d1 = dist(jaccardsim)
# d2 = as.dist(1-jaccardsim)
# head(d1)
# head(d2)
# cor.test(d1,d2)
# plot(d1/max(d1), d2)
######################################

#
signifGO = sgo@result[,1:10] %>%
  filter(qvalues < 0.1 & NES < 0) %>%
  select(ID, NES)
write.table(signifGO, file = 'results/figure4/signifGO_DiCo.csv', row.names = F, quote = F)  



