library(tidyverse)
library(goseq)
# get dico genes
genes = unique(names(which(readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')<0)))
allgo = getgo(genes,"mm9","ensGene")
allgo = allgo[!sapply(allgo,is.null)]
allgo = reshape2::melt(allgo) %>%
  set_names(c('ID','gene'))

GOres = readRDS('./data/processed/raw/dc_gse.rds')
signifGO = GOres@result %>% 
  filter(qvalues<=0.1) %>%
  dplyr::select(1,2,5) 

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
library(pheatmap)
jaccardsim = apply(gogenemat,2,function(x){
  apply(gogenemat,2,function(y){
    sum(x==1 & y==1) / sum(x==1 | y==1)
  })
})

gocatx = signifGO
simmat = jaccardsim[gocatx$ID,gocatx$ID]
k = 25
k=ifelse(k==1,nrow( t(gogenemat)),k)
treex =hclust(dist( t(gogenemat)))
treecl = cutree(treex,k)
reps = sapply(1:k,function(i){
  xx=names(which(treecl==i))
  if(length(xx)>1){
    xx = simmat[xx,xx]
    names(which.max(rowMeans(xx)))
  } else{
    xx
  }
})
repclus = setNames(lapply(1:k,function(i)names(which(treecl==i))),reps)
newdf = reshape2::melt(repclus) %>% 
                set_names(c('ID','rep')) %>%
                arrange(rep) %>% 
                left_join(signifGO) %>%
                unique() 

newdf %>%
  mutate(rep = ID == rep) %>%
  filter(rep) %>% 
  dplyr::select( Description)

