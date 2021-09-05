library(tidyverse)
system('wget https://maayanlab.cloud/static/hdfs/harmonizome/data/transfacpwm/gene_attribute_matrix.txt.gz -P ./data/tfbs')
system('wget https://maayanlab.cloud/static/hdfs/harmonizome/data/transfacpwm/gene_attribute_edges.txt.gz -P ./data/tfbs')
system('wget https://maayanlab.cloud/static/hdfs/harmonizome/data/transfacpwm/gene_list_terms.txt.gz -P ./data/tfbs')
system('wget https://maayanlab.cloud/static/hdfs/harmonizome/data/transfacpwm/attribute_list_entries.txt.gz -P ./data/tfbs')

attr = read_tsv('./data/tfbs/gene_attribute_matrix.txt.gz')
attr = attr[,-2]
colnames(attr)[1] = 'GeneSym'
colnames(attr)[2] = 'GeneID'
tfgeneID = as.character(attr[2,])[-c(1,2)]
attr = attr[-c(1,2),]

attrlist =  read_tsv('./data/tfbs/attribute_list_entries.txt.gz')
attrlist$X2 = NULL
attrlist
identical(as.character(attrlist$GeneID), tfgeneID) #  same order

## hgnc gene id and tfbs target info:
edges = read_tsv('./data/tfbs/gene_attribute_edges.txt.gz')
summary(as.numeric(edges$source_desc))
edges$source_desc = NULL
summary(edges$target_id)
edges$target_desc=NULL

xx = t(as.data.frame(attr[attr$GeneSym%in%'ZNF767P',]))
xx[xx[,1]=="1.000000",]
edges[edges$source%in%'ZNF767P',]$target #  same info

## Hgnc gene list and ids:
genelist = read_tsv('./data/tfbs/gene_list_terms.txt.gz')
genelist$X2 =  NULL
genelist

sum(genelist$GeneSym%in%edges$source) # same list

## convert Hgnc ids to ENS:
martx = biomaRt::useMart(biomart = 'ensembl')
martHs_ = biomaRt::useDataset('hsapiens_gene_ensembl', mart=martx)
genemap = biomaRt::getBM(filters='hgnc_symbol', attributes = c('hgnc_symbol','ensembl_gene_id'),
                         values = genelist$GeneSym, mart = martHs_)
sum(duplicated(genemap$hgnc_symbol)) # 2594
sum(duplicated(genemap$ensembl_gene_id)) # 174
dup1 = unique(genemap$hgnc_symbol[duplicated(genemap$hgnc_symbol)]) #  remove duplicates
genemap = genemap[!genemap$hgnc_symbol%in%dup1,]
sum(duplicated(genemap$ensembl_gene_id)) # 0 # remove duplicates if any
genelist = genelist[genelist$GeneSym%in%genemap$hgnc_symbol,] # remove genes not mapped to ENS id

edges[1,4] = 'target_geneid'
edges[1,3] = 'target_GeneSym'
colnames(edges) = edges[1,]
edges = edges[-1,]
table(edges$weight) #  only target info

# subset  matrix with genes having only ENS ids:
edges = edges %>%  
  filter(GeneSym%in%genelist$GeneSym)
colnames(genemap) = c('GeneSym', 'ENS')
edges = edges %>% 
  left_join(genemap)

martMm_ = biomaRt::useDataset('mmusculus_gene_ensembl', mart=martx)
ensmap = biomaRt::getLDS(attributes = c('ensembl_gene_id'), filters = 'ensembl_gene_id', 
                         values = unique(edges$ENS), mart = martHs_, 
                         attributesL = c('ensembl_gene_id'), martL = martMm_)
colnames(ensmap) = c('ENS_hs', 'ENS_mm')
dup1 = unique(ensmap$ENS_hs[duplicated(ensmap$ENS_hs)]) #  691 duplicate genes
ensmap = ensmap[!ensmap$ENS_hs%in%dup1,]
dup2 = unique(ensmap$ENS_mm[duplicated(ensmap$ENS_mm)]) # 123 duplicate genes
ensmap = ensmap[!ensmap$ENS_mm%in%dup2,] 
# 14150 uniquely matching genes

edges = edges %>% 
  rename(ENS_hs = ENS) %>%
  right_join(ensmap)

saveRDS(edges,
        file = './data/tfbs/attr.rds')

save(list=ls(), file='./data/tfbs/prep.rdata')

