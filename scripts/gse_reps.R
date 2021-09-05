library(tidyverse)
library(clusterProfiler)
ddc_genes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')

cnt=0
library(org.Mm.eg.db)
gse_reps = list()
for(i in 1:20){
  gse_repX = lapply(1:5, function(x){
    cnt <<- cnt+1
    print(cnt)
    dc_gse = gseGO(geneList = sort(ddc_genes, decreasing = T),  ont = "BP",
                   pvalueCutoff = 1, OrgDb = org.Mm.eg.db ,
                   keyType = "ENSEMBL", minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BH',
                   verbose = F)@result
    dc_gse = dc_gse[dc_gse$p.adjust<0.1,c('ID', 'Description', 'NES', 'p.adjust')]
  })
  gse_reps = c(gse_reps, gse_repX)
}





saveRDS(gse_reps,file= './data/processed/raw/gse_reps.rdata')
# dc_gse = readRDS(file="./data/processed/raw/dc_gse.rds")
# dc_gse_genelist =  strsplit(dc_gse@result[,11], split = '/')
# names(dc_gse_genelist) = dc_gse@result[,'ID']
# dc_gse_genelist =reshape2::melt(dc_gse_genelist) %>% 
#   set_names(c('gene_id','GO_ID'))
# dc_gse_table = list(enrichment = dc_gse@result[,1:10],
#                     genelist = dc_gse_genelist)
# write.xlsx(dc_gse_table, './results/SI_tables/TableS11.xlsx')
