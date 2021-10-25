library(tidyverse)
library(openxlsx)
##########

ts.spec = read.csv('./results/source_data/fig4_source/figure4a.csv') # sdata1
ts.spec.dico = read.csv('./results/source_data/fig4_source/figure4b.csv') # sdata2
g1 = read.csv('./results/source_data/fig4_source/figure4c.csv') # sdata3
g2 = read.csv('./results/source_data/fig4_source/figure4d.csv') # sdata4
g3 = read.csv('./results/source_data/fig4_source/figure4e.csv') # sdata5
g4 = read.csv('./results/source_data/fig4_source/figure4f.csv') # sdata6

dc.gse = readRDS('./data/processed/raw/dc_gse.rds') # sdata 7
dc.gse[1:10,1:10]
# dc_gse_genelist =  strsplit(dc.gse@result[,11], split = '/')
# names(dc_gse_genelist) = dc_gse@result[,'ID']
# dc_gse_genelist = reshape2::melt(dc_gse_genelist) %>% 
#   set_names(c('gene_id','GO_ID'))
# 
# dc_gse_table = list(enrichment = dc.gse@result[,1:10],
#                     genelist = dc_gse_genelist)


write.xlsx(list(sdata0 = NULL,
                sdata1 = ts.spec,
                sdata2 = ts.spec.dico,
                sdata3 = g1,
                sdata4 = g2,
                sdata5 = g3,
                sdata6 = g4,
                sdata7 = dc.gse[,1:10]),
           file='./results/source_data/figure_4_source_data_1.xlsx')
