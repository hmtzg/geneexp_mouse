library(tidyverse)
library(openxlsx)
##########
ts.exprch = readRDS('./data/processed/tidy/ts.spec.exprho.rds') # sdata1
ts.es = readRDS('./data/processed/tidy/ts.spec.ES.rds') # sdata1

es.expr = ts.exprch %>%
  left_join(ts.es, by=c('gene id','tissue','spec'))

spec.rev.table = readRDS('./data/processed/raw/tis_spec_rev_enrichment.rds')[[2]]
spec.rev.OR = readRDS('./data/processed/raw/tis_spec_rev_enrichment.rds')[[1]]

spec.rev.table = reshape2::melt(lapply(spec.rev.table, data.frame)) %>% # sdata2
  select(-variable) %>%
  set_names(c('Tis. Spec.','Reversal', 'Count','Tissue'))

spec.rev.OR = sapply(spec.rev.OR, function(x) x) # sdata3
rownames(spec.rev.OR)[2] = 'p.val'

spec.dc.mat = readRDS('./data/processed/raw/ts_spec_dc_enrichment.rds')
spec.dc = data.frame(table(spec.dc.mat)) %>%
  set_names(c('Spec. to a Tiss.','DiCo','Count')) # sdata4

write.xlsx(list(sdata0 = NULL,
                sdata1 = es.expr,
                sdata2 = spec.rev.table,
                sdata3 = spec.rev.OR,
                sdata4 = spec.dc),
           file = 'results/source_data/figure_3_source_data_1.xlsx')
