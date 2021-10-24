library(tidyverse)
library(openxlsx)

source('./scripts/functions.R')
dev = readRDS("./data/htseq/no_blind/development_expression_change.rds")
aging = readRDS("./data/htseq/no_blind/ageing_expression_change.rds")

revg = sapply(names(dev), function(x) revgenes(dev[[x]][,1],aging[[x]][,1], monotonic = T), simplify = F )

saveRDS(revg,'./data/htseq/no_blind/revgenes.tissues.rds')

revgenes = revg %>%
  reshape2::melt(revgenes) %>%
  set_names('gene_id','direction','tissue')

saveRDS(revgenes,'./data/htseq/no_blind/revgenes.tissues_tidy.rds')


rev.common = sapply(names(revg$Cortex), function(y) Reduce(intersect,lapply(revg,function(x) x[[y]] )))
revgenes_common = reshape2::melt(rev.common) %>%
  set_names(c('gene_id','Pattern'))

saveRDS(revgenes_common,'./data/htseq/no_blind/revgenes.shared.rds')

head(revgenes)
head(revgenes_common)

# overlap with original result:
revgenes2 = readRDS('./data/processed/tidy/revgenes.tissues.rds')
head(revgenes2)
revgenes_common2 = readRDS('./data/processed/tidy/revgenes.shared.rds')
head(revgenes_common2)

ntot = revgenes %>% inner_join(revgenes2, by=c('gene_id', 'tissue')) %>% 
  group_by(tissue) %>%
  summarise(n=n())
# tissue     n
# <chr>  <int>
# 1 Cortex 14393
# 2 Liver  14442
# 3 Lung   14474
# 4 Muscle 14465

noverlap = revgenes %>% inner_join(revgenes2, by=c('gene_id', 'tissue')) %>% 
  group_by(gene_id, tissue) %>%
  summarise(sim = direction.x == direction.y) %>%
  group_by(tissue) %>% 
  summarise(n= sum(sim) )
# tissue     n
# <chr>  <int>
# 1 Cortex 10946
# 2 Liver  10153
# 3 Lung   11467
# 4 Muscle 10478

noverlap$n / ntot$n
# 0.7605086 0.7030190 0.7922482 0.7243692 

