library(openxlsx)
library(tidyverse)
## Figure 3

# a
cortex = readRDS('results/source_data/f3/a.rds')

# b
liver = readRDS('results/source_data/f3/b_liver.rds')

# c
lung = readRDS('results/source_data/f3/c_lung.rds')

# d
muscle = readRDS('results/source_data/f3/d_muscle.rds')

tsspec_UD_fisher = readRDS('./results/source_data/f3/tsspec_UD_enrichment.rds')
ud_table = reshape2::melt(tsspec_UD_fisher$table) %>%
  set_names(c('Specific','Reversal','count','Tissue'))

ud_or = reshape2::melt(tsspec_UD_fisher$OR) %>%
  set_names(c('OR','Tissue'))

tsspec_DiCo_fisher = readRDS('./results/source_data/f3/tsspec_dc_enrichment.rds')
dc_table = reshape2::melt(table(tsspec_DiCo_fisher$table)) %>%
  set_names(c('Specific','DiCo','count'))
tst = fisher.test(tsspec_DiCo_fisher$test)
dc_or = c(OR=tst$est,'p-val'=tst$p.value)

#tsspec_DiCo_fisher = unlist(tsspec_DiCo_fisher, recursive = F)

write.xlsx(list(sdata0 = NULL,
                sdata1 = cortex,
                sdata2 = liver,
                sdata3 = lung,
                sdata4 = muscle,
                sdata5 = ud_table,
                sdata6 = ud_or,
                sdata7 = dc_table,
                sdata8 = dc_or
                ),
           file = 'results/source_data/f3_source_data.xlsx')
