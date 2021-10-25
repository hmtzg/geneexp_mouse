library(openxlsx)

## Figure 3

# a
cortex = readRDS('results/source_data/f3/a.rds')

# b
liver = readRDS('results/source_data/f3/b_liver.rds')

# c
lung = readRDS('results/source_data/f3/c_lung.rds')

# d
muscle = readRDS('results/source_data/f3/d_muscle.rds')

write.xlsx(list(sdata0 = NULL,
                sdata1 = cortex,
                sdata2 = liver,
                sdata3 = lung,
                sdata4 = muscle),
           file = 'results/source_data/f3_source_data.xlsx')

