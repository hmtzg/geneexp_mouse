#library(tidyverse)
library(openxlsx)

##########
########## ########## Figure 1- source data 1

## f1a, f1b, f1-fs1, f1-fs2, f1-fs3
pca = readRDS('results/source_data/f1/pca.rds')
pc_cors = readRDS('results/source_data/f1/pc_age_cors.rds')
mdist = readRDS('results/source_data/f1/mdist.rds')
mdistdev = readRDS('results/source_data/f1/mdist_dev.rds')
mdistaging = readRDS('results/source_data/f1/mdist_aging.rds')
pc14anova = readRDS('results/source_data/f1/pc14_anova.rds')
## f1c
cors = readRDS('results/source_data/f1/cors.rds')
corsperm = readRDS('results/source_data/f1/f1pwise_expch_cor_perm_test.rds')
## f1d, f1-fs5, f1-fs7
exp_ch = readRDS('results/source_data/f1/expr_ch.rds')
## f1f
revgenes = readRDS('results/source_data/f1/revgenes.rds')
## f1-fs4
f1fs4 = readRDS('results/source_data/f1/fs4.rds')
## f1-fs6
f1fs6 = readRDS('results/source_data/f1/fs6.rds')
## f1-fs8
f1fs8 = readRDS('results/source_data/f1/fs8.rds')
## f1-fs9
f1fs9 = readRDS('results/source_data/f1/fs9.rds')
## f1-fs10
f1fs10ab = readRDS('results/source_data/f1/fs10pca.rds')
f1fs10c = readRDS('results/source_data/f1/fs10cor.rds')
f1fs10de = readRDS('results/source_data/f1/fs10expr_ch.rds')
f1fs10f = readRDS('results/source_data/f1/fs10revgenes.rds')
fs10_pc_age_cors = readRDS('./results/source_data/f1/fs10_pc_age_cors.rds')
## f1-fs11
f1fs11 = readRDS('results/source_data/f1/fs11.rds')
## f1-fs12
f1fs12 = readRDS('results/source_data/f1/fs12.rds')
f1fs13 = readRDS('results/source_data/f1/fs13.rds')
f1fs14 = readRDS('results/source_data/f1/fs14.rds')
f1fs15 = readRDS('results/source_data/f1/fs15.rds')

write.xlsx(list(sdata0 = NULL,
                sdata1 = pca, 
                sdata1.2 = pc14anova, # added new
                sdata2 = pc_cors, 
                sdata3 = mdist,
                sdata4 = mdistdev, 
                sdata5 = mdistaging, 
                sdata6 = cors,
                sdata7 = corsperm,
                sdata8 = exp_ch,
                sdata9 = revgenes,
                sdata10 = f1fs4,
                sdata11 = f1fs6, 
                sdata12 = f1fs8,
                sdata13 = f1fs9, 
                sdata14 = f1fs10ab,
                sdata14.1 = fs10_pc_age_cors, # added new
                sdata15 = f1fs10c,
                sdata16 = f1fs10de,
                sdata17 = f1fs10f,
                sdata18 = f1fs11,
                sdata19 = f1fs12,
                sdata20 = f1fs13,
                sdata21 = f1fs14,
                sdata22 = f1fs15), 
           file = './results/source_data/f1_source_data.xlsx')
