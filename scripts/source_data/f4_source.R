library(openxlsx)

## Figure 4

# a
all_fisher = readRDS('results/source_data/f4/a.rds')

# b
dc_fisher = readRDS('results/source_data/f4/b.rds')

# c
gr1dat = readRDS( 'results/source_data/f4/c.rds')

# d
gr2dat = readRDS( 'results/source_data/f4/d.rds')

# e
gr3dat = readRDS('results/source_data/f4/e.rds')

# f
gr4dat = readRDS('results/source_data/f4/f.rds')

# g
enricplotdat = readRDS('results/source_data/f4/g.rds')

## fs1
enricplot_ogr_dat = readRDS('results/source_data/f4/fs1.rds')

##
mirtest = readRDS('results/source_data/f4/mirna.rds')

tf = readRDS('./results/source_data/f4/tf.rds')

write.xlsx(list(sdata0 = NULL,
                sdata1 = all_fisher,
                sdata2 = dc_fisher,
                sdata3 = gr1dat,
                sdata4 = gr2dat,
                sdata5 = gr3dat,
                sdata6 = gr4dat, 
                sdata7 = enricplotdat,
                sdata8 = enricplot_ogr_dat,
                sdata9 = mirtest,
                sdata10 = tf
                ), 
           file = './results/source_data/f4_source_data.xlsx' )

