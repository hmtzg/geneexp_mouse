library(openxlsx)

## Figure 5
# a
panel_a_dat = readRDS('results/source_data/f5/a.rds')

# b
maxmin_density = readRDS('results/source_data/f5/b1.rds')
maxminch = readRDS('results/source_data/f5/b2.rds')

# fs1
deconv = readRDS( 'results/source_data/f5/fs1.rds')

# fs2
p1dat = readRDS('results/source_data/f5/fs2.rds')

# fs3
p2dat = readRDS('results/source_data/f5/fs3.rds')

# fs4
p3dat = readRDS('results/source_data/f5/fs4.rds')

# fs5
p4dat = readRDS('results/source_data/f5/fs5.rds')

# fs6
intracovdat = readRDS('results/source_data/f5/fs6.rds')

write.xlsx(list(sdata0 = NULL,
                sdata1 = panel_a_dat,
                sdata2 = maxmin_density,
                sdata3 = maxminch,
                sdata4 = deconv,
                sdata5 = p1dat,
                sdata6 = p2dat,
                sdata7 = p3dat,
                sdata8 = p4dat,
                sdata9 = intracovdat), 
           file = './results/source_data/f5_source_data.xlsx' )
