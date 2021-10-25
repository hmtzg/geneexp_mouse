library(openxlsx)

## Figure 2
cov_dat_sum = readRDS('results/source_data/f2/cov_dat_sum.rds')
top_divcon_cov_dat = readRDS('results/source_data/f2/top_divcon_cov_dat.rds')
top_divcon_gene_dat = readRDS('results/source_data/f2/top_divcon_gene_dat.rds')
cd_props = readRDS('results/source_data/f2/cd_props.rds')
cd_ratio = readRDS('results/source_data/f2/cd_ratio.rds')

# fs1
fs1 = readRDS('results/source_data/f2/fs1.rds') # cov_sum 

# fs2
fs2 = readRDS('results/source_data/f2/fs2.rds') # dccovkmX

# fs3
fs3 =  readRDS('results/source_data/f2/fs3.rds') # expcl

# fs4
fs4 = readRDS('results/source_data/f2/fs4.rds') # cd_count

# fs5
fs5 = readRDS('results/source_data/f2/fs5.rds') # pexpcors

# fs6
fs6a = readRDS('results/source_data/f2/fs6a.rds') # meancors
fs6b = readRDS('results/source_data/f2/fs6b.rds') # meansccors

# fs7
fs7ab = readRDS('results/source_data/f2/fs7_mean_median.rds') # cov_dat_sum
fs7c = readRDS('results/source_data/f2/fs7_pexpcors.rds') # pexpcors

# fs8
fs8abd_g = readRDS('results/source_data/f2/fs8_pca.rds') # pca_data
fs8c = readRDS('results/source_data/f2/fs8_euc.rds') # meanEuc

# fs9
fs9ab = readRDS('results/source_data/f2/fs9_CoV_mean_med.rds') # sumCov
fs9c = readRDS('results/source_data/f2/fs9_pairwisecors.rds') # pairwisedat

# fs10
fs10abd_g = readRDS('results/source_data/f2/fs10_pca.rds') # pca_data
fs10c = readRDS('results/source_data/f2/fs10_euc.rds') # meanEuc

# fs11
fs11ab = readRDS('results/source_data/f2/fs11_cov_mean_med.rds') # sumCov
fs11c = readRDS('results/source_data/f2/fs11_pairwisecors.rds') # pairwiseplotdat

# fs12
fs12 = readRDS('results/source_data/f2/fs12.rds') # dico_prop_permdat

# fs13
fs13 = readRDS('results/source_data/f2/fs13.rds') # indids

# fs14
fs14a = readRDS('results/source_data/f2/fs14_cov_dat_sum.rds') # cov_dat_sum
fs14b = readRDS( 'results/source_data/f2/fs14_top_divcon_cov_dat.rds') # top_divcon_cov_dat
fs14c = readRDS('results/source_data/f2/fs14_top_divcon_gene_dat.rds') # top_divcon_gene_dat
fs14d = readRDS('results/source_data/f2/fs14_cd_props.rds') # cd_props
fs14e = readRDS('results/source_data/f2/fs14_cd_ratio.rds') # cd_ratio

# fs15
fs15a = readRDS('results/source_data/f2/fs15_het.rds') # hetcorplot2dat
fs15b = readRDS('results/source_data/f2/fs15_ncvtest.rds') # hetcorplot3dat

# fs16
fs16ab = readRDS('results/source_data/f2/fs16_cov_mean_med.rds') # sumCovS
fs16cd = readRDS('results/source_data/f2/fs16_pairwisecors.rds') # pairwisedatsex


write.xlsx(list(sdata0 = NULL,
                sdata1 = cov_dat_sum,
                sdata2 = top_divcon_cov_dat,
                sdata3 = top_divcon_gene_dat,
                sdata4 = cd_props,
                sdata5 = cd_ratio,
                sdata6 = fs1, 
                sdata7 = fs2,
                sdata8 = fs3,
                sdata9 = fs4,
                sdata10 = fs5,
                sdata11 = fs6a,
                sdata12 = fs6b,
                sdata13 = fs7ab,
                sdata14 = fs7c,
                sdata15 = fs8abd_g, 
                sdata16 = fs8c, 
                sdata17 = fs9ab, 
                sdata18 = fs9c,
                sdata19 = fs10abd_g,
                sdata20 = fs10c,
                sdata21 = fs11ab,
                sdata22 = fs11c,
                sdata23 = fs12,
                sdata24 = fs13,
                sdata25 = fs14a,
                sdata26 = fs14b,
                sdata27 = fs14c,
                sdata28 = fs14d,
                sdata29 = fs14e,
                sdata30 = fs15a,
                sdata31 = fs15b,
                sdata32 = fs16ab,
                sdata33 = fs16cd), 
           file = './results/source_data/f2_source_data.xlsx' )




