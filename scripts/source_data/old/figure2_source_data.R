library(tidyverse)
library(openxlsx)
##########
##########
cov_dat = readRDS('./data/processed/tidy/CoV.rds')
cov_ch = readRDS('./data/processed//tidy/CoV_change.rds') # sdata 1
expr = readRDS('./data/processed//tidy/expression.rds')
sample_info = readRDS('./data/processed/tidy/sample_info.rds')
cov_dat_wo_cortex = readRDS("./data/processed/tidy/CoV_wo_cortex.rds")

cov_dat_sum = cov_dat %>%
  mutate(ind_id = factor(ind_id)) %>%
  left_join(unique(select(sample_info,-tissue,-sample_id))) %>%
  group_by(ind_id, age) %>%
  summarise(meanCoV = mean(CoV),
            medianCoV = median(CoV)) %>% 
  mutate(period = c('aging','development')[1+(age<=90)]) # sdata 2

cov_sum_wo_cortex = cov_dat_wo_cortex %>%
  mutate(ind_id = factor(ind_id)) %>%
  left_join(unique(select(sample_info,-tissue,-sample_id))) %>%
  group_by(ind_id, age) %>%
  summarise(meanCoV = mean(CoV),
            medianCoV = median(CoV)) %>% 
  mutate(period = c('aging','development')[1+(age<=90)]) # sdata 3

top_divcon_cov_dat = cov_ch %>% # b
  select(-pval,-FDR) %>%
  spread(key = period, value = `CoV_change`) %>%
  mutate(rev = aging*development) %>%
  top_n(n=1, wt=-(rev)) %>%
  left_join(cov_dat) %>%
  mutate(ind_id = factor(ind_id)) %>%
  select(-development,-aging,-rev) %>%
  left_join( unique(select(sample_info,-tissue,-sample_id)) ) %>% 
  mutate(period = c('aging','development')[1+(age<=90)]) %>%
  select(-log2age) # sdata 4

top_divcon_gene_dat = cov_ch %>% # sdata 5
  select(-pval,-FDR) %>%
  spread(key = period, value =`CoV_change`) %>%
  mutate(rev = aging *development) %>%
  top_n(n=1,wt=-(rev)) %>%
  left_join(expr) %>%
  inner_join(sample_info) %>% 
  select(-log2age,-rev,-development,-aging) 

cd_props = cov_ch %>% # sdata 6
  filter(FDR < 0.1) %>%
  mutate(change = ifelse(CoV_change>0, "Diverge","Converge")) %>%
  select(gene_id, change,  period) %>%
  group_by(period,change) %>%
  summarise(n=n()) %>%
  mutate(period =  gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels=c("Development", "Ageing")))

cd_ratio_CI = readRDS('data/processed/tidy/cov_ratio_jk_CI.rds') 

cd_ratio = cd_props %>%
  mutate(period =  gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels=c("Development", "Ageing"))) %>%
  group_by(period) %>%
  summarise(Ratio = n[change=="Converge"]/n[change=="Diverge"]) %>%
  mutate(facet="Conv/Div Ratio") %>%
  mutate(CImin = cd_ratio_CI[1,] ) %>%
  mutate(CImax = cd_ratio_CI[2,] )   # sdata 7

DiCo_cl_genes = readRDS('./data/processed/tidy/DiCo_cl_genes.rds') # sdata 8

pexpcors  = readRDS('./data/processed/tidy/pairwise_tissue_expression_cors.rds') # sdata 9

meancors = pexpcors %>% 
  group_by(id, uage) %>%
  summarise(mean = mean(rho),
            median = median(rho)) %>%
  gather(key='method', value = 'rho', mean, median) %>%
  mutate(period = factor(c('development','ageing')[1 + (uage > 90)], levels = c('development', 'ageing')))
# sdata 10

pcors.sc = pexpcors %>% group_by(pair) %>% mutate(rho = scale(rho)[,1])
meansccors = pcors.sc %>% 
  group_by(id, uage) %>%
  summarise(mean = mean(rho),
            median = median(rho)) %>%
  gather(key='method', value = 'rho', mean, median) %>%
  mutate(period = factor(c('development','ageing')[1 + (uage > 90)], levels = c('development', 'ageing')))
# sdata 11

Jpexpcors = readRDS( './data/other_datasets/jonker/processed/pwise_exp_cors.rds') # sdata 12
Jcov_dat_sum = readRDS(file = './data/other_datasets/jonker/processed/mean_cov.rds') #  sdata 13

Gtex_pca = read.csv('./results/source_data/gtex_sourcedata/fig2ext9/fig2supp9abdefg.csv') # sdata 14
Gtex_mdist = read.csv('./results/source_data/gtex_sourcedata/fig2ext9/fig2supp9c.csv') # sdata 15

Gtex_cov = read.csv('./results/source_data/gtex_sourcedata/fig2ext10/fig2supp10ab.csv') # sdata 16
Gtex_pwisecors = read.csv('./results/source_data/gtex_sourcedata/fig2ext10/fig2supp10c.csv') # sdata 17

Gtex_10ts_pca = read.csv('./results/source_data/gtex_sourcedata/fig2ext11/fig2_supp11abdefg.csv') # sdata 18
Gtex_10ts_mdist = read.csv('./results/source_data/gtex_sourcedata/fig2ext11/fig2_supp11c.csv') # sdata 19

Gtex_10ts_cov = read.csv('./results/source_data/gtex_sourcedata/fig2ext12/fig2_supp12ab.csv') # sdata 20
Gtex_10ts_pwisecors = read.csv('./results/source_data/gtex_sourcedata/fig2ext12/fig2_supp12c.csv') # sdata 21

gtex_heatmap = read.csv('./results/source_data/gtex_sourcedata/fig2extHEATMAP/fig2_suppHEATMAP.csv') # sdata 22

DiCodist = readRDS('./data/processed/tidy/DiCo_sig_perm.rds') # sdata 23

write.xlsx(list(sdata0 = NULL,
                sdata1 = cov_ch,
                sdata2 = cov_dat_sum,
                sdata3 = cov_sum_wo_cortex,
                sdata4 = top_divcon_cov_dat,
                sdata5 = top_divcon_gene_dat,
                sdata6 = cd_props, 
                sdata7 = cd_ratio,
                sdata8 = DiCo_cl_genes,
                sdata9 = pexpcors,
                sdata10 = meancors,
                sdata11 = meansccors,
                sdata12 = Jpexpcors,
                sdata13 = Jcov_dat_sum,
                sdata14 = Gtex_pca,
                sdata15 = Gtex_mdist, 
                sdata16 = Gtex_cov, 
                sdata17 = Gtex_pwisecors, 
                sdata18 = Gtex_10ts_pca,
                sdata19 = Gtex_10ts_mdist,
                sdata20 = Gtex_10ts_cov,
                sdata21 = Gtex_10ts_pwisecors,
                sdata22 = gtex_heatmap,
                sdata23 = DiCodist), 
           file = './results/source_data/figure_2_source_data_1.xlsx' )



