library(tidyverse)
library(openxlsx)
library(ggpubr)

sample_info = readRDS('./data/processed/tidy/sample_info.rds')
expch = readRDS('./data/htseq/blinded/expression_change.rds')
expch_qn = readRDS('./data/processed/tidy/expression_change.rds')
colnames(expch_qn)[3] = 'Exp. Ch. (QN)'
expr = readRDS('./data/htseq/expr_blinded.rds')
expr_qn = readRDS('./data/processed/tidy/expression.rds')
colnames(expr_qn)[3] = 'expression (qn)'
periodcol = c(Development = "#FE6100", Ageing ="#648FFF")

theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)

expch_qn = expch_qn %>% 
  mutate(period = gsub('development','Development', period) ) %>% 
  mutate(period = gsub('aging', 'Ageing', period) )
expch = expch[complete.cases(expch),]

# correlation between age trajectories of QN and VST normalisation methods:
expch_cor = expch %>%
  select(-p, -FDR) %>% 
  inner_join(select(expch_qn, -p, -FDR) ) %>% 
  group_by(tissue, period) %>% 
  summarise(cor = cor(`Expression Change`, `Exp. Ch. (QN)`))

## number of genes :
expch %>%
select(-p, -FDR) %>% 
  inner_join(select(expch_qn, -p, -FDR) ) %>%
  group_by(period, tissue) %>%
  summarise(n = length(gene_id)) 
# period      tissue     n
# <chr>       <chr>  <int>
# 1 Ageing      Cortex 14705
# 2 Ageing      Liver  14689
# 3 Ageing      Lung   14710
# 4 Ageing      Muscle 14708
# 5 Development Cortex 14707
# 6 Development Liver  14705
# 7 Development Lung   14710
# 8 Development Muscle 14710

# same but give range:
expch %>%
  select(-p, -FDR) %>% 
  inner_join(select(expch_qn, -p, -FDR) ) %>%
  group_by(period, tissue) %>%
  summarise(n = length(gene_id)) %>%
  summarise(range = range(n))

# expch_cors = expch_cor %>% 
#   mutate(period=factor(period, levels = c('Development', 'Ageing') )) %>%
#   ggplot(aes(x=tissue, y=cor, fill=period)) +
#   #facet_wrap(~period) +
#   geom_bar(stat='identity', position = 'dodge') +
#   scale_fill_manual(values =periodcol) +
#   scale_y_continuous(limits = c(0,1)) +
#   xlab('') +
#   ylab('Spearman\'s correlation coefficient') +
#   guides(fill = guide_legend('Period'))
# expch_cors

# ggsave('./results/htseq/blinded/expch_cors_qn_vs_vst.pdf', expch_cors, units='cm', height = 5, width = 8,
#        useDingbats=F)
# ggsave('./results/htseq/blinded/expch_cors_qn_vs_vst.png', expch_cors, units='cm', height = 5, width = 8)

expchcorplotdat = expch %>%
  select(-p, -FDR) %>% 
  inner_join(select(expch_qn, -p, -FDR) ) %>%
  mutate(period = factor(period, levels = c('Development', 'Ageing')))

expchcorplot = expchcorplotdat %>%
  ggplot(aes(x = `Exp. Ch. (QN)`, y=`Expression Change`)) +
  facet_grid(period~tissue) +
  geom_point(size = 0.4, alpha=0.4, color= 'lightblue') +
  geom_smooth(method = 'lm',se = F, color='midnightblue' , size=0.4 ) +
  stat_cor(method='spearman', cor.coef.name = 'rho', size=7/pntnorm ) +
  xlab('Expression Change (QN)') +
  ylab('Expression Change (VST)') +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        strip.text = element_text(size=6))

expchcorplot
ggsave('./results/htseq/blinded/expch_cors_plt.pdf', expchcorplot, units='cm', height = 12, width = 16,
       useDingbats=F)
ggsave('./results/htseq/blinded/expch_cors_plt.png', expchcorplot, units='cm', height = 12, width = 16)

ggsave('./results/figure_supplements/f1s/FS11.pdf', expchcorplot, units='cm', height = 12, width = 16,
       useDingbats=F)
ggsave('./results/figure_supplements/f1s/FS11.png', expchcorplot, units='cm', height = 12, width = 16)

saveRDS(expchcorplotdat, 'results/source_data/f1/fs11.rds')

expch_cor %>% arrange(period)
# tissue period        cor
# <chr>  <chr>       <dbl>
# 1 Cortex Ageing      0.825
# 2 Liver  Ageing      0.870
# 3 Lung   Ageing      0.926
# 4 Muscle Ageing      0.845
# 5 Cortex Development 0.937
# 6 Liver  Development 0.810
# 7 Lung   Development 0.897
# 8 Muscle Development 0.811

saveRDS(expch_cor, file='./data/htseq/blinded/expch_cor_qn_vs_vst.rds')

# expch_qn = expch_qn %>% 
#   set_names(c('gene_id', 'tissue', 'QN', 'p_VST', 'FDR_VST','period'))

expch_qn2 = expch_qn %>%
  filter(FDR < 0.1) %>%
  select(-FDR, -p)

head(expch_qn2)

# for significant genes:
xxs = expch %>% 
  filter(FDR < 0.1) %>% 
  select(-FDR, -p) %>% 
  inner_join(expch_qn2) %>%
  group_by(tissue, period) %>%
  summarise(n=length(gene_id) ) %>% arrange(period)

### total sig genes in both datasets:
xxu = expch %>% 
  filter(FDR < 0.1) %>% 
  select(-FDR, -p) %>% 
  full_join(expch_qn2) %>%
  group_by(tissue, period) %>%
  summarise(n=length(gene_id) ) %>% arrange(period)

## proportion of shared sig genes:
xxs %>% full_join(xxu, by=c('tissue', 'period')) %>%
  mutate(prop = n.x/n.y*100)
# tissue period        n.x   n.y  prop
# <chr>  <chr>       <int> <int> <dbl>
#   1 Cortex Ageing         22   125  17.6
# 2 Liver  Ageing         41   241  17.0
# 3 Lung   Ageing       1659  3148  52.7
# 4 Cortex Development  5030  7543  66.7
# 5 Liver  Development  3214  6359  50.5
# 6 Lung   Development  2512  4918  51.1
# 7 Muscle Development  1365  3178  43.0

expch %>% 
  filter(FDR < 0.1) %>% 
  select(-FDR, -p) %>% 
  full_join(expch_qn2) %>% 
  group_by(tissue, period) %>% 
  filter(! (tissue%in%'Muscle' & period =='Ageing') ) %>%
  summarise(cor = cor.test(`Expression Change`, `Exp. Ch. (QN)`, na.rm=T, m='s')$est) %>% arrange(period)
# tissue period        cor
# <chr>  <chr>       <dbl>
#  1 Cortex Ageing      0.805
# 2 Cortex Development 0.901
# 3 Liver  Ageing      0.828
# 4 Liver  Development 0.887
# 5 Lung   Ageing      0.876
# 6 Lung   Development 0.879
# 7 Muscle Development 0.847

#### Compare CoV significant genes

cov_qn = readRDS('./data/processed/tidy/CoV_change.rds') %>% tibble()
cov_vst = readRDS('./data/htseq/blinded/CoV_change.rds') %>% tibble()
cov_qn %>% filter(FDR<0.1) %>%
  mutate(change = ifelse(CoV_change>0, "Diverge","Converge")) %>%
  select(gene_id, change, period) %>%
  group_by(period, change) %>%
  summarise(n=n())
# period      change       n
# <fct>       <chr>    <int>
#   1 development Converge   772
# 2 development Diverge   1809
# 3 aging       Converge    42
# 4 aging       Diverge     20

cov_vst %>% filter(FDR<0.1) %>%
  mutate(change = ifelse(CoV_change>0, "Diverge","Converge")) %>%
  select(gene_id, change, period) %>%
  group_by(period, change) %>%
  summarise(n=n())
# period      change       n
# <chr>       <chr>    <int>
#   1 Ageing      Converge    13
# 2 Ageing      Diverge      6
# 3 Development Converge   398
# 4 Development Diverge   3078

cov_qn %>% filter(period=='development' & FDR < 0.1) %>%
  mutate(change = ifelse(CoV_change>0, "Diverge","Converge")) %>%
  select(gene_id, change) %>%
  rename(change_qn = change) %>%
  inner_join(cov_vst %>% 
              filter(period=='Development' & FDR < 0.1) %>%
              mutate(change = ifelse(CoV_change>0, "Diverge","Converge")) %>%
              select(gene_id, change) %>%
              rename(change_vst = change) ) %>%
  summarise(shared = table(change_qn, change_vst))


