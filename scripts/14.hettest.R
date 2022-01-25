library(tidyverse)
library(car)
library(ggpubr)
source('./scripts/functions.R')
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))

cov =  readRDS("./data/processed/tidy/CoV.rds")
covch = readRDS("./data/processed/tidy/CoV_change.rds")
expch = readRDS('./data/processed/tidy/expression_change.rds')
exp = readRDS('./data/processed/tidy/expression.rds')
sinfo = read_rds('./data/processed/tidy/sample_info.rds')

devdiv =  readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
dico =  names(which(devdiv <0  ))

# heterogeneity using absolute residuals from linear regressin:
het = exp %>% 
  left_join(select(sinfo,-log2age) ) %>%
  filter(age>90) %>%
  mutate(age= log2(age) ) %>%
  group_by(gene_id, tissue) %>%
  mutate( resi = abs( lm(expression ~ age )$resi ) ) %>% ungroup()

# heterogeneity using ncvTest :
het2 = exp %>% 
  left_join(select(sinfo,-log2age) ) %>%
  filter(age>90) %>%
  mutate(age= log2(age)) %>%
  group_by(gene_id, tissue) %>%
  summarise( chi = ncvTest( lm(expression ~ age ) )$ChiSq ) %>% ungroup()

### change in heterogeneity for abs residual data:
hetcor = het %>%
  group_by(gene_id, tissue) %>%
  summarise( rho = cor.test(resi, age, m='s')$est,
             p = cor.test(resi, age, m='s')$p.val) %>% ungroup()

### number of genes showing significant heterogeneity change:
hetcor %>%
  group_by(tissue) %>%
  mutate(FDR = p.adjust(p, method ='BH')) %>%
  mutate(dir = ifelse(rho>0,'inc','dec')) %>%
  group_by(tissue, dir) %>%
  summarise(count = sum(FDR < 0.1, na.rm=T))
# tissue dir   count
# <fct>  <chr> <int>
# 1 Cortex dec      56
# 2 Cortex inc       9
# 3 Liver  dec      15
# 4 Liver  inc       9
# 5 Lung   dec      10
# 6 Lung   inc      70
# 7 Muscle dec      12
# 8 Muscle inc      15

## heterogeneity change for development-divergent genes:
hetcorDi = hetcor %>% 
  filter(gene_id%in%names(devdiv)) %>%
  mutate(type = ifelse(gene_id%in%dico, 'DiCo', 'DiDi') ) %>%
  mutate(type= factor(type, levels=c('DiDi','DiCo')))

hetcorDi %>%
  group_by(tissue) %>%
  mutate(FDR = p.adjust(p, method ='BH')) %>%
  mutate(dir = ifelse(rho>0,'inc','dec')) %>%
  group_by(tissue, dir) %>%
  summarise(count = sum(FDR < 0.1, na.rm=T))
# tissue dir   count
# <fct>  <chr> <int>
# 1 Cortex dec      34
# 2 Cortex inc       5
# 3 Liver  dec      11
# 4 Liver  inc       5
# 5 Lung   dec       2
# 6 Lung   inc      13
# 7 Muscle dec       2
# 8 Muscle inc       3

# violin plot for heterogeneity changes of DiCo and DiDi genes:
# hetcorplot = hetcorDi %>% 
#   ggplot(aes(x=tissue, y= rho, fill=type )) +
#   geom_violin(scale='count', trim=T)  +
#   geom_hline(yintercept=0, linetype='dashed', color='gray50') +
#   geom_boxplot(width=0.1, position = position_dodge(width = .9), show.legend = F, outlier.size = 0.3) +
#   ylab('Heterogeneity change') +
#   xlab('') + 
#   theme(legend.position = 'right')
# hetcorplot
# ggsave('./results/SI_figures/heteroscedasticity/hetchange.pdf', hetcorplot, units='cm', height = 8, 
#        width = 12, useDingbats=F)
# ggsave('./results/SI_figures/heteroscedasticity/hetchange.png', hetcorplot, units='cm', height = 8, 
#        width = 12)

##### ks.test for heterogeneity change  for dico vs didi genes:
annothet = hetcorDi %>% 
  select(-p) %>% 
  spread(key=type, value = rho) %>%
  group_by(tissue) %>%
  summarise(p = ks.test(DiCo, DiDi)$p.val,
            D= ks.test(DiCo, DiDi)$stat)
# tissue      p
# <fct>   <dbl>
#   1 Cortex 0.340 
# 2 Liver  0.534 
# 3 Lung   0.136 
# 4 Muscle 0.0496

# density plots of heterogeneity changes for dico vs didi genes:

hetcorplot2dat = hetcorDi %>% 
  rename('Pattern' = type) %>%
  mutate(Pattern = factor(Pattern, levels=c('DiCo', 'DiDi') ))
  
hetcorplot2 = hetcorDi %>% 
  rename('Pattern' = type) %>%
  mutate(Pattern = factor(Pattern, levels=c('DiCo', 'DiDi') )) %>%
  ggplot(aes(x= rho, fill=Pattern, linetype=Pattern, size=Pattern )) +
  facet_wrap(~tissue) +
  geom_vline(xintercept=0, linetype='dashed', color='gray20') +
  geom_density(alpha=0.3) +
  scale_linetype_manual(values = c('DiCo'='solid', 'DiDi'='dashed') ) + 
  scale_size_manual(values=c(0.4, 0.4)) +
  ylab('Density') +
  xlab('Heterogeneity change') +
  geom_text(data = annothet, inherit.aes = F, parse=T, hjust=0, size = 6/pntnorm,
            mapping = aes(x = -1, y= 0.9, label= paste('p==', round(p,4)  )) ) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))
hetcorplot2

ggsave('./results/SI_figures/heteroscedasticity/hetch_density.pdf', hetcorplot2, units='cm', height = 8, 
       width = 12, useDingbats=F)
ggsave('./results/SI_figures/heteroscedasticity/hetch_density.png', hetcorplot2, units='cm', height = 8, 
       width = 12)

## heterogeneity change of development-divergent genes with ncvTest:
het2Di = het2 %>%
  filter(gene_id%in%names(devdiv)) %>%
  mutate(type = factor(ifelse(gene_id%in%dico, 'DiCo', 'DiDi') ))

# violin plot for heterogeneity changes of DiCo and DiDi genes using ncvTest:
# het2plot = het2Di %>% 
#   group_by(tissue) %>%
#   ggplot(aes(x=tissue, y= chi, fill=type )) +
#   geom_violin(scale='count', trim=T)  +
#   geom_boxplot(width=0.05, position = position_dodge(width = .9), show.legend = F, outlier.size = 0.3) +
#   ylab('ncvTest ChiSquare') +
#   xlab('') + 
#   theme(legend.position = 'right')
# het2plot
# ggsave('./results/ncvTchi.pdf', het2plot, units='cm', height = 8, width = 16, useDingbats=F)
# ggsave('./results/ncvTchi.png', het2plot, units='cm', height = 8, width = 16)


##### ks.test for heterogeneity change  for dico vs didi genes using ncvTest:
annothet2 = het2Di %>% 
  spread(key=type, value = chi) %>%
  group_by(tissue) %>%
  summarise(p = ks.test(DiDi, DiCo)$p.val,
            D= ks.test(DiDi, DiCo)$stat)
# tissue      p      D
# <fct>   <dbl>  <dbl>
#   1 Cortex 0.822  0.0133
# 2 Liver  0.142  0.0242
# 3 Lung   0.0828 0.0266
# 4 Muscle 0.0423 0.0292

# density plots of heterogeneity changes for dico vs didi genes with ncvTest:
hetcorplot3dat = het2Di %>% 
  rename('Pattern' = type)

hetcorplot3 = het2Di %>% 
  rename('Pattern' = type) %>%
  ggplot(aes(x= chi, fill=Pattern, linetype=Pattern, size=Pattern )) +
  facet_wrap(~tissue, scales = 'free_y') +
  geom_density(alpha=0.3) +
  scale_linetype_manual(values=c('DiCo'='solid', 'DiDi'='dashed')) + 
  scale_size_manual(values=c(0.4, 0.4)) +
  ylab('Density') +
  xlab(parse(text='chi^2')) +
  geom_text(data = annothet2, inherit.aes = F, parse=T, hjust=0, size = 6/pntnorm,
            mapping = aes(x = 2, y= c(0.6, 1, 0.4, 0.5), label= paste('p==', round(p,4)  )) ) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))
hetcorplot3

ggsave('./results/SI_figures/heteroscedasticity/chisq_density.pdf', hetcorplot3, units='cm', height = 8, 
       width = 12, useDingbats=F)
ggsave('./results/SI_figures/heteroscedasticity/chisq_density.png', hetcorplot3, units='cm', height = 8, 
       width = 12)

hetplot = ggarrange(hetcorplot2, hetcorplot3, nrow=2, common.legend = T, labels = c('a.','b.'), 
                    vjust = c(0.1,0.1), legend='top', font.label = list(size=8))

ggsave('./results/SI_figures/heteroscedasticity/hetplot.pdf', hetplot, units='cm', height = 14,
       width = 16, useDingbats=F)
ggsave('./results/SI_figures/heteroscedasticity/hetplot.png', hetplot, units='cm', height = 14,
       width = 16, bg='white')

ggsave('./results/figure_supplements/fs2/FS15.pdf', hetplot, units='cm', height = 12,
       width = 15, useDingbats=F)
ggsave('./results/figure_supplements/fs2/FS15.png', hetplot, units='cm', height = 12,
       width = 15, bg='white')

saveRDS(hetcorplot2dat,'results/source_data/f2/fs15_het.rds')
saveRDS(hetcorplot3dat, 'results/source_data/f2/fs15_ncvtest.rds')

# correlation between abs. residual and chiSquare value methods: 
hetcor %>%
  inner_join(het2) %>%
  group_by(tissue) %>%
  summarise(cor = cor(abs(rho), chi, m='s'))
# tissue   cor
# <fct>  <dbl>
#   1 Cortex 0.710
# 2 Liver  0.634
# 3 Lung   0.686
# 4 Muscle 0.700

############

ms = hetcorDi %>% filter(tissue=='Lung')

wilcox.test(ms$rho[ms$type=='DiCo'],ms$rho[ms$type=='DiDi'])

##### wilcox.test:
hetcorDi %>% 
  group_by(tissue) %>%
  summarise(p = wilcox.test(rho~type)$p.val,
            est = wilcox.test(rho~type)$stat)
# tissue       p       est
# <fct>    <dbl>     <dbl>
#   1 Cortex 0.0794  10000788.
# 2 Liver  0.532   10141052.
# 3 Lung   0.00995  9898498.
# 4 Muscle 0.0913  10428368.

##### t.test:
hetcorDi %>% 
  select(-p) %>% 
  ungroup() %>%
  spread(key=type, value = rho) %>%
  group_by(tissue) %>%
  summarise(p = t.test(DiDi, DiCo)$p.val,
            est = t.test(DiDi, DiCo)$stat)
# tissue       p   est
# <fct>    <dbl> <dbl>
#   1 Cortex 0.0442  -2.01
# 2 Liver  0.285   -1.07
# 3 Lung   0.00545 -2.78
# 4 Muscle 0.0806   1.75


#### permutation test for ncvTest chi values:
tis = levels(het2Di$tissue)
obs = c()
perms = list()
for(ts in tis){
  tsX = het2Di %>% filter(tissue==ts )
  Xdi = tsX$chi[tsX$type=='DiDi']
  Xdico = tsX$chi[tsX$type=='DiCo']
  N = length(Xdico)
  pop = c(Xdi, Xdico)
  
  obs[ts] = mean(Xdico)
  permX = sapply(1:1000, function(x){
    mean(sample(pop, size = N))
  })
  perms[[ts]] = permX
}

pvals = sapply(tis, function(x) mean(perms[[x]]>= obs[x]) )
fpr = sapply(tis, function(x) {median(perms[[x]])/ obs[x] }) 

permresult = data.frame(tissue=factor(names(obs)), obs=round(obs,4), p = pvals,
                        row.names = NULL, fpr = round(fpr,2))
permresult

ncvTchipermplot= reshape2::melt(perms) %>% 
  set_names('est', 'tissue') %>% 
  mutate(tissue =  factor(tissue)) %>%
  left_join(data.frame(tissue=names(obs), obs=obs)) %>% 
  ggplot(aes(x=est)) +
  facet_grid(~tissue, scales = 'free') +
  geom_histogram() +
  geom_vline(data = permresult, mapping = aes(xintercept = obs), color='darkred', linetype='dashed') 

ncvTchipermplot
ggsave('./results/ncvTchiperm.pdf', ncvTchipermplot, units='cm', height = 8, width = 16, useDingbats=F)
ggsave('./results/ncvTchiperm.png', ncvTchipermplot, units='cm', height = 8, width = 16)


############ permutation test for het.change(rho) values:
tis = levels(hetcorDi$tissue)
obs = c()
perms = list()
for(ts in tis){
  tsX = hetcorDi %>% filter(tissue==ts )
  Xdi = tsX$rho[tsX$type=='DiDi']
  Xdico = tsX$rho[tsX$type=='DiCo']
  N = length(Xdico)
  pop = c(Xdi, Xdico)
  
  obs[ts] = mean(Xdico)
  permX = sapply(1:10000, function(x){
    mean(sample(pop, size = N, replace = T))
  })
  perms[[ts]] = permX
}
pvals = sapply(tis, function(x) mean(perms[[x]]>= obs[x]) )
fpr = sapply(tis, function(x) {median(perms[[x]])/ obs[x] }) 

permresult = data.frame(tissue=factor(names(obs)), obs=round(obs,4), p = pvals,
                        row.names = NULL, fpr = round(fpr,2))

permresult #(no replacement)
# tissue     obs      p   fpr
# 1 Cortex -0.1306 0.0225  1.06
# 2  Liver  0.0031 0.1443 -0.28
# 3   Lung  0.2516 0.0023  0.96
# 4 Muscle  0.1113 0.9606  1.06

permresult #(with replacement)
# tissue     obs      p   fpr
# 1 Cortex -0.1306 0.0883  1.06
# 2  Liver  0.0031 0.2298 -0.32
# 3   Lung  0.2516 0.0278  0.96
# 4 Muscle  0.1113 0.8866  1.06

hetcorpermplot = reshape2::melt(perms) %>% 
  set_names('est', 'tissue') %>% 
  mutate(tissue =  factor(tissue)) %>%
  left_join(data.frame(tissue=names(obs), obs=obs)) %>% 
  ggplot(aes(x=est)) +
  facet_grid(~tissue, scales = 'free') +
  geom_histogram() +
  geom_vline(data = permresult, mapping = aes(xintercept = obs), color='darkred', linetype='dashed') 

hetcorpermplot
ggsave('./results/hetchperm.pdf', hetcorpermplot, units='cm', height = 8, width = 16, useDingbats=F)
ggsave('./results/hetchperm.png', hetcorpermplot, units='cm', height = 8, width = 16)

save(list=ls(), file = './data/processed/hetsim.rdata')

