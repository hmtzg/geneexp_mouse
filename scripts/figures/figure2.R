library(tidyverse)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(grid)
library(RColorBrewer)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)

tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('Up','Down'))
revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))
periodcol = setNames(c("#FE6100","#648FFF"), c('Development', 'Ageing'))

sample_info = readRDS('./data/processed/tidy/sample_info.rds')
expr = readRDS('./data/processed//tidy/expression.rds')
cov_dat = readRDS('./data/processed/tidy/CoV.rds')
cov_ch = readRDS('./data/processed//tidy/CoV_change.rds')
#cd_ratio_CI = readRDS('data/processed/tidy/cov_ratio_jk_CI.rds')
cd_ratio_pseudo = readRDS('data/processed/raw/cov_ratio_jk_pseudovalues.rds')

######
# expression profile of top gene showing divergence-convergence pattern:
top_divcon_gene_dat = cov_ch %>%
  select(-pval,-FDR) %>%
  spread(key = period, value =`CoV_change`) %>%
  mutate(rev = aging *development) %>%
  top_n(n=1,wt=-(rev)) %>%
  left_join(expr) %>%
  inner_join(sample_info) 

## CoV average for each individual and mean(CoV)-age correlation:
cov_dat_sum = cov_dat %>%
  mutate(ind_id = factor(ind_id)) %>%
  left_join(unique(select(sample_info,-tissue,-sample_id))) %>%
  group_by(ind_id, age) %>%
  summarise(meanCoV = mean(CoV),
            medianCoV = median(CoV)) %>% 
  mutate(period = c('aging','development')[1+(age<=90)])

cov_cordat = group_by(ungroup(cov_dat_sum),period) %>%
  summarise(cor = cor.test(meanCoV, age, method = 's')$est,
            cor.p = cor.test(meanCoV, age, method = 's')$p.val)

cov_cordat_median = group_by(ungroup(cov_dat_sum),period) %>%
  summarise(cor = cor.test(medianCoV, age, method = 's')$est,
            cor.p = cor.test(medianCoV, age, method = 's')$p.val)

# panel a.
szx=6
panel_a = ggplot(cov_dat_sum, aes(x = age, y= meanCoV)) +
  geom_smooth(se=T, color = 'midnightblue', fill = 'lightblue') +
  geom_point(size = 1.5, color='steelblue', alpha= 0.9) +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept = 90, linetype='dashed', color = 'gray25') +
  xlab('Age in days (in log2 scale)') + ylab('Mean CoV') +
  ggtitle('All expressed genes') +
  annotate('text', x=94, y=0.505, label='Ageing', hjust=0, size = szx/pntnorm, fontface='bold') +
  annotate('text', x=94, y=0.500, hjust=0, size = szx/pntnorm,
           label = parse(text = paste0('rho["CoV,age"]==' ,
                                      round(filter(cov_cordat,period=='aging')$cor,2)))) +
  annotate('text', x=94, y=0.496, hjust=0, size = szx/pntnorm, 
           label = parse(text = paste0('p==' ,
                                       (round(filter(cov_cordat,period=='aging')$cor.p,3))))) +
  annotate('text', x=5.5, y=0.485, label='Development', hjust=0, size = szx/pntnorm, fontface='bold') +
  annotate('text', x=6, y=0.480, hjust=0, size = szx/pntnorm,
           label = parse(text = paste0('rho["CoV,age"]==' ,
                                      (round(filter(cov_cordat,period=='development')$cor,2))))) +
  annotate('text', x=6, y=0.476, hjust=0, size = szx/pntnorm,
           label = parse(text = paste0('p==' ,
                                       (round(filter(cov_cordat,period=='development')$cor.p,3))))) +
  theme_bw() +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        plot.title = element_text(size=8))
panel_a
ggsave('./results/figure2/Figure_2a.pdf',panel_a, units = 'cm', width = 8, height = 8, useDingbats = F)
ggsave('./results/figure2/Figure_2a.png',panel_a, units = 'cm', width = 8, height = 8)

saveRDS(cov_dat_sum, 'results/source_data/f2/cov_dat_sum.rds')

# pabel_b:  CoV change of top gene showing divergence-convergence pattern
top_divcon_cov_dat = cov_ch %>%
  select(-pval,-FDR) %>%
  spread(key = period, value = `CoV_change`) %>%
  mutate(rev = aging*development) %>%
  top_n(n=1, wt=-(rev)) %>%
  left_join(cov_dat) %>%
  mutate(ind_id = factor(ind_id)) %>%
  select(-development,-aging,-rev) %>%
  left_join( unique(select(sample_info,-tissue,-sample_id)) ) %>%
  mutate(period = c('aging','development')[1+(age<=90)])

top_divcon_cov_cordat = top_divcon_cov_dat %>% 
  group_by(period) %>%
  summarise(cor = cor.test(CoV, age, method = 's')$est,
            cor.p = cor.test(CoV, age, method = 's')$p.val) %>%
  mutate(cor.p2 = c('p==4*x10^-4', 'p==10^-5'))

panel_b = top_divcon_cov_dat %>%
  ggplot(aes(x = age, y= CoV)) +
  geom_smooth(color='midnightblue', fill = 'lightblue') +
  geom_point(size=1.5, color = 'steelblue', alpha = 0.9) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray25') +
  scale_x_continuous(trans = 'log2') +
  xlab('Age in days (in log2 scale)') + ylab('Coeff. of Variation (CoV)') +
  ggtitle('Cd93') +
  annotate('text', x=94, y=0.70, label='Ageing', hjust=0, size = szx/pntnorm, fontface='bold') +
  annotate('text', x=95, y=0.66, hjust=0, size = szx/pntnorm, 
           label = parse(text = paste0('rho["CoV,age"]==' ,
                                      round(filter(top_divcon_cov_cordat, period=='aging')$cor,2)))) +
  annotate('text', x=95, y=0.62, hjust=0, size = szx/pntnorm,
           label = filter(top_divcon_cov_cordat,period=='aging')$cor.p2, parse=T) +
  annotate('text', x=3, y=0.96, label='Development', hjust=0, size = szx/pntnorm, fontface='bold')  +
  annotate('text', x=3, y=0.91, hjust=0, size = szx/pntnorm,
           label = parse(text = paste0('rho["CoV,age"] ==' ,
                                      round(filter(top_divcon_cov_cordat, period=='development')$cor,2)))) +
  annotate('text', x=3, y=0.87, hjust=0, size = szx/pntnorm,
           label = filter(top_divcon_cov_cordat, period=='development')$cor.p2, parse=T) +
  theme_bw() +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        title = element_text(size=6),
        plot.title = element_text(face='italic'))
panel_b

ggsave('./results/figure2/Figure_2b.pdf',panel_b, units = 'cm', width = 8, height = 8,
       useDingbats = F)
ggsave('./results/figure2/Figure_2b.png',panel_b, units = 'cm', width = 8, height = 8)

saveRDS(top_divcon_cov_dat,'results/source_data/f2/top_divcon_cov_dat.rds')

# panel c.
panel_c = ggplot(top_divcon_gene_dat, aes(x =age, y= expression, color = tissue)) +
  geom_smooth(se=F) +
  geom_point(alpha= 0.5) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  scale_x_continuous(trans = 'log2') +
  scale_color_manual(values = tissuecol) +
  ggtitle('Cd93') +
  guides(color = guide_legend('Tissue', 
                              override.aes = list(size = 2))) +
  # theme(legend.position = 'top',
  #       legend.background = element_rect(fill = 'gray85',color = 'gray25')) +
  theme_bw() +
  xlab('Age in days (in log2 scale)') + ylab('Gene Expression') +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        title = element_text(size=6),
        legend.title = element_text(size=10),
        plot.title = element_text(face = 'italic'))
panel_c
leg = get_legend(panel_c)
ggsave('./results/figure2/Figure_2c.pdf', panel_c, units = 'cm', width = 8, height = 8, useDingbats = F)
ggsave('./results/figure2/Figure_2c.png', panel_c, units = 'cm', width = 8, height = 8)

saveRDS(top_divcon_gene_dat,'results/source_data/f2/top_divcon_gene_dat.rds')

## Conv/Div proportions and ratio:
cd_props = cov_ch %>% 
  filter(FDR < 0.1) %>%
  mutate(change = ifelse(CoV_change>0, "Diverge","Converge")) %>%
  select(gene_id, change,  period) %>%
  group_by(period,change) %>%
  summarise(n=n()) %>%
  mutate(period =  gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels=c("Development", "Ageing")))

# panel d:
panel_d = cd_props %>%
  ggplot(aes(x=change, y=n, fill=change)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~ period, scales="free", ncol=2) +
  scale_fill_manual(values=brewer.pal(3,"Set1")[c(2,1)]) +
  theme_bw() + theme(legend.position = "none") +
  xlab(" ") + ylab("Count") +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6))
panel_d
ggsave('./results/figure2/Figure_2d.pdf',panel_d, units = 'cm', width = 8, height = 8, useDingbats = F)
ggsave('./results/figure2/Figure_2d.png',panel_d, units = 'cm', width = 8, height = 8)

saveRDS(cd_props,'results/source_data/f2/cd_props.rds')

# panel e:
pseudorange = cbind( round(range(log2(cd_ratio_pseudo$dev[-1]), na.rm=T),2),
                     round(range(log2(cd_ratio_pseudo$aging[-1]), na.rm=T),2))
rownames(pseudorange) = c('Lmin', 'Lmax')

cd_ratio = cd_props %>%
  mutate(period =  gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels=c("Development", "Ageing"))) %>%
  group_by(period) %>%
  summarise(Ratio = n[change=="Converge"]/n[change=="Diverge"]) %>%
  mutate(facet="Conv/Div Ratio") %>%
  cbind( t(pseudorange) ) %>%
  mutate(LRatio = log2(Ratio))

panel_e = cd_ratio %>%
  ggplot(aes(x=period, y=LRatio, fill=period)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~ facet) +
  scale_fill_manual(values=c("#FE6100","#648FFF")) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin=Lmin, ymax=Lmax), width=0.1) +
  xlab(" ") + 
  ylab("Log2 Conv/Div Ratio") +
  geom_hline(yintercept = 0, linetype='dashed') +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6))
panel_e
ggsave('./results/figure2/Figure_2e.pdf',panel_e, units = 'cm', width = 8, height = 8, useDingbats = F)
ggsave('./results/figure2/Figure_2e.png',panel_e, units = 'cm', width = 8, height = 8)

saveRDS(cd_ratio,'results/source_data/f2/cd_ratio.rds')
# figure 2
# blank= grid.rect(gp=gpar(col='white'))
# figure_2 = ggarrange(
#   ggarrange(panel_a, panel_b, panel_c, ncol=3, labels = c('a.','b.','c.'), align='h',
#             widths = c(0.8, 0.9, 1.1), font.label = list(size=10) ),
#   ggarrange(panel_d, panel_e, blank, ncol= 3,labels = c('d.','e.',NA), widths = c(1,0.5,0.3),
#             font.label = list(size=10)),
#   heights = c(1,0.8),nrow= 2)
# figure_2

figure_2 = ggarrange(
  ggarrange(panel_a, panel_b, panel_c, ncol=3, labels = c('a.','b.','c.'), align='h',
            widths = c(1, 1, 0.8), font.label = list(size=10), legend = F ),
  ggarrange(panel_d, panel_e, legend.grob = leg, legend = 'right',
            ncol= 3,labels = c('d.','e.',NA), widths = c(2, 1, 0.2), aling='hv',
            font.label = list(size=10)),
  heights = c(1,0.8),nrow= 2)
figure_2

ggsave("results/figure2/Figure_2.pdf", figure_2, units='cm', width = 16, height = 9.5, useDingbats=F)
ggsave("results/figure2/Figure_2.png", figure_2, units='cm', width = 16, height = 9.5)


