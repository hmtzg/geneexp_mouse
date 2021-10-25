library(tidyverse)
library(ggpubr)
library(RColorBrewer)

tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))

celltype_cors = readRDS('./data/other_datasets/scRNA-seq/processed/celltype_pairwise_cors.rds')[['allgenes']]
ct_cors = readRDS('./data/other_datasets/scRNA-seq/processed/celltype_3m_pairwise_cors_allgenes.rds')

dc = readRDS('./data/other_datasets/scRNA-seq/processed/deconvolution_dc_genes.rds')
nondc = readRDS('./data/other_datasets/scRNA-seq/processed/deconvolution_nondc_genes.rds')
age = readRDS('./data/processed/raw/ages.rds')

coeffs = reshape2::melt(list(DiCo=dc,`Non-DiCo`=nondc)) %>%
  set_names(c('cell type','id','proportion','tissue','Geneset')) %>%
  left_join(data.frame(id=names(age), age= age))

######### choose cell types having highest relative contribution to bulk tissues:

top_celltype = coeffs %>%
  group_by(tissue, `cell type`, Geneset) %>%
  summarise(mprop = mean(proportion)) %>%
  group_by(tissue, Geneset) %>%
  top_n(n=1, wt=mprop) %>%
  left_join(coeffs, by=c('tissue','cell type','Geneset')) %>%
  mutate(period= ifelse(age<90,'development','ageing')) %>%
  group_by(tissue,period, Geneset,`cell type`) %>%
  summarise(rho= round(cor.test(age, proportion, m='s')$est,2),
            p.val= round(cor.test(age, proportion, m='s')$p.val,3))

annt = coeffs  %>%
  group_by(tissue, `cell type`, Geneset) %>%
  summarise(mprop = mean(proportion)) %>%
  group_by(tissue, Geneset) %>%
  top_n(n=1, wt=mprop) %>%
  mutate(`cell type` = gsub('bronchial smooth muscle cell','bronchial smooth \nmuscle cell',`cell type`) ) %>%
  mutate(`cell type` = gsub('skeletal muscle satellite cell','skeletal muscle \nsatellite cell',`cell type`) )

#tscol= tissuecol[c(1,3,2,4)]
#names(tscol) = unique(paste(top_celltype$tissue,top_celltype$`cell type`, sep=':\n'))

########### panel a
panel_a_dat = coeffs %>%
  right_join(top_celltype, by=c('tissue','cell type', 'Geneset'))
panel_a = coeffs %>%
  right_join(top_celltype, by=c('tissue','cell type', 'Geneset')) %>%
  #mutate(v1 = paste(tissue,`cell type`,sep=':\n')) %>% 
  #select(-tissue, -`cell type`) %>%
  ggplot(aes(x=age, y=proportion, color=tissue )) +
  facet_grid(tissue~Geneset,scales = 'free') +
  geom_point(alpha = 0.8, size = .7) +
  scale_x_continuous(trans= 'log2') +
  #scale_color_manual(values = tscol, breaks = rev(names(tscol)[c(2,3,4,1)]) ) + 
  scale_color_manual(values = tissuecol) + 
  geom_vline(xintercept = 90, linetype = 'dashed', alpha=0.3) +
  geom_smooth(aes(fill=tissue), method = 'loess', se = T, size = .7, alpha=0.2) +
  xlab('Age in days (in log2 scale)') +
  ylab('Relative Contribution') +
  ggtitle('') +
  guides(color=guide_legend(title = 'Tissue', reverse = T), fill=F) +
  theme_bw(base_size = 12) +
  theme(text = element_text(size=8),
        axis.text.y = element_text(size=4),
        legend.position = 'none',
        legend.key.size = unit(8,'pt'),
        strip.text = element_text(size=8, margin=margin(b=0.5, l=0.5))) +
  geom_text(data=annt,
            mapping = aes(x=c(5,NA,7,NA,14,NA,10,NA), y=rep(c(1.35,4,1.6,1.8), each=2), label = `cell type`), 
            size= 2 ) 
panel_a
ggsave('./results/figure5/Figure_5a.pdf', panel_a, units='cm', height=7, width = 10, useDingbats=F)
ggsave('./results/figure5/Figure_5a.png', panel_a, units='cm', height=7, width = 10)

saveRDS(panel_a_dat,'results/source_data/f5/a.rds')
########################
######################## plot max and min cors change with age
########################

maxcors = celltype_cors %>%
  filter(age==3) %>%
  group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  top_n(n=1, wt=rho) %>% ungroup() %>%
  select(-rho, -`age gr`,-age) %>%
  left_join(celltype_cors, by =c('2nd cell type', '2nd tissue','1st cell type', '1st tissue')  ) %>%
  #select(-`age gr`) %>%
  group_by(`1st tissue`,`1st cell type`) %>% 
  summarise(`rho ch` = round(cor(rho, age, m='s'),2)) %>%
  ungroup() %>%
  arrange(-`rho ch`)

mincors = celltype_cors %>%
  filter(age==3) %>%
  group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  top_n(n=1, wt = -rho) %>% ungroup() %>%
  select(-rho, -`age gr`,-age) %>%
  left_join(celltype_cors, by =c('2nd cell type', '2nd tissue','1st cell type', '1st tissue')  ) %>%
  #select(-`age gr`) %>%
  group_by(`1st tissue`,`1st cell type`) %>% 
  summarise(`rho ch` = round(cor(rho, age, m='s'),2)) %>%
  ungroup() %>%
  arrange(-`rho ch`)


############## panel b
maxminch = maxcors %>% left_join(mincors, by = c('1st tissue','1st cell type')) %>%
  mutate(`1st tissue`= str_to_title(`1st tissue`)) %>%
  reshape2::melt() %>% 
  set_names(c('1st tissue','1st cell type', 'minx','Correlation') ) %>%
  mutate(minx = factor(ifelse(minx=='rho ch.x','Maximum','Minimum'), levels = c('Minimum','Maximum'))) %>%
  mutate(var = 'Tissue Similarity Change')
saveRDS(maxminch, './data/other_datasets/scRNA-seq/processed/max_min_cors_ch.rds')

panel_b2 = maxminch %>%
  ggplot(aes(x= minx, y=Correlation)) +
  #facet_wrap(Correlation~minx) +
  facet_grid(var~minx, scales="free_x") +
  geom_hline(yintercept = 0, linetype ='dashed', alpha=0.5) +
  geom_violin(aes(fill=minx), alpha=0.8) +
  geom_boxplot(width=0.1, alpha = 0.4) +
  scale_fill_manual(values = brewer.pal(10, "Paired")[c(8,10)]) +
  xlab('') +
  ylab('Spearman\'s Correlation coefficient ') +
  guides(fill=guide_legend(title = 'Tissue')) +
  theme_bw(base_size = 12) +
  theme(text = element_text(size=8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.text = element_text(size=8, margin=margin(b=0.5,l=0.5)))
panel_b2

ggsave('./results/figure5/Figure_5b2.pdf', panel_b2, units='cm', height=7, width = 7, useDingbats=F)
ggsave('./results/figure5/Figure_5b2.png', panel_b2, units='cm', height=7, width = 7)

saveRDS(maxminch,'results/source_data/f5/b2.rds')
########### panel b1
celltype_cors = celltype_cors %>%
  mutate(`1st tissue` = str_to_title(`1st tissue`))
  
maxrho = celltype_cors %>%
  filter(age==3) %>%
  group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  top_n(n=1, wt=rho) %>%
  ungroup() %>%
  select(-`age gr`,-age)

maxrho24 = maxrho %>%
  select(-rho) %>% 
  left_join(celltype_cors, by=c('1st cell type','1st tissue','2nd cell type','2nd tissue')) %>% 
  filter(age==24) %>%
  select(-`age gr`,-age)

minrho = celltype_cors %>%
  filter(age==3) %>%
  group_by(`1st tissue`, `1st cell type`,`2nd tissue`) %>%
  top_n(n=1, wt=-rho) %>%
  ungroup() %>%
  select(-`age gr`,-age)

minrho24 = minrho %>%
  select(-rho) %>% 
  left_join(celltype_cors, by=c('1st cell type','1st tissue','2nd cell type','2nd tissue')) %>% 
  filter(age==24) %>%
  select(-`age gr`,-age)

maxmin_density = list(Maximum = list(m3=maxrho, m24=maxrho24), Minimum = list(m3=minrho, m24=minrho24)) %>%
  reshape2::melt() %>%
  rename(Similarity=L1, Correlation=value, age=L2) %>%
  select(-variable) %>%
  mutate(Similarity = factor(Similarity, levels=c('Minimum','Maximum'))) %>%
  mutate(var='Similarity \nDistribution')

saveRDS(maxmin_density, './data/other_datasets/scRNA-seq/processed/max_min_cors_density.rds')

panel_b1 = maxmin_density %>% 
  ggplot(aes(x = Correlation, fill = interaction(Similarity,age), linetype=age, size=age  )) +
  facet_grid(var~Similarity, scales = 'free_x') +
  geom_density( alpha=0.6) +
  scale_fill_manual(values = brewer.pal(10, "Paired")[c(7,9,8,10)]) +
  scale_linetype_manual(values=c('m3'='solid','m24'='dashed') ) +
  scale_size_manual(values=c('m3'=0.05,'m24'=0.1)) +
  xlab('') +
  ylab('Density') +
  theme_bw(base_size = 12) +
  theme(text = element_text(size=8),
        axis.text.y = element_blank(),
        legend.position = 'none',
        axis.title.y = element_text(margin = margin(l=10) ),
        strip.text = element_text(size=8, margin=margin(l=1,b=1,t=1)))
panel_b1

saveRDS(maxmin_density,'results/source_data/f5/b1.rds')

ggsave('./results/figure5/Figure_5b1.pdf', panel_b1, units='cm', height=7, width = 7, useDingbats=F)
ggsave('./results/figure5/Figure_5b1.png', panel_b1, units='cm', height=7, width = 7)

panel_b = ggarrange(panel_b1, panel_b2, nrow=2, heights = c(0.6,1), align = 'v')
panel_b
figure_5 = ggarrange(panel_a, panel_b, ncol=2, labels = c('a.','b.'), font.label = list(size=8))

ggsave('./results/figure5/Figure5.pdf', figure_5, units='cm', height=10, width = 15, useDingbats=F)
ggsave('./results/figure5/Figure5.png', figure_5, units='cm', height=10, width = 15)


# minmaxd = list(Max = maxrho, Min = minrho) %>%
#   reshape2::melt() %>%
#   rename(Similarity=L1, Correlation=value) %>% 
#   mutate(Similarity = as.factor(Similarity)) %>%
#   select(-variable) %>%
#   ggplot(aes(x = Correlation, fill =  Similarity)) +
#   geom_density(aes(group=Similarity), alpha=0.7 ) +
#   facet_wrap(~`1st tissue`, scales = 'free') +
#   xlab('Correlation') +
#   ylab('Density') +
#   theme(text = element_text(size=6),
#       axis.text.y = element_text(size=4),
#       legend.position = 'right',
#       legend.key.size = unit(10,'pt'),
#       strip.text = element_text(size=6, margin=margin(t = 1.5, b=2)))
# 
# ggsave('./results/figure5/figure_5bx.pdf', minmaxd, units='cm', height=7, width = 7, useDingbats=F)
# ggsave('./results/figure5/figure_5bx.png', minmaxd, units='cm', height=7, width = 7)

############### each tissue separate
# tissuecol2 = tissuecol
# names(tissuecol2)[1] = c('Brain')
# p = maxcors %>% left_join(mincors, by = c('1st tissue','1st cell type')) %>%
#   mutate(`1st tissue`= str_to_title(`1st tissue`)) %>%
#   reshape2::melt() %>% 
#   set_names(c('1st tissue','1st cell type', 'minx','Correlation') ) %>%
#   mutate(minx = factor(ifelse(minx=='rho ch.x','Maximum','Minumum'), levels = c('Minumum','Maximum'))) %>%
#   ggplot(aes(x=minx, y= Correlation)) +
#   facet_wrap(~`1st tissue`) +
#   geom_hline(yintercept = 0, linetype ='dashed', alpha=0.5) +
#   geom_violin(aes(fill=`1st tissue`), alpha=0.8) +
#   geom_boxplot(width=0.1, alpha = 0.4) +
#   scale_fill_manual(values = tissuecol2) +
#   xlab('') +
#   guides(fill=guide_legend(title = 'Tissue')) +
#   theme_bw(base_size = 12) +
#   theme(text = element_text(size=t),
#         axis.text.y = element_text(size=4),
#         legend.position = 'none',
#         legend.key.size = unit(2,'pt'),
#         strip.text = element_text(size=8, margin=margin(t = 1.5, b=2)))
# 
# p
# ggsave('./results/figure5/figure_5x.pdf', p, units='cm', height=7, width = 7, useDingbats=F)
# ggsave('./results/figure5/figure_5x.png', p, units='cm', height=7, width = 7)



