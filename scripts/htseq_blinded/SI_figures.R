library(tidyverse)
library(ggpubr)
library(ggrepel)
library(cowplot)
#library(magick)
library(RColorBrewer)
library(corrplot)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('up','down'))
revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))

#############################
sample_info = readRDS('./data/processed/tidy/sample_info.rds')
expr = readRDS('./data/htseq/expr_mat_blinded.rds')
expr_ch = readRDS('./data/htseq/blinded/expression_change.rds')
pca_dat = readRDS('./data/htseq/blinded/pca_data.rds')
# perm_overlaps=readRDS('./data/processed/raw/tissue_gene_overlaps_perm.rds')
# perm_overlaps_fdr=readRDS('./data/processed/raw/tissue_siggene_fdr01_overlaps_perm.rds')
# dev.perm = readRDS("./data/processed/raw/permutations_dev.rds")
# aging.perm = readRDS("./data/processed/raw/permutations_ageing.rds")
cov_dat = readRDS('./data/htseq/blinded/CoV.rds')
cov_dat_wo_cortex = readRDS("./data/htseq/blinded/CoV_wo_cortex.rds")
cov_ch = readRDS('./data/htseq/blinded/CoV_change.rds')
pexpcors  = readRDS('./data/htseq/blinded/pairwise_tissue_expression_cors.rds')
#intra_ts_cov = readRDS('./data/other_datasets/scRNA-seq/processed/intra_tissue_mean_cov.rds')
#deconv = readRDS('./data/other_datasets/scRNA-seq/processed/deconvolutions_combined.rds')

# cov_gsea = readRDS('./data/processed/figures/tidy/CoV_GSEA.rds')
# divcon_gsea = readRDS('./data/processed/figures/tidy/divcon_GSEA.rds')
# revgenes = readRDS('./data/processed/figures/tidy/revgenes.tissues.rds')

# calculate variations explained in all kinds of pca results:
all_raw_var = (pca_dat %>%
                 filter(period == 'all', type == 'raw') %>%
                 select(varExp, PC) %>%
                 unique())$varExp

aging_raw_var= (pca_dat %>%
                  filter(period == 'aging', type == 'raw') %>%
                  select(varExp, PC) %>%
                  unique())$varExp

dev_raw_var = (pca_dat %>%
                 filter(period == 'development', type == 'raw') %>%
                 select(varExp, PC) %>%
                 unique())$varExp

all_notissue_var = (pca_dat %>%
                      filter(period == 'all', type == 'notissue') %>%
                      select(varExp, PC) %>%
                      unique())$varExp

aging_notissue_var = (pca_dat %>%
                        filter(period == 'aging', type == 'notissue') %>%
                        select(varExp, PC) %>%
                        unique())$varExp

dev_notissue_var = (pca_dat %>%
                      filter(period == 'development', type == 'notissue') %>%
                      select(varExp, PC) %>%
                      unique())$varExp

############
############
############  Figure S2, PCA with all samples and tissue effect removed -----------------------
############

all_nt_pca12 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'notissue') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2') +
  coord_fixed(ratio = all_notissue_var[2]/all_notissue_var[1], clip = 'off') +
  xlab(paste('PC1 (', round(all_notissue_var[1]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(all_notissue_var[2]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) +
  theme(legend.position = c(0.7,0.9),
        legend.direction = c('vertical'),
        legend.box = 'horizontal',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) 

# ggsave('./results/SI_figures/SI_panels/Figure_S2a.pdf',all_nt_pca12, units = 'cm', width = 8, height = 8,
#        useDingbats =F)
# ggsave('./results/SI_figures/SI_panels/Figure_S2a.png',all_nt_pca12, units = 'cm', width = 8, height = 8)

all_nt_pc1age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'notissue') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC1, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age in days (in log2)')

all_nt_pc2age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'notissue') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC2, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age in days (in log2)')

all_nt_pc12 = ggarrange(all_nt_pca12, ggarrange(all_nt_pc1age, all_nt_pc2age, ncol = 1, nrow= 2), ncol = 2,
                        nrow= 1, widths = c(1,0.5), common.legend = F, align = 'v')

ggsave('./results/htseq/blinded/Figure_S2.pdf',all_nt_pc12, units = 'cm', width = 16, height = 9,
       useDingbats = F)
ggsave('./results/htseq/blinded/Figure_S2.png',all_nt_pc12, units = 'cm', width = 16, height = 9)

############
############
############  Figure S3, PCA with development samples only -----------------------
############

# panel a
dev_raw_pca12 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2',breaks = c(2,8,30,60)) +
  coord_fixed(ratio = dev_raw_var[2]/dev_raw_var[1], clip = 'off') +
  xlab(paste('PC1 (', round(dev_raw_var[1]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(dev_raw_var[2]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) +
  theme(legend.position = c(0.06,0.65),
        legend.justification=c(0,0),
        legend.direction = 'vertical',
        legend.box = 'horizontal',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) 

# ggsave('./results/SI_figures/SI_panels/Figure_S3a.pdf',dev_raw_pca12, units = 'cm', width = 8, height = 5,
#        useDingbats =F)
# ggsave('./results/SI_figures/SI_panels/Figure_S3a.png',dev_raw_pca12, units = 'cm', width = 8, height = 5)

# panel b
dev_raw_pca34 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2',breaks = c(2,8,30,60)) +
  coord_fixed(ratio = dev_raw_var[4]/dev_raw_var[3], clip = 'off') +
  xlab(paste('PC3 (', round(dev_raw_var[3]*100),'%)',sep='')) +
  ylab(paste('PC4 (', round(dev_raw_var[4]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'),
         size = guide_legend('Age')) +
  theme(legend.position = c(0.15,0.05),
        legend.justification=c(0,0),
        legend.direction = 'horizontal',
        legend.box = 'vertical',
        legend.background = element_rect(fill = 'gray85',color = 'gray25'))

# ggsave('./results/SI_figures/SI_panels/Figure_S3b.pdf',dev_raw_pca34, units = 'cm', width = 8, height = 5,
#        useDingbats =F)
# ggsave('./results/SI_figures/SI_panels/Figure_S3b.png',dev_raw_pca34, units = 'cm', width = 8, height = 5)

dev_raw_pc1age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC1, color = tissue)) +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age in days (in log2)')

dev_raw_pc2age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC2, color = tissue)) +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age in days (in log2)')

dev_raw_pc3age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC3, color = tissue)) +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age in days (in log2)')

dev_raw_pc4age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC4, color = tissue)) +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age in days (in log2)')

dev_raw_pc1234 = ggarrange(ggarrange(dev_raw_pca12, dev_raw_pca34 + theme(legend.position = 'none'),
                                     ncol = 2, align = 'h', labels = c('a.', 'b.'), vjust = 2.7),
                           ggarrange(dev_raw_pc1age, dev_raw_pc2age,dev_raw_pc3age, dev_raw_pc4age,
                                     align = 'hv', ncol =4), nrow=2, heights = c(1, 0.5), labels = c(NA, 'c.'))

ggsave('./results/htseq/blinded/Figure_S3.pdf', dev_raw_pc1234, units = 'cm', width = 16, height = 12,
       useDingbats = F)
ggsave('./results/htseq/blinded/Figure_S3.png', dev_raw_pc1234, units = 'cm', width = 16, height = 12)

############
############
############  Figure S4, permutation test for shared up/down genes across tissues w/o sig. cutoff ---------------------
############

# obs_overlap = expr_ch %>%
#   filter(`Expression Change` !=0 ) %>% # drop expression - age rho = 0 values
#   mutate(direction = `Expression Change` > 0) %>%
#   mutate(direction = ifelse( direction == TRUE, 'Up', 'Down')) %>%
#   mutate(period = gsub("aging", "ageing", period),
#          period = str_to_title(period),
#          period = factor(period, levels = c('Development','Ageing'))) %>%
#   group_by(gene_id,period, direction) %>%
#   summarise(`N Tissue` = length(tissue)) %>%
#   group_by(`N Tissue`, period, direction) %>% 
#   summarise(Obs = n()) %>%
#   ungroup() %>%
#   slice(-c(1:4)) %>%
#   mutate(`N Tissue` = paste(`N Tissue`, c('Tissues')) )
# 
# overlap_test = perm_overlaps %>%
#   left_join(obs_overlap) %>%
#   group_by(direction, period, `N Tissue`) %>%
#   summarise(p_val = mean(`N Overlap` >= unique(Obs)),
#             FPR =  round(median(`N Overlap`)/ unique(Obs),2)) %>%
#   left_join(obs_overlap)
# 
# sizex = 2.5
# plot_overlaps = perm_overlaps %>%
#   ggplot(aes(x=`N Overlap`)) +
#   facet_grid(period+direction ~  `N Tissue`, scales = 'free') +
#   geom_histogram() +
#   xlab('Number of Overlap Genes') +
#   ylab('Frequency')  +
#   geom_vline(data = obs_overlap, mapping = aes(xintercept= Obs), linetype ='dashed', color = 'darkred' ) +
#   geom_text(data= overlap_test, size=  sizex,
#             mapping = aes(x = rep(c(4500, 3900, 1900),4) , y=110, label = paste('Obs =', Obs) )) +
#   geom_text(data= overlap_test,size=  sizex,
#             mapping = aes(x = rep(c(4500, 3900, 1900),4) , y=100, label = paste('FPR =', FPR) )) +
#   geom_text(data= overlap_test,size=  sizex,
#             mapping = aes(x = rep(c(4500, 3900, 1900),4) , y=90, label = paste('p.val =', p_val) )) +
#   theme_bw() +
#   theme(axis.text = element_text(size =6))
# 
# ggsave('./results/SI_figures/Figure_S4.pdf', plot_overlaps, width = 16, height = 15, units='cm', useDingbats = F )
# ggsave('./results/SI_figures/Figure_S4.png', plot_overlaps, width = 16, height = 15, units='cm' )  

############
############
############  Figure S5, # of sig. overlap across tissues, dev-ageing magnitude comparison -----------------------
############

# significant ones:
# star = data.frame(direction = rep(c('down','up'),3),
#                   y_pos = c(461, 285, 139, 46, 34, 24),
#                   period= rep(c('Development','Ageing'),c(4,2)),
#                   n= c(3,3,4,4,2,2))
# 
# siggenes_overlap = expr_ch %>%
#   filter(FDR < 0.1) %>%
#   mutate(direction = `Expression Change` > 0) %>%
#   mutate(direction = ifelse( direction == TRUE, 'up', 'down')) %>%
#   mutate(period = str_to_title(period)) %>%
#   mutate(period = gsub('Aging','Ageing',period)) %>%
#   mutate(period = factor(period, levels = c('Development','Ageing'))) %>%
#   group_by(gene_id,period, direction) %>%
#   summarise(n = length(tissue)) %>%
#   group_by(n, period, direction) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
#   slice(-c(1:4)) %>%
#   ggplot(aes(x = n, y = count, fill = direction)) +
#   facet_wrap(~period) +
#   geom_bar(stat='identity', position = 'dodge') +
#   scale_fill_manual(values = regcol, drop =F) +
#   scale_y_continuous(trans = 'log10') +
#   geom_text(aes(label = count), color = 'gray15', position = position_dodge(width = 1), angle = 90,
#             hjust = 1.1, size = 6/pntnorm) +
#   guides(fill = guide_legend('Direction Of Expression Change',
#                              override.aes = list(size=2))) +
#   geom_text(data=star, aes(y=y_pos), label = '*', position=position_dodge(width = 1)) +
#   theme(legend.position = 'bottom',
#         legend.background = element_rect(fill = 'gray85', color = 'gray25')) +
#   xlab(NULL) + ylab(NULL) 
# 
# ggsave('./results/SI_figures/SI_panels/Figure_S5a.pdf', siggenes_overlap, width = 15, height = 8, units='cm',
#        useDingbats = F )
# ggsave('./results/SI_figures/SI_panels/Figure_S5a.png', siggenes_overlap, width = 15, height = 8, units='cm' )  

# magnitude change comparison between development and ageing with paired wilcox test:
# expr_ch %>%
#   dplyr::select(-p, -FDR) %>%
#   spread(key = period, value = `Expression Change`) %>%
#   group_by(tissue, gene_id) %>%
#   summarise(diff =  abs(development) - abs(aging)) %>%
#   summarise(p_val = wilcox.test(diff, mu = 0)$p.val) %>%
#   mutate(BH = p.adjust(p_val, method='BH'))
# 
# label_b = data.frame(tissue = c('Cortex','Liver','Lung','Muscle'),
#                      diff = c(rep(1,4)))
# 
# magnitude_comparison = expr_ch %>%
#   dplyr::select(-p, -FDR) %>%
#   spread(key = period, value = `Expression Change`) %>%
#   group_by(tissue, gene_id) %>%
#   summarise(diff =  abs(development) - abs(aging)) %>%
#   ggplot(aes(x=tissue, y= diff,  fill= tissue)) +
#   geom_boxplot() +
#   scale_fill_manual(values = tissuecol) +
#   geom_hline(yintercept = 0, linetype = 'dashed', col = 'darkred') +
#   xlab('') +
#   ylab( bquote('Abs('*rho[dev]~')-Abs('*rho[ageing]~')' ) ) +
#   theme(legend.position = 'none') +
#   geom_text(data = label_b, label = "*")
# 
# ggsave('./results/SI_figures/SI_panels/Figure_S5b.pdf',magnitude_comparison, width = 15, height = 8, units='cm',
#        useDingbats = F )
# ggsave('./results/SI_figures/SI_panels/Figure_S5b.png', magnitude_comparison, width = 15, height = 8, units='cm' )  
# 
# fig_s5 = ggarrange(siggenes_overlap, magnitude_comparison, ncol= 2, labels = c('a.','b.'), widths = c(1, 0.7),
#                    hjust = -0.1)
# 
# ggsave('./results/SI_figures/Figure_S5.pdf', fig_s5, width = 15, height = 8, units='cm', useDingbats = F )
# ggsave('./results/SI_figures/Figure_S5.png', fig_s5, width = 15, height = 8, units='cm' )  

############
############
############  Figure S6, permutation test for sig. overlap genes across tissues -----------------------
############

# obs_overlap_fdr = expr_ch %>%
#   #filter(`Expression Change` !=0 ) %>% # drop expression - age rho = 0 values
#   filter(FDR < 0.1) %>%
#   mutate(direction = `Expression Change` > 0) %>%
#   mutate(direction = ifelse( direction == TRUE, 'Up', 'Down')) %>%
#   mutate(period = gsub("aging", "ageing", period),
#          period = str_to_title(period),
#          period = factor(period, levels = c('Development','Ageing'))) %>%
#   group_by(gene_id,period, direction) %>%
#   summarise(`N Tissue` = length(tissue)) %>%
#   group_by(`N Tissue`, period, direction) %>% 
#   summarise(Obs = n()) %>%
#   ungroup() %>%
#   slice(-c(1:4))  %>%
#   mutate(`N Tissue` = paste(`N Tissue`, c('Tissues')) )
# 
# overlap_fdr_test = perm_overlaps_fdr %>%
#   mutate(`N Tissue` = paste(`N Tissue`, c('Tissues')) ) %>% 
#   right_join(obs_overlap_fdr, by = c('N Tissue', 'period', 'direction')) %>% 
#   group_by(direction, period, `N Tissue`) %>%
#   summarise(p_val = mean(`N Overlap` >= unique(Obs)),
#             FPR =  round(median(`N Overlap`)/ unique(Obs),2)) %>%
#   ungroup() %>%
#   mutate(p_val = ifelse(p_val == 0, '< 0.001', p_val)) %>%
#   left_join(obs_overlap_fdr)
# 
# sizex = 2
# fig_s6_a = perm_overlaps_fdr %>% 
#   mutate(`N Tissue` = paste(`N Tissue`, c('Tissues')) ) %>% 
#   right_join(obs_overlap_fdr, by = c('N Tissue', 'period', 'direction')) %>%
#   filter(period == 'Development') %>%
#   ggplot(aes(x=`N Overlap`)) +
#   facet_grid(direction ~  `N Tissue`, scales = 'free') +
#   geom_histogram() +
#   xlab('Number of Overlap Genes') +
#   ylab('Frequency')  +
#   geom_vline(data = filter(overlap_fdr_test, period == 'Development'),
#              mapping = aes(xintercept= Obs), linetype ='dashed', color = 'darkred' ) +
#   geom_text(data= filter(overlap_fdr_test, period == 'Development'), size = sizex,
#             mapping = aes(x = c(1350, 340, 90, 1320,390,80) , y=200, label = paste('Obs =', Obs) )) +
#   geom_text(data= filter(overlap_fdr_test, period == 'Development'), size=  sizex,
#             mapping = aes(x = c(1350, 340, 90, 1320,390,80) , y=180, label = paste('FPR =', FPR) )) +
#   geom_text(data= filter(overlap_fdr_test, period == 'Development'), size=  sizex,
#             mapping = aes(x = c(1350, 340, 90, 1320,390,80) , y=160, label = paste('p.val =', p_val) )) +
#   theme_bw() +
#   theme(axis.text = element_text(size =6))
# 
# fig_s6_b = perm_overlaps_fdr %>% 
#   mutate(`N Tissue` = paste(`N Tissue`, c('Tissues')) ) %>% 
#   right_join(obs_overlap_fdr, by = c('N Tissue', 'period', 'direction')) %>%
#   filter(period == 'Ageing' & `N Tissue`=='2 Tissues') %>% 
#   mutate(`N Tissue`= factor(`N Tissue`, levels = '2 Tissues')) %>% 
#   ggplot(aes(x=`N Overlap`)) +
#   #facet_wrap(~direction, scales = 'free',nrow = 2, strip.position = 'right') +
#   facet_grid(direction ~  `N Tissue`, scales = 'free') +
#   geom_histogram() +
#   xlab('Number of Overlap Genes') +
#   ylab('Frequency')  +
#   geom_vline(data = filter(overlap_fdr_test, period == 'Ageing'),
#              mapping = aes(xintercept= Obs), linetype ='dashed', color = 'darkred' ) +
#   geom_text(data= filter(overlap_fdr_test, period == 'Ageing'), size = sizex,
#             mapping = aes(x = c(41, 30), y=100, label = paste('Obs =', Obs) )) +
#   geom_text(data= filter(overlap_fdr_test, period == 'Ageing'), size=  sizex,
#             mapping = aes(x = c(41, 31) , y=88, label = paste('FPR =', FPR) )) +
#   geom_text(data= filter(overlap_fdr_test, period == 'Ageing'), size=  sizex,
#             mapping = aes(x = c(42, 32) , y=76, label = paste('p.val =', p_val) )) +
#   theme_bw() +
#   theme(axis.text = element_text(size =6))
# 
# fig_s6 = ggarrange(fig_s6_a,fig_s6_b, ncol = 2, widths = c(1, 0.5), labels = c('a.','b.'))
# 
# ggsave('./results/SI_figures/Figure_S6.pdf', fig_s6 , width = 18, height = 10, units='cm', useDingbats = F )
# ggsave('./results/SI_figures/Figure_S6.png', fig_s6 , width = 18, height = 10, units='cm' )  

############
############
############  Figure S7, correlation plot: tissue similarity of expression changes abs(rho)>0.6 -----------------------
############

cors.06 = expr_ch %>%
  filter(abs(`Expression Change`)> 0.6 ) %>%
  select(-p, -FDR) %>%
  mutate(tissue = factor(tissue, levels= c('Cortex','Liver', 'Lung','Muscle'))) %>%
  pivot_wider(names_from = c(tissue, period), values_from=`Expression Change`) %>% 
  select(-gene_id) %>% 
  as.matrix() %>%
  cor(method='s', use ='pairwise.complete.obs')
colnames(cors.06) = sub('\\_.*','',colnames(cors.06))
rownames(cors.06) = sub('\\_.*','',rownames(cors.06))
cors.06 = cors.06[c(1,2,3,4,7,6,5,8),c(1,2,3,4,7,6,5,8)]
periodcode = c(Development = "#FE6100", Ageing ="#648FFF")
diag(cors.06) = 0

#pdf("./results/SI_figures/Fig_S7.pdf")
png("./results/htseq/SI/Fig_S7.png")
corrplot(cors.06, order = "ori", tl.pos = "lt", diag=T, tl.col = rep(periodcode, each=4), font =2,
         tl.cex = 1.2, method="square", outline=T, type="upper",
         col = colorRampPalette(rev(brewer.pal(n=11, name='RdBu')[c(2:4,6,9:11)]))(200) )
corrplot(cors.06, order="ori",tl.pos = "n", diag=F, tl.col = rep(periodcode, each=4), font=2,
         tl.cex = 1.2, method="shade", outline=T, addCoef.col = "black", add=T, type="lower",
         col = colorRampPalette(rev(brewer.pal(n=11, name='RdBu')[c(2:3,6,9:11)]))(101) )
dev.off()

############
############
############  Figure S8, permutation test for reversals in each tissue: -----------------------
############ (plot in 06.reversal_analysis.R script)

############
############
############  Figure S9, permutation test for shared reversals among tissue: -----------------------
############ (plot in 06.reversal_analysis.R script)

############
############
############ Figure S10, CoV change with 16 samples excluding cortex tissue -----------------------
############ 

cov_sum_wo_cortex = cov_dat_wo_cortex %>%
  mutate(ind_id = factor(ind_id)) %>%
  left_join(unique(select(sample_info,-tissue,-sample_id))) %>%
  group_by(ind_id, age) %>%
  summarise(meanCoV = mean(CoV),
            medianCoV = median(CoV)) %>% 
  mutate(period = c('aging','development')[1+(age<=90)])

cov_sumch1_wo = cov_sum_wo_cortex %>%
  ungroup() %>%
  group_by(period) %>%
  summarise(cor = cor.test(meanCoV, age, method = 's')$est,
            cor.p = cor.test(meanCoV, age, method = 's')$p.val)

cov_sumch2_wo = cov_sum_wo_cortex %>%
  ungroup() %>%
  group_by(period) %>%
  summarise(cor = cor.test(medianCoV, age, method = 's')$est,
            cor.p = cor.test(medianCoV, age, method = 's')$p.val)

## CoV average for each individual and mean(CoV)-age correlation (without cortex):
cov_mean_wo_cortex = ggplot(cov_sum_wo_cortex, aes(x = age, y= meanCoV)) +
  geom_point(size=1.5, color="steelblue", alpha=0.9) +
  geom_smooth(se=T,color = 'midnightblue', fill='lightblue') +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept = 90, linetype='dashed',color = 'gray35') +
  xlab('Age in days (in log2 scale)') + 
  ylab('Mean CoV') +
  annotate('text',x=95,y=0.13,label='Aging',hjust=0,size = 8/pntnorm, fontface='bold') +
  annotate('text',x = 95, y=0.125,
           label = parse(text = paste('rho["CoV,age"] ==' ,
                                      (round(filter(cov_sumch1_wo, period=='aging')$cor,3)))),
           hjust=0,size = 8/pntnorm) +
  annotate('text',x = 95, y=0.120,
           label = parse(text = paste0('p ==' ,(round(filter(cov_sumch1_wo, period=='aging')$cor.p,3)))),
           hjust=0,size = 8/pntnorm) +
  annotate('text',x=9,y=0.11,label='Development',hjust=0,size = 8/pntnorm, fontface='bold') +
  annotate('text',x = 9, y=0.105, 
           label = parse(text = paste('rho["CoV,age"] ==' ,
                                      (round(filter(cov_sumch1_wo, period=='development')$cor,3)))),
           hjust=0,size = 8/pntnorm) +
  annotate('text',x = 9, y=0.10, 
           label = parse(text = paste0('p ==' ,
                                       (round(filter(cov_sumch1_wo, period=='development')$cor.p,2)))),
           hjust=0,size = 8/pntnorm)

cov_median_wo_cortex = ggplot(cov_sum_wo_cortex, aes(x = age, y= medianCoV)) +
  geom_point(size=1.5, color="steelblue", alpha=0.9) +
  geom_smooth(se=T,color = 'midnightblue', fill='lightblue') +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept = 90, linetype='dashed',color = 'gray35') +
  xlab('Age in days (in log2 scale)') +
  ylab('Median CoV') +
  annotate('text',x=95,y=0.09,label='Aging',hjust=0,size = 8/pntnorm, fontface='bold') +
  annotate('text',x = 95, y=0.088, 
           label = parse(text = paste('rho["CoV,age"] ==' ,
                                      (round(filter(cov_sumch2_wo, period=='aging')$cor,3)))),
           hjust=0,size = 8/pntnorm) +
  annotate('text',x = 95, y=0.086, 
           label = parse(text = paste0('p ==' ,(round(filter(cov_sumch2_wo, period=='aging')$cor.p,2)))),
           hjust=0,size = 8/pntnorm) +
  annotate('text',x=10,y=0.07,label='Development',hjust=0,size = 8/pntnorm, fontface='bold') +
  annotate('text',x=10, y=0.068, 
           label = parse(text = paste('rho["CoV,age"] ==' ,
                                      (round(filter(cov_sumch2_wo, period=='development')$cor,3)))),
           hjust=0,size = 8/pntnorm) +
  annotate('text',x=10, y=0.066, 
           label = parse(text = paste0('p ==' ,(round(filter(cov_sumch2_wo, period=='development')$cor.p,3)))),
           hjust=0,size = 8/pntnorm)

fig_s10 = ggarrange(cov_mean_wo_cortex , cov_median_wo_cortex, labels = c('a.','b.'),
                    ncol = 2, nrow = 1, align = 'hv', hjust = c(-0.2,-0.2))

ggsave('./results/htseq/SI/Figure_S10.pdf',fig_s10, units = 'cm', width = 16, height = 8, useDingbats = F)
ggsave('./results/htseq/SI/Figure_S10.png',fig_s10, units = 'cm', width = 16, height = 8)

############
############
############ Figure S11, median CoV change with all samples -----------------------
############ 

cov_sum = cov_dat %>%
  mutate(ind_id = factor(ind_id)) %>%
  left_join(unique(select(sample_info,-tissue,-sample_id))) %>%
  group_by(ind_id, age) %>%
  summarise(meanCoV = mean(CoV),
            medianCoV = median(CoV)) %>% 
  mutate(period = c('aging','development')[1+(age<=90)])

cov_sumch = cov_sum %>%
  ungroup() %>%
  group_by(period) %>%
  summarise(cor = cor.test(medianCoV, age, method = 's')$est,
            cor.p = cor.test(medianCoV, age, method = 's')$p.val)


cov_median = ggplot(cov_sum, aes(x = age, y= medianCoV)) +
  geom_point(size=1.5, color="steelblue", alpha=0.9) +
  geom_smooth(se=T,color = 'midnightblue', fill='lightblue') +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept = 90, linetype='dashed',color = 'gray35') +
  xlab('Age in days (in log2 scale)') +
  ylab('Median CoV') +
  annotate('text',x=95,y=0.115,label='Aging',hjust=0,size = 8/pntnorm, fontface='bold') +
  annotate('text',x = 95, y=0.113, 
           label = parse(text = paste('rho["CoV,age"] ==' ,(round(filter(cov_sumch, period=='aging')$cor,3)))),
           hjust=0,size = 8/pntnorm) +
  annotate('text',x = 95, y=0.111, 
           label = parse(text = paste0('p ==' ,(round(filter(cov_sumch, period=='aging')$cor.p,2)))),
           hjust=0,size = 8/pntnorm) +
  annotate('text',x=8,y=0.1,label='Development',hjust=0,size = 8/pntnorm, fontface='bold') +
  annotate('text',x = 8, y=0.098, 
           label = parse(text = paste('rho["CoV,age"] ==' ,
                                      (round(filter(cov_sumch, period=='development')$cor,3)))),
           hjust=0,size = 8/pntnorm) +
  annotate('text',x = 8, y=0.096, 
           label = parse(text = paste0('p ==' ,(round(filter(cov_sumch, period=='development')$cor.p,2)))),
           hjust=0,size = 8/pntnorm)

ggsave('./results/htseq/SI/Figure_S11.pdf',cov_median, units = 'cm', width = 8, height = 8, useDingbats = F)
ggsave('./results/htseq/SI/Figure_S11.png',cov_median, units = 'cm', width = 8, height = 8)

############
############
############ Figure S12, CoV change direction without significance cutoff -----------------------
############ 

cd_count = cov_ch %>% 
  filter(`CoV_change` !=0 ) %>% 
  #filter(FDR < 0.1) %>%
  mutate(change = ifelse(CoV_change>0, "Diverge","Converge")) %>%
  select(gene_id, change,  period) %>%
  group_by(period,change) %>%
  summarise(n=n()) %>%
  mutate(period =  gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels=c("Development", "Ageing")))

cd_count_p = cd_count %>%
  ggplot(aes(x=change, y=n, fill=change)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~ period, ncol=2) +
  scale_fill_manual(values=brewer.pal(3,"Set1")[c(2,1)]) +
  theme_bw() + theme(legend.position = "none") +
  xlab("") + ylab("Count") +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        strip.text = element_text(size=6))

ggsave("results/htseq/SI/Figure_S12.pdf", cd_count_p, units='cm', width = 8, height = 7, useDingbats=F)
ggsave("results/htseq/SI/Figure_S12.png", cd_count_p, units='cm', width = 8, height = 7)

############
############
############ Figure S13, age-related change in pairwise expression correlations among tissues -----------------------
############ 

annottext = pexpcors %>%
  mutate(period = ifelse(uage < 90,'development','aging')) %>%
  group_by(pair, period) %>%
  summarise(p.val = round(cor.test(uage, rho, m='s')$p.val,2),
            rho = round(cor.test(uage, rho, m='s')$est,3)) %>%
  pivot_longer(cols=c(p.val, rho), names_to='stat')

szx = 7
pwisecors = ggplot(pexpcors) +
  aes(x = uage, y = rho) +
  facet_wrap(~pair,scales = 'free', ncol=2) +
  geom_point(size=1.5, color="steelblue", alpha=0.9) +
  geom_smooth(se=T,color = 'midnightblue', fill='lightblue') +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray30') +
  ylab('Sperman correlation coefficient') + xlab('Age in days (in log2 scale)') +
  geom_text(data = filter(annottext, period =='development' & stat == 'rho'), vjust = 2, 
            mapping = aes(x = 26, y= Inf, label = paste('rho[dev]==',value) ), parse =  T, size = szx/pntnorm) +
  geom_text(data = filter(annottext, period =='development' & stat == 'p.val'), vjust = 4,
            mapping = aes(x = 28, y= Inf, label = paste('p==',value) ), parse =  T,size = szx/pntnorm) +
  geom_text(data = filter(annottext, period =='aging' & stat == 'rho'), vjust = 2,
            mapping = aes(x = 330, y= Inf, label = paste('rho[ageing]==',value) ), parse =  T,size = szx/pntnorm) +
  geom_text(data = filter(annottext, period =='aging' & stat == 'p.val'), vjust = 4,
            mapping = aes(x = 330, y= Inf, label = paste('p==',value) ), parse =  T,size = szx/pntnorm) 

ggsave('./results/htseq/SI/Figure_S13.pdf', pwisecors, units = 'cm', width = 10, height = 10, useDingbats =F)
ggsave('./results/htseq/SI/Figure_S13.png', pwisecors, units = 'cm', width = 10, height = 10)

############
############
############ Figure S14, mean or median expression correlation changes with age -----------------------
############ 

meancors = pexpcors %>% 
  group_by(id, uage) %>%
  summarise(mean = mean(rho),
            median = median(rho)) %>%
  gather(key='method', value = 'rho', mean, median) %>%
  mutate(period = factor(c('development','ageing')[1 + (uage > 90)], levels = c('development', 'ageing')))

annotm = meancors %>%
  group_by(method, period) %>%
  summarise(p.val = round(cor.test(uage, rho, m='s')$p.val,3),
            rho = round(cor.test(uage, rho, m='s')$est,2)) %>%
  pivot_longer(cols=c(p.val, rho), names_to='stat')

szx = 8
mcorsplot = meancors %>%
  ggplot(aes(x=uage,y=rho)) +
  facet_wrap(~method, strip.position = 'left',
             labeller = as_labeller(c(mean = 'Mean P.wise Expr. rho',
                                      median = 'Median P.wise Expr. rho'))) +
  geom_point(size=1.5, color="steelblue", alpha=0.9) +
  geom_smooth(se=T,color = 'midnightblue', fill='lightblue') +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  xlab('Age in days (in log2 scale)') + ylab(NULL) +
  geom_text(data = filter(annotm, stat == 'rho' & period == 'development'),
            mapping = aes(x = 20, y = 0.67, label = paste('rho[dev]==', value)), parse = T, size=szx/pntnorm) +
  geom_text(data = filter(annotm, stat == 'p.val' & period == 'development'),
            mapping = aes(x = 20, y = 0.66, label = paste('p==', value)), parse = T, size=szx/pntnorm) +
  geom_text(data = filter(annotm, stat == 'rho' & period == 'ageing'),
            mapping = aes(x = 270, y = 0.67, label = paste('rho[aging]==', value)), parse = T, size=szx/pntnorm) +
  geom_text(data = filter(annotm, stat == 'p.val' & period == 'ageing'),
            mapping = aes(x = 270, y = .66, label = paste('p==', value)), parse = T, size=szx/pntnorm) +
  theme(strip.background = element_blank(),
        strip.placement = 'outside', strip.text = element_text(size = 8)) 

#### scale correlations and take mean
#pcors.sc = sapply(pexpcors, function(x) scale(x))
psccors = pexpcors %>% group_by(pair) %>% mutate(rho = scale(rho)[,1])

# psccors = reshape2::melt(pcors.sc) %>%
#   select(-Var2) %>%
#   set_names('id' ,'rho', 'pair') %>%
#   mutate(id = factor(id)) %>%
#   left_join(data.frame(uage, id = names(uage)) ) 

meansccors = psccors %>% 
  group_by(id, uage) %>%
  summarise(mean = mean(rho),
            median = median(rho)) %>%
  gather(key='method', value = 'rho', mean, median) %>%
  mutate(period = factor(c('development','ageing')[1 + (uage > 90)], levels = c('development', 'ageing')))

annotscm = meansccors %>%
  group_by(method, period) %>%
  summarise(p.val = round(cor.test(uage, rho, m='s')$p.val,3),
            rho = round(cor.test(uage, rho, m='s')$est,2)) %>%
  pivot_longer(cols=c(p.val, rho), names_to='stat')

mcorsplotsc = meansccors %>%
  ggplot(aes(x=uage,y=rho)) +
  facet_wrap(~method, strip.position = 'left',
             labeller = as_labeller(c(mean = 'Mean  Scaled P.wise Expr. rho',
                                      median = 'Median Scaled P.wise Expr. rho'))) +
  geom_point(size=1.5, color="steelblue", alpha=0.9) +
  geom_smooth(se=T,color = 'midnightblue', fill='lightblue') +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_smooth(se = T, color ='midnightblue') +
  xlab('Age in days (in log2 scale)') + ylab(NULL) +
  geom_text(data = filter(annotscm, stat == 'rho' & period =='development'),
            mapping = aes(x = 30, y = 1.5, label = paste('rho[dev]==', value)), parse = T, size=szx/pntnorm) +
  geom_text(data = filter(annotscm, stat == 'p.val' & period =='development'),
            mapping = aes(x = 30, y = 1.3, label = paste('p==', value)), parse = T, size=szx/pntnorm) +
  geom_text(data = filter(annotscm, stat == 'rho' & period =='ageing'),
            mapping = aes(x = 270, y = 1, label = paste('rho[aging]==', value)), parse = T, size=szx/pntnorm) +
  geom_text(data = filter(annotscm, stat == 'p.val' & period =='ageing'),
            mapping = aes(x = 270, y = .7, label = paste('p==', value)), parse = T, size=szx/pntnorm) +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',strip.text = element_text(size = 8)) 

# ggsave('./results/SI_figures/SI_panels/Figure_S14b.pdf', mcorsplotsc, units = 'cm', width = 10, height = 6,
#        useDingbats =F)
# ggsave('./results/SI_figures/SI_panels/Figure_S14b.png', mcorsplotsc, units = 'cm', width = 10, height = 6)

figure_S14  = ggarrange(mcorsplot, mcorsplotsc, nrow = 2, labels = c('a.', 'b.') )

ggsave('./results/htseq/SI/Figure_S14.pdf', figure_S14, units = 'cm', width = 15, height = 12, useDingbats =F)
ggsave('./results/htseq/SI/Figure_S14.png', figure_S14, units = 'cm', width = 15, height = 12)

########################### Figure S22
##################### plot all, dc, non-dc genes:


############################
################## Figure S23 ----- intra tissue cov
# intracov = intra_ts_cov %>%
#   mutate(`age group` = factor(`age group`, levels = c('m3', 'm18', 'm24')) ) %>%
#   mutate(age = as.numeric( gsub('[a-z]','',`age group`) ) ) %>%
#   ggplot(aes(x = `age group`, y = `mean CoV`)) +
#   facet_wrap(~tissue, scales = 'free') +
#   geom_point(color='indianred4')  +
#   ylab('Mean CoV of Genes Among Cell Types') +
#   xlab('Age (in months)')
# 
# addline = intra_ts_cov %>%
#   mutate(age = as.numeric( gsub('[a-z]','',`age group`) ) ) %>%
#   group_by(tissue,`age group`) %>%
#   summarise(med = median(`mean CoV`))
# 
# intracov1 = intracov +
#   geom_text(data = cell.types,
#             size = 6/pntnorm, mapping = aes(x=0, y= Inf, label= paste('n[celltype]==',n)),
#             parse=T, vjust = 2, hjust=-0.1) +
#   geom_segment(data=addline, mapping = aes(x=rep(c(2,3,1)-0.1,4), y=med, xend = rep(c(2,3,1)+0.1,4), yend = med),
#                color=adjustcolor('gray10',alpha.f = 0.8))
# 
# intracov1
# ggsave('./results/SI_figures/Figure_S23.pdf', intracov1, units = 'cm', width = 12, height = 8, useDingbats =F)
# ggsave('./results/SI_figures/Figure_S23.png', intracov1, units = 'cm', width = 12, height = 8)
















