library(tidyverse)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(magick)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('up','down'))
revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))

sample_info = readRDS('./data/processed/figures/tidy/sample_info.rds')
expr = readRDS('./data/processed/figures/tidy/expression.rds')
expr_ch = readRDS('./data/processed/figures/tidy/expression_change.rds')
pca_dat = readRDS('./data/processed/figures/tidy/pca_data.rds')
cov_dat = readRDS('./data/processed/figures/tidy/CoV.rds')
cov_dat_wo_cortex = readRDS('./data/processed/figures/tidy/CoV_wo_cortex.rds')
cov_ch = readRDS('./data/processed/figures/tidy/CoV_change.rds')
cov_gsea = readRDS('./data/processed/figures/tidy/CoV_GSEA.rds')
divcon_gsea = readRDS('./data/processed/figures/tidy/divcon_GSEA.rds')
revgenes = readRDS('./data/processed/figures/tidy/revgenes.tissues.rds')

## plot age distribution in tissues:
# ages_log2 = sample_info %>%
#   ggplot(aes(y = age, x= tissue, color = tissue)) +
#   geom_hline(yintercept = 90, linetype = 'dashed', color = 'gray35') +
#   geom_jitter(width = 0.1, size = 1.5) + 
#   coord_flip() + 
#   scale_color_manual(values = tissuecol) +
#   scale_y_continuous(trans = 'log2', breaks = c(2,10,30,90,300, 900)) +
#   annotate('text', y = 100, label = 'Aging', x = 4.5, hjust = 0, size = 8/pntnorm) +
#   annotate('text', y = 80, label = 'Development', x = 4.5, hjust = 1, size = 8/pntnorm) +
#   ylab('Age in days (in log2 scale)') +
#   xlab(NULL) +
#   guides(color = F)  
# 
# # system('mkdir -p results')
# ggsave('./results/SI_figures/Figure_S1.pdf',ages_log2, units = 'cm', width = 8, height = 6, useDingbats = F)
# ggsave('./results/SI_figures/Figure_S1.png',ages_log2, units = 'cm', width = 8, height = 6)
# 
# ages = sample_info %>%
#   ggplot(aes(y = age, x= tissue, color = tissue)) +
#   geom_hline(yintercept = 90, linetype = 'dashed', color = 'gray35') +
#   geom_jitter(width = 0.1, size = 1.5) + 
#   coord_flip() + 
#   scale_color_manual(values = tissuecol) +
#   # scale_y_continuous(trans = 'log2', breaks = c(2,10,30,90,300, 900)) +
#   annotate('text', y = 100, label = 'Aging', x = 4.5, hjust = 0, size = 4/pntnorm) +
#   annotate('text', y = 80, label = 'Development', x = 4.5, hjust = 1, size = 4/pntnorm) +
#   ylab('Age') +
#   xlab(NULL) +
#   guides(color = F)  +
#   theme_pubr(base_size = 8)
# 
# 
# ggsave('./results/ages.pdf',ages, units = 'cm', width = 8, height = 6, useDingbats = F)
# ggsave('./results/ages.png',ages, units = 'cm', width = 8, height = 6)

# # calculate variations explained in all kinds of pca results:
# all_raw_var = (pca_dat %>%
#   filter(period == 'all', type == 'raw') %>%
#   select(varExp, PC) %>%
#   unique())$varExp
# 
# aging_raw_var= (pca_dat %>%
#                 filter(period == 'aging', type == 'raw') %>%
#                 select(varExp, PC) %>%
#                 unique())$varExp
# 
# dev_raw_var = (pca_dat %>%
#                   filter(period == 'development', type == 'raw') %>%
#                   select(varExp, PC) %>%
#                   unique())$varExp
# 
# all_notissue_var = (pca_dat %>%
#                 filter(period == 'all', type == 'notissue') %>%
#                 select(varExp, PC) %>%
#                 unique())$varExp
# 
# aging_notissue_var = (pca_dat %>%
#                   filter(period == 'aging', type == 'notissue') %>%
#                   select(varExp, PC) %>%
#                   unique())$varExp
# 
# dev_notissue_var = (pca_dat %>%
#                 filter(period == 'development', type == 'notissue') %>%
#                 select(varExp, PC) %>%
#                 unique())$varExp

###################
##### plot pca results:
## raw pca(without removing tissue effect):
# all_raw_pca12 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2') +
#   coord_fixed(ratio = all_raw_var[2]/all_raw_var[1], clip = 'off') +
#   xlab(paste('PC1 (', round(all_raw_var[1]*100),'%)',sep='')) +
#   ylab(paste('PC2 (', round(all_raw_var[2]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) +
#   theme(legend.position = c(0.05,0.6),
#         legend.justification=c(0,0),
#         legend.direction = 'vertical',
#         legend.box = 'horizontal',
#         legend.background = element_rect(fill = 'gray85',color = 'gray25')) 
# 
# ggsave('./results/pca/pca_all_raw12.pdf',all_raw_pca12, units = 'cm', width = 8, height = 5, useDingbats =F)
# ggsave('./results/pca/pca_all_raw12.png',all_raw_pca12, units = 'cm', width = 8, height = 5)
# 
# all_raw_pc1age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC1, color = tissue)) +
#   geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# all_raw_pc2age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC2, color = tissue)) +
#   geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# all_raw_pc12 = ggarrange(all_raw_pca12, ggarrange(all_raw_pc1age, all_raw_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')
# 
# ggsave('./results/pca/pca_all_raw12_wage.pdf',all_raw_pc12, units = 'cm', width = 16, height = 6,
#        useDingbats = F)
# ggsave('./results/pca/pca_all_raw12_wage.png',all_raw_pc12, units = 'cm', width = 16, height = 6)
# 
# all_raw_pca34 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2') +
#   coord_fixed(ratio = all_raw_var[4]/all_raw_var[3], clip = 'off') +
#   xlab(paste('PC3 (', round(all_raw_var[3]*100),'%)',sep='')) +
#   ylab(paste('PC4 (', round(all_raw_var[4]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) +
#   theme(legend.position = 'top',
#         legend.direction = 'horizontal',
#         legend.box = 'vertical',
#         legend.background = element_rect(fill = 'gray85',color = 'gray25'))
# 
# ggsave('./results/pca/pca_all_raw34.pdf',all_raw_pca34, units = 'cm', width = 8, height = 5, useDingbats =F)
# ggsave('./results/pca/pca_all_raw34.png',all_raw_pca34, units = 'cm', width = 8, height = 5)
# 
# all_raw_pc3age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC3, color = tissue)) +
#   geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# all_raw_pc4age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC4, color = tissue)) +
#   geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# ### pc4-age correlations:
# pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   mutate(`age period`= c('development', 'ageing')[1 + (age>90) ] ) %>%
#   group_by(tissue, `age period`) %>%
#   summarise(rho = cor.test(PC4, age, m='s')$est,
#             pval = cor.test(PC4, age, m='s')$p.val)
# 
# all_raw_pc34 = ggarrange(all_raw_pca34, ggarrange(all_raw_pc3age, all_raw_pc4age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), legend = T, labels = c('a',NA))
# 
# 
# ggsave('./results/pca/pca_all_raw34_wage.pdf', all_raw_pc34, units = 'cm', width = 16, height = 6,
#        useDingbats = F)
# ggsave('./results/pca/pca_all_raw34_wage.png', all_raw_pc34, units = 'cm', width = 16, height = 6)
# 
## scaled pca (removing tissue effect):
#pc1,2 age correlations:

# pca_dat %>%
#   select(-varExp) %>%
#   filter(period =='all', type =='notissue') %>%
#   spread(key ='PC', value = 'value') %>%
#   left_join(sample_info) %>%
#   mutate(`age period`= c('development', 'ageing')[1 + (age>90) ] ) %>%
#   group_by(tissue, `age period`) %>%
#   summarise(rho = cor.test(PC2, age, m='s')$est,
#             pval = cor.test(PC2, age, m='s')$p.val)
# 
# all_nt_pca12 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2') +
#   coord_fixed(ratio = all_notissue_var[2]/all_notissue_var[1], clip = 'off') +
#   xlab(paste('PC1 (', round(all_notissue_var[1]*100),'%)',sep='')) +
#   ylab(paste('PC2 (', round(all_notissue_var[2]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) +
#   theme(legend.position = c(0.7,0.9),
#         legend.direction = c('vertical'),
#         legend.box = 'horizontal',
#         legend.background = element_rect(fill = 'gray85',color = 'gray25')) 
# 
# ggsave('./results/pca/pca_all_nt12.pdf',all_nt_pca12, units = 'cm', width = 8, height = 8, useDingbats =F)
# ggsave('./results/pca/pca_all_nt12.png',all_nt_pca12, units = 'cm', width = 8, height = 8)
# 
# all_nt_pc1age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC1, color = tissue)) +
#   geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# all_nt_pc2age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC2, color = tissue)) +
#   geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# all_nt_pc12 = ggarrange(all_nt_pca12, ggarrange(all_nt_pc1age, all_nt_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F, align = 'v')
# 
# ggsave('./results/pca/pca_all_nt12_wage.pdf',all_nt_pc12, units = 'cm', width = 16, height = 9,
#        useDingbats = F)
# ggsave('./results/pca/pca_all_nt12_wage.png',all_nt_pc12, units = 'cm', width = 16, height = 9)

all_nt_pca34 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'notissue') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2') +
  coord_fixed(ratio = all_notissue_var[2]/all_notissue_var[1], clip = 'off') +
  xlab(paste('PC3 (', round(all_notissue_var[3]*100),'%)',sep='')) +
  ylab(paste('PC4 (', round(all_notissue_var[4]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) +
  theme(legend.position = c(0.25,0.8),
        legend.direction = c('vertical'),
        legend.box = 'horizontal',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) 

ggsave('./results/pca/pca_all_nt34.pdf',all_nt_pca34, units = 'cm', width = 8, height = 8, useDingbats =F)
ggsave('./results/pca/pca_all_nt34.png',all_nt_pca34, units = 'cm', width = 8, height = 8)

## development and ageing separately: 
# dev_raw_pca12 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2',breaks = c(2,8,30,60)) +
#   coord_fixed(ratio = dev_raw_var[2]/dev_raw_var[1], clip = 'off') +
#   xlab(paste('PC1 (', round(dev_raw_var[1]*100),'%)',sep='')) +
#   ylab(paste('PC2 (', round(dev_raw_var[2]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) +
#   theme(legend.position = c(0.06,0.65),
#         legend.justification=c(0,0),
#         legend.direction = 'vertical',
#         legend.box = 'horizontal',
#         legend.background = element_rect(fill = 'gray85',color = 'gray25')) 
# 
# ggsave('./results/pca/pca_dev_raw12.pdf',dev_raw_pca12, units = 'cm', width = 8, height = 5, useDingbats =F)
# ggsave('./results/pca/pca_dev_raw12.png',dev_raw_pca12, units = 'cm', width = 8, height = 5)
# 
# dev_raw_pc1age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC1, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# dev_raw_pc2age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC2, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# dev_raw_pc12 = ggarrange(dev_raw_pca12, ggarrange(dev_raw_pc1age, dev_raw_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')
# 
# ggsave('./results/pca/pca_dev_raw12_wage.pdf',dev_raw_pc12, units = 'cm', width = 16, height = 6,
#        useDingbats = F)
# ggsave('./results/pca/pca_dev_raw12_wage.png',dev_raw_pc12, units = 'cm', width = 16, height = 6)

# dev_raw_pca34 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2',breaks = c(2,8,30,60)) +
#   coord_fixed(ratio = dev_raw_var[4]/dev_raw_var[3], clip = 'off') +
#   xlab(paste('PC3 (', round(dev_raw_var[3]*100),'%)',sep='')) +
#   ylab(paste('PC4 (', round(dev_raw_var[4]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) +
#   theme(legend.position = c(0.15,0.05),
#         legend.justification=c(0,0),
#         legend.direction = 'horizontal',
#         legend.box = 'vertical',
#         legend.background = element_rect(fill = 'gray85',color = 'gray25')) 
# 
# ggsave('./results/pca/pca_dev_raw34.pdf',dev_raw_pca34, units = 'cm', width = 8, height = 5, useDingbats =F)
# ggsave('./results/pca/pca_dev_raw34.png',dev_raw_pca34, units = 'cm', width = 8, height = 5)
# 
# dev_raw_pc3age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC3, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# dev_raw_pc4age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC4, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# dev_raw_pc34 = ggarrange(dev_raw_pca34, ggarrange(dev_raw_pc3age, dev_raw_pc4age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F, align = 'v')
# 
# ggsave('./results/pca/pca_dev_raw34_wage.pdf',dev_raw_pc34, units = 'cm', width = 16, height = 6,
#        useDingbats = F)
# ggsave('./results/pca/pcget_legend(dev_raw_pca34)a_dev_raw34_wage.png',dev_raw_pc34, units = 'cm', width = 16, height = 6)

### dev pc1,2,3,4 with ages

# dev_raw_pc1234 = ggarrange(ggarrange(dev_raw_pca12, dev_raw_pca34 + theme(legend.position = 'none'),
#                                      ncol = 2, align = 'h', labels = c('a', 'b'), vjust = 2.7),
#                            ggarrange(dev_raw_pc1age, dev_raw_pc2age,dev_raw_pc3age, dev_raw_pc4age,
#                                      align = 'hv', ncol =4), nrow=2, heights = c(1, 0.5), labels = c(NA, 'c'))
# 
# ggsave('./results/pca/pca_dev_raw1234_wage.pdf', dev_raw_pc1234, units = 'cm', width = 16, height = 12,
#        useDingbats = F)
# ggsave('./results/pca/pca_dev_raw1234_wage.png', dev_raw_pc1234, units = 'cm', width = 16, height = 12)
# 
# ## pc-age cors:
# pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'raw') %>%
#   left_join(sample_info) %>%
#   select(-period, -type, -log2age) %>%
#   group_by(PC, tissue) %>%
#   summarise(rho = cor.test(value, age, m='s')$est,
#             pval = cor.test(value, age, m='s')$p.val)



######

# dev_nt_pca12 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2',breaks = c(2,8,30,60)) +
#   coord_fixed(ratio = dev_notissue_var[2]/dev_notissue_var[1], clip = 'off') +
#   xlab(paste('PC1 (', round(dev_notissue_var[1]*100),'%)',sep='')) +
#   ylab(paste('PC2 (', round(dev_notissue_var[2]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) +
#   theme(legend.position = c(0.26,0.65),
#         legend.justification=c(0,0),
#         legend.direction = 'vertical',
#         legend.box = 'horizontal',
#         legend.background = element_rect(fill = 'gray85',color = 'gray25')) 
# 
# ggsave('./results/pca/pca_dev_nt12.pdf', dev_nt_pca12, units = 'cm', width = 8, height = 5, useDingbats =F)
# ggsave('./results/pca/pca_dev_nt12.png', dev_nt_pca12, units = 'cm', width = 8, height = 5)
# 
# dev_nt_pc1age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC1, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# dev_nt_pc2age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC2, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# dev_nt_pc12 = ggarrange(dev_nt_pca12, ggarrange(dev_nt_pc1age, dev_nt_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')
# 
# ggsave('./results/pca/pca_dev_nt12_wage.pdf',dev_nt_pc12, units = 'cm', width = 16, height = 6,
#        useDingbats = F)
# ggsave('./results/pca/pca_dev_nt12_wage.png',dev_nt_pc12, units = 'cm', width = 16, height = 6)

# dev_nt_pca34 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2',breaks = c(2,8,30,60)) +
#   coord_fixed(ratio = dev_notissue_var[4]/dev_notissue_var[3], clip = 'off') +
#   xlab(paste('PC1 (', round(dev_notissue_var[3]*100),'%)',sep='')) +
#   ylab(paste('PC2 (', round(dev_notissue_var[4]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) 
#   # theme(legend.position = c(0.35,0.05)
#   #       legend.justification=c(0,0),
#   #       legend.direction = 'horizontal',
#   #       legend.box = 'vertical',
#   #       legend.background = element_rect(fill = 'gray85',color = 'gray25'))
# 
# ggsave('./results/pca/pca_dev_nt34.pdf',dev_nt_pca34, units = 'cm', width = 8, height = 5, useDingbats =F)
# ggsave('./results/pca/pca_dev_nt34.png',dev_nt_pca34, units = 'cm', width = 8, height = 5)
# 
# dev_nt_pc3age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC3, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# dev_nt_pc4age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'development', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC4, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# dev_nt_pc34 = ggarrange(dev_nt_pca34, ggarrange(dev_nt_pc3age, dev_nt_pc4age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA))
# 
# ggsave('./results/pca/pca_dev_nt34_wage.pdf',dev_nt_pc34, units = 'cm', width = 16, height = 6,
#        useDingbats = F)
# ggsave('./results/pca/pca_dev_nt34_wage.png',dev_nt_pc34, units = 'cm', width = 16, height = 6)

# aging_raw_pca12 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2') +
#   coord_fixed(ratio = aging_raw_var[2]/aging_raw_var[1], clip = 'off') +
#   xlab(paste('PC1 (', round(aging_raw_var[1]*100),'%)',sep='')) +
#   ylab(paste('PC2 (', round(aging_raw_var[2]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) +
#   theme(legend.position = c(0.06,0.65),
#         legend.justification=c(0,0),
#         legend.direction = 'vertical',
#         legend.box = 'horizontal',
#         legend.background = element_rect(fill = 'gray85',color = 'gray25')) 
# 
# ggsave('./results/pca/pca_aging_raw12.pdf',aging_raw_pca12, units = 'cm', width = 8, height = 5,
#        useDingbats =F)
# ggsave('./results/pca/pca_aging_raw12.png',aging_raw_pca12, units = 'cm', width = 8, height = 5)
# 
# aging_raw_pc1age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC1, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# aging_raw_pc2age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC2, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# aging_raw_pc12 = ggarrange(aging_raw_pca12, ggarrange(aging_raw_pc1age, aging_raw_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')
# 
# ggsave('./results/pca/pca_aging_raw12_wage.pdf',aging_raw_pc12, units = 'cm', width = 16, height = 6,
#        useDingbats = F)
# ggsave('./results/pca/pca_aging_raw12_wage.png',aging_raw_pc12, units = 'cm', width = 16, height = 6)
# 
# aging_raw_pca34 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2') +
#   coord_fixed(ratio = aging_raw_var[4]/aging_raw_var[3], clip = 'off') +
#   xlab(paste('PC3 (', round(aging_raw_var[3]*100),'%)',sep='')) +
#   ylab(paste('PC4 (', round(aging_raw_var[4]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) +
#   theme(legend.position = c(0.16,0.65),
#         legend.justification=c(0,0),
#         legend.direction = 'vertical',
#         legend.box = 'horizontal',
#         legend.background = element_rect(fill = 'gray85',color = 'gray25')) 
# 
# ggsave('./results/pca/pca_aging_raw34.pdf',aging_raw_pca34, units = 'cm', width = 8, height = 5,
#        useDingbats =F)
# ggsave('./results/pca/pca_aging_raw34.png',aging_raw_pca34, units = 'cm', width = 8, height = 5)
# 
# aging_raw_pc3age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC3, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# aging_raw_pc4age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC4, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# aging_raw_pc34 = ggarrange(aging_raw_pca34, ggarrange(aging_raw_pc3age, aging_raw_pc4age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')
# 
# ggsave('./results/pca/pca_aging_raw34_wage.pdf',aging_raw_pc34, units = 'cm', width = 16, height = 6,
#        useDingbats = F)
# ggsave('./results/pca/pca_aging_raw34_wage.png',aging_raw_pc34, units = 'cm', width = 16, height = 6)

# aging_nt_pca12 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2') +
#   coord_fixed(ratio = aging_notissue_var[2]/aging_notissue_var[1], clip = 'off') +
#   xlab(paste('PC1 (', round(aging_notissue_var[1]*100),'%)',sep='')) +
#   ylab(paste('PC2 (', round(aging_notissue_var[2]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) +
#   theme(legend.position = c(0.36,0.1),
#         legend.justification=c(0,0),
#         legend.direction = 'vertical',
#         legend.box = 'horizontal',
#         legend.background = element_rect(fill = 'gray85',color = 'gray25')) 
# 
# ggsave('./results/pca/pca_aging_nt12.pdf',aging_nt_pca12, units = 'cm', width = 8, height = 5,
#        useDingbats =F)
# ggsave('./results/pca/pca_aging_nt12.png',aging_nt_pca12, units = 'cm', width = 8, height = 5)
# 
# aging_nt_pc1age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC1, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# aging_nt_pc2age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC2, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# aging_nt_pc12 = ggarrange(aging_nt_pca12, ggarrange(aging_nt_pc1age, aging_nt_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')
# 
# ggsave('./results/pca/pca_aging_nt12_wage.pdf',aging_nt_pc12, units = 'cm', width = 16, height = 6,
#        useDingbats = F)
# ggsave('./results/pca/pca_aging_nt12_wage.png',aging_nt_pc12, units = 'cm', width = 16, height = 6)
# 
# aging_nt_pca34 = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2') +
#   coord_fixed(ratio = aging_notissue_var[4]/aging_notissue_var[3], clip = 'off') +
#   xlab(paste('PC1 (', round(aging_notissue_var[3]*100),'%)',sep='')) +
#   ylab(paste('PC2 (', round(aging_notissue_var[4]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'), 
#          size = guide_legend('Age')) +
#   theme(legend.position = c(0.16,0.05),
#         legend.justification=c(0,0),
#         legend.direction = 'vertical',
#         legend.box = 'horizontal',
#         legend.background = element_rect(fill = 'gray85',color = 'gray25')) 
# 
# ggsave('./results/pca/pca_aging_nt34.pdf',aging_nt_pca34, units = 'cm', width = 8, height = 5,
#        useDingbats =F)
# ggsave('./results/pca/pca_aging_nt34.png',aging_nt_pca34, units = 'cm', width = 8, height = 5)
# 
# aging_nt_pc3age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC3, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# aging_nt_pc4age = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'aging', type == 'notissue') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = age, y = PC4, color = tissue)) +
#   geom_smooth(alpha = 0.1) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = tissuecol) +
#   scale_x_continuous(trans = 'log2') +
#   guides(color = F) +
#   xlab('Age (log2)')
# 
# aging_nt_pc34 = ggarrange(aging_nt_pca34, ggarrange(aging_nt_pc3age, aging_nt_pc4age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')
# 
# ggsave('./results/pca/pca_aging_nt34_wage.pdf',aging_nt_pc34, units = 'cm', width = 16, height = 6,
#        useDingbats = F)
# ggsave('./results/pca/pca_aging_nt34_wage.png',aging_nt_pc34, units = 'cm', width = 16, height = 6)
# 
#### expression change distrubiton compare between development and ageing:

## top reversal gene:
top_rev_gene_dat = expr_ch %>%
  select(-p,-FDR) %>%
  spread(key = period, value =`Expression Change`) %>%
  mutate(rev = aging *development) %>%
  top_n(n=1, wt=-rev) %>%
  left_join(expr) %>%
  inner_join(sample_info)
  
reversalgene_most = ggplot(top_rev_gene_dat,aes(x =age, y= expression, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_smooth() +
  geom_point() +
  scale_x_continuous(trans = 'log2', breaks = c(2,10,30,90,300, 900)) +
  scale_color_manual(values = tissuecol) +
  ggtitle(unique(top_rev_gene_dat$gene_id)) +
  xlab('Age (log2)') +
  ylab('Gene Expression') +
  guides(color = F)

ggsave('./results/reversalgene.pdf',reversalgene_most, units = 'cm', width = 8, height = 6, useDingbats = F)
ggsave('./results/reversalgene.png',reversalgene_most, units = 'cm', width = 8, height = 6)

# expression change proportions of all genes in development and ageing:
allgenes_exp_ch_p = expr_ch %>%
  group_by(tissue,period) %>%
  summarise( up = sum(`Expression Change`>0,na.rm=T),
             down = sum(`Expression Change`<0,na.rm=T)) %>%
  gather(key= 'direction',value ='n',-tissue,-period) %>%
  mutate(period = factor(period, levels = c('development','aging'))) %>%
  ggplot(aes(x = period, y= n, fill = direction)) +
  facet_wrap(~tissue, ncol=4) +
  geom_bar(stat= 'identity', position ='fill') +
  scale_fill_manual(values = regcol) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray35') +
  xlab(NULL) +ylab(NULL) +
  guides(fill = guide_legend('Direction of Expression Change', 
                             override.aes = list(size = 2))) +
  theme(legend.position = 'top',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) 

ggsave('./results/allgenes_exp_ch.pdf',allgenes_exp_ch_p, units = 'cm', width = 10, height = 6,
       useDingbats = F)
ggsave('./results/allgenes_exp_ch.png',allgenes_exp_ch_p, units = 'cm', width = 10, height = 6)

# expression changes of significant genes in development and ageing:
siggenes_exp_ch_p = expr_ch %>%
  filter(FDR<=0.1) %>%
  group_by(tissue,period) %>%
  summarise( up = sum(`Expression Change`>0,na.rm=T),
             down = sum(`Expression Change`<0,na.rm=T)) %>%
  gather(key= 'direction',value ='n',-tissue,-period) %>%
  mutate(period = factor(period, levels = c('development','aging'))) %>%
  mutate(n = ifelse(n==0,NA,n)) %>%
  ggplot(aes(x = period, y= n, fill = direction)) +
  facet_wrap(~tissue, ncol=4) +
  geom_bar(stat= 'identity', position = 'dodge') +
  scale_fill_manual(values = regcol, drop=FALSE) +
  scale_y_continuous(trans = 'log10') +
  geom_text(aes(label = n), color='gray15', position = position_dodge(width = 1), angle=90, hjust = 1.1,
            vjust = 0.5, size = 6/pntnorm) +
  guides(fill = guide_legend('Direction of Expression Change', 
                             override.aes = list(size = 2))) +
  theme(legend.position = 'bottom',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) +
  xlab(NULL)  + ylab(NULL) +
  ggtitle('Number of Genes with a Significant Change (FDR≤0.1)')

ggsave('./results/siggenes_exp_ch.pdf',siggenes_exp_ch_p, units = 'cm', width = 10, height = 6,
       useDingbats = F)
ggsave('./results/siggenes_exp_ch.png',siggenes_exp_ch_p, units = 'cm', width = 10, height = 6)


## significantly changing genes (fdr<=0.2) overlap across tissues:
# siggenes_overlap_fdr0.2 = expr_ch %>%
#   filter(FDR <= 0.2) %>%
#   mutate(direction = `Expression Change` > 0) %>%
#   mutate(direction = ifelse(direction == TRUE, 'up', 'down')) %>%
#   mutate(period = factor(period, levels = c('development', 'aging'))) %>%
#   group_by(gene_id, period, direction) %>%
#   summarise(n = length(tissue)) %>%
#   group_by(n, period, direction) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
#   slice(-c(1:4)) %>%
#   ggplot(aes(x = n, y = count, fill = direction)) +
#   facet_wrap(~period) +
#   geom_bar(stat='identity', position = 'dodge') +
#   scale_fill_manual(values = regcol, drop = F) +
#   scale_y_continuous(trans = 'log10') +
#   geom_text(aes(label = count), color = 'gray15', position = position_dodge(width = 1), angle = 90,
#             hjust = 1.1, size = 6/pntnorm) +
#   guides(fill = guide_legend('Direction of Expression Change'),
#          override.aes = list(size =2)) +
#   theme(legend.position = 'bottom',
#         legend.background = element_rect(fill = 'gray85', color = 'gray25')) +
#   xlab(NULL) + ylab(NULL) +
#   ggtitle('Overlap of Significantly Changing Genes Among Tissues (FDR<=0.2)')
# 
# ggsave('./results/siggenes_exp_ch_overlap_fdr0.2.pdf', siggenes_overlap_fdr0.2, units = 'cm',
#        width = 10, height = 6, useDingbats = F)  
# ggsave('./results/siggenes_exp_ch_overlap_fdr0.2.png', siggenes_overlap_fdr0.2, units = 'cm',
#        width = 10, height = 6)

## all genes (no cutoff) overlap across tissues:
# siggenes_overlap_nocutoff = expr_ch %>%
#   filter(!is.na(`Expression Change`)) %>%
#   filter(`Expression Change` != 0) %>%
#   mutate(direction = `Expression Change` > 0) %>%
#   mutate(direction = ifelse(direction == TRUE, 'up', 'down')) %>%
#   mutate(period = factor(period, levels = c('development', 'aging'))) %>%
#   group_by(gene_id, period, direction) %>%
#   summarise(n = length(tissue)) %>%
#   group_by(n, period, direction) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
#   slice(-c(1:4)) %>%
#   ggplot(aes(x = n, y = count, fill = direction)) +
#   facet_wrap(~period) +
#   geom_bar(stat='identity', position = 'dodge') +
#   scale_fill_manual(values = regcol, drop = F) +
#   geom_text(aes(label = count), color = 'gray15', position = position_dodge(width = 1), angle = 90,
#             hjust = 1.1, size = 6/pntnorm) +
#   guides(fill = guide_legend('Direction of Expression Change'),
#          override.aes = list(size =2)) +
#   theme(legend.position = 'bottom',
#         legend.background = element_rect(fill = 'gray85', color = 'gray25')) +
#   xlab(NULL) + ylab(NULL) +
#   ggtitle('Overlap of Expression Change Among Tissues (no cutoff)')

# ggsave('./results/exp_ch_overlap_nocutoff.pdf', siggenes_overlap_nocutoff, units = 'cm',
#        width = 10, height = 6, useDingbats = F)  
# ggsave('./results/exp_ch_overlap_nocutoff.png', siggenes_overlap_nocutoff, units = 'cm',
#        width = 10, height = 6)
# 
## maybe combine them to one plot:

# siggenes_overlap_two_cutoff = ggarrange(siggenes_overlap, siggenes_overlap_fdr0.2, ncol =2, nrow = 1,
#                                         common.legend = T, labels = c('a','b'), align = 'hv',
#                                         legend = 'bottom')
# ggsave('./results/siggenes_exp_ch_overlap_two_cutoff.pdf', siggenes_overlap_two_cutoff, units = 'cm',
#        width = 16, height = 6, useDingbats = F)
# ggsave('./results/siggenes_exp_ch_overlap_two_cutoff.png', siggenes_overlap_two_cutoff, units = 'cm',
#        width = 16, height = 6)

###buradasin

## choose more interesting gene (top third) showing div-conv pattern:
# top3rd_divcon_gene_dat = cov_ch %>%
#   select(-p,-FDR) %>%
#   spread(key = period, value =`CoV_change`) %>%
#   mutate(rev = aging *development) %>%
#   top_n(n=3, wt=-(rev)) %>%
#   top_n(n=1, wt=rev) %>%
#   left_join(expr) %>%
#   inner_join(sample_info) 
# 
# top3rd_divcon_gene = ggplot(top3rd_divcon_gene_dat, aes(x =age, y= expression, color = tissue)) +
#   geom_smooth(se=F,size = 0.7) +
#   geom_point(size=0.5) +
#   geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
#   scale_x_continuous(trans = 'log2') +
#   scale_color_manual(values = tissuecol) +
#   ggtitle(unique(top3rd_divcon_gene_dat$gene_id)) +
#   guides(color = guide_legend('Tissue', 
#                               override.aes = list(size = 2))) +
#   theme(legend.position = 'top', 
#         legend.background = element_rect(fill = 'gray85',color = 'gray25')) +
#   xlab('Age (days)') + ylab('Gene Expression')
# 
# ggsave('./results/top3rd_divcon_gene.pdf',top3rd_divcon_gene, units = 'cm', width = 8, height = 8,
#        useDingbats = F)
# ggsave('./results/top3rd_divcon_gene.png',top3rd_divcon_gene, units = 'cm', width = 8, height = 8)
# 
# CoV change of top 3rd gene showing divergence-convergence pattern:
# top3rd_divcon_gene_cov = top3rd_divcon_gene_dat %>%
#   group_by(gene_id, ind_id, age) %>%
#   summarise(sd = sd(expression),
#             mean = mean(expression),
#             cov = sd/mean) %>%
#   ggplot(aes(x = age, y= cov)) +
#   geom_smooth(size = 0.7) +
#   geom_point(size=0.5) +
#   geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
#   scale_x_continuous(trans = 'log2') +
#   ggtitle(unique(top3rd_divcon_gene_dat$gene_id)) +
#   xlab('Age (days)') + ylab('Inter-Tissue Coefficient of Variation(CoV)')
# 
# ggsave('./results/top3rd_divcon_gene_cov.pdf',top3rd_divcon_gene_cov, units = 'cm', width = 8, height = 8,
#        useDingbats = F)
# ggsave('./results/top3rd_divcon_gene_cov.png',top3rd_divcon_gene_cov, units = 'cm', width = 8, height = 8)
# 
# top3rd_divcon_gene_p = ggarrange(top3rd_divcon_gene + ggtitle(NULL) +
#                                    theme(legend.title = element_text(size = 5),
#                                          legend.text = element_text(size = 5)),
#                                  top3rd_divcon_gene_cov , labels = 'auto', font.label=list(size = 10),
#                                ncol=2, nrow=1, align = 'v')
# 
# ggsave('./results/top3rd_divcon_gene_expcov.pdf',top3rd_divcon_gene_p, units = 'cm', width = 10, height = 5,
#        useDingbats = F)
# ggsave('./results/top3rd_divcon_gene_expcov.png',top3rd_divcon_gene_p, units = 'cm', width = 10, height = 5)
# 
# 
# cov_median = ggplot(cov_dat_sum, aes(x = age, y= medianCoV)) +
#   geom_point() +
#   geom_smooth(se=T,color = 'midnightblue') +
#   scale_x_continuous(trans = 'log2') +
#   geom_vline(xintercept = 90, linetype='dashed',color = 'gray35') +
#   xlab('Age (days)') + ylab('Median CoV') +
#   annotate('text',x=95,y=0.34,label='Aging',hjust=0,size = 8/pntnorm) +
#   annotate('text',x = 95, y=0.3375, label = parse(text = paste('rho["CoV,age"] ==' ,(round(filter(cov_cordat_median, period=='aging')$cor,3)))),hjust=0,size = 8/pntnorm) +
#   annotate('text',x = 95, y=0.335, label = parse(text = paste0('p ==' ,(round(filter(cov_cordat_median, period=='aging')$cor.p,2)))),hjust=0,size = 8/pntnorm) +
#   annotate('text',x=8,y=0.315,label='Development',hjust=0,size = 8/pntnorm) +
#   annotate('text',x = 8, y=0.3125, label = parse(text = paste('rho["CoV,age"] ==' ,(round(filter(cov_cordat_median, period=='development')$cor,3)))),hjust=0,size = 8/pntnorm) +
#   annotate('text',x = 8, y=0.31, label = parse(text = paste0('p ==' ,(round(filter(cov_cordat_median, period=='development')$cor.p,2)))),hjust=0,size = 8/pntnorm)
# 
# ggsave('./results/cov_median.pdf',cov_median, units = 'cm', width = 8, height = 8, useDingbats = F)
# ggsave('./results/cov_median.png',cov_median, units = 'cm', width = 8, height = 8)
# 
# 
# ## CoV average for each individual and mean(CoV)-age correlation (without cortex):
# cov_dat_sum_wo_cortex = cov_dat_wo_cortex %>%
#   mutate(ind_id = factor(ind_id)) %>%
#   left_join(unique(select(sample_info,-tissue,-sample_id))) %>%
#   group_by(ind_id, age) %>%
#   summarise(meanCoV = mean(CoV),
#             medianCoV = median(CoV)) %>% 
#   mutate(period = c('aging','development')[1+(age<=90)])
# 
# cov_cordat_wo_cortex = group_by(ungroup(cov_dat_sum_wo_cortex),period) %>%
#   summarise(cor = cor.test(meanCoV, age, method = 's')$est,
#             cor.p = cor.test(meanCoV, age, method = 's')$p.val)
# 
# cov_cordat_median_wo_cortex = group_by(ungroup(cov_dat_sum_wo_cortex),period) %>%
#   summarise(cor = cor.test(medianCoV, age, method = 's')$est,
#             cor.p = cor.test(medianCoV, age, method = 's')$p.val)
# 
# cov_mean_wo_cortex = ggplot(cov_dat_sum_wo_cortex, aes(x = age, y= meanCoV)) +
#   geom_point() +
#   geom_smooth(se=T,color = 'midnightblue') +
#   scale_x_continuous(trans = 'log2') +
#   geom_vline(xintercept = 90, linetype='dashed',color = 'gray35') +
#   xlab('Age (days)') + ylab('Mean CoV') +
#   annotate('text',x=95,y=0.42,label='Aging',hjust=0,size = 8/pntnorm) +
#   annotate('text',x = 95, y=0.416, label = parse(text = paste('rho["CoV,age"] ==' ,(round(filter(cov_cordat_wo_cortex, period=='aging')$cor,3)))),hjust=0,size = 8/pntnorm) +
#   annotate('text',x = 95, y=0.412, label = parse(text = paste0('p ==' ,(round(filter(cov_cordat_wo_cortex, period=='aging')$cor.p,3)))),hjust=0,size = 8/pntnorm) +
#   annotate('text',x=9,y=0.385,label='Development',hjust=0,size = 8/pntnorm) +
#   annotate('text',x = 9, y=0.381, label = parse(text = paste('rho["CoV,age"] ==' ,(round(filter(cov_cordat_wo_cortex, period=='development')$cor,3)))),hjust=0,size = 8/pntnorm) +
#   annotate('text',x = 9, y=0.377, label = parse(text = paste0('p ==' ,(round(filter(cov_cordat_wo_cortex, period=='development')$cor.p,2)))),hjust=0,size = 8/pntnorm)
# 
# ggsave('./results/cov_mean_wo_cortex.pdf',cov_mean_wo_cortex, units = 'cm', width = 8, height = 8,
#        useDingbats = F)
# ggsave('./results/cov_mean_wo_cortex.png',cov_mean_wo_cortex, units = 'cm', width = 8, height = 8)
# 
# cov_median_wo_cortex = ggplot(cov_dat_sum_wo_cortex, aes(x = age, y= medianCoV)) +
#   geom_point() +
#   geom_smooth(se=T,color = 'midnightblue') +
#   scale_x_continuous(trans = 'log2') +
#   geom_vline(xintercept = 90, linetype='dashed',color = 'gray35') +
#   xlab('Age (days)') + ylab('Median CoV') +
#   annotate('text',x=95,y=0.275,label='Aging',hjust=0,size = 8/pntnorm) +
#   annotate('text',x = 95, y=0.27, label = parse(text = paste('rho["CoV,age"] ==' ,(round(filter(cov_cordat_median_wo_cortex, period=='aging')$cor,3)))),hjust=0,size = 8/pntnorm) +
#   annotate('text',x = 95, y=0.265, label = parse(text = paste0('p ==' ,(round(filter(cov_cordat_median_wo_cortex, period=='aging')$cor.p,2)))),hjust=0,size = 8/pntnorm) +
#   annotate('text',x=10,y=0.225,label='Development',hjust=0,size = 8/pntnorm) +
#   annotate('text',x=10, y=0.22, label = parse(text = paste('rho["CoV,age"] ==' ,(round(filter(cov_cordat_median_wo_cortex, period=='development')$cor,3)))),hjust=0,size = 8/pntnorm) +
#   annotate('text',x=10, y=0.215, label = parse(text = paste0('p ==' ,(round(filter(cov_cordat_median_wo_cortex, period=='development')$cor.p,3)))),hjust=0,size = 8/pntnorm)
# 
# ggsave('./results/cov_median_wo_cortex.pdf',cov_median_wo_cortex, units = 'cm', width = 8, height = 8,
#        useDingbats = F)
# ggsave('./results/cov_median_wo_cortex.png',cov_median_wo_cortex, units = 'cm', width = 8, height = 8)
# 
# # combine mean and median CoV (without cortex):
# cov_mean_median_wo_cortex = ggarrange(cov_mean_wo_cortex , cov_median_wo_cortex,
#                                       ncol = 2, nrow = 1, align = 'hv')
# 
# ggsave('./results/cov_mean_median_wo_cortex.pdf',cov_mean_median_wo_cortex, units = 'cm',
#        width = 16, height = 8, useDingbats = F)
# ggsave('./results/cov_mean_median_wo_cortex.png',cov_mean_median_wo_cortex, units = 'cm',
#        width = 16, height = 8)

###

### revgene proportions:
# revprops = revgenes %>%
#   group_by(direction,tissue) %>%
#   summarise(n = length(gene_id)) %>%
#   group_by(tissue) %>%
#   #mutate(freq = n / sum(n)) %>%
#   #select(-n) %>%
#   mutate(direction = factor(direction,levels = c('UpDown','DownUp','UpUp','DownDown'))) %>%
#   ggplot(aes(x=tissue, y=n, fill = direction)) +
#   geom_bar(stat='identity', position = position_fill(reverse=T)) + 
#   scale_fill_manual(values = revcol) + 
#   geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray35') +
#   xlab(NULL) + ylab('Proportion of All Genes') +
#   guides(fill = guide_legend('Direction',
#                              override.aes = list(size=2))) +
#   theme(legend.position = 'top',
#         legend.background = element_rect(fill= 'gray85', color = 'gray25'))
# 
# ggsave('./results/reversal_props.pdf', revprops, units='cm', width = 8,height = 8, useDingbats = F)
# ggsave('./results/reversal_props.png', revprops, units='cm', width = 8,height = 8)

######### FIGURES

######## Fig1.
# fig1a_pca = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
#   geom_point(alpha = 0.7, show.legend = T) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.5,3), trans= 'log2') +
#   coord_fixed(ratio = all_raw_var[4]/all_raw_var[3], clip = 'off') +
#   xlab(paste('PC3 (', round(all_raw_var[3]*100),'%)',sep='')) +
#   ylab(paste('PC4 (', round(all_raw_var[4]*100),'%)',sep='')) +
#   guides(color = guide_legend('Tissue'),
#          size = guide_legend('Age')) +
#   theme(legend.position = c(0.7, 0.28),
#         legend.direction = 'vertical',
#         legend.box = 'horizontal',
#         legend.text = element_text(size = 3),
#         legend.title = element_text(size = 3),
#         legend.background = element_rect(fill = 'gray85',color = 'gray25'))
# 
# inset = pca_dat %>%
#   select(-varExp) %>%
#   filter(period == 'all', type == 'raw') %>%
#   spread(key = 'PC',value = 'value') %>%
#   left_join(sample_info) %>%
#   ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
#   geom_point(alpha = 0.7, show.legend = F) +
#   scale_color_manual(values = tissuecol) +
#   scale_size_continuous(range = c(0.2,1), trans= 'log2') +
#   coord_fixed(ratio = all_raw_var[2]/all_raw_var[1], clip = 'off') +
#   xlab(paste('PC1 (', round(all_raw_var[1]*100),'%)',sep='')) +
#   ylab(paste('PC2 (', round(all_raw_var[2]*100),'%)',sep='')) +
#   theme(axis.text = element_text(size = 3),
#         plot.background = element_rect(fill = "transparent",colour ="NA"),
#         axis.title = element_text(size = 3.5)) 
# 
# corfig = image_read("pairwise_tissue_cors.png")
# corimage = ggdraw() + draw_image(corfig, x = 1, hjust = 1, width = 1, height = 1)
# 
# # fig1_inset = ggdraw(fig1a_pca) +
# #   theme_half_open(font_size = 8) +
# #   draw_plot(inset, 0.83, 0.42, 0.5,0.5, width = 0.7, height = 0.7) 
# 
# ab = ggarrange(fig1a_pca, ggarrange(all_raw_pc3age, all_raw_pc4age, nrow = 2, ncol = 1), nrow = 1, ncol = 2,
#                widths = c(1, 0.25), labels = c("a","b"), font.label = list(size = 7))
# fig1 = ggarrange(ggarrange(ab, siggenes_exp_ch_p, nrow = 2, ncol = 1, labels = c("a","d"),
#                            font.label = list(size = 7)), corimage, nrow = 1, ncol = 2, widths = c(1, 0.5),
#                  labels = c("a","c"),font.label = list(size = 7)) +
#   theme_half_open(font_size = 12) +
#   draw_plot(inset, 0.46, 0.89, 0.5,0.5, width = 1, height = 1, scale = 0.22) +
#   theme_void()
# 
# ggsave('./results/fig1.pdf',fig1, units = 'cm', width = 16, height = 7, useDingbats = F)
# ggsave('./results/fig1.png',fig1, units = 'cm', width = 16, height = 7)


####################
#################### Compare magnitude of expression changes between development and ageing:
####################
####################

## significantly changing genes (fdr<0.1) overlap across tissues:
# siggenes_overlap = expr_ch %>%
#   filter(FDR < 0.1) %>%
#   mutate(direction = `Expression Change` > 0) %>%
#   mutate(direction = ifelse( direction == TRUE, 'up', 'down')) %>%
#   mutate(period = str_to_title(period)) %>%
#   mutate(period = factor(period, levels = c('Development','Aging'))) %>%
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
#   theme(legend.position = 'bottom',
#         legend.background = element_rect(fill = 'gray85', color = 'gray25')) +
#   xlab(NULL) + ylab(NULL) +
#   ggtitle('Overlap of Significantly Changing Genes Among Tissues (FDR<0.1)')


####### Compare magnitude of expression changes between development and ageing:
#######
# expr_ch %>%
#   dplyr::select(-p, -FDR) %>%
#   spread(key = period, value = `Expression Change`) %>%
#   group_by(tissue, gene_id) %>%
#   summarise(diff =  abs(development) - abs(aging)) %>%
#   summarise(p_val = wilcox.test(diff, mu = 0)$p.val) %>%
#   mutate(BH = p.adjust(p_val, method='BH'))
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
#   #ylab( expression( paste('Abs(',rho, '^["Dev"]' ) ) )
#   ylab('Abs. Rho Difference Between Dev. and Ageing')  +
#   theme(legend.position = 'none')
# 
# fig_s5 = ggarrange(siggenes_overlap, magnitude_comparison, ncol= 2, labels = c('a.','b.'), widths = c(1, 0.7),
#                     hjust = -0.1, vjust = c(1.5, 1.9), align = 'hv')
# 
# ggsave('./results/SI_figures/fig_SI5.pdf', fig_s5, width = 15, height = 8, units='cm', useDingbats = F )
# ggsave('./results/SI_figures/fig_SI5.png', fig_s5, width = 15, height = 8, units='cm' )  
# 
# 
####################
####################
####################
####################
























