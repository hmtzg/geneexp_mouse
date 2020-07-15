library(tidyverse)
library(ggpubr)
library(ggrepel)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)
tissuecol = setNames(c("#233789", "#f49e92", "#801008","#dbb32e"),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('up','down'))

sample_info = readRDS('./data/processed/figures/tidy/sample_info.rds')
expr = readRDS('./data/processed/figures/tidy/expression.rds')
expr_ch = readRDS('./data/processed/figures/tidy/expression_change.rds')
pca_dat = readRDS('./data/processed/figures/tidy/pca_data.rds')
cov_dat = readRDS('./data/processed/figures/tidy/CoV.rds')
cov_ch = readRDS('./data/processed/figures/tidy/CoV_change.rds')
cov_gsea = readRDS('./data/processed/figures/tidy/CoV_GSEA.rds')
divcon_gsea = readRDS('./data/processed/figures/tidy/divcon_GSEA.rds')

ages_log2 = sample_info %>%
  ggplot(aes(y = age, x= tissue, color = tissue)) +
  geom_hline(yintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_jitter(width = 0.1, size = 1.5) + 
  coord_flip() + 
  scale_color_manual(values = tissuecol) +
  scale_y_continuous(trans = 'log2', breaks = c(2,10,30,90,300, 900)) +
  annotate('text', y = 100, label = 'Aging', x = 4.5, hjust = 0, size = 8/pntnorm) +
  annotate('text', y = 80, label = 'Development', x = 4.5, hjust = 1, size = 8/pntnorm) +
  ylab('Age (log2)') +
  xlab(NULL) +
  guides(color = F)  
ages_log2
# system('mkdir -p results')
ggsave('./results/ages_log2.pdf',ages_log2, units = 'cm', width = 8, height = 6, useDingbats = F)
ggsave('./results/ages_log2.png',ages_log2, units = 'cm', width = 8, height = 6)

ages = sample_info %>%
  ggplot(aes(y = age, x= tissue, color = tissue)) +
  geom_hline(yintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_jitter(width = 0.1, size = 1.5) + 
  coord_flip() + 
  scale_color_manual(values = tissuecol) +
  # scale_y_continuous(trans = 'log2', breaks = c(2,10,30,90,300, 900)) +
  annotate('text', y = 100, label = 'Aging', x = 4.5, hjust = 0, size = 4/pntnorm) +
  annotate('text', y = 80, label = 'Development', x = 4.5, hjust = 1, size = 4/pntnorm) +
  ylab('Age') +
  xlab(NULL) +
  guides(color = F)  +
  theme_pubr(base_size = 8)
ages
# system('mkdir -p results')
ggsave('./results/ages.pdf',ages, units = 'cm', width = 8, height = 6, useDingbats = F)
ggsave('./results/ages.png',ages, units = 'cm', width = 8, height = 6)


all_raw_var= (pca_dat %>%
  filter(period == 'all', type == 'raw') %>%
  select(varExp, PC) %>%
  unique())$varExp

aging_raw_var= (pca_dat %>%
                filter(period == 'aging', type == 'raw') %>%
                select(varExp, PC) %>%
                unique())$varExp

dev_raw_var= (pca_dat %>%
                  filter(period == 'development', type == 'raw') %>%
                  select(varExp, PC) %>%
                  unique())$varExp

all_notissue_var= (pca_dat %>%
                filter(period == 'all', type == 'notissue') %>%
                select(varExp, PC) %>%
                unique())$varExp

aging_notissue_var= (pca_dat %>%
                  filter(period == 'aging', type == 'notissue') %>%
                  select(varExp, PC) %>%
                  unique())$varExp

dev_notissue_var= (pca_dat %>%
                filter(period == 'development', type == 'notissue') %>%
                select(varExp, PC) %>%
                unique())$varExp

all_raw_pca12 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2') +
  coord_fixed(ratio = all_raw_var[2]/all_raw_var[1], clip = 'off') +
  xlab(paste('PC1 (', round(all_raw_var[1]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(all_raw_var[2]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) +
  theme(legend.position = c(0.05,0.6),
        legend.justification=c(0,0),
        legend.direction = 'vertical',
        legend.box = 'horizontal',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) 

ggsave('./results/pca_all_raw12.pdf',all_raw_pca12, units = 'cm', width = 8, height = 5, useDingbats =F)
ggsave('./results/pca_all_raw12.png',all_raw_pca12, units = 'cm', width = 8, height = 5)

all_raw_pc1age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC1, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age (log2)')

all_raw_pc2age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC2, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age (log2)')

all_raw_pc12 = ggarrange(all_raw_pca12, ggarrange(all_raw_pc1age, all_raw_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')

ggsave('./results/pca_all_raw12_wage.pdf',all_raw_pc12, units = 'cm', width = 16, height = 6, useDingbats = F)
ggsave('./results/pca_all_raw12_wage.png',all_raw_pc12, units = 'cm', width = 16, height = 6)

all_raw_pca34 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2') +
  coord_fixed(ratio = all_raw_var[4]/all_raw_var[3], clip = 'off') +
  xlab(paste('PC3 (', round(all_raw_var[3]*100),'%)',sep='')) +
  ylab(paste('PC4 (', round(all_raw_var[4]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        legend.box = 'vertical',
        legend.background = element_rect(fill = 'gray85',color = 'gray25'))

ggsave('./results/pca_all_raw34.pdf',all_raw_pca34, units = 'cm', width = 8, height = 5, useDingbats =F)
ggsave('./results/pca_all_raw34.png',all_raw_pca34, units = 'cm', width = 8, height = 5)

all_raw_pc3age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC3, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age (log2)')

all_raw_pc4age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC4, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age (log2)')

all_raw_pc34 = ggarrange(all_raw_pca34, ggarrange(all_raw_pc3age, all_raw_pc4age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), legend = F, labels = c('a',NA))

ggsave('./results/pca_all_raw34_wage.pdf', all_raw_pc34, units = 'cm', width = 16, height = 6, useDingbats = F)
ggsave('./results/pca_all_raw34_wage.png', all_raw_pc34, units = 'cm', width = 16, height = 6)

all_raw_pc34_2 = ggarrange(all_raw_pca34, ggarrange(all_raw_pc3age, all_raw_pc4age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), legend = F, labels = c('b',NA), align = 'v')

all_raw_pca = ggarrange(all_raw_pc12,all_raw_pc34_2, align = 'hv', ncol = 1, nrow= 2)

ggsave('./results/all_raw_pca1234.pdf',all_raw_pca, units = 'cm',width = 13, height = 10,useDingbats = F)
ggsave('./results/all_raw_pca1234.png',all_raw_pca, units = 'cm',width = 13, height = 10)

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

ggsave('./results/pca_all_nt12.pdf',all_nt_pca12, units = 'cm', width = 8, height = 8, useDingbats =F)
ggsave('./results/pca_all_nt12.png',all_nt_pca12, units = 'cm', width = 8, height = 8)

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
  xlab('Age (log2)')

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
  xlab('Age (log2)')

all_nt_pc12 = ggarrange(all_nt_pca12, ggarrange(all_nt_pc1age, all_nt_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')

ggsave('./results/pca_all_nt12_wage.pdf',all_nt_pc12, units = 'cm', width = 16, height = 9,
       useDingbats = F)
ggsave('./results/pca_all_nt12_wage.png',all_nt_pc12, units = 'cm', width = 16, height = 9)

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

ggsave('./results/pca_dev_raw12.pdf',dev_raw_pca12, units = 'cm', width = 8, height = 5, useDingbats =F)
ggsave('./results/pca_dev_raw12.png',dev_raw_pca12, units = 'cm', width = 8, height = 5)

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
  xlab('Age (log2)')

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
  xlab('Age (log2)')

dev_raw_pc12 = ggarrange(dev_raw_pca12, ggarrange(dev_raw_pc1age, dev_raw_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')

ggsave('./results/pca_dev_raw12_wage.pdf',dev_raw_pc12, units = 'cm', width = 16, height = 6, useDingbats = F)
ggsave('./results/pca_dev_raw12_wage.png',dev_raw_pc12, units = 'cm', width = 16, height = 6)

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
  xlab(paste('PC1 (', round(dev_raw_var[3]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(dev_raw_var[4]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) +
  theme(legend.position = c(0.15,0.05),
        legend.justification=c(0,0),
        legend.direction = 'horizontal',
        legend.box = 'vertical',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) 

ggsave('./results/pca_dev_raw34.pdf',dev_raw_pca34, units = 'cm', width = 8, height = 5, useDingbats =F)
ggsave('./results/pca_dev_raw34.png',dev_raw_pca34, units = 'cm', width = 8, height = 5)

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
  xlab('Age (log2)')

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
  xlab('Age (log2)')

dev_raw_pc34 = ggarrange(dev_raw_pca34, ggarrange(dev_raw_pc3age, dev_raw_pc4age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')

ggsave('./results/pca_dev_raw34_wage.pdf',dev_raw_pc34, units = 'cm', width = 16, height = 6, useDingbats = F)
ggsave('./results/pca_dev_raw34_wage.png',dev_raw_pc34, units = 'cm', width = 16, height = 6)

######

dev_nt_pca12 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'notissue') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2',breaks = c(2,8,30,60)) +
  coord_fixed(ratio = dev_notissue_var[2]/dev_notissue_var[1], clip = 'off') +
  xlab(paste('PC1 (', round(dev_notissue_var[1]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(dev_notissue_var[2]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) +
  theme(legend.position = c(0.26,0.65),
        legend.justification=c(0,0),
        legend.direction = 'vertical',
        legend.box = 'horizontal',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) 

ggsave('./results/pca_dev_nt12.pdf', dev_nt_pca12, units = 'cm', width = 8, height = 5, useDingbats =F)
ggsave('./results/pca_dev_nt12.png', dev_nt_pca12, units = 'cm', width = 8, height = 5)

dev_nt_pc1age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'notissue') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC1, color = tissue)) +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age (log2)')

dev_nt_pc2age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'notissue') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC2, color = tissue)) +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age (log2)')

dev_nt_pc12 = ggarrange(dev_nt_pca12, ggarrange(dev_nt_pc1age, dev_nt_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')

ggsave('./results/pca_dev_nt12_wage.pdf',dev_nt_pc12, units = 'cm', width = 16, height = 6, useDingbats = F)
ggsave('./results/pca_dev_nt12_wage.png',dev_nt_pc12, units = 'cm', width = 16, height = 6)

dev_nt_pca34 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'notissue') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2',breaks = c(2,8,30,60)) +
  coord_fixed(ratio = dev_notissue_var[4]/dev_notissue_var[3], clip = 'off') +
  xlab(paste('PC1 (', round(dev_notissue_var[3]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(dev_notissue_var[4]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) 
  # theme(legend.position = c(0.15,0.05),
  #       legend.justification=c(0,0),
  #       legend.direction = 'horizontal',
  #       legend.box = 'vertical',
  #       legend.background = element_rect(fill = 'gray85',color = 'gray25')) 

ggsave('./results/pca_dev_nt34.pdf',dev_nt_pca34, units = 'cm', width = 8, height = 5, useDingbats =F)
ggsave('./results/pca_dev_nt34.png',dev_nt_pca34, units = 'cm', width = 8, height = 5)

dev_nt_pc3age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'notissue') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC3, color = tissue)) +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age (log2)')

dev_nt_pc4age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'development', type == 'notissue') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC4, color = tissue)) +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age (log2)')

dev_nt_pc34 = ggarrange(dev_nt_pca34, ggarrange(dev_nt_pc3age, dev_nt_pc4age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA))

ggsave('./results/pca_dev_nt34_wage.pdf',dev_nt_pc34, units = 'cm', width = 16, height = 6, useDingbats = F)
ggsave('./results/pca_dev_nt34_wage.png',dev_nt_pc34, units = 'cm', width = 16, height = 6)

aging_raw_pca12 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'aging', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2') +
  coord_fixed(ratio = aging_raw_var[2]/aging_raw_var[1], clip = 'off') +
  xlab(paste('PC1 (', round(aging_raw_var[1]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(aging_raw_var[2]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) +
  theme(legend.position = c(0.06,0.65),
        legend.justification=c(0,0),
        legend.direction = 'vertical',
        legend.box = 'horizontal',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) 

ggsave('./results/pca_aging_raw12.pdf',aging_raw_pca12, units = 'cm', width = 8, height = 5, useDingbats =F)
ggsave('./results/pca_aging_raw12.png',aging_raw_pca12, units = 'cm', width = 8, height = 5)

aging_raw_pc1age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'aging', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC1, color = tissue)) +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age (log2)')

aging_raw_pc2age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'aging', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC2, color = tissue)) +
  geom_smooth(alpha = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  guides(color = F) +
  xlab('Age (log2)')

aging_raw_pc12 = ggarrange(aging_raw_pca12, ggarrange(aging_raw_pc1age, aging_raw_pc2age, ncol = 1, nrow= 2), ncol = 2, nrow= 1, widths = c(1,0.5), common.legend = F,labels = c('a',NA), align = 'v')

ggsave('./results/pca_aging_raw12_wage.pdf',aging_raw_pc12, units = 'cm', width = 16, height = 6,
       useDingbats = F)
ggsave('./results/pca_aging_raw12_wage.png',aging_raw_pc12, units = 'cm', width = 16, height = 6)






####

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

top_divcon_gene_dat = cov_ch %>%
  select(-p,-FDR) %>%
  spread(key = period, value =`CoV_change`) %>%
  mutate(rev = aging *development) %>%
  top_n(n=1,wt=-(rev)) %>%
  left_join(expr) %>%
  inner_join(sample_info) 
top_divcon_gene = ggplot(top_divcon_gene_dat, aes(x =age, y= expression, color = tissue)) +
  geom_smooth(se=F,size = 0.7) +
  geom_point(size=0.5) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  scale_x_continuous(trans = 'log2') +
  scale_color_manual(values = tissuecol) +
  ggtitle(unique(top_divcon_gene_dat$gene_id)) +
  guides(color = guide_legend('Tissue', 
                             override.aes = list(size = 2))) +
  theme(legend.position = 'top',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) +
  xlab('Age (log2)') + ylab('Gene Expression')

ggsave('./results/top_divcon_gene.pdf',top_divcon_gene, units = 'cm', width = 8, height = 8, useDingbats = F)
ggsave('./results/top_divcon_gene.png',top_divcon_gene, units = 'cm', width = 8, height = 8)

top_divcon_gene_cov = top_divcon_gene_dat %>%
  group_by(gene_id, ind_id, age) %>%
  summarise(sd = sd(expression),
            mean = mean(expression),
            cov = sd/mean) %>%
  ggplot(aes(x = age, y= cov)) +
  geom_smooth(size = 0.7) +
  geom_point(size=0.5) +
  geom_vline(xintercept = 90, linetype = 'dashed', color = 'gray35') +
  scale_x_continuous(trans = 'log2') +
  ggtitle(unique(top_divcon_gene_dat$gene_id)) +
  xlab('Age (log2)') + ylab('CoV')

ggsave('./results/top_divcon_gene_cov.pdf',top_divcon_gene_cov, units = 'cm', width = 8, height = 8,
       useDingbats = F)
ggsave('./results/top_divcon_gene_cov.png',top_divcon_gene_cov, units = 'cm', width = 8, height = 8)

top_divcon_gene_p = ggarrange(top_divcon_gene + ggtitle(NULL), top_divcon_gene_cov + ggtitle(NULL),
                              labels = 'auto', ncol=2, nrow=1, align = 'hv')

ggsave('./results/top_divcon_gene_expcov.pdf',top_divcon_gene_p, units = 'cm', width = 10, height = 5,
       useDingbats = F)
ggsave('./results/top_divcon_gene_expcov.png',top_divcon_gene_p, units = 'cm', width = 10, height = 5)

cov_dat_sum = cov_dat %>%
  mutate(ind_id = factor(ind_id)) %>%
  left_join(unique(select(sample_info,-tissue,-sample_id))) %>%
  group_by(ind_id, age) %>%
  summarise(meanCoV = mean(CoV)) %>% 
  mutate(period = c('aging','development')[1+(age<=90)])

cov_cordat = group_by(ungroup(cov_dat_sum),period) %>%
  summarise(cor = cor.test(meanCoV, age, method = 's')$est,
            cor.p = cor.test(meanCoV, age, method = 's')$p.val)

cov_mean = ggplot(cov_dat_sum, aes(x = age, y= meanCoV)) +
  geom_point() +
  geom_smooth(se=T,color = 'midnightblue') +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept = 90, linetype='dashed',color = 'gray35') +
  xlab('Age (log2)') + ylab('Mean CoV') +
  annotate('text',x=95,y=0.5075,label='Aging',hjust=0,size = 8/pntnorm) +
  annotate('text',x = 95, y=0.505, label = parse(text = paste('rho["CoV,age"] ==' ,(round(filter(cov_cordat,period=='aging')$cor,3)))),hjust=0,size = 8/pntnorm) +
  annotate('text',x = 95, y=0.5025, label = parse(text = paste0('p ==' ,(round(filter(cov_cordat,period=='aging')$cor.p,2)))),hjust=0,size = 8/pntnorm) +
  annotate('text',x=8,y=0.4825,label='Development',hjust=0,size = 8/pntnorm) +
  annotate('text',x = 8, y=0.48, label = parse(text = paste('rho["CoV,age"] ==' ,(round(filter(cov_cordat,period=='development')$cor,3)))),hjust=0,size = 8/pntnorm) +
  annotate('text',x = 8, y=0.4775, label = parse(text = paste0('p ==' ,(round(filter(cov_cordat,period=='development')$cor.p,2)))),hjust=0,size = 8/pntnorm)

ggsave('./results/cov_mean.pdf',cov_mean, units = 'cm', width = 8, height = 8, useDingbats = F)
ggsave('./results/cov_mean.png',cov_mean, units = 'cm', width = 8, height = 8)
