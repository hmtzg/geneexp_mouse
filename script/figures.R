library(tidyverse)
library(ggpubr)
library(ggrepel)
theme_set(theme_pubr(base_size = 10, legend = 'bottom'))
pntnorm <- (1/0.352777778)
tissuecol = setNames(c("#233789", "#f49e92", "#801008","#dbb32e"),c('Cortex','Lung','Liver','Muscle'))
regcol = setNames(c('#8a6578','#7497b5'),c('up','down'))
varcol = setNames(c('#52233B','#234361'),c('div','con'))

sample_info = readRDS('./data/processed/figures/tidy/sample_info.rds')
expr = readRDS('./data/processed/figures/tidy/expression.rds')
expr_ch = readRDS('./data/processed/figures/tidy/expression_change.rds')
pca_dat = readRDS('./data/processed/figures/tidy/pca_data.rds')
cov_dat = readRDS('./data/processed/figures/tidy/CoV.rds')
cov_ch = readRDS('./data/processed/figures/tidy/CoV_change.rds')
cov_gsea = readRDS('./data/processed/figures/tidy/CoV_GSEA.rds')
divcon_gsea = readRDS('./data/processed/figures/tidy/divcon_GSEA.rds')

sample_info %>%
  ggplot(aes(y = age, x= tissue, color = tissue)) +
  geom_jitter(width = 0.1, size = 3) +
  coord_flip() + 
  scale_color_manual(values = tissuecol) +
  scale_y_continuous(trans = 'log2')


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
  guides(color = F)  +
  theme_pubr(base_size = 8)
ages_log2
# system('mkdir -p results')
ggsave('./results/ages_log2.pdf',ages_log2, units = 'cm', width = 8, height = 6, useDingbats = F)
ggsave('./results/ages_log2.png',ages_log2, units = 'cm', width = 8, height = 6)

ages= sample_info %>%
  ggplot(aes(y = age, x= tissue, color = tissue)) +
  geom_hline(yintercept = 90, linetype = 'dashed', color = 'gray35') +
  geom_jitter(width = 0.1, size = 1.5) + 
  coord_flip() + 
  scale_color_manual(values = tissuecol) +
  # scale_y_continuous(trans = 'log2', breaks = c(2,10,30,90,300, 900)) +
  annotate('text', y = 100, label = 'Aging', x = 4.5, hjust = 0, size = 8/pntnorm) +
  annotate('text', y = 80, label = 'Development', x = 4.5, hjust = 1, size = 8/pntnorm) +
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

all_raw_pca = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,4), trans= 'log2') +
  coord_fixed(ratio = all_raw_var[2]/all_raw_var[1]) +
  xlab(paste('PC1 (', round(all_raw_var[1]*100),'%)',sep=''))+
  ylab(paste('PC2 (', round(all_raw_var[2]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) +
  theme(legend.position = 'right')

all_raw_pc1age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC1, color = tissue)) +
  geom_point(size = 0.5)+
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  geom_smooth(se=F) +guides(color = F)

all_raw_pc2age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC2, color = tissue)) +
  geom_point(size = 0.5)+
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  geom_smooth(se=F) +guides(color = F) 

all_raw_pc12 = ggarrange(all_raw_pca, ggarrange(all_raw_pc1age, all_raw_pc2age), ncol = 1, nrow= 2, heights = c(1,0.5), common.legend = T)

all_raw_pca34 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,4), trans= 'log2') +
  coord_fixed(ratio = all_raw_var[4]/all_raw_var[3]) +
  xlab(paste('PC3 (', round(all_raw_var[3]*100),'%)',sep=''))+
  ylab(paste('PC4 (', round(all_raw_var[4]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'), 
         size = guide_legend('Age')) +
  theme(legend.position = 'right')

all_raw_pc3age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC3, color = tissue)) +
  geom_point(size = 0.5)+
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  geom_smooth(se=F) +guides(color = F)

all_raw_pc4age = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  ggplot(aes(x = age, y = PC4, color = tissue)) +
  geom_point(size = 0.5)+
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  geom_smooth(se=F) +guides(color = F) 

all_raw_pc34 = ggarrange(all_raw_pca34, ggarrange(all_raw_pc3age, all_raw_pc4age), ncol = 1, nrow= 2, heights = c(1,0.5), common.legend = T)

all_raw_pca = ggarrange(all_raw_pc12,all_raw_pc34, align = 'hv', common.legend = T, legend = 'top')

ggsave('./results/all_raw_pca1234.pdf',all_raw_pca, units = 'cm',width = 14, height = 8,useDingbats = F)

expr_ch %>%
  select(-p,-FDR) %>%
  spread(key = period, value =`Expression Change`) %>%
  mutate(rev = aging *development) %>%
  top_n(n=1,wt=-rev) %>%
  left_join(expr) %>%
  inner_join(sample_info) %>%
  ggplot(aes(x =age, y= expression, color = tissue)) +
  geom_smooth() +
  geom_point() +
  scale_x_continuous(trans = 'log2') +
  scale_color_manual(values = tissuecol)


expr_ch %>%
  group_by(tissue,period) %>%
  summarise( up = sum(`Expression Change`>0,na.rm=T),
             down = sum(`Expression Change`<0,na.rm=T)) %>%
  gather(key= 'direction',value ='n',-tissue,-period) %>%
  mutate(period = factor(period, levels = c('development','aging'))) %>%
  ggplot(aes(x = period, y= n, fill = direction)) +
  facet_wrap(~tissue, ncol=4) +
  geom_bar(stat= 'identity', position ='fill') +
  scale_fill_manual(values = regcol) +
  geom_hline(yintercept = 0.5)

expr_ch %>%
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
  geom_text(aes(label = n), color='white', position = position_dodge(width = 1), angle=90, hjust = 1, vjust = 0.5)

cov_ch %>%
  select(-p,-FDR) %>%
  spread(key = period, value =`CoV_change`) %>%
  mutate(rev = aging *development) %>%
  top_n(n=10,wt=-rev) %>%
  left_join(expr) %>%
  inner_join(sample_info) %>%
  ggplot(aes(x =age, y= expression, color = tissue)) +
  geom_smooth(se=F) +
  # geom_point() +
  geom_vline(xintercept = 90)+
  scale_x_continuous(trans = 'log2') +
  scale_color_manual(values = tissuecol) +
  facet_wrap(~gene_id, scales = 'free')


cov_dat %>%
  mutate(ind_id = factor(ind_id)) %>%
  left_join(unique(select(sample_info,-tissue,-sample_id))) %>%
  group_by(ind_id, age) %>%
  summarise(meanCoV = mean(CoV),
            medCoV = median(CoV),
            q3CoV = quantile(CoV,0.75),
            q1CoV = quantile(CoV,0.25)) %>%
  ggplot(aes(x = age, y= meanCoV)) +
  # geom_pointrange(aes(ymin = q1CoV, ymax = q3CoV)) +
  geom_point() +
  geom_smooth(se=F) +
  scale_x_continuous(trans = 'log2') +
  geom_vline(xintercept = 90)

cov_gsea %>%
  filter(p.adj <= 0.1) %>%
  select(Description) %>%
  left_join(cov_gsea) %>%
  select(Description, NES, p.adj, period) %>%
  head()
  ggplot(aes(x = period, y= NES, fill = p.adj<=0.1)) +
  geom_boxplot()


