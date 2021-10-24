library(tidyverse)
library(ggpubr)
library(RColorBrewer)
expr = readRDS('./data/processed/tidy/expression.rds')
expch = readRDS('./data/processed/tidy/expression_change.rds')
sinfo = readRDS('./data/processed/tidy/sample_info.rds')
expr = expr %>% left_join(sinfo)
sinfo
alph= 0.4
szx = 0.3
age2cor = expr %>%
  filter(sample_id%in%c('s16','s48')) %>%
  select(-ind_id, -log2age, -sample_id, -age ) %>%
  spread(key = 'tissue', value = 'expression') %>%
  ggplot(aes(x= Lung, y=Muscle)) +
  geom_point(size=szx, alpha=alph, color='lightsalmon4') +
  geom_smooth(method = 'lm', se=F, color='midnightblue', size= 0.81) +
  stat_cor(cor.coef.name = 'rho', method = 'spearman', ) + theme_bw() +
  ggtitle('Age: 2 days') + xlab('')
age93cor = expr %>%
  filter(sample_id%in%c('s24','s56')) %>%
  select(-ind_id, -log2age, -sample_id, -age ) %>%
  spread(key = 'tissue', value = 'expression') %>%
  ggplot(aes(x= Lung, y=Muscle)) +
  geom_point(size=szx, alpha=alph, color='lightsalmon4') +
  geom_smooth(method = 'lm', se=F, color='midnightblue', size= 0.81) +
  stat_cor(cor.coef.name = 'rho', method = 'spearman')+ theme_bw() +
  ggtitle('Age: 93 days') + ylab('')
age904cor = expr %>%
  filter(sample_id%in%c('s31','s63')) %>%
  select(-ind_id, -log2age, -sample_id, -age ) %>%
  spread(key = 'tissue', value = 'expression') %>%
  ggplot(aes(x= Lung, y=Muscle)) +
  geom_point(size=szx, alpha=alph, color='lightsalmon4') +
  geom_smooth(method = 'lm', se=F, color='midnightblue', size= 0.81) +
  stat_cor(cor.coef.name = 'rho', method = 'spearman')+ theme_bw() +
  ggtitle('Age: 904 days')  + ylab('') + xlab('')

exprcors = ggarrange(age2cor, age93cor, age904cor, ncol = 3)
exprcors
ggsave('/home/hmt/Dropbox/conferences/hibit2021/expcors.pdf', exprcors, units='cm', height = 10, width = 20,
       useDingbats=F)
ggsave('/home/hmt/Dropbox/conferences/hibit2021/expcors.png', exprcors, units='cm', height = 10, width = 20)


age2corX = expr %>%
  filter(sample_id%in%c('s16','s32')) %>%
  select(-ind_id, -log2age, -sample_id, -age ) %>%
  spread(key = 'tissue', value = 'expression') %>%
  ggplot(aes(x= Lung, y=Liver)) +
  geom_point(size=szx, alpha=alph, color='lightsalmon4') +
  geom_smooth(method = 'lm', se=F, color='midnightblue', size= 0.81) +
  stat_cor(cor.coef.name = 'rho', method = 'spearman', ) + theme_bw() +
  ggtitle('Age: 2 days') + xlab('')
age93corX = expr %>%
  filter(sample_id%in%c('s24','s40')) %>%
  select(-ind_id, -log2age, -sample_id, -age ) %>%
  spread(key = 'tissue', value = 'expression') %>%
  ggplot(aes(x= Lung, y=Liver)) +
  geom_point(size=szx, alpha=alph, color='lightsalmon4') +
  geom_smooth(method = 'lm', se=F, color='midnightblue', size= 0.81) +
  stat_cor(cor.coef.name = 'rho', method = 'spearman')+ theme_bw() +
  ggtitle('Age: 93 days') + ylab('')
age904corX = expr %>%
  filter(sample_id%in%c('s30','s46')) %>%
  select(-ind_id, -log2age, -sample_id, -age ) %>%
  spread(key = 'tissue', value = 'expression') %>%
  ggplot(aes(x= Lung, y=Liver)) +
  geom_point(size=szx, alpha=alph, color='lightsalmon4') +
  geom_smooth(method = 'lm', se=F, color='midnightblue', size= 0.81) +
  stat_cor(cor.coef.name = 'rho', method = 'spearman')+ theme_bw() +
  ggtitle('Age: 904 days')  + ylab('') + xlab('')

exprcorsX = ggarrange(age2corX, age93corX, age904corX, ncol = 3)
exprcorsX
ggsave('/home/hmt/Dropbox/conferences/hibit2021/expcors2.pdf', exprcorsX, units='cm', height = 10, width = 20,
       useDingbats=F)
ggsave('/home/hmt/Dropbox/conferences/hibit2021/expcors2.png', exprcorsX, units='cm', height = 10, width = 20)

#
revg = readRDS('./data/processed/tidy/revgenes.tissues.rds')
du = revg %>% filter(direction=='DownUp' & tissue=='Cortex')
expch %>%
  filter(tissue=='Cortex') %>%
  filter(gene_id%in%du$gene_id) %>% 
  group_by(tissue) %>%
  top_n(`Expression Change`)

expr %>%
  filter(gene_id%in%'ENSMUSG00000025900' & tissue=='Cortex') %>%
  left_join(sinfo) %>%
  ggplot(aes(x=age, y=expression)) +
  geom_point() +
  geom_smooth() +
  scale_x_continuous(trans='log2') +
  geom_vline(xintercept = 90, linetype='dashed')








