library(tidyverse)
#library(ggthemes)
library(ggpubr)
#library(Rvislib)
#theme_set(theme_rvis(base_size = 6))
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
tsspec_es = readRDS('./data/processed/tidy/ts.spec.ES.rds') %>%
  set_names(c('geneid','tissue','ES','spec'))
tsspec_rho = readRDS('./data/processed/tidy/ts.spec.exprho.rds') %>%
  set_names(c('geneid','tissue','rho','period','spec'))
# load('./data/processed/figures/raw/exp.rdata')

tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('Up','Down'))
revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))
# revcol = setNames(c('gray35', 'gray35', 'gray75', 'gray75'), c('UpDown','DownUp','UpUp','DownDown'))

specdat = full_join(tsspec_es,tsspec_rho) %>%
  mutate(period = factor(ifelse(period == 'Development', 'Development', 'Ageing'), 
                         levels = c('Development','Ageing')))

revdat = specdat %>%
  spread(period,rho) %>%
  mutate(dev = c('Down','Up')[1+(Development>0)], age = c('Down','Up')[1+(Ageing>0)]) %>%
  mutate(rev = paste(dev,age,sep=''))

specdat = left_join(specdat,revdat)

cortdat = specdat %>%
  filter(spec == 'Cortex')
n = length(unique(cortdat$geneid))
q1 = function(x)quantile(x,0.25,na.rm = T)
q3 = function(x)quantile(x,0.75,na.rm = T)
cort_sp = cortdat %>%
  select(tissue,spec, ES,geneid)%>%
  unique() %>%
  ggplot(aes(x = tissue , y = ES, color = tissue)) +
  # geom_violin(trim = T) +
  # geom_boxplot(outlier.size = 0.1, size = 0.1, outlier.shape = NA) +
  stat_summary(geom = 'linerange', fun.min = q1, fun.max = q3, size = 0.5) +
  stat_summary(geom = 'point', fun = median, size = 1) +
  scale_color_manual(values = tissuecol) +
  guides(color = F) +
  xlab(NULL) + ylab('Effect Size') +
  ggtitle('Cortex Specific Genes') +
  # ylim(-2,20) +
  annotate(geom='text',label = paste('n = ',n,sep=''), x = 4, y = 7.5, hjust = 1)
ggsave('./results/figure3/cort_sp.pdf',units = 'cm', height = 2.5, width = 6, useDingbats=F)
ggsave('./results/figure3/cort_sp.png',units = 'cm', height = 2.5, width = 6)

cortnums = filter(cortdat,tissue =='Cortex') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
cort_cort = ggplot(filter(cortdat,tissue == 'Cortex'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Cortex'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.1) + 
  geom_line(alpha = 0.2, aes(group=geneid, color = rev)) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #              fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab('Change in Expression')+
  ggtitle(paste('Cortex, DU:',downup,'%',' UD:',updown,'%',sep=''))
ggsave('./results/figure3/cort_cort.pdf',cort_cort,units = 'cm', height = 5, width = 6, useDingbats=F)
ggsave('./results/figure3/cort_cort.png',cort_cort,units = 'cm', height = 5, width = 6)

cortnums = filter(cortdat,tissue =='Lung') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
cort_lung = ggplot(filter(cortdat,tissue == 'Lung'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Lung'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev), size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #              fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Lung, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/cort_lung.pdf',cort_lung, units = 'cm', height = 2.5, width = 3.1, useDingbats=F)
ggsave('./results/figure3/cort_lung.pdf',cort_lung, units = 'cm', height = 2.5, width = 3.1)

cortnums = filter(cortdat,tissue =='Liver') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
cort_liver = ggplot(filter(cortdat,tissue == 'Liver'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Liver'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #              fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Liver, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/cort_liver.pdf',cort_liver,units = 'cm', height = 2.5, width = 3.1, useDingbats=F)
ggsave('./results/figure3/cort_liver.pdf',cort_liver,units = 'cm', height = 2.5, width = 3.1)

cortnums = filter(cortdat,tissue =='Muscle') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
cort_muscle = ggplot(filter(cortdat,tissue == 'Muscle'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Muscle'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #              fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Muscle, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/cort_muscle.pdf',cort_muscle,units = 'cm', height = 2.5, width = 3.1, useDingbats=F)
ggsave('./results/figure3/cort_muscle.pdf',cort_muscle,units = 'cm', height = 2.5, width = 3.1)

cort_c1 = ggarrange(cort_sp,cort_cort,ncol=1,nrow=2,heights = c(1,2))
cort_c2 = ggarrange(cort_liver,cort_lung,cort_muscle,ncol=1,nrow=3)
cort_p = ggarrange(cort_c1,cort_c2,ncol=2,nrow=1,widths = c(2,1))
ggsave('./results/figure3/cort.pdf',cort_p,units = 'cm',height = 7.5,width = 9, useDingbats=F)
ggsave('./results/figure3/cort.png',cort_p,units = 'cm',height = 7.5,width = 9)

###################
cortdat = specdat %>%
  filter(spec == 'Liver')
n = length(unique(cortdat$geneid))
liver_sp = cortdat %>%
  select(tissue,spec, ES,geneid)%>%
  unique() %>%
  ggplot(aes(x = tissue , y = ES, color = tissue)) +
  # geom_violin(trim = T) +
  # geom_boxplot(outlier.size = 0.1, size = 0.1, outlier.shape = NA) +
  stat_summary(geom = 'linerange', fun.min = q1, fun.max = q3, size = 0.5) +
  stat_summary(geom = 'point', fun = median, size = 1) +
  scale_color_manual(values = tissuecol) +
  guides(color = F) +
  xlab(NULL) +ylab('Effect Size') +
  ggtitle('Liver Specific Genes') +
  annotate(geom='text',label = paste('n = ',n,sep=''), x = 4, y = 6, hjust = 1)
ggsave('./results/figure3/liver_sp.pdf',liver_sp,units = 'cm', height = 2.5, width = 6, useDingbats=F)
ggsave('./results/figure3/liver_sp.png',liver_sp,units = 'cm', height = 2.5, width = 6)

cortnums = filter(cortdat,tissue =='Liver') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
liver_liver = ggplot(filter(cortdat,tissue == 'Liver'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Liver'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.1, outlier.shape = NA) + 
  geom_line(alpha = 0.2, aes(group=geneid, color = rev)) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #              fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab('Change in Expression')+
  ggtitle(paste('Liver, DU:',downup,'%',' UD:',updown,'%',sep=''))
ggsave('./results/figure3/liver_liver.pdf',liver_liver,units = 'cm', height = 5, width = 6, useDingbats=F)
ggsave('./results/figure3/liver_liver.png',liver_liver,units = 'cm', height = 5, width = 6)

cortnums = filter(cortdat,tissue =='Cortex') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
liver_cort = ggplot(filter(cortdat,tissue == 'Cortex'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Cortex'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #              fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Cortex, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/liver_cortex.pdf',liver_cort, units = 'cm', height = 2.5, width = 3.1, useDingbats=F)
ggsave('./results/figure3/liver_cortex.png',liver_cort, units = 'cm', height = 2.5, width = 3.1)

cortnums = filter(cortdat,tissue =='Lung') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
liver_lung = ggplot(filter(cortdat,tissue == 'Lung'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Lung'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #              fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Lung, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/liver_lung.pdf',liver_lung, units = 'cm', height = 2.5, width = 3.1, useDingbats=F)
ggsave('./results/figure3/liver_lung.png',liver_lung, units = 'cm', height = 2.5, width = 3.1)

cortnums = filter(cortdat,tissue =='Muscle') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
liver_muscle = ggplot(filter(cortdat,tissue == 'Muscle'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Muscle'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #              fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Muscle, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/liver_muscle.pdf',liver_muscle, units = 'cm', height = 2.5,width = 3.1,useDingbats=F)
ggsave('./results/figure3/liver_muscle.png',liver_muscle, units = 'cm', height = 2.5, width = 3.1)

liver_c1 = ggarrange(liver_sp,liver_liver,ncol=1,nrow=2,heights = c(1,2))
liver_c2 = ggarrange(liver_cort,liver_lung,liver_muscle,ncol=1,nrow=3)
liver_p = ggarrange(liver_c1,liver_c2,ncol=2,nrow=1,widths = c(2,1))
ggsave('./results/figure3/liver.pdf',liver_p,units = 'cm',height = 7.5,width = 9, useDingbats=F)
ggsave('./results/figure3/liver.png',liver_p,units = 'cm',height = 7.5,width = 9)

###################
cortdat = specdat %>%
  filter(spec == 'Lung')
n = length(unique(cortdat$geneid))
lung_sp = cortdat %>%
  select(tissue,spec, ES,geneid)%>%
  unique() %>%
  ggplot(aes(x = tissue , y = ES, color = tissue)) +
  # geom_violin(trim = T) +
  # geom_boxplot(outlier.size = 0.1, size = 0.1, outlier.shape = NA) +
  stat_summary(geom = 'linerange', fun.min = q1, fun.max = q3, size = 0.5) +
  stat_summary(geom = 'point', fun = median, size = 1) +
  scale_color_manual(values = tissuecol) +
  guides(color = F) +
  xlab(NULL) +ylab('Effect Size') +
  ggtitle('Lung Specific Genes') +
  annotate(geom='text',label = paste('n = ',n,sep=''), x = 1, y = 3, hjust = 0)
ggsave('./results/figure3/lung_sp.pdf',lung_sp,units = 'cm', height = 2.5, width = 6, useDingbats=F)
ggsave('./results/figure3/lung_sp.png',lung_sp,units = 'cm', height = 2.5, width = 6)

cortnums = filter(cortdat,tissue =='Lung') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
lung_lung = ggplot(filter(cortdat,tissue == 'Lung'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Lung'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.1, outlier.shape = NA) + 
  geom_line(alpha = 0.2, aes(group=geneid, color = rev)) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab('Change in Expression')+
  ggtitle(paste('Lung, DU:',downup,'%',' UD:',updown,'%',sep=''))
ggsave('./results/figure3/lung_lung.pdf',lung_lung,units = 'cm', height = 5, width = 6, useDingbats=F)
ggsave('./results/figure3/lung_lung.png',lung_lung,units = 'cm', height = 5, width = 6)

cortnums = filter(cortdat,tissue =='Cortex') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
lung_cortex = ggplot(filter(cortdat,tissue == 'Cortex'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Cortex'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Cortex, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/lung_cortex.pdf',lung_cortex, units = 'cm', height = 2.5, width = 3.1, useDingbats=F)
ggsave('./results/figure3/lung_cortex.png',lung_cortex, units = 'cm', height = 2.5, width = 3.1)

cortnums = filter(cortdat,tissue =='Liver') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
lung_liver = ggplot(filter(cortdat,tissue == 'Liver'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Liver'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Liver, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/lung_liver.pdf',lung_liver, units = 'cm', height = 2.5, width = 3.1, useDingbats=F)
ggsave('./results/figure3/lung_liver.png',lung_liver, units = 'cm', height = 2.5, width = 3.1)

cortnums = filter(cortdat,tissue =='Muscle') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
lung_muscle = ggplot(filter(cortdat,tissue == 'Muscle'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Muscle'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Muscle, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/lung_muscle.pdf',lung_muscle, units = 'cm', height = 2.5, width = 3.1, useDingbats=F)
ggsave('./results/figure3/lung_muscle.png',lung_muscle, units = 'cm', height = 2.5, width = 3.1)

lung_c1 = ggarrange(lung_sp,lung_lung,ncol=1,nrow=2,heights = c(1,2))
lung_c2 = ggarrange(lung_cortex,lung_liver,lung_muscle,ncol=1,nrow=3)
lung_p = ggarrange(lung_c1,lung_c2,ncol=2,nrow=1,widths = c(2,1))
ggsave('./results/figure3/lung.pdf',lung_p,units = 'cm',height = 7.5,width = 9, useDingbats=F)
ggsave('./results/figure3/lung.png',lung_p,units = 'cm',height = 7.5,width = 9)

###################
cortdat = specdat %>%
  filter(spec == 'Muscle')
n = length(unique(cortdat$geneid))
muscle_sp = cortdat %>%
  select(tissue,spec, ES,geneid)%>%
  unique() %>%
  ggplot(aes(x = tissue , y = ES, color = tissue)) +
  # geom_violin(trim = T) +
  # geom_boxplot(outlier.size = 0.1, size = 0.1, outlier.shape = NA) +
  stat_summary(geom = 'linerange', fun.min = q1, fun.max = q3, size = 0.5) +
  stat_summary(geom = 'point', fun = median, size = 1) +
  scale_color_manual(values = tissuecol) +
  guides(color = F) +
  xlab(NULL) +ylab('Effect Size') +
  ggtitle('Muscle Specific Genes') +
  annotate(geom='text',label = paste('n = ',n,sep=''), x = 1, y = 2.5, hjust = 0)
ggsave('./results/figure3/muscle_sp.pdf',muscle_sp,units = 'cm', height = 2.5, width = 6, useDingbats=F)
ggsave('./results/figure3/muscle_sp.png',muscle_sp,units = 'cm', height = 2.5, width = 6)

cortnums = filter(cortdat,tissue =='Muscle') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
muscle_muscle = ggplot(filter(cortdat,tissue == 'Muscle'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Muscle'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.1, outlier.shape = NA) + 
  geom_line(alpha = 0.2, aes(group=geneid, color = rev)) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab('Change in Expression')+
  ggtitle(paste('Muscle, DU:',downup,'%',' UD:',updown,'%',sep=''))
ggsave('./results/figure3/muscle_muscle.pdf',muscle_muscle,units = 'cm', height = 5, width = 6, useDingbats=F)
ggsave('./results/figure3/muscle_muscle.png',muscle_muscle,units = 'cm', height = 5, width = 6)

cortnums = filter(cortdat,tissue =='Cortex') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
muscle_cort = ggplot(filter(cortdat,tissue == 'Cortex'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Cortex'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Cortex, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/muscle_cortex.pdf',muscle_cort, units = 'cm',height = 2.5,width = 3.1,useDingbats=F)
ggsave('./results/figure3/muscle_cortex.png',muscle_cort, units = 'cm', height = 2.5, width = 3.1)


cortnums = filter(cortdat,tissue =='Liver') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
muscle_liver = ggplot(filter(cortdat,tissue == 'Liver'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Liver'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Liver, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/muscle_liver.pdf',muscle_liver, units = 'cm',height = 2.5, width = 3.1,useDingbats=F)
ggsave('./results/figure3/muscle_liver.png',muscle_liver, units = 'cm', height = 2.5, width = 3.1)

cortnums = filter(cortdat,tissue =='Lung') %>%
  select(geneid,rev) %>%
  group_by(rev) %>%
  summarise(n = length(unique(geneid)))
cortnums = setNames(cortnums$n, cortnums$rev)
downup = round(100 * cortnums['DownUp'] / (cortnums['DownUp'] + cortnums['DownDown']))
updown = round(100 * cortnums['UpDown'] / (cortnums['UpDown'] + cortnums['UpUp']))
muscle_lung = ggplot(filter(cortdat,tissue == 'Lung'), aes(x = period , y= rho)) +
  geom_hline(yintercept = 0, color = 'darkgray', linetype = 'dashed') +
  geom_violin(fill = tissuecol['Lung'], alpha = 0.5) +
  geom_boxplot(fill = 'white', width = 0.15, outlier.shape = NA) + 
  geom_line(alpha = 0.1, aes(group=geneid, color = rev),size = 0.1) +
  # stat_summary(color = tissuecol['cortex'], fun = median, 
  #fun.min = function(x)quantile(x,0.25), fun.max = function(x)quantile(x,0.75)) +
  scale_color_manual(values = revcol) +
  # guides(color = guide_legend(title = NULL,override.aes = list(alpha = 1, size = 1))) +
  guides(color =F) +
  xlab(NULL) +ylab(NULL)+
  ggtitle(paste('Lung, DU:',downup,'%','\nUD:',updown,'%',sep=''))
ggsave('./results/figure3/muscle_lung.pdf',muscle_lung, units = 'cm', height = 2.5, width = 3.1,useDingbats=F)
ggsave('./results/figure3/muscle_lung.png',muscle_lung, units = 'cm', height = 2.5, width = 3.1)

muscle_c1 = ggarrange(muscle_sp,muscle_muscle,ncol=1,nrow=2,heights = c(1,2))
muscle_c2 = ggarrange(muscle_cort,muscle_liver,muscle_lung,ncol=1,nrow=3)
muscle_p = ggarrange(muscle_c1,muscle_c2,ncol=2,nrow=1,widths = c(2,1))
ggsave('./results/figure3/muscle.pdf',muscle_p,units = 'cm',height = 7.5,width = 9, useDingbats=F)
ggsave('./results/figure3/muscle.png',muscle_p,units = 'cm',height = 7.5,width = 9)

all_p_rev = ggarrange(cort_p,liver_p,lung_p,muscle_p,ncol=2,nrow=2,labels = c('a.','b.','c.','d.'),
                      font.label = list(size = 10))
ggsave('./results/figure3/fig3_rev.pdf',all_p_rev,units = 'cm',width = 16.7,height = 15, useDingbats = F)
ggsave('./results/figure3/fig3_rev.png',all_p_rev,units = 'cm',width = 16.7,height = 15)

