## other datasets
# gtex
library(tidyverse)
pntnorm <- (1/0.352777778)
library(ggpubr)
library(grid)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))

cov = readRDS('/home/black/Dropbox/projects/ageing/proc_data/gtex/cov.rds')
age = readRDS('/home/black/Dropbox/projects/ageing/proc_data/gtex/age.rds')
pcors = readRDS('/home/black/Dropbox/projects/ageing/proc_data/gtex/pwisecors.rds')

meansd = colMeans(cov)
meancor = cor.test(meansd, age, m= "s")
mediansd = apply(cov,2,median)
mediancor = cor.test(mediansd, age, m= "s")

gtex = data.frame(mcov = meansd, mdcov = mediansd, ind_id = names(meansd), row.names = NULL) %>%
  left_join(data.frame(age = age, ind_id=names(age)))

gtexmp = gtex %>%
  ggplot(aes(x=age, y=mcov)) +
  geom_point() +
  xlab('Age in years') +
  ylab('Mean Cov') +
  ggtitle('GTEX')+
  geom_smooth(method='lm', color='midnightblue') +
  annotate('text', x=30, y= 0.28, label=  parse(text = paste('rho["CoV,age"] ==' ,round(meancor$est,3))),
           size = 8/pntnorm) +
  annotate('text', x =30, y=0.275, label = parse(text = paste0('p ==' ,round(meancor$p.val,2))),
           size = 8/pntnorm)

gtexmdp = gtex %>%
  ggplot(aes(x=age, y=mdcov)) +
  geom_point() +
  xlab('Age in years') +
  ylab('Median Cov') +
  geom_smooth(method='lm', color='midnightblue') +
  annotate('text', x=30, y= 0.24, label=  parse(text = paste('rho["CoV,age"] ==' ,round(mediancor$est,3))),
           size = 8/pntnorm) +
  annotate('text', x =30, y=0.235, label = parse(text = paste0('p ==' ,round(mediancor$p.val,2))),
           size = 8/pntnorm)


gtexp = ggarrange(gtexmp, gtexmdp, ncol = 2, align= 'hv')

covJ = readRDS('/home/black/Dropbox/projects/ageing/proc_data/jonker/cov.rds')
ageJ = readRDS('/home/black/Dropbox/projects/ageing/proc_data/jonker/age.rds')
pcorsJ = readRDS('/home/black/Dropbox/projects/ageing/proc_data/gtex/pwisecors.rds')

meansdJ = colMeans(covJ)
meancorJ = cor.test(meansdJ, ageJ, m= "s")
mediansdJ = apply(covJ, 2, median)
mediancorJ = cor.test(mediansdJ, ageJ, m= "s")

jonker = data.frame(mcov = meansdJ, mdcov = mediansdJ, ind_id = names(meansdJ), row.names = NULL) %>%
  left_join(data.frame(age = ageJ, ind_id=names(ageJ)))

jonkermp = jonker %>%
  ggplot(aes(x=age, y=mcov)) +
  geom_point() +
  scale_x_continuous(trans='log2') +
  xlab('Age in days (log2 scale)') +
  ylab('Mean Cov') +
  ggtitle('Jonker')+
  geom_smooth(method='lm', color='midnightblue') +
  annotate('text', x=128, y= 0.152, label=  parse(text = paste('rho["CoV,age"] ==' ,round(meancorJ$est,3))),
           size = 8/pntnorm) +
  annotate('text', x =128, y=0.1515, label = parse(text = paste0('p ==' ,round(meancorJ$p.val,2))),
           size = 8/pntnorm)

jonkermdp = jonker %>%
  ggplot(aes(x=age, y=mdcov)) +
  geom_point() +
  scale_x_continuous(trans='log2') +
  xlab('Age in days (log2 scale)') +
  ylab('Median Cov') +
  geom_smooth(method='lm', color='midnightblue') +
  annotate('text', x=128, y= 0.106, label=  parse(text = paste('rho["CoV,age"] ==' ,round(mediancorJ$est,3))),
           size = 8/pntnorm) +
  annotate('text', x =128, y=0.1055, label = parse(text = paste0('p ==' ,round(mediancorJ$p.val,2))),
           size = 8/pntnorm)

jonkerp = ggarrange(jonkermp, jonkermdp, ncol = 2, align= 'hv')
p = ggarrange(gtexp, jonkerp, nrow =2, labels='auto')

ggsave('results/other_datasets/gtex_jonker_meanCoV.pdf',  p, units= 'cm',  width = 16, height = 12, useDingbats = F)
ggsave('results/other_datasets/gtex_jonker_meanCoV.png', p, units= 'cm',width = 16, height = 12)


##pwise cors GTEX:
#### pairwise correlations:
pcors
agecors = sapply(names(pcors), function(x){
  c(round(cor.test(pcors[[x]], age, m='s')$est,2),
    round(cor.test(pcors[[x]], age, m='s')$p.val,3))
})
rownames(agecors)[2] = 'p.val'

pexpcors = reshape2::melt(pcors) %>%
  set_names('rho', 'pair') %>%
  #mutate(id = rep(sp, 6) ) %>%
  mutate(age = rep(age, 6))

annottext = reshape2::melt(agecors) %>%
  set_names('stat', 'pair', 'value') 

pwisecors = ggplot(pexpcors) +
  aes(x = age, y = rho) +
  facet_wrap(~pair,scales = 'free') +
  geom_point(color = adjustcolor('gray30', alpha.f = 0.8)) +
  scale_x_continuous(trans = 'log2') +
  geom_smooth(method = 'lm', color = 'midnightblue') +
  ylab('Sperman correlation coefficient') + xlab('Age in years (log2 scale)') +
  geom_text(data = filter(annottext, stat == 'rho'), vjust = 2,
            mapping = aes(x = 30, y= Inf, label = paste('rho==',value) ), parse =  T,size = 7/pntnorm) +
  geom_text(data = filter(annottext,  stat == 'p.val'), vjust = 3,
            mapping = aes(x = 30, y= Inf, label = paste('p==',value) ), parse =  T,size = 7/pntnorm) +
  theme(strip.background = element_rect(fill='gray80'))

ggsave('results/gtex/pwisecor_gtex.pdf', pwisecors, units = 'cm', width = 16, height = 10, useDingbats =F)
ggsave('results/gtex/pwisecor_gtex.png', pwisecors, units = 'cm', width = 16, height = 10)


#### jonker

pcorsJ = readRDS('/home/black/Dropbox/projects/ageing/proc_data/jonker/pwisecors.rds')
ind.id = readRDS('/home/black/Dropbox/projects/ageing/proc_data/jonker/ind.id.rds')

agecorsJ = sapply(names(pcorsJ), function(x){
  c(round(cor.test(pcorsJ[[x]], ageJ, m='s')$est,2),
    round(cor.test(pcorsJ[[x]], ageJ, m='s')$p.val,3))
})
rownames(agecorsJ)[2] = 'p.val'

pexpcorsJ = reshape2::melt(pcorsJ) %>%
  set_names('rho', 'pair') %>%
  mutate(id = rep(ind.id,10) ) %>%
  left_join(data.frame(ageJ, id = ind.id) ) 

annottextJ = reshape2::melt(agecorsJ) %>%
  set_names('stat', 'pair', 'value') 

pwisecorsJ = ggplot(pexpcorsJ) +
  aes(x = ageJ, y = rho) +
  facet_wrap(~pair,scales = 'free') +
  geom_point(color = adjustcolor('gray30', alpha.f = 0.8)) +
  scale_x_continuous(trans = 'log2') +
  geom_smooth(method ='lm', color = 'midnightblue') +
  ylab('Sperman correlation coefficient') + xlab('Age in days (log2 scale)') +
  geom_text(data = filter(annottextJ, stat == 'rho'), vjust = 2,
            mapping = aes(x = 240, y= Inf, label = paste('rho==',value) ), parse =  T,size = 7/pntnorm) +
  geom_text(data = filter(annottextJ,  stat == 'p.val'), vjust = 3,
            mapping = aes(x = 240, y= Inf, label = paste('p==',value) ), parse =  T,size = 7/pntnorm) 


ggsave('results/other_datasets/pwisecor_jonker.pdf',pwisecorsJ, units = 'cm', width = 16, height = 12, useDingbats =F)

ggsave('results/other_datasets/pwisecor_jonker.png',pwisecorsJ, units = 'cm', width = 16, height = 12)









