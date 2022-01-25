library(tidyverse)
library(ggpubr)
library(ggforce)
library(RColorBrewer)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)

sexcolors <- setNames(c('#ffc9b5', '#8fb8de'), c('Female', 'Male'))
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('Up','Down'))
# revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))
revcol = setNames(c('gray25', 'gray25', 'gray25', 'gray25'), c('UpDown','DownUp','UpUp','DownDown'))

#### Functions ####

plotsave <- function(ggobj, prefix, width, height, ...){
  path = strsplit(prefix,'/')[[1]]
  path = paste(path[-length(path)], collapse = '/')
  if(!file.exists(path)){
    system(paste('mkdir -p',path))
  }
  saveRDS(object = ggobj, file = paste(prefix,'.rds',sep=''))
  ggsave(file = paste(prefix, '.pdf', sep = ''), plot = ggobj, units = 'cm', 
         width = width, height = height, useDingbats = F, limitsize = F)
  ggsave(file = paste(prefix, '.png', sep = ''), plot = ggobj, units = 'cm', bg='white',
         width = width, height = height, limitsize = F)
}

tablesave <- function(tib, prefix, ...){
  path = strsplit(prefix,'/')[[1]]
  path = paste(path[-length(path)], collapse = '/')
  if(!file.exists(path)){
    system(paste('mkdir -p',path))
  }
  readr::write_csv(tib, path = paste(prefix,'.csv', sep = ''), append = F)
  readr::write_tsv(tib, path = paste(prefix,'.tsv', sep = ''), append = F)
  saveRDS(object = tib, file = paste(prefix,'.rds',sep=''))
}

#### Data ####
sampleinfo = readRDS('./data/processed/tidy/sample_info.rds') %>%
  mutate(period = c('Development','Ageing')[1+(age>=90)]) 
expx = readRDS('./data/processed/tidy/expression.rds') %>%
  left_join(sampleinfo) 
tisspec = readRDS('./data/processed/raw/ts.specQ3.genes.rds')
tisspec = reshape2::melt(tisspec) %>%
  set_names(c('gene_id','native_tissue'))
divgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds') #<0 DC
dico = names(which(divgenes<0))
exp_ch = readRDS('./data/processed/tidy/expression_change.rds') %>%
  filter(period == 'aging')

cortex = as.matrix(filter(expx, tissue == 'Cortex' & period == 'Ageing') %>%
                     select(gene_id, sample_id, expression) %>%
                     reshape2::acast(gene_id~sample_id))
lung = as.matrix(filter(expx, tissue == 'Lung' & period == 'Ageing') %>%
                   select(gene_id, sample_id, expression) %>%
                   reshape2::acast(gene_id~sample_id))
liver = as.matrix(filter(expx, tissue == 'Liver' & period == 'Ageing') %>%
                    select(gene_id, sample_id, expression) %>%
                    reshape2::acast(gene_id~sample_id))
muscle = as.matrix(filter(expx, tissue == 'Muscle' & period == 'Ageing') %>%
                     select(gene_id, sample_id, expression) %>%
                     reshape2::acast(gene_id~sample_id))

sampage = setNames(sampleinfo$log2age,sampleinfo$sample_id)
cortage = sampage[colnames(cortex)]
cortbeta = t(apply(cortex,1,function(x){
  summary(lm(x~cortage))$coef[2,c(1,4)]
}))
lungage = sampage[colnames(lung)]
lungbeta = t(apply(lung,1,function(x){
  summary(lm(x~lungage))$coef[2,c(1,4)]
}))
liverage = sampage[colnames(liver)]
liverbeta = t(apply(liver,1,function(x){
  summary(lm(x~liverage))$coef[2,c(1,4)]
}))
muscleage = sampage[colnames(muscle)]
musclebeta = t(apply(muscle,1,function(x){
  summary(lm(x~muscleage))$coef[2,c(1,4)]
}))

cortbeta = as.data.frame(cortbeta) %>% set_names(c('beta','beta.p')) %>% 
  mutate(gene_id = rownames(cortbeta)) %>% mutate(tissue ='Cortex')
lungbeta = as.data.frame(lungbeta) %>% set_names(c('beta','beta.p')) %>% 
  mutate(gene_id = rownames(lungbeta)) %>% mutate(tissue ='Lung')
liverbeta = as.data.frame(liverbeta) %>% set_names(c('beta','beta.p')) %>% 
  mutate(gene_id = rownames(liverbeta)) %>% mutate(tissue ='Liver')
musclebeta = as.data.frame(musclebeta) %>% set_names(c('beta','beta.p')) %>% 
  mutate(gene_id = rownames(musclebeta)) %>% mutate(tissue ='Muscle')

exp_beta = rbind(cortbeta,lungbeta,liverbeta,musclebeta)

saveRDS(exp_beta,'results/figure4/exp_beta.rds')
#### Loss of tissue specific expression ####

expdat = exp_beta %>%
  group_by(gene_id) %>%
  summarise(max_exp_ch = list(tissue[abs(beta) == max(abs(beta))]) ) %>%
  ungroup() %>%
  right_join(exp_beta) %>%
  left_join(sampleinfo) %>%
  left_join(expx) %>%
  left_join(tisspec) %>%
  mutate(`DevDiv` = gene_id %in% names(divgenes)) %>%
  mutate(DC = gene_id %in% dico) %>%
  left_join(exp_ch)

expdat$exc = sapply(expdat$max_exp_ch,length)>1
expdat = filter(expdat,!exc)
expdat$max_exp_ch = sapply(expdat$max_exp_ch, c)
expdat = expdat %>%
  mutate(sametis = c('Max Expr. Change\nin Other Tissues',
                     'Max Expr. Change\nin the Native Tissue')[1+(max_exp_ch == `native_tissue`)]) %>%
  mutate(expchange = c('Down','Up')[1+(beta>0)])

dc_fisher = expdat %>%
  filter(DC) %>%
  filter(tissue == max_exp_ch) %>%
  group_by(sametis, expchange) %>%
  summarise(n = length(unique(gene_id))) %>%
  filter(!is.na(sametis)) %>%
  ungroup()

fidc = fisher.test(reshape2::acast(dc_fisher, sametis~expchange))
fidcOR = round(1/fidc$estimate,2)
fidc$p.value
fidctitle = expression('DiCo Genes - OR=74.81 p=5.9x10'^'-203') 
dccontplot = dc_fisher %>%
  ggplot(aes( x = sametis, fill = expchange, y = n)) +
  geom_bar(stat = 'identity', position = 'fill', width = 1, color = 'black') +
  geom_text(aes(label = paste('n=',n,sep='')), position = 'fill', 
            size = 7/pntnorm, vjust = 1.2, color = 'black') +
  scale_fill_manual(values = regcol,
                    guide = guide_legend(
                      title = 'Direction of Expression Change During Ageing',
                      title.vjust = 1, title.hjust = 0.5)) +
  xlab(NULL) + ylab('') +
  annotate(geom='text',x=c(1.1, 1.1, 2.1, 2.1)-0.5,y=c(0.95,0.05,0.95,0.05),
           label=paste('GR',c(1,2,3,4), sep=''), size = 6/pntnorm, fontface = 'bold') +
  #theme_rvis(base_size = 6, legend.pos = 'bottom') +
  theme(axis.ticks.length.x = unit(0,'pt'),
        axis.text.x = element_text(vjust = 2),
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(vjust = -0.5),
        legend.background = element_rect(color='black',size = 0.1),
        legend.key.size = unit(2.5, 'pt'),
        legend.spacing.x = unit(2.5, 'pt')) +
  ggtitle(fidctitle  )  
  #ggtitle(paste('DiCo Genes - OR=',round(1/fidc$estimate,2),scales::pvalue(fidc$p.value, add_p = T))) 

saveRDS(dc_fisher, 'results/source_data/f4/b.rds')

all_fisher = expdat %>%
  # filter(DC) %>%
  filter(tissue == max_exp_ch) %>%
  group_by(sametis, expchange) %>%
  summarise(n = length(unique(gene_id))) %>%
  filter(!is.na(sametis)) 

fiall = fisher.test(reshape2::acast(all_fisher, sametis~expchange))

fiall$p.value
1/fiall$estimate
fiall.title = expression('All Genes - OR=5.50 p=2.1x10'^'-129') 
allcontplot = all_fisher %>%
  ggplot(aes( x = sametis, fill = expchange, y = n)) +
  geom_bar(stat = 'identity', position = 'fill', width = 1, color = 'black') +
  geom_text(aes(label = paste('n=',n,sep='')), position = 'fill', size = 7/pntnorm, 
            vjust = 1.2, color = 'black') +
  scale_fill_manual(values = regcol, 
                    guide = guide_legend(title = 'Direction of Expression Change During Ageing',
                    title.vjust = 1, title.hjust = 0.5)) +
  xlab(NULL) + ylab('Proportion') +
  theme(axis.ticks.length.x = unit(0,'pt'),
        axis.text.x = element_text(vjust = 2),
        panel.border = element_blank(), 
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(vjust = -0.5),
        legend.background = element_rect(color='black', size=0.1),
        legend.key.size = unit(2.5, 'pt'),
        legend.spacing.x = unit(2.5, 'pt')) +
  ggtitle(fiall.title)  
  #ggtitle(paste('All Genes - OR=',round(1/fiall$estimate,2),scales::pvalue(fiall$p.value, add_p = T)))

saveRDS(all_fisher,'results/source_data/f4/a.rds')

# p1 = ggarrange(allcontplot, dccontplot, common.legend = T, legend = 'top', labels = c('a.','b.'), 
#                font.label = list(size = 8))

ourdat = readRDS('data/processed/specloss_wdico_fisher.rds')
jonker = readRDS('data/other_datasets/jonker/processed/specloss_fisher_co.rds')
schaum4 = readRDS('data/other_datasets/schaum/4tissue/specloss_fisher_co.rds')
schaum8 = readRDS('data/other_datasets/schaum/8tis/specloss_fisher_co.rds')
gtex4 = readRDS('data/other_datasets/GTEx/specloss_fisher_co.rds')
gtex10 = readRDS('results/GTEx/alltissues/specloss_fisher_co.rds')

dat = rbind(
  Our_data = c(OR = unname(ourdat$fisher$est), p = ourdat$fisher$p.value),
  Jonker = c(OR = unname(jonker$fisher$est), p = jonker$fisher$p.value),
  Schaum4 = c(OR = unname(schaum4$fisher$est), p = schaum4$fisher$p.value),
  Schaum8 = c(OR = unname(schaum8$fisher$est), p = schaum8$fisher$p.value),
  GTEx4 = c(OR = unname(gtex4$fisher$est), p = gtex4$fisher$p.value),
  GTEx10 = c(OR = unname(gtex10$fisher$est), p = gtex10$fisher$p.value) ) %>%
  reshape2::melt() %>%
  spread(key='Var2', value='value') %>%
  set_names(c('dataset', 'OR', 'P value'))
# dataset        OR       P value
# 1 Our_data 74.808000 5.899456e-203
# 2   Jonker  7.518605 6.464742e-109
# 3  Schaum4 58.025610 1.505944e-197
# 4  Schaum8 84.200686  9.733053e-96
# 5    GTEx4  7.208754  7.040190e-87
# 6   GTEx10 13.014020 5.680559e-114


#annot = data.frame(dataset = dat$dataset, yval = log2(dat$OR), )
p.adjust(dat$`P value`, method = 'BH')
ORplot = ggplot(dat, aes(y= log2(OR), x=dataset)) +
  geom_bar(stat='identity', fill = '#9B6A6C') +
  xlab('') + ylab('Log2 (OR)') +
  geom_text(x = 1:6, y = log2(dat$OR)+0.02, label='***') +
  theme(axis.text.x  = element_text(size=4) )
ORplot
ggsave(filename = 'results/figure4/alldatasets.pdf',ORplot, units='cm', height = 8, width = 8, useDingbats=F)
ggsave(filename = 'results/figure4/alldatasets.png',ORplot, units='cm', height = 8, width = 8)

saveRDS(dat, 'results/source_data/f4/new_c.rds')
# p1 = ggarrange(allcontplot, dccontplot, ORplot,  ncol=3, common.legend = T, legend = 'top', 
#           labels = c('a.','b.','c.'),font.label = list(size = 8))

p1 = ggarrange(
  ggarrange(allcontplot, dccontplot, common.legend = T, legend = 'top', labels = c('a.','b.'),
            font.label = list(size = 8) ), widths = c(2,1),
               ORplot, labels=c(NA, 'c.'), font.label = list(size = 8)) 

p1.1 = ggarrange(
  ggarrange(allcontplot, dccontplot, common.legend = T, legend = 'top', labels = c('a.','b.'),
            font.label = list(size = 8) ), widths = c(2,1),
  ggarrange(NULL, ORplot, nrow=2, heights = c(0.1, 0.8)),
  labels=c(NA, 'c.'), font.label = list(size = 8), vjust = 4) 

##
GR1 = unique((expdat %>% filter((max_exp_ch != native_tissue) & (beta<0) & DC) %>% 
                filter(max_exp_ch == tissue) %>% 
                select(beta, gene_id) %>% unique() %>%
                arrange(-abs(beta)) )$gene_id)
GR2 = unique((expdat %>% filter((max_exp_ch != native_tissue) & (beta>0) & DC) %>% 
                filter(max_exp_ch == tissue) %>% 
                select(beta, gene_id) %>% unique() %>%
                arrange(-abs(beta)) )$gene_id)
GR3 = unique((expdat %>% filter((max_exp_ch == native_tissue) & (beta<0) & DC) %>% 
                filter(max_exp_ch == tissue) %>% 
                select(beta, gene_id) %>% unique() %>%
                arrange(-abs(beta)) )$gene_id)
GR4 = unique((expdat %>% filter((max_exp_ch == native_tissue) & (beta>0) & DC) %>% 
                filter(max_exp_ch == tissue) %>% 
                select(beta, gene_id) %>% unique() %>%
                arrange(-abs(beta)) )$gene_id)

genelist = c(GR1,GR2,GR3,GR4)
martx = biomaRt::useMart('ensembl','mmusculus_gene_ensembl')
genex = biomaRt::getBM(attributes = c('ensembl_gene_id', 'mgi_symbol', 'mgi_description'), 
                       filters = 'ensembl_gene_id', values = genelist, mart = martx)
genex = setNames(genex$mgi_symbol,genex$ensembl_gene_id)
saveRDS(genex, 'results/figure4/gene_names.rds')
i=6
gr1dat = expdat %>% filter(gene_id == GR1[i])
gr1plot = ggplot(gr1dat, aes(x = age, y = expression, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed') +
  geom_smooth(se = F, size = 0.5) +
  geom_point(size = 0.2) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  ggtitle(paste('GR1 Example,',genex[GR1[i]])) +
  #theme_rvis(base_size = 6, legend.pos='bottom') +
  theme(plot.title = element_text(face = 'italic', size = 7),
        legend.background = element_rect(color='black', size=0.1),
        panel.border = element_rect(color='black', fil=NA)) +
  xlab('Age (in log2 scale)') + ylab('Gene Expression')

saveRDS(gr1dat, 'results/source_data/f4/c.rds')

i=1
gr2dat = expdat %>% filter(gene_id == GR2[i])
gr2plot = expdat %>%
  filter(gene_id == GR2[i]) %>%
  ggplot(aes(x = age, y = expression, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed') +
  geom_smooth(se = F, size = 0.5) +
  geom_point(size = 0.2) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  ggtitle(paste('GR2 Example,',genex[GR2[i]])) +
  #theme_rvis(base_size = 6, legend.pos='bottom') +
  theme(plot.title = element_text(face = 'italic', size = 7),
        legend.background = element_rect(color='black', size=0.1),
        panel.border = element_rect(color='black', fil=NA)) +
  xlab('Age (in log2 scale)') + ylab('Gene Expression')

saveRDS(gr2dat, 'results/source_data/f4/d.rds')

i=2
gr3dat = expdat %>% filter(gene_id == GR3[i])
gr3plot = expdat %>%
  filter(gene_id == GR3[i]) %>%
  ggplot(aes(x = age, y = expression, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed') +
  geom_smooth(se = F, size = 0.5) +
  geom_point(size = 0.2) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  ggtitle(paste('GR3 Example,',genex[GR3[i]])) +
  #theme_rvis(base_size = 6, legend.pos='bottom') +
  theme(plot.title = element_text(face = 'italic', size = 7),
        legend.background = element_rect(color='black', size=0.1),
        panel.border = element_rect(color='black', fil=NA)) +
  xlab('Age (in log2 scale)') + ylab('Gene Expression')

saveRDS(gr3dat, 'results/source_data/f4/e.rds')

i=1
gr4dat = expdat %>% filter(gene_id == GR4[i])
gr4plot = expdat %>%
  filter(gene_id == GR4[i]) %>%
  ggplot(aes(x = age, y = expression, color = tissue)) +
  geom_vline(xintercept = 90, linetype = 'dashed') +
  geom_smooth(se = F, size = 0.5) +
  geom_point(size = 0.2) +
  scale_color_manual(values = tissuecol) +
  scale_x_continuous(trans = 'log2') +
  ggtitle(paste('GR4 Example,',genex[GR4[i]]))+
  #theme_rvis(base_size = 6, legend.pos='bottom') +
  theme(plot.title = element_text(face = 'italic', size = 7),
        legend.background = element_rect(color='black', size=0.1),
        panel.border = element_rect(color='black', fil=NA)) +
  xlab('Age (in log2 scale)') + ylab('Gene Expression') 

saveRDS(gr4dat, 'results/source_data/f4/f.rds')

p2 = ggarrange(gr1plot, gr2plot, gr3plot,gr4plot,font.label = list(size = 8), 
               labels = c('d.','e.','f.','g.'), ncol = 4, nrow = 1, common.legend = T, legend = 'top')
fig4af = ggarrange(p1.1, p2, ncol = 1 , nrow = 2, heights = c(1,1))

##### Fig 4g (Dico Enrichment plot)
expch = readRDS('./data/processed/tidy/expression_change.rds') %>%
  mutate(period = gsub('aging', 'Ageing', period)) %>%
  mutate(period = str_to_title(period) ) %>%
  mutate(period = factor(period, levels = c('Development', 'Ageing'))) %>%
  dplyr::rename(Period = period)

reprg = readRDS('./results/figure4/gorepresentatives.rds')

enricplotdat = expch %>%
  inner_join(reprg) %>%
  group_by(Period, ID, tissue, Description, repNames, n) %>%
  summarise(mrho = mean(`Expression Change`),
            medrho = median(`Expression Change`))
enricplot = expch %>%
  inner_join(reprg) %>%
  group_by(Period, ID, tissue, Description, repNames, n) %>%
  summarise(mrho = mean(`Expression Change`),
            medrho = median(`Expression Change`)) %>%
  ggplot(aes(fill=Period, y=mrho, x=reorder(repNames, n) ) ) +
  geom_bar(stat='identity', position=position_dodge()) +
  facet_wrap(~tissue, ncol=4) +
  scale_fill_manual(values=brewer.pal(3,"Set1")[c(2,1)]) +
  coord_flip() +
  geom_hline(yintercept = 0, size=0.3, linetype='solid', color='gray30') +
  geom_vline(xintercept = seq(1.5,25, by=1), linetype = 'dashed', size = 0.2, color = 'gray') +
  theme(legend.position = 'top',
        axis.text.x = element_text(size=4, vjust=2), 
        axis.ticks.length.x = unit(0,'pt'),
        axis.text.y = element_text(size=5),
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(vjust = -0.5),
        legend.background = element_rect(color='black', size=0.1),
        legend.key.size = unit(3, 'pt'),
        axis.title.x = element_text(size=6)) +
  xlab('') +
  ylab(bquote('Mean Expression Change ('*rho*')'))
enricplot

saveRDS(enricplotdat,'results/source_data/f4/g.rds')

ggsave('./results/figure4/dicoGO.pdf',enricplot, units='cm', width = 16, height = 8, useDingbats=F)
ggsave('./results/figure4/dicoGO.png',enricplot, units='cm', width = 16, height = 8)  

fig4  = ggarrange(fig4af, enricplot, nrow =2, labels=c(NA, 'h.'), font.label = list(size = 8))
fig4
plotsave(fig4,'./results/figure4/figure4', fig4, units='cm', width = 16, height = 16)


#### enriched GO with expr. changes, FDR<0.1
#########
######### Figure 4- figure supplement 1
#########

# enricdat = expch %>%
#   inner_join(reprgenes) %>%
#   group_by(Period, GOID, tissue, Representative) %>%
#   filter(FDR<0.1, .preserve = T) %>% 
#   summarise(n = n(),
#             mrho = mean(`Expression Change`),
#             medrho = median(`Expression Change`))
# enrichpfdr =  enricdat %>% 
#   ggplot(aes(fill=Period, y=mrho, x=Representative ) ) +
#   geom_bar(stat='identity', position=position_dodge2(preserve = 'single', padding=0.2) ) +
#   facet_wrap(~tissue, ncol=4) +
#   scale_fill_manual(values=brewer.pal(3,"Set1")[c(2,1)]) +
#   coord_flip() +
#   geom_hline(yintercept = 0, size=0.3, linetype='solid', color='gray30') +
#   geom_vline(xintercept = seq(1.5,17, by=1), linetype = 'dashed', size = 0.2, color = 'gray') +
#   geom_text(aes(label=n), color='gray20', position=position_dodge(width = .9),
#            angle=0, hjust=1, vjust=0.5, size=2 ) +
#   #geom_text(data = filter(enricdat, medrho<0, .preserve = T ), aes(label=n),color='gray20', 
#   #          position=position_dodge(width = .9), angle=0, hjust=1, vjust=0.5, size=1.5 )  +
#   #geom_text(data = filter(enricdat, medrho>0 ), aes(label=n),color='gray20', 
#   #          position=position_dodge(width = .9), angle=0, hjust=1, vjust=0.5, size=1.5 )  +
#   theme(legend.position = 'top',
#         axis.text.x = element_text(size=4), 
#         axis.text.y = element_text(size=5),
#         axis.title.x = element_text(size=6)) +
#   xlab('') +
#   ylab(bquote('Median Significant Expression Change ('*rho*')'))
# enrichpfdr
# ggsave('./results/figure4/dicoGOfdr.pdf',enrichpfdr, units='cm', width = 16, height = 12, useDingbats=F)
# ggsave('./results/figure4/dicoGOfdr.png',enrichpfdr, units='cm', width = 16, height = 12)  
