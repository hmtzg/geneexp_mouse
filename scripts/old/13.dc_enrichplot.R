library(clusterProfiler)
#library(enrichplot)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
#library(goseq)

theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)

dc = readRDS('./data/processed/raw/dc_gse.rds')
signifGO = dc@result[,1:10] %>%
  filter(qvalues < 0.1 & NES < 0) %>%
  select(ID, NES)
write.table(signifGO, file = 'results/figure4/signifGO_DiCo.csv', row.names = F, quote = F)  

repr = readRDS('./results/figure4/revigo_representative.rds')

library(GO.db)
library(AnnotationDbi)
library(org.Mm.eg.db)
term = select(GO.db, columns = c("GOID","TERM"), keytype = "GOID", keys = keys(GO.db))
terms = as.character(term[,2])
names(terms) = as.character(term[,1])
#save(terms, file = "./data/goterms_20210831.RData")
saveRDS(terms, './data/goterms_20210831.rds')
x <- c(as.list(GOBPOFFSPRING)) # children and all their children
xx <- sapply(1:length(x), function(i) x[[i]] = c(names(x)[i], x[[i]]) )
xxx <- sapply(1:length(xx), function(i) xx[[i]][ complete.cases(xx[[i]]) ] )
length(x)
length(xx)
length(xxx)
names(xx) <- names(x)
names(xxx) <- names(x)
xxx[[1]]
go2gene = select(org.Mm.eg.db, keys = names(xxx), columns = c("GO","ENSEMBL"), keytype = "GO")
GO2Gene = list()
for(i in 1:length(xxx)){
  print(i/length(xxx))
  myg = go2gene[go2gene[,1]%in%xxx[[i]], 4]
  myg = unique(myg[complete.cases(myg)])
  GO2Gene[[i]] = myg
}
names(GO2Gene) = names(xxx)
head(GO2Gene)
length(GO2Gene)
save(GO2Gene, file="./data/GO2GeneBP_20210831.RData")
reprgenes = GO2Gene[repr$term_ID]

reprgenes = reshape2::melt(reprgenes) %>% 
  set_names(c('gene_id', 'GOID') ) %>%
  left_join( set_names(repr, c('GOID', 'Representative') )  )

reprlevels = reprgenes %>% group_by(GOID, Representative) %>%
  summarise(n=n()) %>%
  arrange(-n)
reprgenes = reprgenes %>%
  mutate(GOID = factor(GOID, levels = reprlevels$GOID),
         Representative = factor(Representative, levels = reprlevels$Representative ))

expch = readRDS('./data/processed/tidy/expression_change.rds') %>%
  mutate(period = gsub('aging', 'Ageing', period)) %>%
  mutate(period = str_to_title(period) ) %>%
  mutate(period = factor(period, levels = c('Development', 'Ageing'))) %>%
  dplyr::rename(Period = period)

enricplot = expch %>%
  inner_join(reprgenes) %>%
  group_by(Period, GOID, tissue, Representative) %>%
  summarise(mrho = mean(`Expression Change`),
            medrho = median(`Expression Change`)) %>% 
  ggplot(aes(fill=Period, y=mrho, x=Representative)) +
  geom_bar(stat='identity', position=position_dodge()) +
  facet_wrap(~tissue, ncol=4) +
  scale_fill_manual(values=brewer.pal(3,"Set1")[c(2,1)]) +
  coord_flip() +
  geom_hline(yintercept = 0, size=0.3, linetype='solid', color='gray30')+
  theme(legend.position = 'top',
        axis.text.x = element_text(size=4), 
        axis.text.y = element_text(size=5),
        axis.title.x = element_text(size=6)) +
  xlab('') +
  ylab(bquote('Mean Expression Change ('*rho*')'))
enricplot
ggsave('./results/figure4/dicoGO.pdf',enricplot, units='cm', width = 16, height = 8, useDingbats=F)
ggsave('./results/figure4/dicoGO.png',enricplot, units='cm', width = 16, height = 8)  

enricplot2 = expch %>%
  inner_join(reprgenes) %>%
  group_by(Period, GOID, tissue, Representative) %>%
  summarise(mrho = mean(`Expression Change`),
            medrho = median(`Expression Change`)) %>% 
  ggplot(aes(fill=Period, y=medrho, x=Representative)) +
  geom_bar(stat='identity', position=position_dodge()) +
  facet_wrap(~tissue, ncol=4) +
  scale_fill_manual(values=brewer.pal(3,"Set1")[c(2,1)]) +
  coord_flip() +
  geom_hline(yintercept = 0, size=0.3, linetype='solid', color='gray30')+
  theme(legend.position = 'top',
        axis.text.x = element_text(size=4), 
        axis.text.y = element_text(size=5),
        axis.title.x = element_text(size=6)) +
  xlab('') +
  ylab(bquote('Median Expression Change ('*rho*')'))
enricplot2
ggsave('./results/figure4/dicoGOmed.pdf',enricplot2, units='cm', width = 16, height = 8, useDingbats=F)
ggsave('./results/figure4/dicoGOmed.png',enricplot2, units='cm', width = 16, height = 8)  

enrichpfdr = expch %>%
  inner_join(reprgenes) %>%
  group_by(Period, GOID, tissue, Representative) %>%
  filter(FDR<0.1) %>% 
  summarise(n = n(),
            mrho = mean(`Expression Change`),
            medrho = median(`Expression Change`)) %>% 
  ggplot(aes(fill=Period, y=mrho, x=Representative, width = (n/max(n))^0.5 ) ) +
  geom_bar(stat='identity', position=position_dodge() ) +
  facet_wrap(~tissue, ncol=4) +
  scale_fill_manual(values=brewer.pal(3,"Set1")[c(2,1)]) +
  coord_flip() +
  geom_hline(yintercept = 0, size=0.3, linetype='solid', color='gray30')+
  theme(legend.position = 'top',
        axis.text.x = element_text(size=4), 
        axis.text.y = element_text(size=5),
        axis.title.x = element_text(size=6)) +
  xlab('') +
  ylab(bquote('Mean Expression Change ('*rho*')'))
enrichpfdr
ggsave('./results/figure4/dicoGOfdr.pdf',enrichpfdr, units='cm', width = 16, height = 8, useDingbats=F)
ggsave('./results/figure4/dicoGOfdr.png',enrichpfdr, units='cm', width = 16, height = 8)  


### for tissue specific genes:
ts.spec = readRDS('./data/processed/raw/ts.specQ3.genes.rds') %>%
  reshape2::melt() %>%
  set_names(c('gene_id', 'Spec'))
reprgenes %>%
  left_join(ts.spec) %>%
  group_by(GOID, Spec, Representative) %>%
  summarise(n = length(gene_id))

reprgenes %>%
  inner_join(ts.spec) %>% 
  inner_join(expch) %>% 
  group_by(period, GOID, tissue, Representative) %>%
  summarise(mrho = mean(`Expression Change`),
            medrho = median(`Expression Change`)) %>% 
  ggplot(aes(fill=period, y=mrho, x=Representative)) +
  geom_bar(stat='identity', position=position_dodge()) +
  facet_wrap(~tissue, ncol=4) +
  scale_fill_manual(values=brewer.pal(3,"Set1")[c(2,1)]) +
  #scale_y_continuous(trans='log2') +
  coord_flip() +
  geom_hline(yintercept = 0, size=0.3, linetype='solid', color='gray30',)+
  theme(legend.position = 'top',
        axis.text.x = element_text(size=4), 
        axis.text.y = element_text(size=7)) +
  xlab('') +
  ylab('Mean Rho')

# dc_rep = dc
# dc_rep@result = dc@result[dc@result$ID%in%repr$term_ID,]
# emapplot(dc_rep, color='NES')
#dc_rep@result = dc@result[dc@result$Description%in%rep, 1:10]
#dc_emap = emapplot(dc, color = 'NES', showCategory = 25, layout = 'nicely')
#dc_emap
#dc_emap$layers =  dc_emap$layers[-3]
#dc_emap = dc_emap + 
#  geom_text_repel(aes(x = x, y=y, label = name), size = 1.5)
#  
# dc_emap$labels$colour = 'NES'
# 
# dc_gse_plot = dc_emap +
#   scale_color_gradient(low = 'indianred4', high = 'palevioletred2') +
#   theme(legend.text = element_text(size = 4),
#         text = element_text(size = 4))
# 
# ggsave('Dropbox/projects/ageing/results.n/sd.method/dc_gse_plot.pdf', dc_gse_plot, units = 'cm',
#        height = 8, width = 12, useDingbats=F)
# ggsave('Dropbox/projects/ageing/results.n/sd.method/dc_gse_plot.png', dc_gse_plot, units = 'cm',
#        height = 8, width = 12)

######
# a=readRDS('Dropbox/projects/repos/geneexp_mouse/data/processed/raw/dev_divergent_genes_dc_rank.rds')
# head(a)
# b= strsplit(dc@result[1,11], '/')[[1]]
# 
# 
# head(dc@result[,1:10])
# head(dc@result[dc@result$NES>0,1:10],20)
# 
# which(dc@result$ID%in%'GO:0099504')
# b= strsplit(dc@result[212,11], '/')[[1]]