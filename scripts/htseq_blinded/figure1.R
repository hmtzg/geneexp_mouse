library(tidyverse)
library(reshape2)
library(ggpubr)
#library(pheatmap)
library(corrplot)
library(grid)
library(stringr)
library(RColorBrewer)

theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)

tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('Up','Down'))
revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))

pca_dat = readRDS('data/htseq/blinded/pca_data.rds')
sample_info = readRDS('data/processed/tidy/sample_info.rds')
expr_ch = readRDS('data/htseq/blinded/expression_change.rds')
revgenes = readRDS('data/htseq/blinded/revgenes.tissues_tidy.rds')

all_raw_var = (pca_dat %>%
                 filter(period == 'all', type == 'raw') %>%
                 select(varExp, PC) %>%
                 unique())$varExp
names(all_raw_var) = paste0('PC',1:4)
all_raw_var = round(all_raw_var,2)

### panel a ####
pa1 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info, by=c("sample_id")) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
  geom_point(alpha = 0.7, show.legend = T) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2') +
  coord_fixed(ratio = all_raw_var[2]/all_raw_var[1], clip = 'off') +
  xlab(paste('PC1 (', round(all_raw_var[1]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(all_raw_var[2]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'),
         size = guide_legend('Age')) + theme_bw(base_size = 6)

pa2 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info, by=c("sample_id")) %>%
  ggplot(aes(x = PC3, y= PC4, color = tissue, size = age))  +
  geom_point(alpha = 0.7, show.legend = T) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2') +
  coord_fixed(ratio = all_raw_var[4]/all_raw_var[3], clip = 'off') +
  xlab(paste('PC3 (', round(all_raw_var[3]*100),'%)',sep='')) +
  ylab(paste('PC4 (', round(all_raw_var[4]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'),
         size = guide_legend('Age'))  + theme_bw(base_size = 6)

### panel b ####
pb1 = pca_dat %>%
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
  xlab('Age in days (in log2 scale)') + theme_bw(base_size = 6)

pb2 = pca_dat %>%
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
  xlab('Age in days (in log2 scale)') + theme_bw(base_size = 6)

pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info) %>%
  mutate(period = ifelse(age<90,'dev', 'aging'))  %>%
  group_by(tissue, period) %>%
  summarise(rho = cor.test(PC4, age, m='s')$est,
            p = cor.test(PC4, age, m='s')$p.val) %>% arrange(period)
# tissue period     rho         p
# <fct>  <chr>    <dbl>     <dbl>
# 1 Cortex aging  -0.0120 0.978    
# 2 Liver  aging   0.798  0.00989  
# 3 Lung   aging   0.908  0.000721 
# 4 Muscle aging   0.395  0.293    
# 5 Cortex dev    -0.955  0.000806 
# 6 Liver  dev    -0.991  0.0000146
# 7 Lung   dev    -0.847  0.0162   
# 8 Muscle dev    -0.577  0.175  
  
## to get legend only:
pb3 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info, by=c("sample_id")) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
  geom_point(alpha = 0.7, show.legend = T) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2') +
  coord_fixed(ratio = all_raw_var[2]/all_raw_var[1], clip = 'off') +
  xlab(paste('PC1 (', round(all_raw_var[1]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(all_raw_var[2]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'),
         size = guide_legend('Age')) +
  theme_bw(base_size = 6) +
  theme(legend.box = 'horizontal',
        legend.key.size = unit(0.05, 'cm'),
        legend.text = element_text(margin = margin(t=1) ),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        legend.justification = c(0.2, 2) )

pb3 = get_legend(pb3)

#pb = ggarrange(pb1, pb2, pb3, ncol =3, widths = c(0.7, 0.7, 0.4))

### panel c - correlation heatmap #####
expdf = expr_ch %>%
  select(gene_id, tissue, period, `Expression Change`) %>%
  dcast(gene_id + period ~ tissue, value.var = "Expression Change")

expa = expdf %>% 
  filter(period == "Ageing") %>%
  select(- c("gene_id", "period")) %>%
  magrittr::set_colnames(paste0("age_", colnames(.)))

cors = expdf %>% 
  filter(period == "Development") %>%
  select(- c("gene_id", "period")) %>%
  magrittr::set_colnames(paste0("dev_", colnames(.))) %>%
  bind_cols(expa) %>%
  as.matrix() %>%
  cor(method="s", use="pairwise.complete.obs")

saveRDS(cors, './data/htseq/blinded/pwise_expch_cors.rds')

#####
periodcode = c(Development = "#FE6100", Ageing ="#648FFF")
# colnames(cors) = gsub('dev_|age_','',colnames(cors))
# rownames(cors) = gsub('dev_|age_','',colnames(cors))
colnames(cors) = rep(c('CTX', 'LV', 'LNG', 'MS'),2)
rownames(cors) = rep(c('CTX', 'LV', 'LNG', 'MS'),2)
diag(cors) = 0

pdf('results/htseq/blinded/figure1/corplot.pdf', height = 8, width = 9)
corrplot(cors, order = "ori", tl.pos = "lt", diag=T, tl.col = rep(periodcode, each=4), font =2,
         tl.cex = 1.2, method="square", outline=T, type="lower",
         col = colorRampPalette(rev(brewer.pal(n=11, name='RdBu')[c(1:3,6,9:11)]))(100) )
corrplot(cors, order="ori",tl.pos = "n", diag=F, tl.col = rep(periodcode, each=4), font=2,
         tl.cex = 1.2, method="shade", outline=T, addCoef.col = "black", add=T, type="upper",
         col = colorRampPalette(rev(brewer.pal(n=11, name='RdBu')[c(1:3,6,9:11)]))(100) )
dev.off()

png('results/htseq/blinded/figure1/corplot.png', width = 500)
corrplot(cors, order = "ori", tl.pos = "lt", diag=T, tl.col = rep(periodcode, each=4), font =2,
         tl.cex = 1.2, method="square", outline=T, type="lower",
         col = colorRampPalette(rev(brewer.pal(n=11, name='RdBu')[c(1:3,6,9:11)]))(100) )
corrplot(cors, order="ori",tl.pos = "n", diag=F, tl.col = rep(periodcode, each=4), font=2,
         tl.cex = 1.2, method="shade", outline=T, addCoef.col = "black", add=T, type="upper",
         col = colorRampPalette(rev(brewer.pal(n=11, name='RdBu')[c(1:3,6,9:11)]))(100) )
dev.off()

#####
# #annotations
# annot = data.frame(id = rownames(cors))%>%
#   mutate(id =  as.character(id),
#          Period = sapply(id, function(x) strsplit(x, "_")[[1]][1]),
#          Period = setNames(c("Development", "Ageing"),c("dev","age"))[Period]) %>%
#   column_to_rownames(var="id")
# my_col = list(Period = c(Development = "#FE6100", Ageing ="#648FFF"))
# 
# corplot = pheatmap(cors, breaks = seq(-1, 1, length.out=101),
#                    color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[c(1:3,6,9:11)]))(101), 
#                    annotation_row = annot, annotation_col = annot, annotation_colors = my_col, 
#                    legend_labels = c(1,0,-1), legend_breaks = c(1,0,-1),  
#                    annotation_names_row = F, annotation_names_col = F, 
#                    width = 4, height = 4, 
#                    cellwidth = 25, cellheight = 25,
#                    cluster_rows = F, cluster_cols = F, border_color = "gray60",
#                    labels_row = setNames(c("CTX", "LV", "LNG", "MS"),
#                                          c("Cortex", "Liver", "Lung", "Muscle"))[gsub("age_|dev_", "",
#                                                                                       rownames(cors))], 
#                    labels_col = setNames(c("CTX", "LV", "LNG", "MS"),
#                                          c("Cortex", "Liver", "Lung", "Muscle"))[gsub("age_|dev_", "",
#                                                                                       rownames(cors))], 
#                    show_rownames = T, show_colnames = T, display_numbers=T, number_format="%.2f",
#                    number_color="black", 
#                    fontsize_number =8, fontsize_col = 6, fontsize_row = 6, fontsize = 8)
# pc=corplot
# corplot
# svg("results/htseq/blinded/figure1/fig1c_corplot.svg", width = 5, height = 7/2.54)
# grid::grid.newpage()
# grid::grid.draw(corplot$gtable)
# dev.off()

### panel d  ####
pd = expr_ch %>%
  filter(FDR<=0.1) %>%
  group_by(tissue,period) %>%
  summarise( up = sum(`Expression Change`>0,na.rm=T),
             down = sum(`Expression Change`<0,na.rm=T)) %>%
  gather(key= 'direction',value ='n',-tissue,-period) %>%
  mutate(direction = str_to_title(direction), 
         period = gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels = c('Development','Ageing'))) %>%
  mutate(n = ifelse(n==0,NA,n)) %>%
  ggplot(aes(x = period, y= n, fill = direction)) +
  facet_wrap(~tissue, ncol=2) +
  geom_bar(stat= 'identity', position = 'dodge') +
  scale_fill_manual(values = regcol, drop=FALSE) +
  scale_y_continuous(trans = 'log10') +
  geom_text(aes(label = n), color='gray15', position = position_dodge(width = 0.9), 
            angle=90, hjust = 1.1, vjust = 0.5, size=6/pntnorm) +
  guides(fill = guide_legend('Direction of\nExpression Change', 
                             override.aes = list(size = 2))) +
  theme(legend.position = 'right',
        legend.background = element_rect(fill = 'gray85',color = 'gray25'),
        axis.text.x = element_text(size=6/pntnorm)) +
  xlab("Period")  + ylab("No. of genes (in log scale)") +
  theme_bw(base_size = 6)
pd
ggsave('./results/htseq/blinded/figure1/fig1d.pdf', pd, units='cm', height = 9, width = 14, useDingbats=F)
ggsave('./results/htseq/blinded/figure1/fig1d.png', pd, units='cm', height = 9, width = 14)

### panel e ####
pe = expr_ch %>%
  drop_na() %>%
  filter(`Expression Change` !=0 ) %>% # drop expression - age rho = 0 values
  mutate(direction = `Expression Change` > 0) %>%
  mutate(direction = ifelse( direction == TRUE, 'Up', 'Down')) %>%
  mutate(period = gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels = c('Development','Ageing'))) %>%
  group_by(gene_id,period, direction) %>%
  summarise(n = length(tissue)) %>%
  group_by(n, period, direction) %>% 
  summarise(count = n()) %>%
  ungroup() %>%
  slice(-c(1:4)) %>% 
  ggplot(aes(x = n, y = count, fill = direction)) +
  facet_wrap(~period) +
  geom_bar(stat='identity', position = 'dodge') +
  scale_fill_manual(values = regcol, drop =F) +
  geom_text(aes(label = count), color = 'gray15', position = position_dodge(width = 1), angle = 90,
            hjust = 1.1, size=6/pntnorm) +
  theme(legend.position = 'none') +
  theme_bw(base_size = 6) +
  guides(fill = guide_legend('Direction of\nExpression Change')) +
  xlab("Overlap") + ylab("No. of genes") 
pe
ggsave('./results/htseq/blinded/figure1/fig1e.pdf', pe, units='cm', height = 8, width = 14, useDingbats=F)
ggsave('./results/htseq/blinded/figure1/fig1e.png', pe, units='cm', height = 8, width = 14)

###  panel f ####
# rev props:
revgenes %>%
  group_by(direction,tissue) %>%
  summarise(n = length(gene_id)) %>%
  group_by(tissue) %>%
  mutate( tot = sum(n)) %>%
  group_by(tissue) %>%
  summarise(prop = (n[direction=='DownUp'] + n[direction=='UpDown']) / sum(n)  )
# tissue  prop
# <chr>  <dbl>
# 1 Cortex 0.419
# 2 Liver  0.588
# 3 Lung   0.517
# 4 Muscle 0.563

pf = revgenes %>%
  group_by(direction,tissue) %>%
  summarise(n = length(gene_id)) %>%
  group_by(tissue) %>%
  mutate(direction = factor(direction,levels = c('UpDown','DownUp','UpUp','DownDown')),
         tissue = str_to_title(tissue)) %>%
  ggplot(aes(x=tissue, y=n, fill = direction)) +
  geom_bar(stat='identity', position = position_fill(reverse=T), alpha=0.9) + 
  scale_fill_manual(values = revcol) + 
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray25') +
  xlab("Tissue") + ylab('Proportion of genes') +
  guides(fill = guide_legend('Direction',
                             override.aes = list(size=2), ncol =2)) +
  theme_bw(base_size = 6) +
  theme(legend.position = "top",
        legend.key.size = unit(0.1, 'cm'),
        legend.spacing.x = unit(0.1, 'cm'))

ggsave('./results/htseq/blinded/figure1/fig1f.pdf', pf, units='cm', height = 10, width = 8, useDingbats=F)
ggsave('./results/htseq/blinded/figure1/fig1f.png', pf, units='cm', height = 10, width = 8)

pdeleg = get_legend(pd)

# pab = ggarrange(ggarrange(pa1, pa2, ncol=2, common.legend=TRUE, widths = c(1, 1.5), legend="none"), 
#                 ggarrange(pb1, pb2, pb3, pdeleg, ncol=4, widths=c(1.1, 1.1, 0.8, 0.7) ),
#                 nrow=2, heights=c(2,1),
#                 common.legend=F)
#pab
# ggsave("results/htseq/blinded/figure1/fig1ab.pdf", pab, units = 'cm',width = 12, height = 8, 
#        useDingbats = F)
# ggsave("results/htseq/blinded/figure1/fig1ab.png", pab, units = 'cm', width = 12, height = 8, bg = 'white')
# fig1 = ggarrange( ggarrange(pab, as_ggplot(pc$gtable), ncol=2, widths = c(1.6,1) ),
#           ggarrange(ggarrange(pd, pe, ncol=2, widths = c(1,1), common.legend = T, legend=F),
#                     pf, ncol = 2, widths = c(2.5,1) ), 
#           nrow=2, heights = c(1.5,1) )


## Figure 1: 
## to get legend only:
pb3 = pca_dat %>%
  select(-varExp) %>%
  filter(period == 'all', type == 'raw') %>%
  spread(key = 'PC',value = 'value') %>%
  left_join(sample_info, by=c("sample_id")) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue, size = age))  +
  geom_point(alpha = 0.7, show.legend = T) +
  scale_color_manual(values = tissuecol) +
  scale_size_continuous(range = c(0.5,3), trans= 'log2') +
  coord_fixed(ratio = all_raw_var[2]/all_raw_var[1], clip = 'off') +
  xlab(paste('PC1 (', round(all_raw_var[1]*100),'%)',sep='')) +
  ylab(paste('PC2 (', round(all_raw_var[2]*100),'%)',sep='')) +
  guides(color = guide_legend('Tissue'),
         size = guide_legend('Age')) +
  theme_bw(base_size = 6) +
  theme(legend.box = 'horizontal',
        legend.key.size = unit(0.3, 'cm'), 
        legend.text = element_text(margin = margin(t=1) ),
        legend.spacing.x = unit(0.01, 'cm'),
        legend.spacing.y = unit(0.01, 'cm'),
        legend.justification = c(0.2, 1) )
pb3 = get_legend(pb3)

pab = ggarrange(ggarrange(pa1, pa2, ncol=2, widths = c(1, 2), legend="none"), 
                ggarrange(pb1, pb2, pb3, ncol=3, widths=c(1, 1, 0.8) ),
                nrow=2, heights=c(1.5, 1), common.legend=F, labels=c('a.', 'b.'), font.label = list(size=8),
                vjust = c(3,1))

pdeleg= get_legend(pd)
fig1 = ggarrange(
  ggarrange(pab, ncol=2, widths = c(1, 0.6)),
  NULL,
  ggarrange(ggarrange(pd, as_ggplot(pdeleg), pe,legend = F, ncol=3, labels=c('d.', NA,'e.'),
                      font.label = list(size=8), widths = c(1.2, 0.5, 1)),
            pf, ncol=2, labels=c(NA,'f.'), font.label = list(size=8), widths = c(1, 0.4)),
  nrow=3, heights = c(1, 0.1, 0.8) )

ggsave('./results/htseq/blinded/figure1/fig1n.pdf', fig1, units='cm', width = 16, height = 12,
       useDingbats=F)
ggsave('./results/htseq/blinded/figure1/fig1n.png', fig1, units='cm', width = 16, height = 16, bg='white')

### add panel c in inkscape
