library(tidyverse)
library(reshape2)
library(ggpubr)
library(pheatmap)
library(grid)
library(stringr)
library(RColorBrewer)

#save pheatmap
save_phmap <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  svg(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('Up','Down'))
revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))

pca_dat = readRDS('data/processed/tidy/pca_data.rds')
sample_info = readRDS('data/processed/tidy/sample_info.rds')
expr_ch = readRDS('data/processed/tidy/expression_change.rds')
revgenes = readRDS('data/processed/tidy/revgenes.tissues.rds')

all_raw_var = (pca_dat %>%
                filter(period == 'all', type == 'raw') %>%
                select(varExp, PC) %>%
                unique())$varExp

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
         size = guide_legend('Age')) +
  theme_bw()

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
         size = guide_legend('Age')) +
  theme_bw()

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
  xlab('Age in days (in log2 scale)') +
  theme_bw()

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
  xlab('Age in days (in log2 scale)') +
  theme_bw()

##
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
  theme_bw() +
  theme(legend.box = 'horizontal')

pb3 = get_legend(pb3)
# pb = ggarrange(pb1, pb2, pb3, ncol =3, widths = c(0.7, 0.7, 0.4))
# pb
pab = ggarrange(ggarrange(pa1,pa2,ncol=2,common.legend=TRUE, widths = c(1, 1.8), legend="none"),
                ggarrange(pb1,pb2,pb3 ,ncol=3, widths=c(1.3,1.3,0.7)), nrow=2, heights=c(1.5,1),
                common.legend=F, labels = c('a.','b.'), font.label = list(size=8))

saveRDS(pca_dat %>% left_join(sample_info),
        file='./results/source_data/f1/pca.rds')

ggsave("results/figure1/fig1ab.pdf", pab, units = 'cm',width = 19, height = 10, useDingbats = F)
ggsave("results/figure1/fig1ab.png", pab, units = 'cm',width = 20, height = 10)

### panel c - correlation heatmap #####
expdf = expr_ch %>%
  select(gene_id, tissue, period, `Expression Change`) %>%
  dcast(gene_id + period ~ tissue, value.var = "Expression Change")

expa = expdf %>% 
  filter(period == "aging") %>%
  select(- c("gene_id", "period")) %>%
  magrittr::set_colnames(paste0("age_", colnames(.)))

cors = expdf %>% 
  filter(period == "development") %>%
  select(- c("gene_id", "period")) %>%
  magrittr::set_colnames(paste0("dev_", colnames(.))) %>%
  bind_cols(expa) %>%
  as.matrix() %>%
  cor(method="s", use="complete.obs")

saveRDS(cors, './data/processed/raw/pwise_expch_cors.rds')
saveRDS(cors,'./results/source_data/f1/cors.rds')

diag(cors) = 0

#annotations
annot = data.frame(id = rownames(cors))%>%
  mutate(id =  as.character(id),
         Period = sapply(id, function(x) strsplit(x, "_")[[1]][1]),
         Period = setNames(c("Development", "Ageing"),c("dev","age"))[Period]) %>%
  column_to_rownames(var="id")

my_col = list(Period = c(Development = "#FE6100",
                         Ageing ="#648FFF"))


corplot = pheatmap(cors, breaks = seq(-1, 1, length.out=101), 
                   color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")[c(1:3,6,9:11)]))(101), 
                 annotation_row = annot, annotation_col = annot, annotation_colors = my_col, 
                 legend_labels = c(1,0,-1), legend_breaks = c(1,0,-1),  fontsize = 18,
                 annotation_names_row = F, annotation_names_col = F,
                 width = 5 , height = 4, 
                 cellwidth = 25, cellheight = 25,
                 cluster_rows = F, cluster_cols = F, border_color = "gray60",
                 labels_row = setNames(c("CTX", "LV", "LNG", "MS"), 
                                       c("Cortex", "Liver", "Lung", "Muscle"))[gsub("age_|dev_", "", 
                                                                                    rownames(cors))], 
                 labels_col = setNames(c("CTX", "LV", "LNG", "MS"), 
                                       c("Cortex", "Liver", "Lung", "Muscle"))[gsub("age_|dev_", "", 
                                                                                    rownames(cors))], 
                 show_rownames = T, show_colnames = T,display_numbers=T, number_format="%.2f", 
                 number_color="black", fontsize_number=6)
pc=corplot
save_phmap(corplot,"results/figure1/fig1c_corplot.svg",
           width=15, height=15)

### panel d  ####
pd_dat = expr_ch %>%
  filter(FDR<=0.1) %>%
  group_by(tissue,period) %>%
  summarise( up = sum(`Expression Change`>0,na.rm=T),
             down = sum(`Expression Change`<0,na.rm=T)) %>%
  gather(key= 'direction',value ='n',-tissue,-period) %>%
  mutate(direction = str_to_title(direction), 
         period = gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels = c('Development','Ageing'))) %>%
  mutate(n = ifelse(n==0,NA,n))

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
            angle=90, hjust = 1.1, vjust = 0.5) +
  guides(fill = guide_legend('Direction of\nExpression Change', 
                             override.aes = list(size = 2))) +
  theme(legend.position = 'bottom',
        legend.background = element_rect(fill = 'gray85',color = 'gray25')) +
  xlab("Period")  + ylab("No. of genes (in log scale)") +
  theme_bw()

ggsave('./results/figure1/fig1d.pdf', pd, units='cm', height = 8, width = 16, useDingbats=F)
ggsave('./results/figure1/fig1d.png', pd, units='cm', height = 8, width = 16)

saveRDS(pd_dat, file='./results/source_data/f1/exprch_allg.rds')

### panel e ####
pe_dat = expr_ch %>%
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
  slice(-c(1:4))

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
            hjust = 1.1) +
  guides(fill = guide_legend('Direction of\nExpression Change')) +
  xlab("Overlap") + ylab("No. of genes") +
  theme_bw()

ggsave('./results/figure1/fig1e.pdf', pe, units='cm', height = 8, width = 16, useDingbats=F)
ggsave('./results/figure1/fig1e.png', pe, units='cm', height = 8, width = 16)

saveRDS(pe_dat, file='./results/source_data/f1/exprch_sigg.rds')

###  panel f ####
pf_dat = revgenes %>%
  group_by(direction,tissue) %>%
  summarise(n = length(gene_id)) %>%
  group_by(tissue) %>%
  mutate(direction = factor(direction,levels = c('UpDown','DownUp','UpUp','DownDown')),
         tissue = str_to_title(tissue))
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
  theme_bw() +
  theme(legend.position = "top")

ggsave('./results/figure1/fig1f.pdf', pf, units='cm', height = 10, width = 8, useDingbats=F)
ggsave('./results/figure1/fig1f.png', pf, units='cm', height = 10, width = 8)

saveRDS(pf_dat, file='./results/source_data/f1/revnumbers.rds')

saveRDS(revgenes, file='./results/source_data/f1/revgenes.rds')
