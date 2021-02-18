library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
library(ggpubr)
library(scales)

library(RColorBrewer)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)

dc = readRDS('/home/black/Dropbox/projects/ageing/proc_data/figure_rds/rev.gse.rds')
dc_gse = emapplot(dc, color = 'NES', layout = 'nicely', showCategory = 25)
dc_gse$layers = dc_gse$layers[-3]
dc_gse = dc_gse +
  geom_text_repel(aes(x = x, y = y, label = name), size = 8/pntnorm, max.overlaps = 15) +
  scale_color_gradient2(low = scales::muted('darkred'), high = 'white') +
  guides(fill = guide_legend(title = c('NES','size') ) )
dc_gse$labels$colour = 'NES'
#dc_gse$labels$size = 'size'

ggsave('results/dc_gsea.pdf', dc_gse,units = 'cm', width = 18, height = 12, useDingbats = F)
ggsave('results/dc_gsea.png', dc_gse,units = 'cm', width = 18, height = 12)

