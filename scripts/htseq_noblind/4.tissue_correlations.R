library(tidyverse)
library(openxlsx)
library(corrplot)
library(ggpubr)
library(RColorBrewer)

pntnorm <- (1/0.352777778)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('Up','Down'))
#source('./scripts/functions.R')

########################################
######################################## 
########################################  pairwise tissue expression correlations to confirm CoV
########################################

expr = readRDS('./data/htseq/expr_mat_no_blind.rds')
ind_id = readRDS('./data/processed/raw/individual_ids.rds')
age = readRDS('./data/processed/raw/ages.rds')
samp_tissue = readRDS('./data/processed/raw/tissue_ids.rds')
colnames(expr) = ind_id
names(age) = ind_id

########################
uage = age[16:31]
dage = uage[uage<90]
aage = uage[uage>90]

# pairwise tissue expression correlations, and change with age
exp= expr
chs = combn(4,2)
pcors = list()
for(i in 1:6){
  tsx = chs[,i]
  chsts = unique(samp_tissue)
  pname = paste(chsts[tsx[1]], chsts[tsx[2]], sep = '-')
  e1 = exp[, samp_tissue == chsts[tsx[1]] ]
  e2 = exp[, samp_tissue == chsts[tsx[2]] ]
  sameind = intersect(colnames(e1), colnames(e2))
  pcors[[pname]] = sapply(sameind, function(x){ cor(e1[,x], e2[,x], m='s')})
}
names(pcors) = toupper(names(pcors))

agecors.dev = sapply(names(pcors), function(x){
  c(round(cor.test(pcors[[x]][names(dage)], dage, m='s')$est,2),
    round(cor.test(pcors[[x]][names(dage)], dage, m='s')$p.val,3))
})
rownames(agecors.dev)[2] = 'p.val'
agecors.aging = sapply(names(pcors), function(x){
  c(round(cor.test(pcors[[x]][names(aage)], aage, m='s')$est,2),
    round(cor.test(pcors[[x]][names(aage)], aage, m='s')$p.val,3))
})
rownames(agecors.aging)[2] = 'p.val'

#####
pexpcors = reshape2::melt(pcors) %>%
  set_names('rho', 'pair') %>%
  mutate(id = unname(unlist((sapply(pcors,names))))) %>%
  left_join(data.frame(uage, id = names(uage)) ) 

saveRDS(pexpcors, './data/htseq/no_blind/pairwise_tissue_expression_cors.rds')

####################
#################### 
#################### 
#################### 
####################

expr_ch = readRDS('./data/htseq/no_blind/expression_change.rds')
expr_chqn = readRDS('./data/processed/tidy/expression_change.rds')

sig_genes = expr_ch %>%
  filter(FDR < 0.1) %>%
  mutate(direction = `Expression Change` > 0) %>%
  mutate(direction = ifelse( direction == TRUE, 'Up', 'Down')) %>%
  mutate(period = gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels = c('Development','Ageing'))) %>%
  dplyr::select(-p, -FDR, -`Expression Change`)

sig_genesqn = expr_chqn %>%
  filter(FDR < 0.1) %>%
  mutate(direction = `Expression Change` > 0) %>%
  mutate(direction = ifelse( direction == TRUE, 'Up', 'Down')) %>%
  mutate(period = gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels = c('Development','Ageing'))) %>%
  dplyr::select(-p, -FDR, -`Expression Change`)


siggenes_overlap = sig_genes %>%
  group_by(gene_id, period, direction) %>%
  summarise(n = length(tissue) ) %>%
  group_by(n, period, direction) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  slice(-c(1:4)) %>%
  ggplot(aes(x = n, y = count, fill = direction)) +
  facet_wrap(~period) +
  geom_bar(stat='identity', position = 'dodge') +
  scale_fill_manual(values = regcol, drop =F) +
  scale_y_continuous(trans = 'log10') +
  geom_text(aes(label = count), color = 'gray15', position = position_dodge(width = 1), angle = 90,
            hjust = 1.1, size = 6/pntnorm) +
  guides(fill = guide_legend('Direction Of Expression Change',
                             override.aes = list(size=2))) +
  theme(legend.position = 'bottom',
        legend.background = element_rect(fill = 'gray85', color = 'gray25')) +
  xlab(NULL) + ylab(NULL) +
  ggtitle('VST')

siggenes_overlap2 = sig_genes %>%
  inner_join(sig_genesqn) %>% 
  group_by(gene_id, period, direction) %>%
  summarise(n = length(tissue) ) %>%
  group_by(n, period, direction) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  slice(-c(1:4)) %>%
  ggplot(aes(x = n, y = count, fill = direction)) +
  facet_wrap(~period) +
  geom_bar(stat='identity', position = 'dodge') +
  scale_fill_manual(values = regcol, drop =F) +
  scale_y_continuous(trans = 'log10') +
  geom_text(aes(label = count), color = 'gray15', position = position_dodge(width = 1), angle = 90,
            hjust = 1.1, size = 6/pntnorm) +
  guides(fill = guide_legend('Direction Of Expression Change',
                             override.aes = list(size=2))) +
  theme(legend.position = 'bottom',
        legend.background = element_rect(fill = 'gray85', color = 'gray25')) +
  xlab(NULL) + ylab(NULL) +
  ggtitle('Overlap btw. VST and QN')

#ggarrange(siggenes_overlap, siggenes_overlap2, ncol=2)
ggsave('./results/htseq/no_blind/Figure_S5a.pdf', siggenes_overlap, width = 15, height = 8, units='cm',
       useDingbats = F )
ggsave('./results/htseq/no_blind/Figure_S5a.png', siggenes_overlap, width = 15, height = 8, units='cm' )  
