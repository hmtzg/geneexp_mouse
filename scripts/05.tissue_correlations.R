library(tidyverse)
library(openxlsx)
library(corrplot)
library(ggpubr)
library(RColorBrewer)
pntnorm <- (1/0.352777778)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
source('./scripts/functions.R')

########################################
######################################## 
########################################  pairwise tissue expression correlations to confirm CoV
########################################

expr = readRDS('./data/processed/raw/expression.rds')
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
exp = expr
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

saveRDS(pexpcors, './data/processed/tidy/pairwise_tissue_expression_cors.rds')

########################################
######################################## Permutation Test of Pairwise Tissue Expression-Age Correlations:
######################################## no sig. cutoff (test significance of fig 1c. for each pair)
########################################

expr_ch = readRDS('./data/processed/tidy/expression_change.rds')
aging.perm = readRDS('./data/processed/raw/permutations_ageing.rds')
dev.perm = readRDS('./data/processed/raw/permutations_dev.rds')
# calculate observed correlations:
cors = expr_ch %>%
  select(-p, -FDR) %>%
  pivot_wider(names_from = c(tissue, period), values_from=`Expression Change`) %>% 
  select(-gene_id) %>% 
  as.matrix() %>%
  cor(method='s', use ='pairwise.complete.obs')

colnames(cors) = sub('\\_.*','',colnames(cors))
rownames(cors) = sub('\\_.*','',rownames(cors))

na.genes = names(which(is.na(aging.perm$cortex[,1])))
aging.perm = lapply(aging.perm,function(x) x[!rownames(x)%in%na.genes,])
dev.perm = lapply(dev.perm,function(x) x[!rownames(x)%in%na.genes,])
names(dev.perm)
# reorder cors as permutation tissue order:
cors = cors[c(1,3,2,4,5,7,6,8),c(1,3,2,4,5,7,6,8)]

# compare development vs ageing:
dcors = cors[1:4,1:4][lower.tri(cors[1:4,1:4])]
acors = cors[5:8,5:8][lower.tri(cors[5:8,5:8])]
wilcox.test(dcors, acors, paired = T)
#V = 16, p-value = 0.3125

###
chs = combn(4,2)
corperm.dev = c()
for(i in 1:6){
  k = chs[,i]
  pair = sapply(1:1000,function(j){
    cor.test( dev.perm[[ k[1] ]][,j], dev.perm[[ k[2] ]][,j], m="s" )$est
  })
  corperm.dev = cbind(corperm.dev,pair)
  colnames(corperm.dev)[i] = paste(colnames(cors)[k[1]],colnames(cors)[k[2]],sep="-" )
}
saveRDS(corperm.dev, file = './data/processed/raw/pairwise_tissue_cors_perm_dev.rds')

corperm.aging = c()
for(i in 1:6){
  k = chs[,i]
  pair = sapply(1:1000,function(j){
    cor.test( aging.perm[[ k[1] ]][,j], aging.perm[[ k[2] ]][,j], m="s" )$est
  })
  corperm.aging = cbind(corperm.aging, pair)
  colnames(corperm.aging)[i] = paste(colnames(cors)[k[1]],colnames(cors)[k[2]],sep="-" )
}
saveRDS(corperm.aging, file = './data/processed/raw/pairwise_tissue_cors_perm_aging.rds')

dev.perm.results = c()
chs = combn(4,2)
for(i in 1:6){
  k = chs[,i]
  p1 = mean( corperm.dev[,i] >= cors[k[2], k[1] ] )
  f1 = median(corperm.dev[,i]) / cors[k[2], k[1] ]
  dev.perm.results = cbind(dev.perm.results,c(p1,f1))
}
colnames(dev.perm.results) = colnames(corperm.dev)
rownames(dev.perm.results) = c("permutation p-value","eFPP")

aging.perm.results = c()
chs = combn(4,2)
for(i in 1:6){
  k = chs[,i]
  p1 = mean( corperm.aging[,i] >= cors[k[2]+4, k[1]+4 ] )
  f1 = median(corperm.aging[,i]) / cors[k[2]+4, k[1]+4 ]
  aging.perm.results = cbind(aging.perm.results,c(p1,f1))
}
colnames(aging.perm.results) = colnames(corperm.dev)
rownames(aging.perm.results) = c("permutation p-value","eFPP")

pwise_cor_perm_test = 
  as.data.frame(rbind( round(t(rbind(rho = unname(c(cors[2:4,1], cors[3:4,2],cors[4,3])), 
                                     dev.perm.results)),3),
                       round(t(rbind(rho = unname(c(cors[6:8,5], cors[7:8,6],cors[8,7])), 
                                     aging.perm.results)),3))) %>% 
  mutate(`Tissue A` = rep(sapply(strsplit(colnames(dev.perm.results), '-'),function(x)x[1]),2)) %>%
  mutate(`Tissue B` = rep(sapply(strsplit(colnames(dev.perm.results), '-'),function(x)x[2]),2)) %>%
  #set_names(c('rho', 'p.val', 'FDR','Tissue A','Tissue B')) %>%
  mutate(Period = rep(c('Development', 'Ageing'), each = 6)) %>%
  relocate(`Tissue A`, `Tissue B`, Period, rho, `permutation p-value`, eFPP)
rownames(pwise_cor_perm_test) = NULL
saveRDS(pwise_cor_perm_test,'./data/processed/tidy/pwise_expch_cor_perm_test.rds')
saveRDS(pwise_cor_perm_test,'results/source_data/f1/f1pwise_expch_cor_perm_test.rds')
#write.xlsx(pwise_cor_perm_test, file='./results/SI_tables/TableS3.xlsx', row.names=T)

########################################
######################################## Permutation test for shared up/down genes (no significance cutoff )
######################################## in dev/aging among tissues: 
######################################## 2,3,4-tissue overlap test separately 

## observed 2,3,4-tissue overlaps in development and ageing
obs_overlap = expr_ch %>%
  filter(`Expression Change` !=0 ) %>% # drop expression - age rho = 0 values
  mutate(direction = `Expression Change` > 0) %>%
  mutate(direction = ifelse( direction == TRUE, 'Up', 'Down')) %>%
  mutate(period = gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels = c('Development','Ageing'))) %>%
  group_by(gene_id,period, direction) %>%
  summarise(`N Tissue` = length(tissue)) %>%
  group_by(`N Tissue`, period, direction) %>% 
  summarise(Obs = n()) %>%
  ungroup() %>%
  slice(-c(1:4)) %>%
  mutate(`N Tissue` = paste(`N Tissue`, c('Tissues')) )

perm_overlaps = list()
n_tissue = setNames(c(2,3,4), c('2 Tissues','3 Tissues','4 Tissues'))
for(i in 1:ncol(dev.perm$cortex)){
  permdevcor_i = sapply(dev.perm,function(x) x[,i])
  permagingcor_i = sapply(aging.perm,function(x)x[,i])
  pname = paste0('p',i)
  
  N_overlap = lapply(n_tissue, function(x){
    cbind( Development = c(Up = sum(rowSums(permdevcor_i > 0) == x),
                   Down = sum(rowSums(permdevcor_i < 0) == x)),
           Ageing = c(Up = sum(rowSums(permagingcor_i > 0) == x),
                     Down = sum(rowSums(permagingcor_i < 0) == x)) ) 
  })
  perm_overlaps[[ pname ]] = N_overlap
}

perm_overlaps = reshape2::melt(perm_overlaps) %>% 
  set_names(c('direction', 'period', 'N Overlap', 'N Tissue','permutation'))

saveRDS(perm_overlaps, file='./data/processed/raw/tissue_gene_overlaps_perm.rds')

overlap_test = perm_overlaps %>%
  left_join(obs_overlap) %>%
  group_by(direction, period, `N Tissue`) %>%
  summarise(Perm_p = mean(`N Overlap` >= unique(Obs)),
            eFPP =  round(median(`N Overlap`)/ unique(Obs),2)) %>%
  left_join(obs_overlap)

########################################
######################################## Permutation test for shared up/down significant genes (FDR<0.1)
######################################## in dev/aging among tissues: 
######################################## 2,3,4-tissue overlap test separately 

obs_overlap_fdr = expr_ch %>%
  #filter(`Expression Change` !=0 ) %>% # drop expression - age rho = 0 values
  filter(FDR < 0.1) %>%
  mutate(direction = `Expression Change` > 0) %>%
  mutate(direction = ifelse( direction == TRUE, 'Up', 'Down')) %>%
  mutate(period = gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels = c('Development','Ageing'))) %>%
  group_by(gene_id,period, direction) %>%
  summarise(`N Tissue` = length(tissue)) %>%
  group_by(`N Tissue`, period, direction) %>% 
  summarise(Obs = n()) %>%
  ungroup() %>%
  slice(-c(1:4))  %>%
  mutate(`N Tissue` = paste(`N Tissue`, c('Tissues')) )

choose_n = expr_ch %>%
  filter(FDR < 0.1) %>%
  mutate(direction = `Expression Change` > 0) %>%
  mutate(direction = ifelse( direction == TRUE, 'Up', 'Down')) %>%
  mutate(period = gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels = c('Development','Ageing'))) %>%
  group_by(tissue, period, direction) %>%
  summarise(`N Significant` = length(gene_id))

#####

names(dev.perm) = str_to_title(names(dev.perm))
names(aging.perm) = str_to_title(names(aging.perm))

overlap_fdr_f = function(perm, Period = 'Development', Direction = 'Up', p = 1){
  # function to calculate 2,3,4 tissue overlaps with the chosen number of top highest or lowest genes
  # requires choose_n data to choose the number of genes
  px_genes = sapply(names(perm),function(x){ # requires choose_n data
    obs_sig = filter(choose_n, period==Period & direction == Direction & tissue == x) %>% pull(`N Significant`)
    if(length(obs_sig)==0) return(NULL) else {
      if(Direction == 'Up') {
        return(names(sort(perm[[x]][,p],decreasing = T)[1:obs_sig ]) )
      } else if(Direction == 'Down'){
        return(names(sort(perm[[x]][,p],decreasing = F)[1:obs_sig ]) )
      }
    }
  })
  px_unique = unique(unname(unlist(px_genes)))
  overlaps = table( factor(sapply(1:length(px_unique), function(pxu) {
    length(unlist(sapply(1:4,function(x){ intersect(px_unique[pxu], px_genes[[x]]) }) ))
  }), levels = c(1,2,3,4)))[-1]
  return(overlaps)
}

perm_overlaps_fdr = list()
for(n in 1:ncol(dev.perm$Cortex)){
  print(n/ncol(dev.perm$Cortex))
  N_overlap = list( Development = 
                      list(Up = overlap_fdr_f(dev.perm, Period = 'Development', Direction = 'Up', p = n),
                           Down = overlap_fdr_f(dev.perm, Period = 'Development', Direction = 'Down', p = n) ),
                    Ageing = 
                      list(Up = overlap_fdr_f(aging.perm, Period = 'Ageing', Direction = 'Up', p = n),
                           Down = overlap_fdr_f(aging.perm, Period = 'Ageing', Direction = 'Down', p = n) ) )
  
  pname = paste0('p',n)
  perm_overlaps_fdr[[ pname ]] = N_overlap
}

perm_overlaps_fdr = reshape2::melt(perm_overlaps_fdr) %>%
  set_names(c('N Tissue', 'N Overlap', 'direction', 'period', 'permutation'))

saveRDS(perm_overlaps_fdr, './data/processed/raw/tissue_siggene_fdr01_overlaps_perm.rds')

overlap_fdr_test = perm_overlaps_fdr %>%
  mutate(`N Tissue` = paste(`N Tissue`, c('Tissues')) ) %>% 
  right_join(obs_overlap_fdr, by = c('N Tissue', 'period', 'direction')) %>% 
  group_by(direction, period, `N Tissue`) %>%
  summarise(Perm_p = mean(`N Overlap` >= unique(Obs)),
            eFPP =  round(median(`N Overlap`)/ unique(Obs),2)) %>%
  ungroup() %>%
  mutate(Perm_p = ifelse(Perm_p == 0, '< 0.001', Perm_p)) %>%
  left_join(obs_overlap_fdr)

####################
#################### GO Over representation  Analysis with topGO package:
#################### foreground: shared significant up/down genes among tissues
#################### background: all significant up/down genes across tissues
####################

sig_genes = expr_ch %>%
  filter(FDR < 0.1) %>%
  mutate(direction = `Expression Change` > 0) %>%
  mutate(direction = ifelse( direction == TRUE, 'Up', 'Down')) %>%
  mutate(period = gsub("aging", "ageing", period),
         period = str_to_title(period),
         period = factor(period, levels = c('Development','Ageing'))) %>%
  dplyr::select(-p, -FDR, -`Expression Change`)

shared_genes_go = sapply(unique(sig_genes$period), function(p){
  sapply(unique(sig_genes$direction), function(d){
    
    foreground = sig_genes %>%
      filter(period==p, direction == d) %>%
      group_by(gene_id) %>%
      summarise(genelist = length(tissue)==4) %>% 
      filter(genelist == T) %>%
      pull(gene_id) %>% as.character()
    
    genelist = sig_genes %>%
      filter(period==p) %>%
      mutate(genelist = 0) %>%
      dplyr::select(gene_id, genelist) %>% 
      distinct() %>%
      pull(genelist, name= gene_id)
    
    genelist[names(genelist)%in%foreground] = 1
    if(sum(genelist)>=10){
      go = go_bp_enrich.test.Mm(genelist = genelist, selection = 1,padj = 'BH')
      return(go) } 
  }, simplify = F, USE.NAMES = T)
}, simplify = F, USE.NAMES = T)

saveRDS(shared_genes_go,
        'data/processed/raw/shared_expch_gora.rds')
names(shared_genes_go) =  unique(sig_genes$period)

table_sX = unlist(shared_genes_go, recursive = F)

write.xlsx(table_sX, 
           file='./results/supplementary_files/Supplementary_File_2.xlsx', row.names=F)


####################
####################
####################
####################
