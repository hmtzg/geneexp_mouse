library(tidyverse)
library(openxlsx)
library(ggpubr)
library(RColorBrewer)

theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
pntnorm <- (1/0.352777778)
tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
varcol = setNames(c('dodgerblue','firebrick3'),c('div','con'))
regcol = setNames(c('rosybrown3','paleturquoise3'),c('up','down'))
revcol = setNames(c('brown4', '#1C7AD9', 'indianred', '#6FADEC'), c('UpDown','DownUp','UpUp','DownDown'))

source('./scripts/functions.R')
dev = readRDS("./data/processed/raw/development_expression_change.rds")
aging = readRDS("./data/processed/raw/ageing_expression_change.rds")

revg = sapply(names(dev), function(x) revgenes(dev[[x]][,1],aging[[x]][,1], monotonic = T), simplify = F )

saveRDS(revg,'./data/processed/raw/revgenes.tissues.rds')

revgenes = revg %>%
  reshape2::melt(revgenes) %>%
  set_names('gene_id','direction','tissue')

saveRDS(revgenes,'./data/processed/tidy/revgenes.tissues.rds')


########################################
######################################## GORA for reversal in each tissue:
######################################## foreground: up-down (or down-up)
######################################## background: up-up (or  down-down)

updown.go = list() # enrichment for UD reversal genes vs UD + UU
upup.go = list() # enrichment for UU genes vs UU+UD (depletion for UD)

for(i in 1:4){
  bg = setNames(rep(0, length( c(revg[[i]]$UpDown, revg[[i]]$UpUp))),
                nm = c(revg[[i]]$UpDown, revg[[i]]$UpUp) )
  bg[names(bg)%in%revg[[i]]$UpDown] = 1
  go = go_bp_enrich.test.Mm(bg, selection = 1, padj = 'BH')
  go2 = go_bp_enrich.test.Mm(bg, selection = 0, padj = 'BH')
  updown.go [[ names(dev)[i] ]] = go
  upup.go[[ names(dev)[i] ]] = go2
}
saveRDS(updown.go, './data/processed/raw/updown_gora_each_tissue.rds')
saveRDS(upup.go, './data/processed/raw/upup_gora_each_tissues.rds')

downup.go = list() # enrichment for DU reversal genes vs DU+DD
downdown.go = list() # enrichment for DD genes vs DU+DD (depletion for DU)

for(i in 1:4){
  bg = setNames(rep(0, length(c(revg[[i]]$DownUp, revg[[i]]$DownDown))),
                nm = c(revg[[i]]$DownUp, revg[[i]]$DownDown))
  bg[names(bg)%in%revg[[i]]$DownUp] = 1
  go = go_bp_enrich.test.Mm(bg, selection = 1, padj = 'BH')
  go2 = go_bp_enrich.test.Mm(bg, selection = 0, padj = 'BH')
  downup.go[[ names(dev)[i] ]] = go
  downdown.go[[ names(dev)[i] ]] = go2
}
saveRDS(downup.go, './data/processed/raw/downup_gora_each_tissue.rds')
saveRDS(downdown.go, './data/processed/raw/downdown_gora_each_tissue.rds')

reversal_go_each_tissue = list(UpDown = updown.go, DownUp = downup.go)
table_sX = unlist(reversal_go_each_tissue, recursive = F)

write.xlsx(table_sX, file='./results/SI_tables/TableS5.xlsx', row.names=F)

########################################
######################################## Test significance of reversal genes in each tissue separately
########################################
########################################

devcors = readRDS("./data/processed/raw/development_expression_change.rds")
agingcors = readRDS("./data/processed/raw/ageing_expression_change.rds")
aging.perm = readRDS("./data/processed/raw/permutations_ageing.rds")
lapply(aging.perm,function(x) sum(is.na(x))) # 3 genes NA in cortex
na.genes = names(which(is.na(aging.perm$cortex[,1])))
aging.perm = lapply(aging.perm,function(x) x[!rownames(x)%in%na.genes,])
devcors = lapply(devcors, function(x) x[!rownames(x)%in%na.genes,])
agingcors = lapply(agingcors, function(x) x[!rownames(x)%in%na.genes,])

devcor = sapply(devcors, function(x) x[,1])

# UD/UU: test each tissue for reversal separately:
UD.tissue = list()
for(ts in 1:4){
  devups.ts = devcor[,ts]>0
  obs.ud.ts = length(revg[[ts]]$UpDown) / (length(revg[[ts]]$UpDown)+length(revg[[ts]]$UpUp))
  UDdist.ts = c()
  for(i in 1:1000){
    agingX = aging.perm[[ts]][devups.ts,i]
    UDX = sum(agingX<0)
    UUX = sum(agingX>0)
    UDdist.ts[i] = UDX / (UDX + UUX)
  }
  pUD.ts = mean(UDdist.ts >= obs.ud.ts)
  efpp.ts = median(UDdist.ts) / obs.ud.ts
  UD.tissue[[ts]] = list(dist = UDdist.ts, obs = round(obs.ud.ts,2), pval = round(pUD.ts,2), 
                         efpp = round(efpp.ts,2))
}
names(UD.tissue) = colnames(devcor)

# DU/DD: test each tissue for reversal separately:
DU.tissue = list()
for(ts in 1:4){
  devdowns.ts = devcor[,ts]<0
  obs.du.ts = length(revg[[ts]]$DownUp) / (length(revg[[ts]]$DownUp)+length(revg[[ts]]$DownDown))
  DUdist.ts = c()
  for(i in 1:1000){
    agingX = aging.perm[[ts]][devdowns.ts,i]
    DUX = sum(agingX>0)
    DDX = sum(agingX<0)
    DUdist.ts[i] = DUX / (DUX + DDX)
  }
  pDU.ts = mean(DUdist.ts >= obs.du.ts)
  efpp.ts = median(DUdist.ts) / obs.du.ts
  DU.tissue[[ts]] = list(dist=DUdist.ts, obs = round(obs.du.ts,2), pval = round(pDU.ts,2),
                         efpp = round(efpp.ts,2))
}
names(DU.tissue) = colnames(devcor)

saveRDS(UD.tissue, file = './data/processed/raw/updown_perm_each.rds')
saveRDS(DU.tissue, file = './data/processed/raw/downup_perm_each.rds')

rev.tissue = list(UpDown = lapply(UD.tissue,function(x) x$dist),
                  DownUp = lapply(DU.tissue, function(x)x$dist))

rev_test =  list(UpDown = lapply(UD.tissue, function(x) x[-1] ),
                 DownUp = lapply(DU.tissue, function(x) x[-1] ) ) 
  
rev_test  = reshape2::melt(rev_test) %>% 
  set_names(c('value', 'params', 'tissue', 'Reversal')) %>% 
  spread(key=params, value=value) %>%
  mutate(Reversal = factor(Reversal, levels= c('UpDown', 'DownUp')) )
  
sizex = 2.5
rev_each_plot = reshape2::melt(rev.tissue) %>% 
  set_names(c('value', 'tissue', 'Reversal')) %>% 
  mutate(Reversal = factor(Reversal, levels= c('UpDown', 'DownUp')) ) %>%
  ggplot(aes(x = value)) +
  facet_grid(tissue~Reversal, scales = 'free_x') +
  geom_histogram() +
  ylab('Frequency') +
  xlab('Reversal Proportion') +
  geom_vline(data = rev_test, mapping = aes(xintercept= obs), linetype ='dashed', color = 'darkred' ) +
  geom_text(data= rev_test, size=  sizex,
            mapping = aes(x = c(0.44, 0.58,0.6, 0.67, 0.5, 0.65, 0.6, 0.65) , y = 90, 
                          label = paste('Obs =', obs) )) +
  geom_text(data= rev_test,size=  sizex,
            mapping = aes(x = c(0.44, 0.58,0.6, 0.67, 0.5, 0.65, 0.6, 0.65) , y=80, 
                          label = paste('eFPP =', efpp) )) +
  geom_text(data= rev_test,size=  sizex,
            mapping = aes(x = c(0.44, 0.58,0.6, 0.67, 0.5, 0.65, 0.6, 0.65) , y=70, 
                          label = paste('p.val =', pval))) +
  theme_bw() +
  theme(axis.text = element_text(size =5))

ggsave('./results/figure_supplements/f1s/FS8.pdf', rev_each_plot, width = 16, height = 15, units='cm', 
       useDingbats = F )
ggsave('./results/figure_supplements/f1s/FS8.png', rev_each_plot, width = 16, height = 15, units='cm' )  

########################################
######################################## Test significance of shared reversal genes among tissues
########################################
########################################

obs_rev_overlap = revgenes %>%
  group_by(gene_id, direction) %>%
  summarise(n = length(tissue)) %>%
  filter(n == 4) %>%
  group_by(direction) %>%
  summarise(n= length(gene_id))

# UD/UU test: test among development UP genes (background), use these genes in ageing for null distribution:
updown_obs = as.numeric(obs_rev_overlap[3,2] / (obs_rev_overlap[3,2] + obs_rev_overlap[4,2]) )
upup_obs = as.numeric(obs_rev_overlap[4,2] / (obs_rev_overlap[3,2] + obs_rev_overlap[4,2]) )

UDdist = c()
UUdist = c()
devup = names(which(rowSums(devcor > 0)==4))
for(i in 1:1000){
  agingX = sapply(aging.perm, function(x) x[devup,i])
  UUX = sum(rowSums(agingX>0) == 4)
  UDX = sum(rowSums(agingX<0) == 4)
  UDdist[i] = UDX / (UUX + UDX)
  UUdist[i] = UUX / (UUX + UDX)
}

# DU/DD test: test among development DOWN genes(background), use these genes in ageing for null:
downup_obs = as.numeric(obs_rev_overlap[2,2] / (obs_rev_overlap[2,2] + obs_rev_overlap[1,2]) )
downdown_obs = as.numeric(obs_rev_overlap[1,2] / (obs_rev_overlap[2,2] + obs_rev_overlap[1,2]) )

DUdist = c()
DDdist = c()
devdown = names(which(rowSums(devcor < 0)==4))
for(i in 1:1000){
  agingX = sapply(aging.perm, function(x) x[devdown,i])
  DUX = sum(rowSums(agingX>0) == 4)
  DDX = sum(rowSums(agingX<0) == 4)
  DUdist[i] = DUX / (DUX + DDX)
  DDdist[i] = DDX / (DUX + DDX)
}

rev_shared_test = reshape2::melt(data.frame( pval = round(c(mean(UDdist >= updown_obs),
                                                            mean(DUdist >= downup_obs)),2),
            efpp = round(c(median(UDdist) / updown_obs, median(DUdist) / downup_obs ),2),
            obs = round(c(updown_obs,downup_obs),2), Reversal = c('UpDown','DownUp') )) %>%
  spread(key=variable, value=value) %>%
  mutate(Reversal = factor(Reversal, levels = c('UpDown', 'DownUp')))

shared_rev_perm_test = reshape2::melt(data.frame( UpDown = UDdist, DownUp = DUdist)) %>%
  set_names(c('Reversal','Proportion'))
saveRDS(shared_rev_perm_test, './data/processed/tidy/shared_rev_perm_test.rds')

rev_shared_plot = reshape2::melt(data.frame( UpDown = UDdist, DownUp = DUdist)) %>%
  set_names(c('Reversal','value')) %>%
  mutate(Reversal = factor(Reversal, levels = c('UpDown', 'DownUp'))) %>%
  ggplot(aes(x = value)) +
  facet_grid(~Reversal, scales = 'free_x') +
  geom_histogram() +
  ylab('Frequency') +
  xlab('Reversal Proportion') +
  geom_vline(data = rev_shared_test, mapping = aes(xintercept= obs), linetype ='dashed', 
             color = 'darkred' ) +
  geom_text(data= rev_shared_test, size=  sizex,
            mapping = aes(x = c(0.6, 0.72) , y = 65, label = paste('Obs =', obs) )) +
  geom_text(data= rev_shared_test,size=  sizex,
            mapping = aes(x = c(0.6, 0.72) , y=60, label = paste('eFPP =', efpp) )) +
  geom_text(data= rev_shared_test, size=  sizex,
            mapping = aes(x = c(0.6, 0.72) , y=56, label = paste('p.val =', pval))) +
  theme_bw() +
  theme(axis.text = element_text(size =5))

ggsave('./results/figure_supplements/f1s/FS9.pdf', rev_shared_plot, width = 12, height = 8, units='cm', 
       useDingbats = F )
ggsave('./results/figure_supplements/f1s/FS9.png', rev_shared_plot, width = 12, height = 8, units='cm' )  


########################################
######################################## GORA enrichment of shared reversal genes across tissues
########################################
########################################

rev.common = sapply(names(revg$Cortex), function(y) Reduce(intersect,lapply(revg,function(x) x[[y]] )))
revgenes_common = reshape2::melt(rev.common) %>% 
  set_names(c('gene_id','Pattern'))

saveRDS(revgenes_common,'./data/processed/tidy/revgenes.shared.rds')

bg = setNames(rep(0, length(c(rev.common$UpDown, rev.common$UpUp)) ),
              nm = c(rev.common$UpDown, rev.common$UpUp) )
bg[names(bg)%in%rev.common$UpDown] = 1

updown_shared_go = go_bp_enrich.test.Mm(bg, selection = 1, padj = 'BH')
upup_shared_go = go_bp_enrich.test.Mm(bg, selection = 0, padj = 'BH')

bg = setNames(rep(0, length(c(rev.common$DownUp,rev.common$DownDown))),
              nm = c(rev.common$DownUp,rev.common$DownDown) )
bg[names(bg)%in%rev.common$DownUp] = 1

downup_shared_go = go_bp_enrich.test.Mm(bg, selection = 1, padj = 'BH')
downdown_shared_go = go_bp_enrich.test.Mm(bg, selection = 0, padj = 'BH')

shared_reversal_go = list(updown_shared_go = updown_shared_go,
     upup_shared_go = upup_shared_go, 
     downup_shared_go = downup_shared_go,
     downdown_shared_go = downdown_shared_go )

saveRDS(shared_reversal_go,
        './data/processed/raw/shared_reversal_gora.rds')

################