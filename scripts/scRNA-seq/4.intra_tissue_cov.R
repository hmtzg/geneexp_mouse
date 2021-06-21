library(tidyverse)
library(ggpubr)
pntnorm <- (1/0.352777778)
library(ggpubr)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
tissuecode = rep(c("#592389", "#F49292", "#FF0600","#E0CE00"), times = c(15,16,16,16))
ct.exp = readRDS(file="./data/other_datasets/scRNA-seq/processed/celltype_expr_per_ind.rds")

# get same cell types from age groups in each tissue:
for(i in 1:4){
  scelltype = Reduce(intersect,list(names(ct.exp$m3[[i]]),
                                    names(ct.exp$m18[[i]]),
                                    names(ct.exp$m24[[i]])))
  ct.exp$m3[[i]] = ct.exp$m3[[i]][scelltype]
  ct.exp$m18[[i]] = ct.exp$m18[[i]][scelltype]
  ct.exp$m24[[i]] = ct.exp$m24[[i]][scelltype]
}

#### 
# get cell types that are found in at least 2 individuals in every time point: (intracov1)
for(i in 1:4){
  # cell types found in at least 3 individuals points:
  chs.ct = sapply(ct.exp$m3[[i]], ncol) > 1 &
    sapply(ct.exp$m18[[i]], ncol) > 1 &
    sapply(ct.exp$m24[[i]], ncol) > 1
  
  ct.exp$m3[[i]] = ct.exp$m3[[i]][ chs.ct ]
  ct.exp$m18[[i]] = ct.exp$m18[[i]][ chs.ct ]
  ct.exp$m24[[i]] = ct.exp$m24[[i]][ chs.ct ]
}

ind.id = sapply(names(ct.exp$m3), function(x){
  list(m3 = unique(unname(unlist(lapply(ct.exp$m3[[x]], colnames)))),
       m18 = unique(unname(unlist(lapply(ct.exp$m18[[x]], colnames)))),
       m24 = unique(unname(unlist(lapply(ct.exp$m24[[x]], colnames)))) )
}, simplify = F)

#ctexps = list(m3 = ct.exp.3m, m18 = ct.exp.18m, m24 = ct.exp.24m)

## convert list of cell type expressions to list of individual expressions:
for(i in 1:4){  # 3 months
  ts = sapply(ind.id[[i]]$m3, function(x){
    sapply(ct.exp$m3[[i]], function(y) {
      if(sum(colnames(y)%in%x)==0 ) {
        NULL 
      } else{
        y[,x]
      }
    }, simplify = F)
  }, simplify = F)
  ct.exp$m3[[i]] = sapply(names(ts), function(a) do.call(cbind, ts[[a]]), simplify = F ) 
}
for(i in 1:4){  # 18 months
  ts = sapply(ind.id[[i]]$m18, function(x){
    sapply(ct.exp$m18[[i]], function(y) {
      if(sum(colnames(y)%in%x)==0 ) {
        NULL 
      } else{
        y[, x]
      }
    }, simplify = F)
  }, simplify = F)
  ct.exp$m18[[i]] = sapply(names(ts), function(a) do.call(cbind, ts[[a]]), simplify = F ) 
}
for(i in 1:4){  # 24 months
  ts = sapply(ind.id[[i]]$m24, function(x){
    sapply(ct.exp$m24[[i]], function(y) {
      if(sum(colnames(y)%in%x)==0 ) {
        NULL 
      } else{
        y[, x]
      }
    }, simplify = F)
  }, simplify = F)
  ct.exp$m24[[i]] = sapply(names(ts), function(a) do.call(cbind, ts[[a]]), simplify = F ) 
}
# combine age groups within tissues:
exp = sapply(names(ct.exp$m3), function(x){
  list(m3 = ct.exp$m3[[x]], m18 = ct.exp$m18[[x]], m24 = ct.exp$m24[[x]] )
}, simplify = F)
names(exp) = str_to_title(names(exp))

## filter individuals having at least 3 cell types:
exp = lapply(exp, function(b) sapply(b, function(a) a[ sapply(a, function(x) ncol(x))>2] ))

###### calculate CoV among cell types for each individual in each tissue:
exp.cov = sapply(exp, function(z){
  sapply(z, function(a){
    sapply(a, function(x){
      apply(x, 1, function(y) sd(y)/ mean(y) )
    })
  })
}, simplify = F)

## calculate mean cov:
#exp.cov = exp.cov[-2]
mean.cov = sapply(exp.cov, function(x){
  sapply(x, function(y){ colMeans(y, na.rm = T) })
}, simplify = F)
med.cov = sapply(exp.cov, function(x){
  sapply(x, function(y){ apply(y,2, function(a) median(a, na.rm = T) ) })
}, simplify = F)
mean.cov = lapply(mean.cov, function(x) lapply(x, function(y) y[!is.na(y)]) )
# med.cov = lapply(med.cov, function(x) lapply(x, function(y) y[!is.na(y)]) )

# intracov1: each present cell type represented  in every age group in at least 2 individuals:
names(mean.cov)[4] = 'Cortex'

intra_tissue_mean_cov  = reshape2::melt(mean.cov) %>%
  set_names('mean CoV', 'age group', 'tissue')

saveRDS(intra_tissue_mean_cov, './data/other_datasets/scRNA-seq/processed/intra_tissue_mean_cov.rds')

####################################
# cell.types = list()
# for(i in 1:4){
#   nm = names(exp)
#   cell.types[[ nm[i] ]] = list( m3 =  Reduce(union, lapply(exp[[i]]$m3, colnames)),
#                                 m18 = Reduce(union, lapply(exp[[i]]$m18, colnames)),
#                                 m24 = Reduce(union, lapply(exp[[i]]$m24, colnames)) )
# }
# cell.types = reshape2::melt(lapply(cell.types, function(a) length(a[[1]]))) %>%
#   set_names('n', 'tissue')
# 
# cell.types[4,2] = 'Cortex'



###########################
# ggsave('Dropbox/projects/ageing/results.n/scRNA/intra_tissue_cov/intracov1.pdf', intracov1, units = 'cm',
#        width = 12, height = 8, useDingbats = F)
# ggsave('Dropbox/projects/ageing/results.n/scRNA/intra_tissue_cov/intracov1.png', intracov1, units = 'cm',
#        width = 12, height = 8)

# ## intra cov median: same result, no significant
# intracov2 = reshape2::melt(med.cov) %>%
#   set_names('median CoV', 'age group', 'tissue') %>%
#   mutate(`age group` = factor(`age group`, levels = c('m3', 'm18', 'm24')) ) %>%
#   mutate(age = as.numeric( gsub('[a-z]','',`age group`) ) ) %>%
#   ggplot(aes(x = `age group`, y = `median CoV`)) +
#   facet_wrap(~tissue, scales = 'free') +
#   geom_point(color='indianred4')  +
#   ylab('Median CoV of Genes Among Cell Types') +
#   xlab('Age (in months)')
# 
# 



# save(list=ls(),
#      file = 'Dropbox/projects/ageing/proc_data/scRNA/intra_tissue_CoV/intra_tissue_CoV2.rdata')



