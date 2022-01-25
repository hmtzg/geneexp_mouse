library(tidyverse)
library(openxlsx)
source("scripts/functions.R")

exp = readRDS('./data/processed/raw/expression.rds')
ind.id = readRDS('./data/processed/raw/individual_ids.rds')
age = readRDS('./data/processed/raw/ages.rds')
samp_tissue = readRDS('./data/processed/raw/tissue_ids.rds')

tissuecol = setNames(c('#233789', '#f49e92', '#801008','#dbb32e'),c('Cortex','Lung','Liver','Muscle'))
colnames(exp) = ind.id
names(age) = ind.id
# remove one individual that lack expression in cortex
# exp = exp[, !colnames(exp) == "465"]
# samp_tissue = samp_tissue[!names(age)=="465"]
# age = age[!names(age)=="465"]

eff.size.dev = sapply(unique(samp_tissue), function(y) {
  sapply(rownames(exp),function(x){
    cohens_d(exp[x, samp_tissue == y & age < 90], exp[x, (!samp_tissue == y) & age < 90])
  })
})
saveRDS(eff.size.dev, './data/processed/raw/eff.size.dev.rds')

expr_ch_dev = readRDS("./data/processed/raw/development_expression_change.rds")
expr_ch_ag = readRDS("./data/processed/raw/ageing_expression_change.rds")
cov_ch = readRDS('./data/processed/tidy/CoV_change.rds')
dgenes = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
dcgenes = names(dgenes[dgenes<0])

expch_dev = sapply(expr_ch_dev,function(x) x[,1] )
expch_ag = sapply(expr_ch_ag,function(x) x[,1])

# get tissue with highest ES for each gene:
ts.spec = colnames(eff.size.dev)[apply(eff.size.dev, 1, which.max)]
names(ts.spec) = rownames(eff.size.dev)

# get tissue specific genes using >3Q of genes assigned to a tissue:
ts.specQ3 = sapply(unique(ts.spec), function(x){
  ts.genes = names(which(ts.spec==x))
  cutoff = summary(eff.size.dev[ts.genes,x])['3rd Qu.']
  q3.genes = names(which(eff.size.dev[ts.genes,x] > cutoff))
  return(q3.genes)
})
saveRDS(ts.specQ3,'./data/processed/raw/ts.specQ3.genes.rds')

# get ES values for tissue specific genes:
ts.spec.ES = sapply(names(ts.specQ3),function(x) {
  eff.size.dev[ts.specQ3[[x]],]
}, simplify = F)

ts.spec.ES2 = reshape2::melt(ts.spec.ES) %>% # save for ggplot
  set_names('gene id', 'tissue', 'ES','spec')
saveRDS(ts.spec.ES2,'./data/processed/tidy/ts.spec.ES.rds')

# get exp-age rho values for tissue specific genes, both in development and ageing:
ts.spec.expch = sapply(names(ts.specQ3), function(x){
  a = list(Development = expch_dev[ts.specQ3[[x]],], Aging = expch_ag[ts.specQ3[[x]],])
},simplify = F)

# save exp-age rho values of tissue specific genes for ggplot:
ts.spec.expch2 = reshape2::melt(ts.spec.expch) %>%
  set_names('gene id', 'tissue', 'expage rho','period', 'spec') %>%
  mutate(period = factor(period, levels = c('Development','Aging')))
saveRDS(ts.spec.expch2, file = './data/processed/tidy/ts.spec.exprho.rds')

##
ts.spec.segments = ts.spec.expch2 %>% 
  pivot_wider(names_from = period, values_from = `expage rho`)

# tis. spec. genes development-ageing rho comparison wilcox test:
ts.spec.dev.agetest = ts.spec.expch2 %>%
  group_by(spec,tissue) %>%
  summarise( w.test = wilcox.test(`expage rho`~period, paired = T)$p.val)

##### in which tissue the highest expression change occurs for each gene: dont use rho
ts.expr.ch = colnames(expch_ag)[apply(expch_ag, 1, function(x) which.max(abs(x)))]
names(ts.expr.ch) = rownames(expch_ag)

### use beta from linear regressin:
cortex = exp[, age>90 & samp_tissue=='Cortex']
lung = exp[, age>90 & samp_tissue=='Lung']
liver = exp[, age>90 & samp_tissue=='Liver']
muscle = exp[, age>90 & samp_tissue=='Muscle']

cortage = log2(age[age>90 & samp_tissue=='Cortex'])
lungage = log2(age[age>90 & samp_tissue=='Lung'])
liverage = log2(age[age>90 & samp_tissue=='Liver'])
muscleage = log2(age[age>90 & samp_tissue=='Muscle'])

cortbeta = t(apply(cortex,1,function(x){
  summary(lm(x~cortage))$coef[2,c(1,4)]
}))
lungbeta = t(apply(lung,1,function(x){
  summary(lm(x~lungage))$coef[2,c(1,4)]
}))
liverbeta = t(apply(liver,1,function(x){
  summary(lm(x~liverage))$coef[2,c(1,4)]
}))
musclebeta = t(apply(muscle,1,function(x){
  summary(lm(x~muscleage))$coef[2,c(1,4)]
}))

expbeta = cbind(cortbeta[,1], lungbeta[,1], liverbeta[,1], musclebeta[,1])
colnames(expbeta) = c('Cortex', 'Lung', 'Liver', 'Muscle')
saveRDS(expbeta, 'data/processed/tidy/expbeta.rds')

ts.beta.ch = colnames(expbeta)[apply(expbeta, 1, function(x) which.max(abs(x)))]
names(ts.beta.ch) = rownames(expbeta)

# direction of expression change for those genes :
# ts.expr.ch.dir = sapply(1:length(ts.expr.ch), function(x){ sign(expch_ag[x, ts.expr.ch[x]]) } )
# names(ts.expr.ch.dir) = names(ts.expr.ch)

ts.beta.ch.dir = sapply(1:length(ts.beta.ch), function(x){ sign(expbeta[x, ts.beta.ch[x]]) } )
names(ts.beta.ch.dir) = names(ts.beta.ch)

########
######## Expr change in native tissue :
########
# for tissue-specific genes (>Q3) :
ts.specQ3.genes = unlist(ts.specQ3)
names(ts.specQ3.genes) = gsub('[0-9]','', names(ts.specQ3.genes))

mat = data.frame(sameness = names(ts.specQ3.genes) == ts.beta.ch[ts.specQ3.genes],
                 exp_dir = ts.beta.ch.dir[ts.specQ3.genes] )
table(mat)[,c(2,1)]
fisher.test(table(mat)[,c(2,1)])
# OR: 5.50, p = 2.2e-16

### using only DiCo genes:
# losing expr in native, gain in other tiss:
specsub = ts.specQ3.genes[ts.specQ3.genes%in%dcgenes]
expchsub = ts.beta.ch[dcgenes]
expchdirsub = ts.beta.ch.dir[dcgenes]
matsub = data.frame(sameness = names(specsub) == expchsub[specsub],
           expdir = expchdirsub[specsub])
table(matsub)[,c(2,1)]
sum(table(matsub)[,c(2,1)]) # 1287 genes
fisher.test(table(matsub)[,c(2,1)])
fisher.test(table(matsub)[,c(2,1)])$p.val
# OR = 74.81, p=5.9e-203
saveRDS(list(tbl=table(matsub)[,c(2,1)],
        fisher=fisher.test(table(matsub)[,c(2,1)])),
        file ='data/processed/specloss_wdico_fisher.rds')

## DiCo vs DiDi and tis spec vs non-tis spec
spec.pat.mat = data.frame(pat = names(dgenes)%in%dcgenes,
           ts.spec = names(dgenes)%in%ts.specQ3.genes)
table(spec.pat.mat )
fisher.test(table(spec.pat.mat ))
fisher.test(table(spec.pat.mat ))$p.val
# OR = 1.56, p=1.349306e-18
saveRDS(table(spec.pat.mat ), 
        'results/source_data/f3/tsspec_dico_enrichment.rds')

##################
##################
##################
##################
##################
### among div(significant)-conv genes:
# divconvg2 = names(which(aging.covage[dev.covage[,1]>0 & dev.covage[,3]<0.1,1 ] < 0)) # 1010 genes
# mat.dc2 = mat[rownames(mat)%in%divconvg2,] # 268 genes in total
# table(mat.dc2)[,c(2,1)]
# fisher.test(table(mat.dc2)[,c(2,1)])
# # OR: 49.45, pval < 1e10-16

# table(ts.spec[divconvg]) # tissue specific gene proportion of div-conv genes

# ########## development (confirmation):

# check expression change direction and in development with tissue specific genes:
# maximum expression change in which tissue:
ts.expr.ch.dev = colnames(expch_dev)[apply(expch_dev, 1, function(x)which.max(abs(x)))]
names(ts.expr.ch.dev) = rownames(expch_dev)
# direction of maximum expression change:
ts.expr.ch.dir.dev = sapply(1:length(ts.expr.ch.dev), function(x)sign(expch_dev[x,ts.expr.ch.dev[x]]))
names(ts.expr.ch.dir.dev) = rownames(expch_dev)

# among tissue assigned (all) genes:
mat.dev = data.frame(sameness = ts.spec == ts.expr.ch.dev, exp_dir = ts.expr.ch.dir.dev)
table(mat.dev)
fisher.test(table(mat.dev))
# OR: 1.32, pval: 3.5x10-15
# we see the opposite trend, genes specific to a tissue, tend to show highest increase in that tissue;
# genes specific to a tissue, tend to show highest decrease in non-specific tissue.

# among tissue specific genes (>Q3): actual result:
mat.dev2 = data.frame(sameness = names(ts.specQ3.genes) == ts.expr.ch.dev[ts.specQ3.genes],
                      exp_dir = ts.expr.ch.dir.dev[ts.specQ3.genes])
table(mat.dev2)
fisher.test(table(mat.dev2)) 
# OR: 1.85, pval: 2x10-16, more significant as expected

# repeat among div-conv genes:
table(mat.dev[rownames(mat.dev2)%in%dcgenes,])
fisher.test(table(mat.dev[rownames(mat.dev2)%in%dcgenes,]))
# OR: 1.28, pval: 4.2x10-5, still significant


###########################
###########################
###########################
###########################
###########################
########################### test tissue spec and reversal genes enrichment
###########################
###########################

# revg = readRDS('./data/processed/tidy/revgenes.tissues.rds')
# 
# ### background all tissue specific genes, test UD reversal
# i=c('Cortex','Lung','Liver','Muscle')[3]
# revi = filter(revg, tissue==i & direction =='UpDown') %>% pull(gene_id)
# 
# mat = data.frame(xspec  = ts.specQ3.genes %in% ts.specQ3.genes[names(ts.specQ3.genes)==i],
#                  rev =  ts.specQ3.genes%in%revi )
# fisher.test(table(mat))
# 
# ### background UD + DU genes tis.spec intersection
# i=c('Cortex','Lung','Liver','Muscle')[1]
# revi = filter(revg, tissue==i & direction%in%c('UpDown','DownUp')) %>% pull(gene_id)
# 
# bg = intersect(ts.specQ3.genes, revi )
# mat = data.frame(xspec  = bg %in% ts.specQ3.genes[names(ts.specQ3.genes)==i],
#                  rev =  bg%in%revi )
# table(mat)
# fisher.test(table(mat))
# 
# ### background UD + DU genes, test enrichment of tis. spec. genes among UD reversal genes
# i=c('Cortex','Lung','Liver','Muscle')[4]
# revi = filter(revg, tissue==i & direction%in%c('UpDown','DownUp')) %>% pull(gene_id)
# revud = filter(revg, tissue==i & direction%in%c('UpDown')) %>% pull(gene_id)
# bg = revi
# mat = data.frame(xspec  = bg %in% ts.specQ3.genes[names(ts.specQ3.genes)==i],
#                  rev =  bg%in%revud )
# fisher.test(table(mat)) # enriched
# mat = data.frame(xspec  = bg %in% ts.specQ3.genes[names(ts.specQ3.genes)==i],
#                  rev =  !bg%in%revud )
# fisher.test(table(mat)) # depleted

#################################
#################################
#################################

revg = readRDS('./data/processed/raw/revgenes.tissues.rds')

## old result, wrong:
# eachtissue_OR = list()
# eacttissue_table = list()
# for(i in names(revg)){
#   bg = intersect(ts.specQ3.genes, c(revg[[i]]$UpDown,revg[[i]]$DownUp) )
#   mat = data.frame(xspec  = ts.specQ3.genes%in%ts.specQ3.genes[names(ts.specQ3.genes)==i],
#                    rev =  ts.specQ3.genes%in%revg[[i]]$UpDown )
#   tst = fisher.test(table(mat))
#   eachtissue_OR[[i]]= c('OR' = tst$est, 'pval' = tst$p.val)
#   eacttissue_table[[i]] =  table(mat)
# }
# eachtissue_OR
# eacttissue_table
# saveRDS(list(OR= eachtissue_OR, table= eacttissue_table),
#         file='./data/processed/raw/tis_spec_rev_enrichment.rds')
# saveRDS(list(OR= eachtissue_OR, table= eacttissue_table),
#         file='./results/source_data/f3/tsspec_UD_enrichment.rds')
#table_s9 = lapply(eacttissue_table, function(x) data.frame(rbind(x)))
#write.xlsx(eacttissue_table, 'results/SI_tables/TableS9.xlsx')

# ORs = sapply(eachtissue_OR,function(x) x) %>% reshape2::melt() %>%
#   set_names(c('stat', 'tissue', 'value')) %>%
#   spread(stat, value) %>%
#   set_names(c('tissue', 'OR', 'pval'))
# 
# ts.specrevfisher = ORs %>%
#   mutate(LgOR = log2(OR)) %>%
#   ggplot(aes(x=tissue, y=LgOR)) +
#   geom_bar(stat='identity', fill='aquamarine4') +
#   ylab('Log2(OR)') + xlab('') +
#   geom_text(aes(x=tissue), label = '*', size=4)
# 
# ggsave('./results/figure3/OR_tsspec_rev.pdf', units='cm', height = 10, width = 8, useDingbats=F)
# ggsave('./results/figure3/OR_tsspec_rev.png', units='cm', height = 10, width = 8)
#####   

## Enrichment of UD (vs UU) and tissue specific genes:
##### choose background as up-up:
eachtissue_OR = list()
eacttissue_table = list()
for(i in names(revg)){
  bg = c(revg[[i]]$UpDown, revg[[i]]$UpUp) # background up-up genes
  spectis = ts.specQ3[[i]] # spec to tissue
  mat = data.frame(xspec = bg%in%spectis,
             rev = bg%in%revg[[i]]$UpDown)
  tst = fisher.test(table(mat))
  eachtissue_OR[[i]]= c('OR' = tst$est, 'pval' = tst$p.val)
  eacttissue_table[[i]] =  table(mat)
}
eachtissue_OR
eacttissue_table

saveRDS(list(OR= eachtissue_OR, table= eacttissue_table),
        file='./data/processed/raw/tis_spec_rev_enrichment.rds')
ORs = sapply(eachtissue_OR,function(x) x) %>% reshape2::melt() %>%
  set_names(c('stat', 'tissue', 'value')) %>%
  spread(stat, value) %>%
  set_names(c('tissue', 'OR', 'pval'))
# tissue       OR         pval
# 1 Cortex 1.649600 1.882387e-09
# 2   Lung 6.518892 1.371290e-60
# 3  Liver 0.871200 9.586655e-02
# 4 Muscle 1.255894 1.482466e-02

saveRDS(list(OR= ORs, table= eacttissue_table),'results/source_data/f3/tsspec_UD_enrichment.rds')

# ts.specrevfisher = ORs %>%
#   mutate(LgOR = log2(OR)) %>%
#   ggplot(aes(x=tissue, y=LgOR)) +
#   geom_bar(stat='identity', fill='aquamarine4') +
#   ylab('Log2(OR)') + xlab('') +
#   geom_text(aes(x=tissue), label = '*', size=4)

# ggsave('./results/figure3/OR_tsspec_rev.pdf', units='cm', height = 10, width = 8, useDingbats=F)
# ggsave('./results/figure3/OR_tsspec_rev.png', units='cm', height = 10, width = 8)

sapply(revg$Liver, length)
5157/(5157+3164)

sapply(sapply(revg$Liver, function(x) x[!x%in%ts.specQ3[['Liver']] ] ), length)
4765/(4765+2891)

## UD proportion of non-tissue specific genes:
sapply(names(revg), function(y){
  xx = sapply(revg[[y]], function(x) x[!x%in%ts.specQ3[[y]] ] )
  ud = length(xx$UpDown)
  uu = length(xx$UpUp)
  return( ud / sum(ud+uu) )
})

## UD proportion of tissue specific genes:
sapply(names(revg), function(y){
  xx = sapply(revg[[y]], function(x) x[x%in%ts.specQ3[[y]] ] )
  ud = length(xx$UpDown)
  uu = length(xx$UpUp)
  return( ud / sum(ud+uu) )
})

# eachtissue_OR2 = list()
# eacttissue_table2 = list()
# for(i in names(revg)){
#   #bg = c(revg[[i]]$UpDown,revg[[i]]$UpUp) # background up-up genes
#   spectis = ts.specQ3[[i]] # spec to tissue
#   revg = c(revg[[i]]$UpDown,revg[[i]]$UpUp)
#   mat = data.frame(xspec = spectis%in%revg[[i]]$UpDown,
#                    rev = spectis%in%revg[[i]]$UpDown)
#   tst = fisher.test(table(mat))
#   eachtissue_OR2[[i]]= c('OR' = tst$est, 'pval' = tst$p.val)
#   eacttissue_table2[[i]] =  table(mat)
# }
# eachtissue_OR2
# eacttissue_table2
