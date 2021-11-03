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
exp = exp[, !colnames(exp) == "465"]
samp_tissue = samp_tissue[!names(age)=="465"]
age = age[!names(age)=="465"]

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

##### in which tissue the highest expression change occurs for each gene:
ts.expr.ch = colnames(expch_ag)[apply(expch_ag, 1, function(x) which.max(abs(x)))]
names(ts.expr.ch) = rownames(expch_ag)

# direction of expression change for those genes :
ts.expr.ch.dir = sapply(1:length(ts.expr.ch), function(x){ sign(expch_ag[x, ts.expr.ch[x]]) } )
names(ts.expr.ch.dir) = names(ts.expr.ch)

### for all genes:
# ts.spec: specificity of genes; genes having max ES assigned to tissue
# ts.specQ3: tissue specific genes; genes having ES >Q3 assigned to a tissue
# ts.expr.ch: maximum abs. expression change occur in which tissue

mat = data.frame(sameness = ts.spec == ts.expr.ch, exp_dir = ts.expr.ch.dir)
table(mat) # FALSE: highest expression change occur in not specific tissue
# TRUE: highest expression change occur in specific tissue
# -1 : expression change direction is negative (decrease with age in ageing period)
# 1 : expression change direction is positive (increase with age in ageing period)
# -1_FALSE: highest expression decrease occur in not specific tissue (not inrestested in this: )
# -1_TRUE: highest expression decrease occur in specific tissue (interesting: identity lose?)
# 1_FALSE: highest expression increase occur in not specific tissue (interesting: tissue resembl each other?)
# 1_TRUE: highest expression increase occur in specific tissue (not interested in this:)

fisher.test(table(mat)[,c(2,1)])
# OR: 2.564, pval < 1e10-16


########
########
######## actual result:
########
# for tissue-specific genes (>Q3) :
ts.specQ3.genes = unlist(ts.specQ3)
names(ts.specQ3.genes) = gsub('[0-9]','', names(ts.specQ3.genes))

mat = data.frame(sameness = names(ts.specQ3.genes) == ts.expr.ch[ts.specQ3.genes],
                 exp_dir = ts.expr.ch.dir[ts.specQ3.genes] )
table(mat)[,c(2,1)]
fisher.test(table(mat)[,c(2,1)])
# OR: 4.298, pval < 1e10-16

# exprch_tisspec = table(mat)[c(2,1),c(2,1)]
# colnames(exprch_tisspec) = c('Increase w/ age (+)', 'Decrease w/ age (-)' )
# rownames(exprch_tisspec) = c('Same tissue', 'Other tissue')
# exprch_tisspec
# fisher.test(exprch_tisspec[,c(2,1)])

#### among div-conv genes:
# mat.dc = mat[rownames(mat)%in%dcgenes,] # 1287 genes in total
# table(mat.dc)[,c(2,1)]
# fisher.test(table(mat.dc)[,c(2,1)])
# OR: 38.52, pval < 1e10-16


########
### tis.spec genes enrichment among div-conv genes:
spec.dc.mat = data.frame(spec = rownames(exp)%in%ts.specQ3.genes ,
                         divconv = rownames(exp)%in%dcgenes )
# spec 1: not ts specific, 2: ts specific
# divconv 1: not dc gene, 2: dc gene

fisher.test(table(spec.dc.mat)) # 
# OR: 1.149, pval = 0.00051
saveRDS(spec.dc.mat, file='./data/processed/raw/ts_spec_dc_enrichment.rds')
spec.dc.mat.table = table(spec.dc.mat)
spec.dc.mat.fisher = fisher.test(table(spec.dc.mat))
saveRDS(list(table = spec.dc.mat, test=spec.dc.mat.table),
        file='./results/source_data/f3/tsspec_dc_enrichment.rds')

dc_vs_tisspec = table(spec.dc.mat)[c(2,1),c(2,1)]
colnames(dc_vs_tisspec) = c('DC', ' Non-DC')
rownames(dc_vs_tisspec) = c('Specific to a tissue', 'Not tissue-specific')
dc_vs_tisspec
#write.xlsx(dc_vs_tisspec, file='results/SI_tables/TableS10.xlsx')

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

eachtissue_OR = list()
eacttissue_table = list()
for(i in names(revg)){
  bg = intersect(ts.specQ3.genes, c(revg[[i]]$UpDown,revg[[i]]$DownUp) )
  mat = data.frame(xspec  = ts.specQ3.genes%in%ts.specQ3.genes[names(ts.specQ3.genes)==i],
                   rev =  ts.specQ3.genes%in%revg[[i]]$UpDown )
  tst = fisher.test(table(mat))
  eachtissue_OR[[i]]= c('OR' = tst$est, 'pval' = tst$p.val)
  eacttissue_table[[i]] =  table(mat)
}
eachtissue_OR
eacttissue_table
saveRDS(list(OR= eachtissue_OR, table= eacttissue_table),
        file='./data/processed/raw/tis_spec_rev_enrichment.rds')
saveRDS(list(OR= eachtissue_OR, table= eacttissue_table),
        file='./results/source_data/f3/tsspec_UD_enrichment.rds')
#table_s9 = lapply(eacttissue_table, function(x) data.frame(rbind(x)))
#write.xlsx(eacttissue_table, 'results/SI_tables/TableS9.xlsx')


ORs = sapply(eachtissue_OR,function(x) x) %>% reshape2::melt() %>%
  set_names(c('stat', 'tissue', 'value')) %>%
  spread(stat, value) %>%
  set_names(c('tissue', 'OR', 'pval'))

ts.specrevfisher = ORs %>%
  mutate(LgOR = log2(OR)) %>%
  ggplot(aes(x=tissue, y=LgOR)) +
  geom_bar(stat='identity', fill='aquamarine4') +
  ylab('Log2(OR)') + xlab('') +
  geom_text(aes(x=tissue), label = '*', size=4)

ggsave('./results/figure3/OR_tsspec_rev.pdf', units='cm', height = 10, width = 8, useDingbats=F)
ggsave('./results/figure3/OR_tsspec_rev.png', units='cm', height = 10, width = 8)
#####   

# bg = ts.specQ3.genes
# cort_tst = data.frame(cort_sp = bg%in%ts.specQ3$Cortex,
#            cort_rev = bg%in%c(revg$Cortex$UpDown))
# table(cort_tst)
# fisher.test(table(cort_tst))

