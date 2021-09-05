
library(tidyverse)
library(ggpubr)
pntnorm <- (1/0.352777778)
theme_set(theme_pubr(base_size = 6, legend = 'top') +
            theme(legend.key.size = unit(2,'pt')))
attr = readRDS('./data/processed/GTEx/attr.rds') 

# filter(!exc) #exclude individuals with more than one sample in a minor tissue

expvals = list(Cortex = readRDS('./data/processed/GTEx/expression/tpm/Brain-Cortex.rds'),
               Lung = readRDS('./data/processed/GTEx/expression/tpm/Lung.rds'),
               Liver = readRDS('./data/processed/GTEx/expression/tpm/Liver.rds'),
               Muscle = readRDS('./data/processed/GTEx/expression/tpm/Muscle-Skeletal.rds'))


samplelist = unique(unlist(lapply(expvals, function(x)colnames(x))))

# there are samples from same individual but not having expression:
attr = attr %>%
  filter(minor_tissue %in% c('Lung','Liver','Brain - Cortex','Muscle - Skeletal')) %>%
  filter(sample_id %in% samplelist) %>%
  group_by(id,minor_tissue) %>%
  summarise(n = length(unique(sample_id))) %>%
  ungroup() %>%
  mutate(exc = n>1) %>%
  right_join(attr) 
table(attr$exc) # no duplication
summary(attr$exc)
genelist = unique(unlist(lapply(expvals,function(expx){
  sapply(strsplit(as.character(rownames(expx)),'[.]'),function(x)x[1])
})))
#martx = biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
martx = biomaRt::useEnsembl('ensembl', 'hsapiens_gene_ensembl', mirror='asia')
biotype = biomaRt::getBM(attributes = c('ensembl_gene_id', 'gene_biotype'), filters = 'ensembl_gene_id', 
                         values = genelist, mart = martx)
codinggenes = unique(filter(biotype,gene_biotype == 'protein_coding')$ensembl_gene_id)

incids = sapply(expvals, colnames) %>% 
                            reshape2::melt() %>%
                            set_names(c('sample_id','Tissue')) %>%
                            left_join(attr)
                          
# idmap = setNames(attr$id, attr$sample_id)[samplelist]
# inc_sampleIDs = names(idmap[idmap%in%incids]) # get the sample ids for the individuals to be included

# subset attr data for only samples common in four tissues:
# attr = attr %>% 
#   filter(sample_id %in% inc_sampleIDs ) 

expvals = lapply(expvals, function(expx){
  rownames(expx) = sapply(strsplit(as.character(rownames(expx)),'[.]'),function(x)x[1]) # remove gene version
  expx = expx[rownames(expx) %in% codinggenes,] #get only the coding genes
  dupgenes = unique(rownames(expx)[duplicated(rownames(expx))]) # to remove duplicated genes
  #expx = expx[!(rownames(expx)%in%dupgenes), colnames(expx) %in% inc_sampleIDs] 
  # exc = apply(expx,1,function(x) (max(x)<1)|(mean(x==0)==1))
  # expx = expx[!exc,] 
  expx
})

sample2id = setNames(attr$id, attr$sample_id)
expx = expvals[[4]]
sum(duplicated(unname(sample2id[colnames(expx)])))

genelist = Reduce('intersect',lapply(expvals,rownames))
expvals = cbind(expvals[[1]][genelist,],expvals[[2]][genelist,],
                expvals[[3]][genelist,],expvals[[4]][genelist,]) #to make sure the gene order is the same
expvals = expvals [ ! (rowMeans(expvals==0)>0.25), ] # exclude the genes if expression value is 0 in more than
# 25% of the samples
dim(expvals)
#16211 - 188 (47 samples)  -> 16197

exp_l2 = log2(expvals+1)
exp_l2_qn = preprocessCore::normalize.quantiles(exp_l2)
dimnames(exp_l2_qn) = dimnames(exp_l2)

expvals = reshape2::melt(exp_l2_qn) %>%
  set_names(c('GeneID','sample_id','Expression')) %>%
  left_join(incids) %>%
  mutate(Tissue = ifelse(major_tissue=='Brain','Cortex',major_tissue)) %>%
  mutate(Age = c(25,35,45,55,65,75)[age]) 

#Expression - Age correlation
agecor = expvals %>%
  group_by(Tissue, GeneID, sex) %>%
  summarise(rho = cor.test(Expression, Age, method = 's')$est,
            p = cor.test(Expression, Age, method = 's')$p.val)

agecor = ungroup(agecor) %>%
  mutate(adjusted_p = p.adjust(p, method = 'fdr'))
agecor

###################
incids %>% group_by(Tissue, sex) %>% summarise(n=length(sample_id))
# 1 Cortex female    32
# 2 Cortex male     125
# 3 Liver  female    22
# 4 Liver  male      84
# 5 Lung   female    42
# 6 Lung   male     140
# 7 Muscle female    51
# 8 Muscle male     183

# Expression Change correlations between female and male:
agecor %>%
  select(-p, -adjusted_p) %>% 
  spread(key= sex, value = rho) %>%
  group_by(Tissue) %>%
  summarise(rho = cor.test(female, male, method='s')$est,
            p = cor.test(female, male, method='s')$p.val)
# 1 Cortex 0.605 0        
# 2 Liver  0.119 1.20e- 51
# 3 Lung   0.218 2.79e-173
# 4 Muscle 0.517 0    

sexcorplot = agecor %>%
  select(-p, -adjusted_p) %>% 
  spread(key= sex, value = rho) %>%
  ggplot(aes(x=female, y=male)) +
  facet_grid(~Tissue) +
  geom_point(size=0.4, alpha=0.3, color='gray30') +
  geom_smooth(method = 'lm', se=F, color='darkred') +
  stat_cor(method = 'spearman', size = 6/pntnorm, cor.coef.name = 'rho') 

sexcorplot
ggsave('results/GTEx/sex_agecor_allsamples.pdf', sexcorplot, units='cm', width = 12, height = 8,useDingbats=F)
ggsave('results/GTEx/sex_agecor_allsamples.png', sexcorplot, units='cm', width = 12, height = 8)

save(list=ls(), file = 'results/GTEx/sexdata.rdata')
