library(data.table)
library(tidyverse)

att = read_tsv('./data/other_datasets/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', 
               guess_max = 25000)
phe = read_tsv('./data/other_datasets/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', 
               guess_max = 1000)
att$SUBJID = sapply(strsplit(att$SAMPID,'-'),function(x)paste(x[1],x[2],sep='-'))
phe$SEX = c('male','female')[phe$SEX]
phe = phe %>%
  set_names(c('id','sex','age','death')) %>%
  mutate(sex = as.factor(sex),
         age = as.factor(age))
allattr = att %>%
  select(SAMPID,SMTS,SMTSD,SMNABTCH,SMGEBTCH,SUBJID,SMRIN,SMTSISCH) %>%
  set_names(c('sample_id','major_tissue','minor_tissue','batch1','batch2','id','rin','ischemic_time')) %>%
  unique() %>%
  full_join(phe)

#Filter out the indv with hardy death score of 1 or 2 
allattr = allattr %>%
  filter(death %in% c(1,2))

# file is a symbolic link, read from the original directory:
# dat=fread('../../../../pseudodropbox/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct',
#           select = c('Name',allattr$sample_id))
dat=fread('./data/other_datasets/GTEx/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct',
          select = c('Name',allattr$sample_id))
samplesx=colnames(dat)[-1]

dat = as.data.frame(dat)
rownames(dat) = dat$Name
dat$Name = NULL
dat = as.matrix(dat)
samplesx = intersect(allattr$sample_id, colnames(dat))
samples_by_tissues = tapply(as.character(allattr$sample_id),
                            INDEX=allattr$minor_tissue, FUN = function(x)unique(c(x)))
dat = lapply(samples_by_tissues,function(samps){
  samps = intersect(samps,samplesx)
  dat[,samps]
})

system('mkdir -p data/processed/GTEx/expression/tpm/')
names(dat) = sapply(strsplit(gsub(' ','',names(dat)),'[(]'),function(x)x[[1]])

sapply(names(dat),function(nm){
  saveRDS(dat[[nm]],paste('data/processed/GTEx/expression/tpm/',nm,'.rds',sep=''))
})
saveRDS(allattr,'./data/processed/GTEx/attr.rds')

rm(list=ls())
