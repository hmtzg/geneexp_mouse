library(tidyverse)

# metadat = read.csv('data/other_datasets/schaum/GSE132040_MACA_Bulk_metadata.csv.gz')
# head(metadat)
# colnames(metadat)
# unique(metadat$Sample.name)
# unique(metadat$title) #  all bulk RNA seq
# unique(metadat$source.name) #  tissue data
# unique(metadat$organism) # C57/BL6
# unique(metadat$characteristics..age) #  age
# unique(metadat$characteristics..developmental.stage) # months postnatal
# unique(metadat$characteristics..sex) # m, f, missing
# unique(metadat$molecule) # total RNA
# unique(metadat$description) #  NA
# unique(metadat$processed.data.file) #
# unique(metadat$raw.file) #  raw srr file name
# unique(metadat$BioSample) #  sample id
# unique(metadat$Instrument.Model) # platform
# 
# ##################
# ids = sapply(strsplit(metadat$Sample.name,'Bulk_'),`[[`,1)
# ids = sapply(strsplit(ids,'_384'),`[[`,1)
# plateid = sapply(strsplit(metadat$Sample.name,'Bulk_'),`[[`,2)
# plate = sapply(strsplit(plateid,'_'),`[[`, 1)
# sid = sapply(strsplit(plateid,'_'),`[[`, 2)
# 
# tissues = strsplit(metadat$source.name,'_')
# table(sapply(tissues, length))
# tissues[sapply(tissues, length)==1] # all NA
# tissues[sapply(tissues, length)==2] # one word name tissues
# tissues[sapply(tissues, length)==3] # limb muscle & small intestine
# 
# ## tissue names:
# tis = sapply(strsplit(metadat$source.name,'_'),`[[`,1 )
# # correct the names of two word tissues:
# tis[sapply(tissues, length)==3] = sapply(tissues[sapply(tissues, length)==3], function(x) paste(x[1],x[2]) )
# 
# # tissue ids:
# tis.id = sapply(tissues,function(x){
#   if(length(x)==1) x
#   else if(length(x)==2) x[2]
#   else x[3]
# })
# 
# # age:
# age = metadat$characteristics..age
# table(age)
# 
# # sex:
# sex = metadat$characteristics..sex
# table(sex)
#####
# ## expression dat:
# expdat = read.csv(paste0('data/other_datasets/schaum/',
#                   'GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv.gz'))
# length(unique(expdat$gene)) == nrow(expdat) #  unique gene ids
# rownames(expdat) = expdat$gene
# expdat = expdat[,-1]
# 
# ## check if the order of expr colnames match to  metadat ids:
# head(colnames(expdat))
# # get sample names from expression matrix columns:
# expcols = sapply(strsplit(colnames(expdat),'.gencode'),`[[`,1)
# head(expcols)
# head(metadat$Sample.name)
# identical(expcols, metadat$Sample.name) # not in the same order
# sample_order = match(metadat$Sample.name,expcols)
# identical(expcols[sample_order], metadat$Sample.name)
# # convert expression matrix to same order with metadat:
# expdat = expdat[,sample_order]
# expcols = expcols[sample_order]
# 
# ## remove NA individuals:
# NAposition = grep('NA', tis)
# 
# expdat = expdat[,-NAposition]
# expcols = expcols[-NAposition]
# sex = sex[-NAposition]
# age = as.numeric(age[-NAposition])
# tis.id = tis.id[-NAposition]
# tis = tis[-NAposition]
# tissues = tissues[-NAposition]
# sid = sid[-NAposition]
# plate = plate[-NAposition]
# ids = ids[-NAposition]
# metadat = metadat[-NAposition,]
# 
# ###
# tissues[age==12 & sex=='f']
# 
# sid[tis=='Lung' & age==12 & sex=='m']
# ids[tis=='Lung' & age==12 & sex=='m']
# metadat[tis=='Lung' & age==12 & sex=='m',]
# metadat[age==12 & sex=='f',]
# 
# tis[age==12 & sex=='m']
# table(tis[age==12 & sex=='m'])
# table(tis.id[age==12 & sex=='m'])
# sid[age==12 & sex=='m']
# ids[age==12 & sex=='m']
# 
# 
# expcols
# 
# ##############
# colnames(expdat) = expcols
# saveRDS(expdat, 'data/other_datasets/tabulabulk/rawexp.rds')
# 
# sinfo = data.frame(id = expcols, age = age, sex = sex, tissue = tis, plateid = plate, letterid = ids, sid = sid,
#            tissueid = tis.id) %>% tibble()
# saveRDS(sinfo,  'data/other_datasets/tabulabulk/sinfo.rds')


###########
########### New metadat
###########
###########
#####
mdat = read.csv('data/other_datasets/schaum/MACA_Bulk_metadata.csv')

# metadat_sampname = sapply(strsplit(metadat$Sample.name, '_S'),`[[`, 1)
# length(intersect(metadat_sampname, mdat$Sample_Name))
# setdiff(metadat_sampname, mdat$Sample_Name)
# setdiff(mdat$Sample_Name,metadat_sampname)

## expression dat:
expdat = read.csv(paste0('data/other_datasets/schaum/',
                         'GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv.gz'))
length(unique(expdat$gene)) == nrow(expdat) #  unique gene ids
rownames(expdat) = expdat$gene
expdat = expdat[,-1]

## check if the order of expr colnames match to  metadat ids:
head(colnames(expdat))
# get sample names from expression matrix columns:
table(sapply(strsplit(colnames(expdat),'_S'),length))
expcols = sapply(strsplit(colnames(expdat),'_S'),`[[`,1)
colnames(expdat) = expcols
rm(expcols)
# get same samples:
ss = intersect(colnames(expdat), mdat$Sample_Name)
expdat = expdat[,ss]
# also for metadat:
mdat = mdat %>% filter(Sample_Name%in%ss)

setdiff(colnames(expdat), mdat$Sample_Name)
identical(colnames(expdat), mdat$Sample_Name) # not in the same order
sample_order = match(mdat$Sample_Name, colnames(expdat))
identical(colnames(expdat)[sample_order], mdat$Sample_Name)
# convert expression matrix to same order with metadat:
expdat = expdat[,sample_order]
#expcols = expcols[sample_order]

##

##############
identical(colnames(expdat), mdat$Sample_Name)
table(colSums(expdat)>4e6)

# remove 17 samples
#filt_samples = names(which(colSums(expdat)>4e6))

#expdat = expdat[,colSums(expdat)>4e6]

saveRDS(expdat, 'data/other_datasets/schaum/rawexp.rds')

#mdat = mdat[mdat$Sample_Name%in%filt_samples,]
head(mdat)
tis = mdat$Organ
sex = mdat$Sex
mouseID =  mdat$Mouse_ID
mouseNum = mdat$Mouse_Num
age = mdat$Age
samp_name = mdat$Sample_Name

sinfo = mdat %>%
  select(Sample_Name, Number, Mouse_ID, Sex, Organ, Experiment_ID, Age, Mouse_Num) %>% tibble()

saveRDS(sinfo,  'data/other_datasets/schaum/sinfo.rds')

