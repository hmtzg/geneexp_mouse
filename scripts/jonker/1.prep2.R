getwd()
library(oligo)
library(GEOquery)
library(tidyverse)

s.matrix = getGEO("GSE34378",getGPL = F)
features = pData(s.matrix[[1]])
colnames(features)
as.character(features[,"title"])
gsmall = as.character(features[,"geo_accession"])
tissueall = features[,"tissue:ch1"]
strainall = features[,"strain:ch1"]
ageall = as.numeric(gsub(" wks","",features[,"age:ch1"]))
genderall = "no gender info"
sample_id = names(features[,'characteristics_ch1.2'])
sample_id= paste0('s',1:90)

unique(features[,'tissue:ch1'])
ind_id = sapply(strsplit(as.character(features[,1]), split = '\ '), function(x) 
  paste((x[c(2,5)]),collapse = '-') )

###
path = "/home/hmt/pseudodropbox/raw_exp_data/Jonker2013_GSE34378/"
files = list.celfiles(path,listGzipped = T,full.names=F)
files = files[grep("GSM",files)]
flist = list.celfiles(path,listGzipped = T,full.names=T)

# required package pd.mouse430.2 for read.cellfiles function :
oligo = read.celfiles(filenames=flist)
rma_expset = rma(oligo)
exp = exprs(rma_expset)

identical(substr(colnames(exp),1,9), gsmall)
colnames(exp) = substr(colnames(exp),1,9)
exp = 2^exp
#BiocManager::install("mouse4302.db")
library("mouse4302.db")
x = mouse4302ENSEMBL
mapped_genes = mappedkeys(x)
xx = as.list(x[mapped_genes])
######## exclude probes annotated more than one gene
exclude=c()
for(i in 1:length(xx)){
  if(length(xx[[i]])>1){
    exclude=c(exclude,names(xx[i]))
  }
}
exp=exp[!(rownames(exp)%in%exclude),]
ids=unlist(xx[!names(xx)%in%exclude])
#### get only annotated probes
exp=exp[rownames(exp)%in%names(ids),] 
rownames(exp) = ids[rownames(exp)]
#if multiple probesets annotated to same gene, take average of probesets
exp2 = t(sapply(unique(rownames(exp)),function(x){
  if(sum(rownames(exp)%in%x)>1){
    apply(exp[rownames(exp)%in%x,],2,mean)
  }
  else{
    exp[x,]
  }
}))

############## normalise samples across tissues
library(preprocessCore)
expn=normalize.quantiles(exp2)
dimnames(expn)=dimnames(exp2)

ind_id = rep(paste0('ind',1:length(unique(ind_id))),5)
colnames(expn) = sample_id
names(ageall) = sample_id
names(ind_id) = sample_id
names(tissueall) = sample_id
saveRDS(expn, './data/other_datasets/jonker/raw/exp_qn2.rds')
saveRDS(ageall, './data/other_datasets/jonker/raw/age2.rds')
saveRDS(sample_id, './data/other_datasets/jonker/raw/sample_id2.rds')
saveRDS(tissueall, './data/other_datasets/jonker/raw/tissue_id2.rds')
saveRDS(ind_id, './data/other_datasets/jonker/raw/ind_id2.rds')

sample_info = data.frame(age=ageall, tissue= tissueall, ind_id = ind_id, sample_id= sample_id,
                         log2age= log2(ageall))

expr = reshape2::melt(expn) %>%
  set_names(c('gene_id', 'sample_id', 'expression'))

saveRDS(expr,'./data/other_datasets/jonker/raw/expression2.rds')
saveRDS(sample_info,'./data/other_datasets/jonker/raw/sample_info2.rds')

save(list=ls(), file='./data/other_datasets/jonker/prep_analysis.rdata')
