library(GEOquery)
#library(oligo)
s.matrix = getGEO("GSE34378")
features = pData(s.matrix[[1]])
as.character(features[,"title"])
gsmall = as.character(features[,"geo_accession"])
tissueall = features[,"tissue:ch1"]
# strainall = features[,"strain:ch1"]getwd()
# source("http://bioconductor.org/biocLite.R")
# biocLite("GEOquery") # installed: 26.03.2019
ageall = as.numeric(gsub(" wks","",features[,"age:ch1"]))
#genderall = "no gender info"
sample_id = names(features[,'characteristics_ch1.2'])
sample_id= paste0('s',1:90)

unique(features[,'tissue:ch1'])
ind_id = sapply(strsplit(as.character(features[,1]), split = '\ '), function(x) paste((x[c(2,5)]),collapse = '-') )

kidney = readRDS('./data/other_datasets/jonker/raw/kidney.rds')
liver = readRDS('./data/other_datasets/jonker/raw/liver.rds')
spleen = readRDS('./data/other_datasets/jonker/raw/spleen.rds')
lung = readRDS('./data/other_datasets/jonker/raw/lung.rds')
brain = readRDS('./data/other_datasets/jonker/raw/brain.rds')
identical(rownames(brain),rownames(spleen))

exp = cbind(kidney, liver,spleen, lung,brain)
identical(substr(colnames(exp),1,9), gsmall)
colnames(exp) = substr(colnames(exp),1,9)
exp = 2^exp
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

##################### get protein coding genes
allgenetype = readRDS("./data/genetype.rds")
pgenes = allgenetype[allgenetype[,2]=="protein_coding",1]
shrgenes = rownames(exp2)[(rownames(exp2)%in%pgenes)] 
exp2 = exp2[shrgenes,]

# remove the genes which are not detected in at least 25% of the samples (23).
exp2 = exp[!apply(exp2,1,function(x){ sum(x==0)>23 } ), ] # 6684 genes removed
exp2 = log2(exp2+1)


############## normalise samples across tissues
library(preprocessCore)
expn=normalize.quantiles(exp2)
dimnames(expn)=dimnames(exp2)

ind_id = rep(paste0('ind',1:length(unique(ind_id))),5)
colnames(expn) = sample_id
names(ageall) = sample_id
names(ind_id) = sample_id
names(tissueall) = sample_id
saveRDS(expn, './data/other_datasets/jonker/raw/expression.rds')
saveRDS(ageall, './data/other_datasets/jonker/raw/age.rds')
saveRDS(sample_id, './data/other_datasets/jonker/raw/sample_id.rds')
saveRDS(tissueall, './data/other_datasets/jonker/raw/tissue_id.rds')
saveRDS(ind_id, './data/other_datasets/jonker/raw/ind_id.rds')

library(tidyverse)
sample_info = data.frame(age=ageall, tissue= tissueall, ind_id = ind_id, sample_id= sample_id, log2age= log2(ageall))
expr = reshape2::melt(expn) %>%
  set_names('gene_id','sample_id','expression')

saveRDS(expr,'./data/other_datasets/jonker/processed/expression.rds')
saveRDS(sample_info,'./data/other_datasets/jonker/processed/sample_info.rds')

