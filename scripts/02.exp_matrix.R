getwd()
load("./data/preprocess/sample_features.rdata")

nm = c(ctx_name, lng_name, lv_name, ms_name)
flist = paste0("./raw_data/",nm,".txt")
nms = gsub("-","_",c(ctx_samples, lng_samples, lv_samples, ms_samples))
age = c(ctx_age, lng_age, lv_age, ms_age)
samp_id = paste('s',1:length(age),sep='')
samp_tissue = rep(c("cortex","lung","liver","muscle"),times=c(15,16,16,16))
names(age) = samp_id
names(ind_id) = samp_id
names(samp_tissue) = samp_id

# prep for deseq2:
data.frame(fname = nms, ind_id = names(nms),
           tissue = samp_tissue, 
           age = age, sample_id = samp_id) %>%
  saveRDS(. , file = './data/processed/raw/fname_samples.rds')

# load each FPKM data and assign it to lane number; MA1_3, MA1_5 etc:
for(i in 1:length(flist)) assign(nms[i],read.table(flist[i],sep="\t",header = T))

# check if there is any duplicated genes:
y=sapply(1:length(nms),function(i){
  dups = get(nms[i])[duplicated(get(nms[i])[,1]),1]
  temp = unique(as.character(get(nms[i])[get(nms[i])[,1]%in%dups,1] ))
  return(temp)
})

dim(y) # 50 genes are duplicated for all individuals
for(i in 1:(ncol(y)-1)) print(length( intersect(y[,i],y[,i+1]) ))
# all of them are the same genes

# take the sum of the duplicated genes :
duplist=list() # save duplicated gene expressions
for(i in 1:length(nms)){
  k = get(nms[i])[,2]
  names(k) = get(nms[i])[,1]
  dups = unique(names(k)[duplicated(names(k))]) # duplicates, 66, unique 50
  duplist[[ nms[i] ]] = sapply(dups,function(x){ k[names(k)%in%x] })
  sumdup = sapply(dups,function(x){ sum(k[names(k)%in%x]) })
  k=k[!names(k)%in%dups]
  k=c(k,sumdup)
  assign(nms[i], k)
}

# take common genes: 
genes = Reduce(intersect, lapply(1:63, function(i) names(get(nms[i])) ) )
length(genes)

mat = matrix(NA,nrow = length(genes), ncol = length(nms), dimnames = list(genes,nms))

# assign sample FPKM  values to corresponding cell of the expression matrix.
for(i in 1:length(nms)){
  mt = rownames(mat)[rownames(mat)%in%names(get(nms[i])) ]
  mat[rownames(mat)%in%names(get(nms[i])), i] = get(nms[i])[mt]
}

# remove individual sample matrices:
rm(list=nms)

# save raw matrix:
colnames(mat) = samp_id
rawmat = mat
colnames(rawmat) = samp_id
save(rawmat, age, samp_tissue, ind_id,  file="./data/preprocess/FPKM_rawmatrix.rdata")

# take only protein coding genes:
allgenetype = readRDS("./data/genetype.rds")
pgenes = allgenetype[allgenetype[,2]=="protein_coding",1]
length(pgenes)  
length(unique(pgenes))
shrgenes = rownames(mat)[(rownames(mat)%in%pgenes)] 
mat = mat[shrgenes,]

saveRDS(mat, file="./data/preprocess/FPKM_pcodinggenes_rawmatrix.rds")
saveRDS(age, file="./data/preprocess/ages.rds")
saveRDS(ind_id, file="./data/preprocess/individual_ids.rds")
samp_tissue = str_to_title(samp_tissue)
names(samp_tissue) = samp_id
saveRDS(samp_tissue, file="./data/preprocess/tissue_ids.rds")

# Filtration:
# remove the genes which are not detected in at least 25% of the samples (15).
mat2 = mat[!apply(mat,1,function(x){ sum(x==0)>15 } ), ] # 6684 genes removed
mat3 = log2(mat2+1)

saveRDS(mat3, file='./data/preprocess/log2_pcodinggenes_rawmatrix.rds')

library(preprocessCore) 
expr = normalize.quantiles(mat3)
dimnames(expr) = dimnames(mat3)

saveRDS(expr, file= './data/processed/raw/expression.rds')
saveRDS(age, file="./data/processed/raw/ages.rds")
saveRDS(ind_id, file="./data/processed/raw/individual_ids.rds")
saveRDS(samp_tissue, file="./data/processed/raw/tissue_ids.rds")

library(tidyverse)
sample_info = data.frame(ind_id = ind_id, age = unname(age), tissue = samp_tissue) %>%
  mutate(log2age = log2(age)) %>%
  mutate(sample_id = samp_id) 

saveRDS(sample_info, './data/processed/tidy/sample_info.rds')

expression = reshape2::melt(expr) %>%
  set_names(c('gene_id','sample_id','expression')) 

saveRDS(expression,'./data/processed/tidy/expression.rds')
##########
##########
##########