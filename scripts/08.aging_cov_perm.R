#source("Dropbox/projects/Rfunctions.R")
load("./data/processed/raw/ageperm1000.RData")
#exp = readRDS('Dropbox/projects/repos/geneexp_mouse/data/processed/raw/expression.rds')
#load("/home/hmt/Dropbox/projects/ageing/proc_data/preproc/exp.rdata")
#ind.age = readRDS("Dropbox/projects/ageing/proc_data/preproc/ind.age.rds")
exp = readRDS('./data/processed/raw/expression.rds')
ind.age = readRDS('./data/preprocess/ages.rds')
ind.id  = readRDS('./data/preprocess/individual_ids.rds')
names(ind.age)= ind.id
samp_tissue = readRDS('./data/preprocess/tissue_ids.rds')
colnames(exp) = ind.id

## remove individual who lacks expression in cortex
exp = exp[, !colnames(exp) =='465']
names(samp_tissue) = names(ind.age)
samp_tissue = samp_tissue[!names(ind.age)=="465"]
#age = age[!names(ind.age)=="465"]

ageu = ind.age[unique(names(ind.age))]
ageu = ageu[!names(ageu)%in%'465']

genecov = sapply(unique(colnames(exp)),function(x){
  sapply(rownames(exp),function(y){
    sd(exp[y, colnames(exp) == x]) / mean(exp[y, colnames(exp) == x])
  })
})
covA = genecov[,ageu>90]

#####
names(permageA) = paste0('p',1:length(permageA))
count = 0
randcorA = sapply(permageA, function(perm) {
  count <<- count + 1;
  print(count)
  covchperm =  sapply(rownames(covA), function(x){
    cor(perm[[1]], covA[x, ], m='s')
  })
})
colnames(randcorA) = paste0('p',1:length(permageA))

saveRDS(randcorA, file="./data/processed/raw/aging.cov.perm.rds")
