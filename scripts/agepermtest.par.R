timestart = proc.time()
print(timestart)
arg = commandArgs(trailingOnly = T)
k = as.numeric(arg[1])
#k=10
cl = as.numeric(arg[2])
#cl=4
load("./data/processed/raw/ageperm1000.RData")
exp = readRDS('./data/processed/raw/expression.rds')
ind.age = readRDS('./data/preprocess/ages.rds')
ind.id  = readRDS('./data/preprocess/individual_ids.rds')
samp_tissue = readRDS('./data/preprocess/tissue_ids.rds')
names(ind.age) = ind.id
#load("Dropbox/projects/ageing/proc_data/preproc/exp.rdata")
#ind.age = readRDS("/home/hmt/Dropbox/projects/ageing/proc_data/preproc/ind.age.rds")
colnames(exp) = names(ind.age)
###
library(foreach)
library(doParallel)
registerDoParallel(cl)

print(paste("running permutation test for development with", cl, "clusters and", k, "iterations", sep=" "))
randcordev = list()
for(i in 1:4){
  tissue = samp_tissue == unique(samp_tissue)[i]
  print(paste("tissue:", unique(samp_tissue)[i]))
  exp.ts = exp[,tissue]
  perm.ts = foreach(perm = 1:k, .combine = cbind, .export = c("perm","exp.ts","permageD"), .errorhandling = "stop") %dopar%{
    if(perm%%5==0)print(paste("perm:",perm))
    apply(exp.ts,1,function(x){
      cor.test(permageD[[perm]][[i]], x[names(permageD[[perm]][[i]])],m="s")$est
    })
  }
  randcordev[[ unique(samp_tissue)[i] ]] = perm.ts
}
lapply(randcordev,function(x)x[1:5,1:5])

###
randcoraging = list()
for(i in 1:4){
  tissue = samp_tissue == unique(samp_tissue)[i]
  print(paste("tissue:", unique(samp_tissue)[i]))
  exp.ts = exp[,tissue]
  perm.ts = foreach(perm = 1:k, .combine = cbind, .export = c("perm","exp.ts","permageA"), .errorhandling = "stop") %dopar%{
    if(perm%%5==0)print(paste("perm:",perm))
    apply(exp.ts,1,function(x){
      cor.test(permageA[[perm]][[i]], x[names(permageA[[perm]][[i]])],m="s")$est
    })
  }
  randcoraging[[ unique(samp_tissue)[i] ]] = perm.ts
}

lapply(randcoraging,function(x)x[1:5,1:5])

stopImplicitCluster()

saveRDS(randcordev, file="./data/processed/raw/dev.perm.rds")
saveRDS(randcoraging, file="./data/processed/raw/aging.perm.rds")

elapsed = proc.time()-timestart
print(elapsed)
print("bitti")
