getwd()
library(openxlsx)
id = read.xlsx("./data/preprocess/mm_sample_id.xlsx",sheet = 1)
head(id)
id[,1] = as.character(id[,1])
region = gsub("[0-9]","",id[,1])
mm_id = gsub("[A-Z]","",id[,1])
age = as.numeric(id[,5])

# get sample run ids for each tissue
ctx_samples = as.character(id[,4][region=="CTX"])
lv_samples = as.character(id[,4][region=="LV"])
ms_samples = as.character(id[,4][region=="MS"])
lng_samples = as.character(id[,4][region=="LNG"])

# get mouse ids for each tissue
ctx_id = mm_id[region=="CTX"]
lv_id = mm_id[region=="LV"]
ms_id = mm_id[region=="MS"]
lng_id = mm_id[region=="LNG"]

# get ages for each tissue
ctx_age = as.numeric(id[,5][region=="CTX"])
lv_age = as.numeric(id[,5][region=="LV"])
ms_age = as.numeric(id[,5][region=="MS"])
lng_age = as.numeric(id[,5][region=="LNG"])

# get unique run ids for each tissue from 3 lane run: 
ctx_name=c()
for(i in 1:length(ctx_samples)){
  if (substr(ctx_samples[i],3,3)==1 ) 
    ctx_name[i]=paste0("Sample_WGC088207-",ctx_samples[i])
  else if (substr(ctx_samples[i],3,3)==2) 
    ctx_name[i]=paste0("Sample_WGC088208-",ctx_samples[i])
  else ctx_name[i]=paste0("Sample_WGC088209-",ctx_samples[i])
}
ctx_name 

lv_name=c()
for(i in 1:length(lv_samples)){
  if (substr(lv_samples[i],3,3)==1) 
    lv_name[i]=paste0("Sample_WGC088207-",lv_samples[i])
  else if (substr(lv_samples[i],3,3)==2) 
    lv_name[i]=paste0("Sample_WGC088208-",lv_samples[i])
  else lv_name[i]=paste0("Sample_WGC088209-",lv_samples[i])
}
lv_name

ms_name=c()
for(i in 1:length(ms_samples)){
  if (substr(ms_samples[i],3,3)==1) 
    ms_name[i]=paste0("Sample_WGC088207-",ms_samples[i])
  else if (substr(ms_samples[i],3,3)==2) 
    ms_name[i]=paste0("Sample_WGC088208-",ms_samples[i])
  else ms_name[i]=paste0("Sample_WGC088209-",ms_samples[i])
}
ms_name

lng_name=c()
for(i in 1:length(lng_samples)){
  if (substr(lng_samples[i],3,3)==1) 
    lng_name[i]=paste0("Sample_WGC088207-",lng_samples[i])
  else if (substr(lng_samples[i],3,3)==2) 
    lng_name[i]=paste0("Sample_WGC088208-",lng_samples[i])
  else lng_name[i]=paste0("Sample_WGC088209-",lng_samples[i])
}
lng_name

# change order of samples as age increment
names(ctx_age) = ctx_id
names(ctx_samples) = ctx_id
ord = order(ctx_age)
ord[c(1,2,8,9)] = c(9,1,12,11) 
ctx_age = ctx_age[ord]
ctx_samples = ctx_samples[ord]
ctx_name = ctx_name[ord]
ctx_id = ctx_id[ord]

names(lng_age) = lng_id
names(lng_samples) = lng_id
ord = order(lng_age)
ord[c(1,2,15,16)] = c(11,6,14,1)
lng_age = lng_age[ord]
lng_samples = lng_samples[ord]
lng_name = lng_name[ord]
lng_id = lng_id[ord]

names(lv_age) = lv_id
names(lv_samples) = lv_id
ord = order(lv_age)
ord[c(1,2)] = c(4,2)
lv_age = lv_age[ord]
lv_samples = lv_samples[ord]
lv_name = lv_name[ord]
lv_id = lv_id[ord]

names(ms_age) = ms_id
names(ms_samples) = ms_id
ord = order(ms_age)
ord[c(1,2,8,9,15,16)] = c(10,3,16,8,14,5)
ms_age = ms_age[ord]
ms_samples = ms_samples[ord]
ms_name = ms_name[ord]
ms_id = ms_id[ord]

ind_id = c(ctx_id, lng_id, lv_id, ms_id)
names(ind_id) = rep(c('cortex','lung','liver','muscle'),times = c(15,16,16,16))

save(ctx_age, ctx_samples, ctx_id, ctx_name,
     lng_age, lng_samples, lng_id, lng_name,
     lv_age, lv_samples, lv_id, lv_name,
     ms_age, ms_samples, ms_id, ms_name,
     ind_id,
     file = "./data/preprocess/sample_features.rdata")

