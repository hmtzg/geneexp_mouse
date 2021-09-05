#source("Dropbox/projects/Rfunctions.R")
load('./data/preprocess/sample_features.rdata')
age = c(ctx_age,lng_age, lv_age, ms_age)
uage = lng_age
#l(uage) # 16

uageD = uage[uage < 90]
uageA = uage[uage > 90]
permageD=list()
permageA=list()
for(k in 1:1000){
  sampageD = sample(uageD)
  names(sampageD) = names(uageD)
  sampageA = sample(uageA)
  names(sampageA) = names(uageA)
  
  ctx.sampD = sampageD[ctx_id[ctx_age<90]]
  ctx.sampA = sampageA[ctx_id[ctx_age>90]]
  lng.sampD = sampageD[lng_id[lng_age<90]]
  lng.sampA = sampageA[lng_id[lng_age>90]]
  lv.sampD = sampageD[lv_id[lv_age<90]]
  lv.sampA = sampageA[lv_id[lv_age>90]]
  ms.sampD = sampageD[ms_id[ms_age<90]]
  ms.sampA = sampageA[ms_id[ms_age>90]]
  
  permageD[[k]] = list("ctx" = ctx.sampD, "lng" = lng.sampD, "lv" = lv.sampD, "ms" = ms.sampD)
  permageA[[k]] = list("ctx" = ctx.sampA, "lng" = lng.sampA, "lv" = lv.sampA, "ms" = ms.sampA)
  print(k/1000)
}

save(permageD,permageA,file='./data/processed/raw/ageperm1000.RData')
