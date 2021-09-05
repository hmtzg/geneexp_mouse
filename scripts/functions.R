

########################
######################## Effect size using cohen's D formula:
########################

cohens_d=function(x,y){ 
  #calculates effect size of vectors a and b, using cohen's d by 
  #http://en.wikipedia.org/wiki/Effect_size#Cohen.27s_d
  lx=length(x)
  ly=length(y)
  ( mean(x)-mean(y) ) / ( ( ( (lx-1)*var(x) + (ly-1)*var(y) ) / ( lx+ly-2 ) )^0.5 )
}

########################
######################## Gene Overrepresentation Analysis (GORA) with topgo:
########################
go_bp_enrich.test.Mm = function(genelist, selection, description="",ID="Ensembl", nSizemin = 10, padj="BH",
                                nSizemax = 500){ 
  # selection 1 or 0
  require(topGO)
  require("org.Mm.eg.db")
  GOdata = new("topGOdata",ontology="BP", allGenes= genelist,geneSel=function(p){p==selection},
               description=description,annot=annFUN.org, nodeSize = nSizemin, mapping="org.Mm.eg.db",ID=ID)
  resultFisher = runTest(GOdata,algorithm="classic",statistic="fisher")
  sigGO = resultFisher@geneData["SigTerms"]
  genetable = GenTable(GOdata,classicFisher=resultFisher,topNodes=sigGO)
  genetable = genetable[genetable$Annotated <= nSizemax,]
  genetable_padj=cbind(genetable,p.adjust((genetable$classicFisher),method = padj) )
  names(genetable_padj)[7] = padj
  return(genetable_padj)
}


########################
######################## Calculate proportion of genes in two groups as continuous vs reversal:
########################

revgenes = function(dev,age,revgenes=T,monotonic=F){
  # provide a vector of rho values whose names are gene names.
  upD = names(which( dev>0 ))
  downD = names(which( dev<0 ))
  upA = names(which( age>0 ))
  downA = names(which( age<0 ))
  ud = intersect(upD,downA)
  du = intersect(downD,upA)
  uu = intersect(upD,upA)
  dd = intersect(downD,downA)
  
  if(revgenes & monotonic){
    return(list(UpDown=ud,DownUp=du,UpUp=uu,DownDown=dd))
  } else if(revgenes){
    return(list(UpDown=ud,DownUp=du))
  } else if(!revgenes & monotonic){
    return(list(UpUp=uu,DownDown=dd))
  }
}


########################
######################## leave one out CoDi ratio pseudovalues:
########################

LOO_CoDiR = function(mat, age, method='spearman', qval = 0.1, qmethod='BH', estimate.all=T){
  # mat: cov matrix with columns are individual names
  # age: ages of the individuals with the same order of mat
  # method: cov-age test method
  # qval: multiple testing correction threshold for age-cov test
  # qmethod: which multiple testing correction method to use
  # estimate.all: also calculate for all data together?
  
  if(estimate.all){
    ss = c('S', colnames(mat))
  } else{
    ss =  colnames(mat)
  }
  # exclude each column
  count=0
  codi = sapply(ss, function(ex){
    exmat = mat[, !colnames(mat)%in%ex]
    exage = age[!names(age)%in%ex]
    count <<- count + 1
    print(paste('%', round(count/length(ss),2)*100 ))
    sub_cov = t(sapply(rownames(exmat), function(x){
      a = cor.test(exmat[x,], exage, m=method)
      c(a$est, a$p.val)
    }))
    sub_cov = cbind(sub_cov, p.adjust(sub_cov[,2], m=qmethod))
    colnames(sub_cov) = c('rho', 'p', qmethod)
    sub_cov_sig = sub_cov[sub_cov[,3] < qval, ]
    codiR = sum(sub_cov_sig[,1] < 0) / sum(sub_cov_sig[,1] > 0)
    # Yb = sum(sub_cov_sig[,1] < 0)
    # Xb = sum(sub_cov_sig[,1] > 0)
    # codiR = c(Yb,Xb)
    return(codiR)
  })
  #names(codi)[-1] = paste0('-', names(codi)[-1])
  return(codi)
}

########################
######################## leave one out Confidence interval for pseudovalues:
########################

LOO_JK_CI = function(ps, alpha=0.05){
  #Si = na.omit(ps[-1])
  Si = log(na.omit(ps[-1]))
  #S = ps[1]
  S = log(ps[1])
  N = length(Si)
  
  psi = N*S - (N-1)*Si
  ps_bar = mean(psi)
  var = sum( (psi-ps_bar)^2) / (N-1) # usual variance estimate
  
  t_val =  qt(p = 1-alpha, df = N-1)
  CI_95 = c(ps_bar - t_val*sqrt(var/N),
            ps_bar + t_val*sqrt(var/N))
  return( exp(1)^CI_95 )
}
