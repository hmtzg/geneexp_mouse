


########################
######################## Gene Overrepresentation Analysis (GORA) with topgo:
########################
go_bp_enrich.test.Mm = function(genelist, selection, description="",ID="Ensembl", nSizemin = 10, padj="BY",
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
######################## robust kmeans
########################

robust.kmeans.1.f = function( 
  Mat,  #normalized expression matrix
  k2try,  #which K to try, just one integer
  REPEAT=1000   #number of repetitions
) {
  #this wraps up k2try.2.f and postk2try.2.f in one function, but allows only one K to be tested
  km_l1 = lapply(k2try, function(k) { print(k)
    tab = lapply(1:REPEAT, function(i) {
      km1 = kmeans(Mat, k)$cluster
      x = sapply(1:length(table(km1)), function(i) sort(names(km1)[km1==i])) #order genes in cluster based on gene names
      xo = order(sapply(x, function(x) x[[1]])) #order cluster based on first gene names
      return(x[xo])
    } )
  } )
  names(km_l1) = k2try
  tabx = sapply(km_l1, function(x) sort(sapply(x, function(xx) paste(unlist(xx), collapse="_"))))
  Tab = apply(tabx, 2, function(x) as.numeric(rev(sort(table(x)))[1:3])); print(Tab)
  km_best = names(rev(sort(table(tabx[,as.character(k2try)])))[1])
  km1 = {
    i = 0; while (i != 1) {
      cat(i); km1 = kmeans(Mat, k2try)$cluster;
      x = sapply(1:length(table(km1)), function(i) sort(names(km1)[km1==i])) #order genes in cluster based on gene names
      xo = order(sapply(x, function(x) x[[1]])) #order cluster based on first gene names
      kmx = paste(unlist(x[xo]), collapse="_")
      if (km_best == kmx) { KmX = km1; i = 1 } }
    Kmc = Kmc2 = km1;
    for (J in 1:k2try) { Kmc2[Kmc %in% names(rev(sort(table(Kmc))))[J]] = J }
    Kmc2 }
  return(km1)
}



