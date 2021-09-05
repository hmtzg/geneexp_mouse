library(goseq)
allgo=getgo(allgenes,"hg19","ensGene")
allgo = allgo[!sapply(allgo,is.null)]
allgo = reshape2::melt(allgo) %>%
  set_names(c('category','gene'))
signifGO_genes = left_join(signifGO,allgo) #signifGO GO enrichment table i
### burasi benim analize ozel kisim ###
signifGO_genes$gene_cluster = sapply(signifGO_genes$gene,function(x){
  nm = unique(sapply(strsplit(names(which(sapply(genegr,function(y)x%in%y))),'_'),function(x)x[1]))
  ifelse(length(nm)==0,'-',gsub('cl','',nm))
})
######################################
gogenemat = signifGO_genes %>%
#  filter(gene_cluster!='-') %>% # benim analize ozel, sil
  dplyr::select(category, gene) %>%
  unique() %>%
  mutate(value = 1) %>%
  spread(category, value, fill = 0)

# gen - category matrixi olusturup jaccard hesapla
gogenemat = as.data.frame(gogenemat)
rownames(gogenemat) = gogenemat$gene
gogenemat$gene = NULL
gogenemat = as.matrix(gogenemat)
library(pheatmap)
jaccardsim = apply(gogenemat,2,function(x){
  apply(gogenemat,2,function(y){
    sum(x==1 & y==1) / sum(x==1 | y==1)
  })
})

# ben analizi BP, MF, CC hepsi icin yapiyordum, biz sanirim BP kullaniyoruz sadece, onunla limitleyebiliriz. cluster kismi da benim analize ozeldi o da olmayacak yani for loop olmadan kod bir kez calisacak, kod icindeki cl iceren kisimlar silinecek.
newdf = data.frame()
for(ont in c('BP','MF','CC')){
  for(cl in c('1','2','3','1-2','1-3','2-3','1-2-3')) {
    print(c(ont,cl))
    gocatx = (signifGO %>% filter(ontology == ont & cluster == cl))
    gocat= unique(gocatx$category)
    if(length(gocat)==0){

    } else if( length(gocat)==1){
      newdf = rbind(newdf, mutate(gocatx,rep=gocat))
    } else {
      simmat = jaccardsim[gocat,gocat]
      k = min(which(sapply(1:(nrow(simmat)),function(i){
        mycl = cutree(hclust(dist(simmat)),i)
        all(sapply(1:i,function(k){
          xx = names(which(mycl==k))
          if(length(xx)==1){
            return(1)
          } else {
            xx = simmat[xx,xx]
            median(xx[upper.tri(xx)])
          }
        })>=0.5) # clusterlari median jaccardi 0.5 uzeri olan gruplar olarak tanimlamisim, burayi degistirip daha az ya da daha cok rep. alabiliriz.
      })))
      k=ifelse(k==1,nrow(simmat),k)
      treex =hclust(dist(simmat))
      treecl = cutree(treex,k)
      reps = sapply(1:k,function(i){
        xx=names(which(treecl==i))
        if(length(xx)>1){
          xx = simmat[xx,xx]
          names(which.max(rowMeans(xx)))
        } else{
          xx
        }
      })
      repclus=setNames(lapply(1:k,function(i)names(which(treecl==i))),reps)
      newdf = rbind(newdf, reshape2::melt(repclus) %>%
                      set_names(c('category','rep')) %>%
                      arrange(rep) %>%
                      mutate(ontology = ont, cluster = cl) %>%
                      left_join(signifGO) %>%
                      unique() %>%
                      dplyr::select(1,3:23,2))
    }}}


# tum dosyayi table olarak verip, sadece representativeleri iceren bir listeyi de gorsellestirmistim (bar plot, x ekseni enrichment, y ekseni go)
newdf %>%
  filter(rep==category) %>%
#  filter(cluster=='1') %>% # bu benim analiz icin gerekli olan kisimdi, silinecek
  mutate(termname = ifelse(nchar(term)>40,paste(substr(term,1,37),'...',sep=''),term)) %>% #burada uzun GO termleri kisaltip yazmistim.
  unique()
