library(tidyverse)
library(topGO)
#revgenes = readRDS('./data/processed/raw/revgenes.tissues.rds')
#head(revgenes$Cortex)
revg = readRDS('./data/processed/tidy/revgenes.tissues.rds')
revgsh = readRDS('./data/processed/tidy/revgenes.shared.rds')
head(revg)

revg.cortex = revg %>%
  filter(tissue == 'Cortex')
head(revg.cortex)

xx <- annFUN.org("BP", mapping = "org.Mm.eg.db", ID = "ensembl")
irsp = xx$`GO:0008286`
length(irsp) #  57 gene
revg.cortex %>% 
  filter(gene_id%in%irsp) %>%
  group_by(direction) %>%
  summarise(n=n())
# # A tibble: 4 x 2
# direction     n
# <chr>     <int>
# 1 DownDown     21
# 2 DownUp        7
# 3 UpDown       10
# 4 UpUp         13

#####################

term = select(GO.db, columns = c('GOID', 'TERM'), keytype = 'GOID', keys = keys(GO.db) )
terms = as.character(term[,2])
names(terms) = as.character(term[,1])

term[grep('insulin', term$TERM),]
x <- c(as.list(GOBPOFFSPRING)) # children and all their children
xx <- sapply(1:length(x), function(i) x[[i]] = c(names(x)[i], x[[i]]) )
xxx <- sapply(1:length(xx), function(i) xx[[i]][ complete.cases(xx[[i]]) ] )
length(x)
length(xx)
length(xxx)
names(xx) <- names(x)
names(xxx) <- names(x)
xxx[[1]]
go2gene = select(org.Mm.eg.db, keys = names(xxx), columns = c("GO","ENSEMBL"), keytype = "GO")
GO2Gene = list()
for(i in 1:length(xxx)){
  print(i/length(xxx))
  myg = go2gene[go2gene[,1]%in%xxx[[i]], 4]
  myg = unique(myg[complete.cases(myg)])
  GO2Gene[[i]] = myg
}
names(GO2Gene) = names(xxx)
head(GO2Gene)
length(GO2Gene)

terms['GO:0008286']
irsp = GO2Gene$`GO:0008286`
revg %>%
  filter(gene_id%in%irsp) %>%
  mutate(rev = ifelse(direction =='UpDown' | direction =='DownUp','rev', 'notrev') ) %>%
  group_by(rev, tissue) %>%
  summarise(n=n()) %>%
  spread(key = rev, value = n) %>%
  summarise(tissue, prop = rev/(rev+notrev))
# # A tibble: 4 x 2
# tissue  prop
# <chr>  <dbl>
# 1 Cortex 0.333
# 2 Liver  0.596
# 3 Lung   0.538
# 4 Muscle 0.569


# GO:1901142 
# "insulin metabolic process" 
irsp2 = GO2Gene$`GO:1901142`

revg %>%
  filter(gene_id%in%irsp2) %>%
  mutate(rev = ifelse(direction =='UpDown' | direction =='DownUp','rev', 'notrev') ) %>%
  group_by(rev, tissue) %>%
  summarise(n=n()) %>%
  spread(key = rev, value = n) %>%
  summarise(tissue, prop = rev/(rev+notrev))


devdi = readRDS('./data/processed/raw/dev_divergent_genes_dc_rank.rds')
dico = names(devdi[devdi<0])

sum(irsp%in%dico) # 31
sum(irsp2%in%dico) # 3


