library(biomaRt)
mart = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
convert = getBM(filters = "uniprot_gn_symbol", attributes = c("uniprot_gn_symbol","ensembl_gene_id"),
                mart = mart, values = genepool)
dup1 = unique(convert[duplicated(convert[,1]),1])
dup2 = unique(convert[duplicated(convert[,2]),2])
uconvert = convert[!convert[,1]%in%dup1,]
uuconvert = uconvert[!uconvert[,2]%in%dup2,]
saveRDS(uuconvert, './data/uniprot_to_ENS.rds')
