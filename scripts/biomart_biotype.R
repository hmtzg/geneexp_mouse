library(biomaRt)
listMarts()
mus.mart = useMart(biomart = 'ensembl')
mus.ensembl = useDataset('mmusculus_gene_ensembl', mart = mus.mart )
genetype = getBM(attributes = c('ensembl_gene_id', 'gene_biotype'), filters = 'biotype',
                 values = '',mart = mus.ensembl)

saveRDS(genetype, './data/biomart_genetype.rds')
saveRDS(mus.ensembl, './data/biomart_ensemblmart.rds')

'''
ensemblmart and genetype objects created in 21.03.2019
'''