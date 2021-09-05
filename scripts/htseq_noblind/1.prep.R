getwd()
library(DESeq2)
library(tidyverse)
library(ggpubr)
library("RColorBrewer")

# theme_set(theme(legend.key.size = unit(2,'pt'),
#                 plot.title = element_text(size=6)))

path = './raw_data/htseq/'
nms = readRDS('./data/processed/raw/fname_samples.rds')

# generate sampleTable
fnames = sapply(strsplit(list.files(path), split = '-|_'), function(x) paste(x[3], x[4], sep='_' ))
stable = data.frame(fname = fnames,
                    file = list.files(path))
stable = stable %>%
  left_join(nms) %>%
  mutate(ind_id =  as.factor(ind_id),
         tissue =  as.factor(tissue)) %>%
  mutate(log2age = log2(age))

# generate DESeq object from counts and sampletable with full design:
dsdat = DESeqDataSetFromHTSeqCount(sampleTable = stable,
                           directory = path,
                           design = ~ tissue + age + tissue:age)
# dsdat = DESeqDataSetFromHTSeqCount(sampleTable = stable,
#                                    directory = path,
#                                    design = ~ age )
# dsdat2 = DESeqDataSetFromHTSeqCount(sampleTable = stable,
#                                    directory = path,
#                                    design = ~ 1)
dsdat # 51826 genes

# keep at least 10 reads per gene:
keep = rowSums(counts(dsdat)) >= 10
dsdat = dsdat[keep,]
#dsdat2 = dsdat2[keep,]
dsdat # 30984 genes

# take only protein coding genes:
allgenetype = readRDS("./data/genetype.rds")
pgenes = allgenetype[allgenetype[,2]=="protein_coding",1]
dsdat = dsdat[rownames(counts(dsdat))%in%pgenes,]
#dsdat2 = dsdat2[rownames(counts(dsdat2))%in%pgenes,]
dsdat # 18552 genes

# Filtration:
# remove the genes which are not detected in at least 25% of the samples (15).
remov = apply(counts(dsdat),1,function(x) {sum(x==0)>15} )
dsdat = dsdat[!remov,]
#dsdat2 = dsdat2[!remov,]

dsdat #  14973 genes

# vst parametric normalisation:
vstdat = vst(dsdat, blind = T) # blinded: no design 
vstdat2 = vst(dsdat, blind = F) # no blind: take tissue and age effect into account

# vst non-parametric normalisation:
# vstdat_nonparam = vst(dsdat, blind=T, fitType = 'local')
# vstdat_nonparam2 = vst(dsdat, blind=F, fitType = 'local')

# rlog transformation: dont run this, takes too long. use vst instead
#rld = rlog(dsdat, blind = T)

## effect of transformations on variance
ntd1 = normTransform(dsdat)
library(vsn)
msd1 = meanSdPlot(assay(ntd1)) 
msd21 = meanSdPlot(assay(vstdat)) # no design
msd22 = meanSdPlot(assay(vstdat2)) 
# msd31 = meanSdPlot(assay(vstdat_nonparam))# no design
# msd32 = meanSdPlot(assay(vstdat_nonparam2)) 

norm_effect = 
  ggarrange(msd1$gg + ggtitle('log2+1 (vsn package)'), ggplot() ,
            msd21$gg + ggtitle('VST (parametric) - Blinded'),
            msd22$gg + ggtitle('VST (parametric) - No Blinding (tissue+age+Interaction)'),
            msd31$gg + ggtitle('VST (local fit) - Blinded'),
            msd32$gg + ggtitle('VST (local fit) - No Blinding (tissue+age+Interaction)'), nrow=3, ncol =2)

norm_effect
ggsave('./results/htseq/prep/vst_norms_on_var.pdf', norm_effect, units = 'cm',
       width = 25, height = 20, useDingbats=F)
ggsave('./results/htseq/prep/vst_norms_on_var.png', norm_effect, units = 'cm',
       width = 25, height = 20)

# median ratio method
dsdat_sf = estimateSizeFactors(dsdat)
dsdat_dips = estimateDispersions(dsdat_sf)
dsdat_dips
plotDispEsts(dsdat_dips)

### quality assesment
library(pheatmap)

#select = order(rowMeans(counts(dsdat, normalized=F)), decreasing = T)[1:50]
select = order(rowMeans(counts(dsdat_dips, normalized=T)), decreasing = T)[1:50]
df =  as.data.frame(colData(dsdat)[,c('log2age', 'tissue')])
ord = order(as.numeric(gsub('s','',stable$sample_id)))

p1 = pheatmap(assay(ntd1)[select,ord], cluster_rows = F, show_rownames = F,
              cluster_cols = T,  annotation_col = df)

p2 = pheatmap(assay(vstdat)[select,ord], cluster_rows = F, show_rownames = F,
              cluster_cols = T,  annotation_col = df)

p3 = pheatmap(assay(vstdat2)[select,ord], cluster_rows = F, show_rownames = F,
         cluster_cols = T,  annotation_col = df)

# p4 = pheatmap(assay(vstdat_nonparam)[select,ord], cluster_rows = F, show_rownames = F,
#               cluster_cols = T,  annotation_col = df)

ggsave('./results/htseq/prep/top50Count_vst_blinded.pdf' , p2, units='cm',
       height=15, width = 15, useDingbats=F)
ggsave('./results/htseq/prep/top50Count_vst_blinded.png' , p2, units='cm',
       height=15, width = 15)

ggsave('./results/htseq/prep/top50Count_vst_no_blind.pdf' , p3, units='cm',
       height=15, width = 15, useDingbats=F)

ggsave('./results/htseq/prep/top50Count_vst_no_blind.png' , p3, units='cm',
       height=15, width = 15)

# clustering

sampledists = dist(t(assay(vstdat)))
sampledistMat = as.matrix(sampledists)
rownames(sampledistMat) =  paste(vstdat$tissue, vstdat$age, sep='-')
colnames(sampledistMat) = NULL
colsx = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pcl1 = pheatmap(sampledistMat,
         clustering_distance_rows = sampledists,
         clustering_distance_cols = sampledists,
         col=colsx)

sampledists = dist(t(assay(vstdat2)))
sampledistMat = as.matrix(sampledists)
rownames(sampledistMat) =  paste(vstdat$tissue, vstdat$age, sep='-')
colnames(sampledistMat) = NULL
colsx = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pcl2 = pheatmap(sampledistMat,
                clustering_distance_rows = sampledists,
                clustering_distance_cols = sampledists,
                col=colsx)

ggsave('./results/htseq/prep/cluster_vst_blinded.pdf' , pcl1, units='cm',
       height=20, width = 25, useDingbats=F)
ggsave('./results/htseq/prep/cluster_vst_blinded.png' , pcl1, units='cm',
       height=20, width = 25)

ggsave('./results/htseq/prep/cluster_vst_no_blinded.pdf' , pcl2, units='cm',
       height=20, width = 25, useDingbats=F)
ggsave('./results/htseq/prep/cluster_vst_no_blind.png' , pcl2, units='cm',
       height=20, width = 25)


# pca

pcdat = plotPCA(vstdat, intgroup = c('tissue', 'age'), returnData=T)
expVar = round(100*attr(pcdat, 'percentVar'))
p_pc1 = ggplot(pcdat, aes(PC1, PC2, color= tissue, size=age)) +
  geom_point() + 
  scale_size_continuous(range=c(0.5,3), trans = 'log10') +
  xlab(paste0('PC1: ', expVar[1],'%' )) +
  ylab(paste0('PC1: ', expVar[2],'%' ))

pcdat = plotPCA(vstdat2, intgroup = c('tissue', 'age'), returnData=T)
expVar = round(100*attr(pcdat, 'percentVar'))
p_pc2 = ggplot(pcdat, aes(PC1, PC2, color= tissue, size=age)) +
  geom_point() + 
  scale_size_continuous(range=c(0.5,3), trans = 'log10') +
  xlab(paste0('PC1: ', expVar[1],'%' )) +
  ylab(paste0('PC1: ', expVar[2],'%' ))

ggsave('./results/htseq/prep/pca_vst_blinded.pdf' , p_pc1, units='cm',
       height=12, width = 18, useDingbats=F)
ggsave('./results/htseq/prep/pca_vst_blinded.png' , p_pc1, units='cm',
       height=12, width = 18)
ggsave('./results/htseq/prep/pca_vst_no_blind.pdf' , p_pc2, units='cm',
       height=12, width = 18, useDingbats=F)
ggsave('./results/htseq/prep/pca_vst_no_blind.png' , p_pc2, units='cm',
       height=12, width = 18)

#### save data for downstream analysis:

mat1 = assay(vstdat) # blinded
mat2 = assay(vstdat2) #  no blind

# expr = readRDS('./data/processed/raw/expression.rds')
# head(expr)
# dim(expr)
# head(mat)
# df1 =  setdiff(rownames(mat2), rownames(expr)) # 263 genes
# df2 =  setdiff(rownames(expr), rownames(mat2) ) # 354 genes

identical(stable$fname, colnames(mat1))
identical(stable$fname, colnames(mat2))

colnames(mat1) = stable$sample_id
colnames(mat2) = stable$sample_id
ord = order(as.numeric( gsub('s','', stable$sample_id) ))
mat1 = mat1[,ord]
mat2 = mat2[,ord]

saveRDS(mat1, './data/htseq/expr_mat_blinded.rds')
saveRDS(mat2, './data/htseq/expr_mat_no_blind.rds')

mat1t = reshape2::melt(mat1) %>%
  set_names(c('gene_id','sample_id','expression'))

mat2t = reshape2::melt(mat2) %>%
  set_names(c('gene_id','sample_id','expression'))

saveRDS(mat1t, './data/htseq/expr_blinded.rds')
saveRDS(mat2t, './data/htseq/expr_no_blind.rds')

