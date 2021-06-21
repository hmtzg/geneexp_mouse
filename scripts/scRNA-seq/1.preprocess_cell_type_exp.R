library(Seurat)
# data downloaded from: https://figshare.com/articles/dataset/Tabula_Muris_Senis_Data_Objects/12654728
# metadata downloaded from: 
# https://s3.console.aws.amazon.com/s3/buckets/czb-tabula-muris-senis/Metadata/?region=us-west-2&tab=overview

####################
#################### Load annotated Seurat data and metadata
####################
####################
lung = ReadH5AD(
  "./data/other_datasets/scRNA-seq/raw/tabula-muris-senis-facs-processed-official-annotations-Lung.h5ad")
liver = ReadH5AD(
  "./data/other_datasets/scRNA-seq/raw/tabula-muris-senis-facs-processed-official-annotations-Liver.h5ad")
muscle = ReadH5AD(
  "./data/other_datasets/scRNA-seq/raw/tabula-muris-senis-facs-processed-official-annotations-Limb_Muscle.h5ad")
brain = ReadH5AD(
    "./data/other_datasets/scRNA-seq/raw/tabula-muris-senis-facs-processed-official-annotations-Brain_Non-Myeloid.h5ad")
metadata = read.csv(
  "./data/other_datasets/scRNA-seq/raw/tabula-muris-senis-facs-official-raw-obj__cell-metadata__cleaned_ids.csv")

#################### combine cell type expression matrices of tissues
tis.exp = list(lung = Seurat::GetAssayData(lung, slot="data"),
               liver = Seurat::GetAssayData(liver, slot="data"),
               muscle = Seurat::GetAssayData(muscle, slot="data"),
               brain = Seurat::GetAssayData(brain, slot="data"))

#################### combine metadata of tissues
metadata = list(lung = metadata[match(rownames(lung@meta.data),metadata$obs_names),],
                liver = metadata[match(rownames(liver@meta.data),metadata$obs_names),],
                muscle = metadata[match(rownames(muscle@meta.data),metadata$obs_names),],
                brain = metadata[match(rownames(brain@meta.data),metadata$obs_names),])

rm(lung,liver,muscle,brain)

for(ts in names(metadata)){
  metadata[[ts]]$mouse.id = as.character(metadata[[ts]]$mouse.id)
  metadata[[ts]]$obs_names = as.character(metadata[[ts]]$obs_names)
  metadata[[ts]]$cell_ontology_class = as.character(metadata[[ts]]$cell_ontology_class)
}


####################
#################### remove cell types that are identified in less than 15 cells in all time points. 
#################### (apply only for lung, other organs have high number of cells for each cell type)
####################

range(rowSums(table(metadata$brain$cell_ontology_class,metadata$brain$mouse.id))) # min 17 cells
range(rowSums(table(metadata$muscle$cell_ontology_class,metadata$muscle$mouse.id))) # min 136 cells
range(rowSums(table(metadata$liver$cell_ontology_class,metadata$liver$mouse.id))) # min 30 cells
sort(rowSums(table(metadata$lung$cell_ontology_class,metadata$lung$mouse.id))) # min 4

rm.ct = rowSums(table(metadata$lung$cell_ontology_class,metadata$lung$mouse.id))
rm.ct = names(rm.ct[rm.ct<15]) # remove cell types with < 15 cells: 7 cell types

#################### remove cell types from metadata and expression matrix of lung
metadata$lung = metadata$lung[!metadata$lung$cell_ontology_class%in%rm.ct,]
tis.exp$lung = tis.exp$lung[,colnames(tis.exp$lung)%in%metadata$lung$obs_names]

#################### check if expression matrix columns and metadata observation names are in same order:
for(i in 1:4) print(identical(metadata[[i]]$obs_names, colnames(tis.exp[[i]])))


####################
#################### calculate average expression for each cell types for each tissue and age
#################### 
####################


#################### separate expression matrices into  into different age groups in a list
expr = sapply(c('3m', '18m','24m'), function(ag){
  sapply(names(tis.exp), function(x){
    tis.exp[[x]][, metadata[[x]]$age == ag]
  })
}, simplify = F)
names(expr) = c('m3','m18','m24')

#################### remove genes not expressed in any of the cells for each age group:
expr = lapply(expr,function(ag){
  sapply(names(ag), function(x){
    ag[[x]][!apply(ag[[x]]==0, 1, function(y) sum(y) == ncol(ag[[x]])) , ]
  })
})

rm(tis.exp)

#################### separate metadata into different lists:
metadata = sapply(c('3m', '18m','24m'), function(ag){
  sapply(names(metadata), function(x){
    metadata[[x]][metadata[[x]]$age == ag,]
  }, simplify = F)
}, simplify = F)
names(metadata) = c('m3','m18','m24')

####################  check if expression matrix columns and metadata cell rows are in same order:
sapply(names(metadata), function(ag){
  sapply(names(metadata[[ag]]), function(x){
    identical(metadata[[ag]][[x]]$obs_names, colnames(expr[[ag]][[x]]) )  
  })
})

#################### calculate mean expression for each cell type:
# first take mean of cells, then take mean of individuals for each cell type for each age group,
# taking into account zero expressions when calculating average:

ct_expr = sapply(names(expr), function(ag){
  ctexp = sapply(names(expr[[ag]]), function(ts){
    cell.class = unique(metadata[[ag]][[ts]]$cell_ontology_class)
    mouse.id = unique(metadata[[ag]][[ts]]$mouse.id)
    cell.class.exp = sapply(1:length(cell.class), function(x){
      f = sapply(1:length(mouse.id), function(y){
        cell.names = metadata[[ag]][[ts]][metadata[[ag]][[ts]]$cell_ontology_class == cell.class[x] &
                                            metadata[[ag]][[ts]]$mouse.id == mouse.id[y], 1]
        if(length(cell.names)!=0){
          rowMeans(as.matrix(expr[[ag]][[ts]][,cell.names]))
        }
      },simplify = F)
      f = f[!sapply(f,is.null)]
      rowMeans(do.call(cbind,f))
    })
    colnames(cell.class.exp) = cell.class
    return(cell.class.exp)
  })
  return(ctexp)
}, simplify = F)

saveRDS(ct_expr,'./data/other_datasets/scRNA-seq/processed/celltype_expr.rds')

# saveRDS(cell.type.exp.3m,file="Dropbox/projects/ageing/proc_data/scRNA/3m.celltype.rds")
# saveRDS(cell.type.exp.18m,file="Dropbox/projects/ageing/proc_data/scRNA/18m.celltype.rds")
# saveRDS(cell.type.exp.24m,file="Dropbox/projects/ageing/proc_data/scRNA/24m.celltype.rds")

####################
#################### calculate cell type expression per individual:
#################### 
####################

ct_expr_per = sapply(names(expr), function(ag){
  ctexp = sapply(names(expr[[ag]]),function(ts){
    cell.class = unique(metadata[[ag]][[ts]]$cell_ontology_class)
    mouse.id = unique(metadata[[ag]][[ts]]$mouse.id)
    cell.class.exp = sapply(1:length(cell.class), function(x){
      f = sapply(1:length(mouse.id), function(y){
        cell.names = metadata[[ag]][[ts]][metadata[[ag]][[ts]]$cell_ontology_class == cell.class[x] &
                                            metadata[[ag]][[ts]]$mouse.id == mouse.id[y], 1]
        if(length(cell.names)!=0){
          rowMeans(as.matrix(expr[[ag]][[ts]][,cell.names]))
        }
      },simplify = F)
      rm_null_ind = sapply(f,is.null)
      f = f[!rm_null_ind]
      indX.exp = do.call(cbind,f)
      colnames(indX.exp) = mouse.id[!rm_null_ind]
      return(indX.exp)
    },simplify = F)
    names(cell.class.exp) = cell.class
    return(cell.class.exp)
  })
  return(ctexp)
}, simplify = F)

saveRDS(ct_expr_per,'./data/other_datasets/scRNA-seq/processed/celltype_expr_per_ind.rds')

# saveRDS(cell.type.ind.exp.3m, file="Dropbox/projects/ageing/proc_data/scRNA/3m.celltype.per.ind.rds")
# saveRDS(cell.type.ind.exp.18m, file="Dropbox/projects/ageing/proc_data/scRNA/18m.celltype.per.ind.rds")
# saveRDS(cell.type.ind.exp.24m, file="Dropbox/projects/ageing/proc_data/scRNA/24m.celltype.per.ind.rds")

####################
#################### calculate number of cells per cell type in each age group:
#################### 
####################

ct_prop = sapply(names(expr), function(ag){
  celltype_prop = sapply(names(expr[[ag]]),function(ts){
    cell.class = unique(metadata[[ag]][[ts]]$cell_ontology_class)
    mouse.id = unique(metadata[[ag]][[ts]]$mouse.id)
    cell.class.cc = sapply(1:length(cell.class), function(x){ # cell class Cell Count
      f = sapply(1:length(mouse.id), function(y){
        cell.names = metadata[[ag]][[ts]][metadata[[ag]][[ts]]$cell_ontology_class == cell.class[x] &
                                            metadata[[ag]][[ts]]$mouse.id == mouse.id[y], 1]
        return(length(cell.names))
      })
      names(f) = mouse.id
      return(f)
    }) 
    cell.class.cc = t(scale(t(cell.class.cc),scale=rowSums(cell.class.cc), center = F))
    colnames(cell.class.cc) = cell.class
    return(cell.class.cc)
  })
  return(celltype_prop)
}, simplify = F)


saveRDS(ct_prop, './data/other_datasets/scRNA-seq/processed/celltype_proportions.rds')

# saveRDS(cell.type.prop.3m, file='Dropbox/projects/ageing/proc_data/scRNA/3m.cell.prop.rds')
# saveRDS(cell.type.prop.18m, file='Dropbox/projects/ageing/proc_data/scRNA/18m.cell.prop.rds')
# saveRDS(cell.type.prop.24m, file='Dropbox/projects/ageing/proc_data/scRNA/24m.cell.prop.rds')

