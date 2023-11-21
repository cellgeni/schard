# schard
This package allows one to load scanpy h5ad into R as list, SingleCellExperiment or Seurat object. For now it only loads `X`, `obs`, `var`, `obsm` (as reduced dimensions) if requested and `images` for visium data. 
The package is based on rhdf5 for h5ad manipulation and is pure R (that is reticulate-free).

# Installation
schard depends on some Bioconductor packages that should be installed manually:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("rhdf5","Seurat", "SingleCellExperiment"))
```

Then install schard from github
```
devtools::install_github("cellgeni/schard")
```

# Usage
Download some public h5ad and load them to R:
```
download.file('https://covid19.cog.sanger.ac.uk/baron16.processed.h5ad','ba16.h5ad') # old data from 2016
# visium and single nuclei datasets from cziscience: https://cellxgene.cziscience.com/collections/3116d060-0a8e-4767-99bb-e866badea1ed
download.file('https://datasets.cellxgene.cziscience.com/c5ac5c36-f60c-4680-8018-2d6cb65c0a37.h5ad','vis.heart.h5ad')
download.file('https://datasets.cellxgene.cziscience.com/8cc521c8-c4ff-4cba-a07b-cae67a9dcba9.h5ad','sn.heart.h5ad')

# load h5ad as Single Cell Experiment
ba16.sce = schard::h5ad2sce('ba16.h5ad')
# load h5ad as Seurat
snhx = schard::h5ad2seurat('sn.heart.h5ad')
# load all visium samples as single Seurat object
visx = schard::h5ad2seurat_spatial('vis.heart.h5ad')
# or load as list of visium objects
visl = schard::h5ad2seurat_spatial('vis.heart.h5ad',simplify = FALSE)
# or load raw counts
snhr = schard::h5ad2seurat('sn.heart.h5ad',use.raw = TRUE)
# raw counts for visium
visr = schard::h5ad2seurat_spatial('vis.heart.h5ad',use.raw = TRUE)


# check that it works
Seurat::SpatialPlot(visx,features = 'total_counts')
Seurat::SpatialPlot(visx,features = 'total_counts',images = 'HCAHeartST11702009')
Seurat::SpatialPlot(visl$HCAHeartST11702010,features = 'total_counts')
plot(colSums(visx),colSums(visr),pch=16) # raw counts are different from normolized ones
Seurat::DimPlot(snhx,group.by = 'cell_state') # the name of reduction is 'Xumap_' (autotranslated from scanpy to Seurat), somehow DimPlot manages to find it, but probably safier to specify it manually with reduction = 'Xumap_'
```

There is not need to load whole object if you only need cell metadata:
```
obs = schard::h5ad2data.frame('sn.heart.h5ad','obs')
# one can load umap in the same way.
# not sure about the name? lets see into h5ad:
ls = rhdf5::h5ls('sn.heart.h5ad')
ls[ls$group=='/obsm',] # so we have 'X_umap' here
umap = t(schard::h5ad2Matrix('sn.heart.h5ad','/obsm/X_umap')) # I like it more transposed
plot(umap[,1:2],pch=16,cex=0.4,col=factor(obs$cell_state))
```

# Alternatives
There are two known alternatives:
1. [sceasy](https://github.com/cellgeni/sceasy) uses reticulate and thus depends on python environment. Proved to be unstable and hard to use.
2. [SeuratDisk](https://github.com/mojaveazure/seurat-disk) also uses rhdf5, but uses h5-based Seurat format as an intermediate that looks like overcomplication. Additionally,  SeuratDisk seems to be almost not supported and it fails even on examples from its own tutorial.

Despite all problems of both packages above they have clear advantage over schard: they allow not only to read h5ad into R but also to write it.
