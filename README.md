# schard
This package is intended to be pure R, reticulate-free replacement of sceasy. It is based on rhdf5 for h5ad manipulation. Just two functions are available right now (load visium h5ad into list of Seurat objects and load h5ad into SingleCellExperiment) and their functionality is limited, that is, they do not load full content of input h5ad.
# Installation
devtools::install_github("cellgeni/schard")
