#' Loads scanpy h5ad file with multiple Visium data into list of Seurat object
#'
#' @param filename character, path to h5ad file
#' @param library_id_field name of adata.obs that specifies visium library. Function trys to guess it if set to NULL (default).
#'
#' @return list of Seurat object
#' @import Matrix
#' @export
h5ad2seurat_l_spatial = function (filename,library_id_field=NULL){
  loadRequiredPackages('Seurat')
  obsattr = rhdf5::h5readAttributes(filename,'obs')
  varattr = rhdf5::h5readAttributes(filename,'var')

  a = rhdf5::H5Fopen(filename,flags = 'H5F_ACC_RDONLY')
  tryCatch({
    obs = parseH5ADdataframe(a$obs,obsattr)
    var = parseH5ADdataframe(a$var,varattr)

    m = a$X
    mtx = sparseMatrix(i=m$indices+1, p=m$indptr,x = as.numeric(m$data),dims = c(nrow(var),nrow(obs)))
    rownames(mtx) = rownames(var)
    colnames(mtx) = rownames(obs)

    # find column in obs all values of each has spatial data in a$uns$spatial
    if(is.null(library_id_field)){
      imglibs = names(a$uns$spatial)
      for(obsf in colnames(obs)){
        if(all(obs[,obsf] %in% imglibs)){
          library_id_field = obsf
          break
        }
      }
    }

    res = list()
    for(lid in unique(obs[[library_id_field]])){
      # subset obs
      f = obs[[library_id_field]] == lid
      obs_ = obs[f,]
      if(!is.null(obs_$barcode)){
        rownames(obs_) = obs_$barcode
      }
      # subset counts
      mtx_ = mtx[,which(f)]
      colnames(mtx_) = rownames(obs_)

      object <- CreateSeuratObject(counts = mtx_, assay = 'Spatial')

      # extract image
      image <- aperm(a$uns$spatial[[lid]]$images$hires,3:1)

      scale.factors <- a$uns$spatial[[lid]]$scalefactors
      # looks like in scanpy both images are hires actually (at least they were identical for example I tried)
      scale.factors$tissue_lowres_scalef = scale.factors$tissue_hires_scalef
      tissue.positions = cbind(obs[f,c('in_tissue','array_row','array_col')], t(a$obsm$spatial[,f])[,2:1])
      colnames(tissue.positions) = c("tissue", "row", "col", "imagerow", "imagecol")

      rownames(tissue.positions) = rownames(obs_)

      unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
      spot.radius <- unnormalized.radius/max(dim(x = image))
      image = new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                                                                  fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef,
                                                                                  scale.factors$tissue_lowres_scalef),
                  coordinates = tissue.positions, spot.radius = spot.radius)

      image <- image[Cells(x = object)]
      DefaultAssay(object = image) = 'Spatial'
      object[['slice1']] = image

      object@meta.data = cbind(object@meta.data,obs_)
      rownames(object@meta.data) = rownames(obs_)

      object@assays$Spatial@meta.features = var
      res[[lid]] = object
    }
  },error=function(e){
    rhdf5::H5Fclose(a)
    stop(e)
  })

  rhdf5::H5Fclose(a)
  return(res)
}

loadRequiredPackages = function(pkgs){
  for(pkg in pkgs){
    if(!require(pkg,character.only =TRUE))
      stop(paste0('Please install package "',pkg,'"'))
  }
}

#' Extracts dataframe from H5AD list
#' function is used in conjunction with rhdf5::H5Fopen from myLoadH5AD_Spatials function
#'
#' @param collist list read from H5AD object
#' @param attr
#'
#' @return a data.frame
parseH5ADdataframe = function(collist,attr){
  # slashes in names leads to nested structure, lets fix it
  ll = sapply(collist,length)
  for(i in which(ll==1)){
    n = names(collist)[[i]]
    repeat{
      if(length(collist[[i]])>1)
        break
      n = paste0(n,'/',names(collist[[i]]))
      collist[[i]] = collist[[i]][[1]]
    }
    names(collist)[i] = n
  }

  # factors are stored as list of categories and names, other types are stored as vectors
  # take them first
  ll = sapply(collist,length)
  res = as.data.frame(collist[ll==max(ll)],check.names=FALSE)

  # first way to store factors
  for(fn in names(collist[['__categories']])){
    res[[fn]] = as.vector(collist[['__categories']][[fn]][res[[fn]]+1])
  }
  # another way to store factors
  for(fn in names(ll)[ll==2]){
    if(all(names(collist[[fn]]) %in% c("categories","codes"))){
      res[[fn]] = as.vector(collist[[fn]]$categories[collist[[fn]]$codes+1])
    }
  }
  rownames(res) = collist[[attr$`_index`]]
  res[,c(attr$`_index`,attr$`column-order`),drop=FALSE]
}

#' Loads h5ad scanpy single cell file as SingleCellExperiment
#'
#' functions exports X,obs,var and obsm
#' it is very experimental and comes with no warranty
#'
#' @param filename path to h5ad file
#'
#' @return SingleCellExperiment object
#' @export
#'
#' @import Matrix
#' @examples
#' sce = h5ad2sce('adata.h5ad')
h5ad2sce = function(filename){
  loadRequiredPackages('SingleCellExperiment')
  obsattr = rhdf5::h5readAttributes(filename,'obs')
  varattr = rhdf5::h5readAttributes(filename,'var')
  xattr = rhdf5::h5readAttributes(filename,'X')

  a = rhdf5::H5Fopen(filename,flags = 'H5F_ACC_RDONLY')

  tryCatch({
    obs = parseH5ADdataframe(a$obs,obsattr)
    var = parseH5ADdataframe(a$var,varattr)

    m = a$X
    if(is.array(m)){
      mtx = m
    }else{
      if(xattr$`encoding-type`=='csr_matrix'){
        mtx = sparseMatrix(i=m$indices+1, p=m$indptr,x = as.numeric(m$data),dims = c(nrow(var),nrow(obs)))
      }else{
        mtx = sparseMatrix(i=m$indices+1, p=m$indptr,x = as.numeric(m$data),dims = c(nrow(obs),nrow(var)))
        mtx = Matrix::t(mtx)
      }
    }
    rownames(mtx) = rownames(var)
    colnames(mtx) = rownames(obs)

    sce = SingleCellExperiment(list(X=mtx),
                               colData=obs,
                               rowData=var
    )
    # reduced dims
    for(n in names(a$obsm)){
      reducedDim(sce,n) = t(a$obsm[[n]])
    }
    rhdf5::H5Fclose(a)
  },error=function(e){
    rhdf5::H5Fclose(a)
    stop(e)
  })
  sce
}


#' Title
#'
#' @param obj
#' @param outFile
#'
#' @return
#'
#' @examples
# cnt = Seurat::Read10X('~/nfs/projects/2307.omer/data.nfs/visium_hd/Visium_HD_hprostate_mbrain/mBrain/square_050um/filtered_feature_bc_matrix',gene.column=1)
# obj = CreateSeuratObject(cnt)
# genes = read.table('~/nfs/projects/2307.omer/data.nfs/visium_hd/Visium_HD_hprostate_mbrain/mBrain/square_050um/filtered_feature_bc_matrix/features.tsv.gz',row.names = 1)
# obj@assays[[assay]]@meta.features = genes[rownames(obj),]
seurat2h5ad = function(obj,outFile,assay = "RNA", main_layer = "data"){
  stop("function is not finished yet")
  if (!requireNamespace("Seurat")) {
    stop("This function requires the 'Seurat' package.")
  }

  X = Seurat::GetAssayData(object = obj, assay = assay, slot = main_layer)
  X = Matrix::t(X)
  obs = obj@meta.data
  var = Seurat::GetAssay(obj, assay = assay)@meta.features
  SeuratDisk::WriteH5Group()

}

