#' Loads h5ad file as ordinary R list
#'
#' @param filename path to h5ad file
#' @param load.raw character, either 'no' (do not load raw), 'both' (load both adata.raw nad adata.X), or 'replace' (load raw only if exists)
#'
#' @return list with following elements:
#'  obs - data.frame
#'  var - data.frame
#'  X - matrix (dense or sparse - as it was in X)
#'  reducedDim - list of matrices with reduced dimmensions
#'  raw - list, eitehr empty (if load.raw is false or if raw is not present)
#' @export
h5ad2list = function(filename,load.raw='no'){
  a = rhdf5::H5Fopen(filename,flags = 'H5F_ACC_RDONLY')

  tryCatch({
    obs = parseH5ADdataframe(a,'obs')
    var = parseH5ADdataframe(a,'var')
    raw = list()
    X = NULL

    if(load.raw %in% c('no','both') || !rhdf5::H5Lexists(a,'raw/X')){
      X = pasreH5ADMatrix(a,'X')
      rownames(X) = rownames(var)
      colnames(X) = rownames(obs)
    }
    if(load.raw %in% c('both','replace')){
      if(rhdf5::H5Lexists(a,'raw/X')){
        raw$X = pasreH5ADMatrix(a,'raw/X')
        raw$var = parseH5ADdataframe(a,'raw/var')
        rownames(raw$X) = rownames(raw$var)
        colnames(raw$X) = rownames(obs)
      }else{
        warning('adata.raw is not found, will use adata.X instead')
      }
    }
    if(load.raw == 'replace'){
      X = raw$X
      var = raw$var
      raw = list()
    }
    # reduced dims
    reducedDim = list()
    for(n in names(a$obsm)){
      if(is.array(a$obsm[[n]]))
        reducedDim[[n]] = t(a$obsm[[n]])
    }
    rhdf5::H5Fclose(a)
  },error=function(e){
    rhdf5::H5Fclose(a)
    stop(e)
  })
  list(obs=obs,var=var,X=X,reducedDim = reducedDim,
       raw=raw)
}


#' Loads h5ad scanpy single cell file as SingleCellExperiment
#'
#' functions exports X,obs,var and obsm
#' it is very experimental and comes with no warranty
#'
#' @param filename path to h5ad file
#' @param use.raw logical, whether to attempt to get data from adata.raw. Throws warning (not exception!) if adata.raw is not present and proceeds with adata.X.
#'
#' @return SingleCellExperiment object
#' @export
#'
#' @import Matrix
#' @examples
#' sce = h5ad2sce('adata.h5ad')
h5ad2sce = function(filename,use.raw=FALSE){
  loadRequiredPackages('SingleCellExperiment')
  data = h5ad2list(filename,load.raw = ifelse(use.raw,'replace','no'))

  sce = SingleCellExperiment(list(X=data$X),
                             colData=data$obs,
                             rowData=data$var
  )
  # reduced dims
  for(n in names(data$reducedDim)){
    reducedDim(sce,n) = data$reducedDim[[n]]
  }
  sce
}

#' Loads scanpy h5ad file with Visium data into list of Seurat objects
#'
#' simplifies to single Seurat object if there is just one sample (see simplify parameter)
#'
#' @param filename character, path to h5ad file
#' @param library_id_field name of adata.obs that specifies visium library. Function trys to guess it if set to NULL (default).
#' @param simplify logical, whether return just seurat object (instead of list) if there is only one sample
#' @param use.raw logical, whether to attempt to get data from adata.raw. Throws warning (not exception!) if adata.raw is not present and proceeds with adata.X.
#'
#' @return list of Seurat object
#' @import Matrix
#' @export
h5ad2seurat_l_spatial = function (filename,library_id_field=NULL,simplify=TRUE,use.raw=FALSE){
  loadRequiredPackages('Seurat')
  data = h5ad2list(filename,load.raw = ifelse(use.raw,'replace','no'))

  a = rhdf5::H5Fopen(filename,flags = 'H5F_ACC_RDONLY')
  tryCatch({
    obs = data$obs
    var = data$var
    X = data$X

    # find column in obs all values of each has spatial data in a$uns$spatial
    if(is.null(library_id_field)){
      imglibs = names(a$uns$spatial)
      for(obsf in colnames(obs)){
        if(all(obs[,obsf] %in% imglibs)){
          library_id_field = obsf
          break
        }
      }
      # if there is not library field in obs we will only be able to proceed if the h5ad contains single samples:
      if(is.null(library_id_field)){
        if(length(imglibs) == 1){
          # generate filed name that is not already in use
          library_id_field = 'library_id'
          repeat{
            if(!(library_id_field %in% colnames(obs))) break
            library_id_field = paste0(library_id_field,'_')
          }
          obs[,library_id_field] = imglibs
        }else{
          stop("The h5ad seems to contain multiple visium samples but library is not specified correctly in adata.obs. There should be a column that matches keys is adata.uns.spatial.")
        }
      }
    }
    # load images and assemble Seurat objects
    # coordinates can be either in spatial or in X_spatial
    coor.names = 'spatial'
    if(!(coor.names %in% names(a$obsm)))
      coor.names = 'X_spatial'
    res = list()
    for(lid in unique(obs[[library_id_field]])){
      # subset obs
      f = obs[[library_id_field]] == lid
      obs_ = obs[f,]
      if(!is.null(obs_$barcode)){
        rownames(obs_) = obs_$barcode
      }
      # subset counts
      X_ = X[,which(f)]
      colnames(X_) = rownames(obs_)
      object = CreateSeuratObject(counts = X_, assay = 'Spatial')

      # extract image
      image = loadImageFromH5ad(a,lid,obs_,t(a$obsm[[coor.names]][,f])[,2:1])

      image = image[Cells(x = object)]
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

  if(simplify & length(res)==1)
    res = res[[1]]

  rhdf5::H5Fclose(a)
  return(res)
}


#' Loads Seurat VisiumV1 object from visium h5ad
#'
#' there is some magic in this function, I'm not sure why and how it works (but it works).
#' the function is not supposed to be used for its own, only from h5ad2seurat_l_spatial
#'
#' @param a H5IdComponent (for whole h5ad file)
#' @param lid library_id, the key in a.uns.spatia
#' @param obs_ subset of adata.obs for this particular object
#'
#' @return VisiumV1 object
loadImageFromH5ad = function(a,lid,obs_,obsm_){
  # extract image
  ires = 'hires'
  if(!(ires %in% names(a$uns$spatial[[lid]]$images)))
    ires = 'lowres'
  image <- aperm(a$uns$spatial[[lid]]$images[[ires]],3:1)

  scale.factors <- a$uns$spatial[[lid]]$scalefactors
  # looks like in scanpy both images are hires actually (at least they were identical for example I tried)
  scale.factors$tissue_lowres_scalef = scale.factors$tissue_hires_scalef

  tissue.positions = cbind(obs_[,c('in_tissue','array_row','array_col')], obsm_)
  colnames(tissue.positions) = c("tissue", "row", "col", "imagerow", "imagecol")

  rownames(tissue.positions) = rownames(obs_)

  unnormalized.radius = scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius = unnormalized.radius/max(dim(x = image))
  image = new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                                                              fiducial = scale.factors$fiducial_diameter_fullres,
                                                                              hires = scale.factors$tissue_hires_scalef,
                                                                              scale.factors$tissue_lowres_scalef),
              coordinates = tissue.positions, spot.radius = spot.radius)
  image
}
