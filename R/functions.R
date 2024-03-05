#' Loads main parts of h5ad file as ordinary R list
#'
#' @param filename path to h5ad file
#' @param load.raw logical, whether to load adata.raw.X instead of adata.X
#' @param load.obsm logical, whether to load adata.obsm.
#' @param load.X logical, whether to load expression. Set all expressions to zeroes if FALSE. It can be much faster not to load expression
#'
#' @return list with following elements:
#'  obs - data.frame
#'  var - data.frame
#'  X - matrix (dense or sparse - as it was in X)
#'  obsm - list of matrices from obsm (presumably reduced dimensions)
#' @export
h5ad2list = function(filename,use.raw=FALSE,load.obsm=FALSE,load.X = TRUE){
  h5struct = rhdf5::h5ls(filename)
  res = list()
  res$obs = h5ad2data.frame(filename,'obs')

  if(use.raw){
    if(!any(h5struct$group=='/raw/X')){
      stop(paste0('There is not raw slot in "',filename,'" please set "use.raw" to FALSE'))
    }
    res$var = h5ad2data.frame(filename,'raw/var')
    expname = 'raw/X'
  }else{
    res$var = h5ad2data.frame(filename,'var')
    expname = 'X'
  }
  if(load.X){
    res$X = h5ad2Matrix(filename,expname)
  }else{
    res$X = Matrix::sparseMatrix(i=integer(0),p=0L,x=numeric(0),dims=c(nrow(res$var),nrow(res$obs)))
  }
  rownames(res$X) = rownames(res$var)
  colnames(res$X) = rownames(res$obs)

  res$obsm = list()
  if(load.obsm){
    for(n in h5struct$name[h5struct$group=='/obsm']){
      obsm = h5ad2Matrix(filename,paste0('obsm/',n)) # they are transposed for some reason...
      if(is.array(obsm)){
        if(nrow(obsm) != nrow(res$obs))
          obsm = t(obsm)
        res$obsm[[n]] = obsm
      }
    }
  }
  res
}


#' Loads h5ad as SingleCellExperiment object
#'
#' functions exports X,obs,var and obsm
#' it is very experimental and comes with no warranty
#'
#' @param filename path to h5ad file
#' @param use.raw logical, whether to attempt to get data from adata.raw. Throws warning (not exception!) if adata.raw is not present and proceeds with adata.X.
#' @param load.obsm logical, whether to load adata.obsm. All matrix-like objects from there will be added as reducedDim to the output object.
#' @param load.X logical, whether to load expression. Set all expressions to zeroes if FALSE. It can be much faster not to load expression
#'
#' @return SingleCellExperiment object
#' @export
#'
#' @import Matrix
#' @examples
#' sce = h5ad2sce('adata.h5ad')
h5ad2sce = function(filename,use.raw=FALSE,load.obsm=TRUE,load.X=TRUE){
  loadRequiredPackages('SingleCellExperiment')
  data = h5ad2list(filename,use.raw = use.raw,load.obsm=load.obsm,load.X=load.X)

  sce = SingleCellExperiment(list(X=data$X),
                             colData=data$obs,
                             rowData=data$var
  )
  # reduced dims
  for(n in names(data$obsm)){
    reducedDim(sce,n) = data$obsm[[n]]
  }
  sce
}

#' Loads h5ad as Seurat object
#'
#' For non-spatial data only, use h5ad2seurat_spatial for Visium
#'
#' @param filename path to h5ad file
#' @param use.raw logical, whether to use adata.raw instead of adata.X
#' @param load.obsm logical, whether to load adata.obsm. All matrix-like objects from there will be added as DimReduc to the output object with names coerced to Seurat style (that is not underscores in the middle, single underscore at the end).
#' @param assay what assay to put data it (RNA by default)
#' @param load.X logical, whether to load expression. Set all expressions to zeroes if FALSE. It can be much faster not to load expression
#'
#' @return Seurat object
#' @export
#'
#' @import Matrix
#' @examples
#' seu = h5ad2seurat('adata.h5ad')
h5ad2seurat = function(filename,use.raw=FALSE,load.obsm=TRUE,assay='RNA',load.X=TRUE){
  loadRequiredPackages('Seurat')
  data = h5ad2list(filename,use.raw = use.raw,load.obsm = load.obsm,load.X=load.X)

  # guess whether data is scaled
  scaled = FALSE
  # scaled data cannot be sparse (hopefully!)
  if(is.array(data$X)){
    n0 = sum(data$X!=0,na.rm=TRUE)
    density = n0/sum(!is.na(data$X))
    # if it dense or has no chance to fit into memory
    scaled = density > 0.8 | n0 > 2^31-1
  }
  # if scaled then supply it as data to prevent Seurat making it sparse
  if(scaled){
    counts = CreateAssayObject(data = data$X)
  }else{
    counts = CreateAssayObject(counts = data$X)
  }

  seu = CreateSeuratObject(counts = counts,assay = assay)
  seu = AddMetaData(seu,data$obs)
  seu[[assay]] = AddMetaData(seu[[assay]],metadata = data$var)

  # reduced dims
  # move spatial to X_spatial (that is Xspatial_ in Seurat) to a) do not interact with Seurat[['spatial_']] b) consistency between adata versions
  names(data$obsm)[names(data$obsm) == 'spatial'] = 'X_spatial'

  for(n in names(data$obsm)){
    if(ncol(data$obsm[[n]]) == 0) next
    nn = paste0(gsub('_','',n),'_') # Seurat naming requirements
    colnames(data$obsm[[n]]) = paste0(nn,seq_len((data$obsm[[n]]))) # Seurat wants to have colnames
    rownames(data$obsm[[n]]) = rownames(data$obs)
    seu[[nn]] <- CreateDimReducObject(embeddings = data$obsm[[n]], key = nn, assay = assay)
  }
  seu
}


#' Loads scanpy h5ad file with Visium data into Seurat object
#'
#'
#' @param filename character, path to h5ad file
#' @param use.raw logical, whether to attempt to get data from adata.raw. Throws warning (not exception!) if adata.raw is not present and proceeds with adata.X.
#' @param load.obsm logical, whether to load adata.obsm. All matrix-like objects from there will be added as DimReduc to the output object with names coerced to Seurat style (that is not underscores in the middle, single underscore at the end).
#' @param simplify logical, whether to merge loaded samples into single object, return list otherwise.
#' @param img.res which of lowres' or 'hires' image to be loaded
#' @param load.X logical, whether to load expression. Set all expressions to zeroes if FALSE. It can be much faster not to load expression
#'
#' @return list of Seurat object
#' @import Matrix
#' @export
#' @examples
#' vs = h5ad2seurat_spatial('adata.h5ad')
h5ad2seurat_spatial = function(filename,use.raw=FALSE,load.obsm=TRUE,simplify=TRUE,img.res = 'lowres',load.X=TRUE){
  loadRequiredPackages('Seurat')
  data = h5ad2list(filename,use.raw = use.raw,load.obsm = TRUE,load.X=load.X) # load obsm in any case - we need spatial info
  images = h5ad2images(filename)
  results = list()

  # find column all values of each has images
  library_id_field = NULL
  for(col in colnames(data$obs)){
    if(all(data$obs[,col] %in% names(images))){
      library_id_field = col
      break
    }
  }
  # if there is not library field in obs we will only be able to proceed if the h5ad contains single samples:
  if(is.null(library_id_field)){
    if(length(images) == 1){
      # generate filed name that is not already in use
      library_id_field = 'library_id'
      repeat{
        if(!(library_id_field %in% colnames(data$obs))) break
        library_id_field = paste0(library_id_field,'_')
      }
      data$obs[,library_id_field] = names(images)
    }else{
      stop("The h5ad seems to contain multiple visium samples but library is not specified correctly in adata.obs. There should be a column that matches keys is adata.uns.spatial.")
    }
  }
  # rename spatial to X_spatial
  names(data$obsm)[names(data$obsm)=='spatial'] = 'X_spatial'

  # create per sample Seurats
  for(lid in unique(data$obs[,library_id_field])){
    f = data$obs[,library_id_field] == lid
    obs_ = data$obs[f,]
    x_spatial_ = data$obsm$X_spatial[f,2:1] # manually discovered that coordinates are transposed

    seu = CreateSeuratObject(counts = data$X[,rownames(obs_)], assay = 'Spatial')
    seu = AddMetaData(seu,obs_)
    seu[['Spatial']] = AddMetaData(seu[['Spatial']],data$var)

    # add dim reductions (if any)
    if(load.obsm){
      for(n in setdiff(names(data$obsm),'X_spatial')){
        nn = paste0(gsub('_','',n),'_') # Seurat naming requirements
        obsm_ = data$obsm[[n]][f,,drop=FALSE]
        if(ncol(obsm_) == 0) next
        colnames(obsm_) = paste0(rep(nn,ncol(obsm_)),seq_len(ncol(obsm_))) # Seurat wants to have colnames
        rownames(obsm_) = rownames(obs_)
        seu[[nn]] = CreateDimReducObject(embeddings = obsm_, key = nn, assay = 'Spatial')
      }
    }

    # prepare image
    # set to mock if mesh coordinaates are not in obs
    tissue.positions = obs_
    if(is.null(tissue.positions$in_tissue)) tissue.positions$in_tissue = 1
    if(is.null(tissue.positions$array_row)) tissue.positions$array_row = 0
    if(is.null(tissue.positions$array_col)) tissue.positions$array_col = 0

    tissue.positions = cbind(tissue.positions[,c('in_tissue','array_row','array_col')], x_spatial_)
    colnames(tissue.positions) = c("tissue", "row", "col", "imagerow", "imagecol")
    rownames(tissue.positions) = rownames(obs_)

    scale.factors = images[[lid]]$scale.factors
    # if specified resolution is not in h5ad, try another
    if(!(img.res %in% names(images[[lid]]))){
      warning(paste0("'",img.res,"' is not avaliable, trying another resolution",))
      img.res = setdiff(c('hires','lowres'),img.res)
    }
    if(!(img.res %in% names(images[[lid]]))){
      stop('No image avaliable in h5ad')
    }

    if(img.res == 'hires'){
      scale.factors$tissue_lowres_scalef = scale.factors$tissue_hires_scalef
    }

    image = images[[lid]][[img.res]]

    # just copied from Seurat::Read10X_Image
    unnormalized.radius = scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
    spot.radius = unnormalized.radius/max(dim(x = image))
    image = new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                                                                fiducial = scale.factors$fiducial_diameter_fullres,
                                                                                hires = scale.factors$tissue_hires_scalef,
                                                                                scale.factors$tissue_lowres_scalef),
                coordinates = tissue.positions, spot.radius = spot.radius)

    image = image[Cells(x = seu)]
    DefaultAssay(object = image) = 'Spatial'
    seu[[lid]] = image
    results[[lid]] = seu
  }
  if(simplify){
    if(length(results)==1)
      results = results[[1]]
    else{
      results = merge(results[[1]],results[-1])
    }
  }
  results
}

