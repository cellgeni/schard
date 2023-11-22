loadRequiredPackages = function(pkgs){
  for(pkg in pkgs){
    if(!require(pkg,character.only =TRUE))
      stop(paste0('Please install package "',pkg,'"'))
  }
}

#' Extracts data.frame from h5ad file
#'
#' @param filename file name or H5IdComponent to read data.frame from
#' @param name name of group that contains data.frame to be loaded
#'
#' @return a data.frame
#' @export
#' @examples
#' obs = h5ad2data.frame('adata.h5ad','obs')
h5ad2data.frame = function(filename,name){
  collist = rhdf5::h5read(filename,name,read.attributes = TRUE)
  attr = attributes(collist)
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
    codes = res[[fn]]+1
    codes[codes==0] = NA
    res[[fn]] = as.vector(collist[['__categories']][[fn]][codes])
  }
  # another way to store factors
  for(fn in names(ll)[ll==2 & ll != max(ll)]){
    if(all(names(collist[[fn]]) %in% c("categories","codes"))){
      codes = collist[[fn]]$codes+1
      codes[codes==0] = NA
      res[[fn]] = as.vector(collist[[fn]]$categories[codes])
    }
  }
  if('_index' %in% names(attr))
    rownames(res) = collist[[attr$`_index`]]
  if('index' %in% colnames(res))
    rownames(res) = res$index
  if(all(c('_index',attr$`column-order`) %in% colnames(res)))
    res = res[,c(attr$`_index`,attr$`column-order`),drop=FALSE]
  res
}




#' Parse matrix from h5ad file
#'
#' @param filename file name or H5IdComponent to read data.frame from
#' @param name name of group that contains matrix to be loaded
#'
#' @return R matrix, dense or sparse - in dependence on input
#' @export
#' @examples
#' obs = h5ad2data.frame('adata.h5ad','X')
h5ad2Matrix = function(filename,name){
  m = rhdf5::h5read(filename,name,read.attributes = TRUE)
  attr = attributes(m)
  format = attr$`encoding-type`
  if(is.null(format))
    format = attr$h5sparse_format
  shape = attr$shape
  if(is.null(shape))
    shape = attr$h5sparse_shape
  if(is.array(m)){
    mtx = m
  }else{
    if(startsWith(format,'csr')){
      mtx = Matrix::sparseMatrix(i=m$indices+1, p=m$indptr,x = as.numeric(m$data),dims = rev(shape))
    }else{
      mtx = Matrix::sparseMatrix(i=m$indices+1, p=m$indptr,x = as.numeric(m$data),dims = shape)
      mtx = Matrix::t(mtx)
    }
  }
  mtx
}



#' Load visium images from h5ad
#'
#'
#' @param filename path 2 visium h5ad
#'
#' @return list of slides. Each slide is list with scale.fctors, hires and lowres images
#'
#' @examples
#' imgs = h5ad2images('adata.h5ad')
h5ad2images = function(filename){
  h5struct = rhdf5::h5ls(filename)
  library_ids = h5struct$name[h5struct$group=='/uns/spatial']
  result = list()

  for(lid in library_ids){
    result[[lid]] = list()
    result[[lid]]$scale.factors = rhdf5::h5read(filename,paste0('/uns/spatial/',lid,'/scalefactors'))
    for(res in h5struct$name[h5struct$group==paste0('/uns/spatial/',lid,'/images')])
      result[[lid]][[res]] = aperm(rhdf5::h5read(filename,paste0('/uns/spatial/',lid,'/images/',res)),3:1)
  }
  result
}
