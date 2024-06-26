loadRequiredPackages = function(pkgs){
  for(pkg in pkgs){
    if(!require(pkg,character.only =TRUE))
      stop(paste0('In order to use this functionality please install package "',pkg,'"'))
  }
}

#' Extracts data.frame from h5ad file
#'
#' @param filename file name or H5IdComponent to read data.frame from
#' @param name name of group that contains data.frame to be loaded
#' @param keep.rownames.as.column whether to keep rownames (index) as column in output. TRUE by default
#'
#' @return a data.frame
#' @export
#' @examples
#' obs = h5ad2data.frame('adata.h5ad','obs')
h5ad2data.frame = function(filename,name,keep.rownames.as.column=TRUE){
  collist = rhdf5::h5read(filename,name,read.attributes = TRUE)
  attr = attributes(collist)

  # slashes in names may leads to nested structure, lets fix it
  repeat{
    ll = sapply(collist,length)
    # candidates
    inx_ = which(ll!=max(ll) & names(collist) != '__categories')
    # remove factors
    inx = c()
    for(i in inx_){
      if(length(collist[[i]]) != 2 | !all(names(collist[[i]]) %in% c("categories","codes")))
        inx = c(inx,i)
    }
    if(length(inx) == 0)
      break
    # fix first
    inx = inx[1]
    for(i in seq_len(length(collist[[inx]]))){
      n = paste0(names(collist)[inx],'/',names(collist[[inx]])[i])
      collist[[n]] = collist[[inx]][[i]]
    }
    collist[[inx]] = NULL
  }

  # sometime vector columns are stored as arrays, it'll cause problems latter on, lets coerce them to vectors
  for(i in which(sapply(collist,is.array)))
    collist[[i]] = as.vector(collist[[i]])

  # factors are stored as list of categories and names, other types are stored as vectors
    # first way to store factors (all levels are in collist[['__categories']])
  for(fn in names(collist[['__categories']])){
    codes = collist[[fn]]+1
    codes[codes==0] = NA
    collist[[fn]] = as.vector(collist[['__categories']][[fn]][codes])
  }
  collist[['__categories']] = NULL

  # another way to store factors (each factor is 2-item list)
  for(fn in names(ll)[ll==2]){
    if(all(names(collist[[fn]]) %in% c("categories","codes"))){
      codes = collist[[fn]]$codes+1
      codes[codes==0] = NA
      collist[[fn]] = as.vector(collist[[fn]]$categories[codes])
    }
  }
  ll = sapply(collist,length)
  if(any(ll != max(ll)))
    warning("unexpected data.frame format, some columns can be missed")
  res = as.data.frame(collist[ll==max(ll)],check.names=FALSE)

  index.col = NULL
  if('index' %in% colnames(res))
    index.col = 'index'
  if('_index' %in% names(attr) && attr$`_index` %in% colnames(res))
    index.col = attr$`_index`

  if(!is.null(index.col))
    rownames(res) = res[,index.col]

  ord =NULL
  if(!is.null(attr$`column-order`) & all(attr$`column-order` %in% colnames(res))){
    ord = attr$`column-order`
    if(keep.rownames.as.column && !is.null(index.col))
      ord = c(index.col,ord)
    res = res[,ord,drop=FALSE]
  }
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
  attr = rhdf5::h5readAttributes(filename,name)
  # load as dataframe if it appers to be a dataframe, but then convert to matrix as the matrix was requested
  if(!is.null(attr$`encoding-type`) && attr$`encoding-type` =='dataframe'){
    mtx = h5ad2data.frame(filename,name,keep.rownames.as.column=FALSE)
    mtx = as.matrix(mtx)
    return(mtx)
  }
  m = rhdf5::h5read(filename,name)
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
