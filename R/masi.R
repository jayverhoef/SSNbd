#-------------------------------------------------------------------------------
#
#           subSampIndxSSN
#
#-------------------------------------------------------------------------------

#' make a subsample index for a SpatialStreamNetwork object
#'
#' make a subsample index for a SpatialStreamNetwork object so that 
#'   it can be used with fast methods for big data
#'
#' @param ssn.object a SpatialStreamNetwork object, usually imported 
#'   with SSN package
#' @param nIndx the number of indexes (number of groups for subsampling)
#'
#' @return an object of class "SpatialStreamNetwork", which is the same 
#'   as the input ssn.object, except the points data.frame within now
#'   has a new field named subSampIndx which indexes the groups for 
#'   subsampling. 
#'
#' @author Jay Ver Hoef
#' @export

masi = function(ssn.object, nIndx)
{
  xy = ssn.object@obspoints@SSNPoints[[1]]@point.coords
  xychar = paste0(xy[,1],xy[,2])
  unik = !duplicated(xychar)  
  xyunik = xy[seq_along(xychar)[unik],]
  indx = rep(c(1:nIndx), times = ceiling(dim(xyunik)[1]/nIndx))[1:dim(xyunik)[1]]
  indx = indx[order(runif(length(indx)))]
  xyindx = merge(data.frame(xychar = xychar), 
    data.frame(xychar = xychar[seq_along(xychar)[unik]], indx = indx),
    by.x = "xychar", all.x = TRUE, sort = FALSE)
  DF = getSSNdata.frame(ssn.object)
  DF$subSampIndx = xyindx$indx
  ssn.object = putSSNdata.frame(DF, ssn.object)
  ssn.object
}
