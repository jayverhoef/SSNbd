#-------------------------------------------------------------------------------
#
#           mapi
#
#-------------------------------------------------------------------------------

#' make a partition index for a SpatialStreamNetwork object
#'
#' make a partition index for a SpatialStreamNetwork object so that 
#'   it can be used with fast methods for big data
#'
#' @param ssn.object a SpatialStreamNetwork object, usually imported 
#'   	with SSN package
#' @param partMeth method for subsampling.  Options are 'rand' for 
#'		random, 'comp' for compact blocks using kmeans on coordinates,
#'		and 'zimm' for compact blocks with 10% re-assigned at random.
#'		Default is 'zimm'.
#' @param nIndx the number of indexes (number of groups for partitioning)
#' @param set.seed.no set the random number seed so results are repeatable.
#'
#' @return an object of class "SpatialStreamNetwork", which is the same 
#'   as the input ssn.object, except the points data.frame within now
#'   has a new field named partIndx which indexes the groups for 
#'   partitioning. 
#'
#' @author Jay Ver Hoef
#' @export

mapi = function(ssn.object, nIndx, partMeth = 'comp',
	set.seed.no = 1001)
{
  xy = ssn.object@obspoints@SSNPoints[[1]]@point.coords
  xychar = paste0(xy[,1],xy[,2])
  unik = !duplicated(xychar)  
  xyunik = xy[seq_along(xychar)[unik],]
  indx = rep(c(1:nIndx), 
		times = ceiling(dim(xyunik)[1]/nIndx))[1:dim(xyunik)[1]]
	#set the random number seed to results are 
	set.seed(set.seed.no)
  if(partMeth == 'rand') {
		indx = indx[order(runif(length(indx)))]
	}
	if(partMeth == 'comp') {
		indx = kmeans(as.data.frame(xyunik), nIndx, iter.max = 500)$cluster
	}
	if(partMeth == 'zimm') {
		indx = kmeans(as.data.frame(xyunik), nIndx, iter.max = 500)$cluster
		#randomly re-assign 10% of each group at random
		gsamp = sample(1:length(indx), round(length(indx)/10))
		indx[gsamp] = indx[gsamp][order(runif(round(length(indx)/10)))]
	}
  xyindx = merge(data.frame(xychar = xychar), 
    data.frame(xychar = xychar[seq_along(xychar)[unik]], indx = indx),
    by.x = "xychar", all.x = TRUE, sort = FALSE)
  DF = getSSNdata.frame(ssn.object)
  DF$partIndx = xyindx$indx
  ssn.object = putSSNdata.frame(DF, ssn.object)
  ssn.object
}
