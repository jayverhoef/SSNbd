#-------------------------------------------------------------------------------
#
#           cope
#
#-------------------------------------------------------------------------------

#' covariance parameter estimates with fast methods for big data with SSN
#'
#' covariance parameter estimates with fast methods for big data with SSN
#'
#' @param formula an object of class "formula" (or one that can be 
#'   coerced to that class): a symbolic description of the model to be 
#    fitted. The details of model specification are given under 'Details'.
#' @param ssn.object an object of class
#'   \code{SpatialStreamNetwork}, 
#'   representing a spatial stream network. This contains the variables 
#'   used in the model.
#' @param CorModels a vector of spatial autocorrelation models for 
#'   stream networks.
#' @param use.nugget add a nugget effect, default is TRUE. This can be 
#'   thought of as a variance component for independent errors, adding 
#'   a variance component only along the diagonal of the covariance matrix.
#' @param addfunccol the name of the variable in the 
#'   SpatialStreamNetwork object that is used to define spatial 
#'   weights. For the tailup models, weights need to be used for 
#'   branching. This is an additive function and is described in 
#'   Ver Hoef and Peterson (2010). See example below.
#' @param EstMeth Estimation method; either "ML" for maximum 
#'   likelihood, or "REML" for restricted maximum likelihood (default).
#' @param subSampIndxCol the column in the \code{points data.frame}
#'   within the \code{SpatialStreamNetwork} object that indexes the
#'   grouping to be used when subsampling.
#'
#' @return an object of class "estCovParSSNbd", which is a list, where 
#'   \code{estCovPar} is a vector of estimated covariance parameters, 
#'   \code{optimOut} is the output from \code{optim} used to estimate
#'   the parameters. 
#'
#' @author Jay Ver Hoef
#' @export
#' @importFrom filematrix fm.open fm.load
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach registerDoSEQ 

cope = function(formula, ssn.object, 
	CorModels = c("Exponential.tailup", "Exponential.taildown", 
  "Exponential.Euclid"), use.nugget = TRUE, addfunccol = NULL,
	EstMeth = "REML", partIndxCol, parallel = FALSE)
{
    DF = getSSNdata.frame(ssn.object)
    ng = max(DF[,partIndxCol])
    ordi = order(DF[,partIndxCol], 
      as.integer(as.character(DF$netID)), DF$pid)
    DF = DF[ordi,]
    xy = ssn.object@obspoints@SSNPoints[[1]]@point.coords
    xy = xy[ordi,]
    REind <- which(names(DF) %in% CorModels)
    if(length(REind) > 0 & sum(apply(is.na(DF[REind]),2,any)) > 0) 
			stop("Missing values for random effects are not allowed.")
    cl = match.call()
    mf = match.call(expand.dots = FALSE)
    mfcall = mf
    m = match("ssn.object", names(mf), 0L)
    names(mf)[m] = "data"
    mf[[m]] = quote(DF)
    m = match(c("formula", "data"), names(mf), 0L)
    mf = mf[c(1L, m)]
    mf$drop.unused.levels = TRUE
    mf[[1L]] = quote(stats::model.frame)
    mf = eval(mf, DF)
    mt = attr(mf, "terms")
    y = model.response(mf, "numeric")
    ind = rownames(DF) %in% names(y)
    offset = as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
		#design matrix
		X <- model.matrix(mf, DF)
    #design matrices for random effects
    Zs <- NULL
      REind <- which(names(DF) %in% CorModels)
      if(length(REind)) {
        Zs <- list()
        REnames <- sort(names(DF)[REind])
        ## model matrix for a RE factor
        for(ii in 1:length(REind)) {
          DF[,REnames[ii]] = as.factor(as.character(DF[,REnames[ii]]))
          Zs[[REnames[ii]]] = model.matrix(~DF[,REnames[ii]] - 1)
          rownames(Zs[[REnames[ii]]]) = DF$pid
        }
        names(Zs) = REnames
      }
    #X = model.matrix(mt, mf, contrasts)
    distLi = distList(ssn.object, DF, y, X, Zs, xy, CorModels = CorModels, 
			addfunccol = addfunccol, subSampIndxCol = partIndxCol, 
			distPath = ssn.object@path, parallel = parallel)
    #initial estimate of theta
    theta = thetainibd(distLi, ng, CorModels)
		TH.scale <- attributes(theta)$scale
		TH.type <- attributes(theta)$type
		TH.terms <- attributes(theta)$terms

		n = length(y)
		p = length(X[1,])

		#estimate covariance parameters
		optimOut = optim(par = theta, fn = m2LLstrbd, distLi = distLi,
			CorModels = CorModels, use.nugget = use.nugget, 
			use.anisotropy = FALSE, useTailDownWeight = FALSE,
			EstMeth = EstMeth, n = n, p = p, scale = TH.scale, maxrang = NULL)

		estCovPar <- SSN:::untrans.theta(theta = optimOut$par, scale = TH.scale)

    outpt = list(estCovPar = estCovPar, optimOut = optimOut, 
      distanceList = distLi, xy = xy[ind,], mfcall = mfcall, DF = DF[ind,], ind = ind,
      ssnr = ssn.object)
    class(outpt) <- "estCovParSSNbd"
    
#    stopCluster(cl)
#		registerDoSEQ()
    return(outpt)
}
