#-------------------------------------------------------------------------------
#
#           estCovParSSNbd
#
#-------------------------------------------------------------------------------

#' estimates covariance parameters with fast methods for big data with SSN
#'
#' estimates covariance parameters with fast methods for big data with SSN
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

estCovParSSNbd = function(formula, ssn.object, 
	CorModels = c("Exponential.tailup", "Exponential.taildown", 
  "Exponential.Euclid"), use.nugget = TRUE, addfunccol = NULL,
	EstMeth = "REML", subSampIndxCol)
{
    DF = getSSNdata.frame(ssn.object)
    ng = max(DF[,subSampIndxCol])
    ordi = order(DF[,subSampIndxCol], 
      as.integer(as.character(DF$netID)), DF$pid)
    DF = DF[ordi,]
    xy = ssn.object@obspoints@SSNPoints[[1]]@point.coords
    xy = xy[ordi,]
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
    offset = as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    X = model.matrix(mt, mf, contrasts)
    distLi = distList(ssn.object, DF, y, X, xy, CorModels = CorModels, 
			addfunccol = addfunccol, subSampIndxCol = subSampIndxCol, 
			distPath = ssn.object@path)
    #initial estimate of theta
    theta = NULL
    for(i in 1:ng) {
      yi = distLi[[i]]$y
      Xi = distLi[[i]]$X
      netDi = distLi[[i]]$netD
      Zsi = distLi[[i]]$Zs
      xci = distLi[[i]]$xc
      yci = distLi[[i]]$yc
	    theta = cbind(theta, SSN:::theta.ini(z = yi, X = Xi,
	      CorModels= CorModels,
	      use.nugget = TRUE, use.anisotropy = FALSE,
	      dist.hydro.data = netDi, x.dat = xci,
	      y.dat = yci, REs = Zsi))
    }
    theta = apply(theta,1,mean)
		attributes(theta) = attributes(SSN:::theta.ini(z = yi, X = Xi,
				CorModels= CorModels,
				use.nugget = TRUE, use.anisotropy = FALSE,
				dist.hydro.data = netDi, x.dat = xci,
				y.dat = yci, REs = Zsi))
		attributes(theta) ## scale, type, terms
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
      distanceList = distLi, xy = xy, mfcall = mfcall, DF = DF, 
      ssnr = ssn.object)
    class(outpt) <- "estCovParSSNbd"
    return(outpt)
}
