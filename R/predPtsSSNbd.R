#-------------------------------------------------------------------------------
#
#           predSSNbdNN
#
#-------------------------------------------------------------------------------

#' predictions with fast methods for big data with SSN
#'
#' Predictions fixed effects with fast methods for big data with SSN
#'
#' @param ecp object from estCovParSSNbd function.
#' @param efe object from estFixEffSSNbd function
#' @param predsID name of prediction data set in ssn object (passed with ecp)
#' @param nNN number of nearest neighbors for predictions 
#'
#' @return a data.frame with predictions in first column, and prediction standard errors in the second column.
#'
#' @author Jay Ver Hoef
#' @export
#' @importFrom nabor knn
predSSNbdNN <- function(ecp, efe, predsID, nNN)
{
	d1 <- ecp$ssnr@obspoints@SSNPoints[[1]]@point.data
	coords <- ecp$ssnr@obspoints@SSNPoints[[1]]@point.coords
	theta <- ecp$estCovPar
	nobs <- dim(d1)[1]
	CorModels <- as.character(ecp$mfcall$CorModels)
	CorModels = CorModels[2:length(CorModels)]
#	useTailDownWeight <-object$args$useTailDownWeight
#	  REind <- which(names(data) %in% CorModels)
	REs <- names(d1)[which(names(d1) %in% CorModels)]
  dp = ecp$ssnr@predpoints@SSNPoints[[
		which(ecp$ssnr@predpoints@ID == predsID)]]@point.data
	np = dim(dp)[1]
  coordsp = ecp$ssnr@predpoints@SSNPoints[[
		which(ecp$ssnr@predpoints@ID == predsID)]]@point.coords
	nearxy = knn(data = coords, query = coordsp, k = nNN)
	storePreds = matrix(NA, nrow = np, ncol = 2)
	i = 2
		for(i in 1:np) {
		DF1 = d1[nearxy$nn.idx[i,],]
		xy1 = coords[nearxy$nn.idx[i,],]
		distord <- order(as.integer(as.character(DF1[,"netID"])), DF1[,"pid"])
		DF1 = DF1[distord,]
		xy1 = xy1[distord,]
		DF2 = dp[i,]
		xy2 = coordsp[i,, drop = FALSE]
		#undebug(SSNbd:::dMatsEtc)
		dmts = dMatsEtc(ecp$ssn, CorModels, 'Obs', DF1, xy1, ecp$mfcall$addfunccol)
		V <- SSN:::makeCovMat(theta = theta, dist.hydro = dmts$dist.hydro,
				a.mat = dmts$a.mat, b.mat = dmts$b.mat, w.matrix = dmts$w.matrix,
				net.zero = dmts$net.zero, x.row = xy1[,1], y.row = xy1[,2],
				x.col = xy1[,1], y.col = xy1[,2],
				CorModels, useTailDownWeight = FALSE,
				use.nugget = ecp$mfcall$use.nugget,
				use.anisotropy = FALSE, dmts$REs)
		dmts = dMatsEtc(ecp$ssn, CorModels, 'Obs', DF1, xy1, ecp$mfcall$addfunccol, 
			predsID, DF2, xy2)	
		Vpred <- SSN:::makeCovMat(theta = theta, dist.hydro = dmts$dist.hydro,
				a.mat = dmts$a.mat, b.mat = dmts$b.mat, w.matrix = dmts$w.matrix,
				net.zero = dmts$net.zero, x.row = xy1[,1], y.row = xy1[,2],
				x.col = xy2[,1], y.col = xy2[,2],
				CorModels, useTailDownWeight = FALSE,
				use.nugget = FALSE,
				use.anisotropy = FALSE, dmts$REPs)	
		
		# get a list of response and covariate names
		response.col <- as.character(ecp$mfcall$formula[[2]])
		DF2[,response.col] <- -1
		formula <- ecp$mfcall$formula
		# create design matrix for prediction data set
		Xobs <- model.matrix(as.formula(ecp$mfcall$formula), DF1)
		Xp <- model.matrix(as.formula(ecp$mfcall$formula), DF2)

		# get the sum of partial sills
		sumparsil <- sum(theta[attr(theta,"type") == "parsill"])
		Vi = solve(V)
		covb <- efe$covB
		z <- DF1[, response.col]
		n <- length(z)
		p <- dim(Xp)[2]

		parsilvec <- rep(sumparsil, times = length(Vpred[1,]))


			M <- rbind(Vpred, t(Xp), parsilvec)
			XXSiXi <- Xobs %*% covb
			XSi <- t(Xobs) %*% Vi
			pred.out <- t(apply(M, 2, SSN:::UK4Apply, covb = covb,
				XXSiXi = XXSiXi, XSi = XSi, Vi = Vi, z = z, n = n, p = p))
			storePreds[i,] = pred.out
		}
		storePreds
}

