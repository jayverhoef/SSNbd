#-------------------------------------------------------------------------------
#
#           pure
#
#-------------------------------------------------------------------------------

#' prediction of unobserved responses with fast methods for big data with SSN
#'
#' Prediction of unobserved responses with fast methods for big data with SSN
#'
#' @param ecp object from cope function.
#' @param efe object from fefe function
#' @param predsID name of prediction data set in ssn object (passed with ecp)
#' @param nNN number of nearest neighbors for predictions 
#'
#' @return a data.frame with predictions in first column, and prediction standard errors in the second column.
#'
#' @author Jay Ver Hoef
#' @export
#' @importFrom nabor knn
pure <- function(ecp, efe, predsID, nNN)
{
	d1 <- ecp$ssnr@obspoints@SSNPoints[[1]]@point.data
	coords <- ecp$ssnr@obspoints@SSNPoints[[1]]@point.coords
	theta <- ecp$estCovPar
	nobs <- dim(d1)[1]
	CorModels <- as.character(ecp$mfcall$CorModels)
	CorModels = CorModels[2:length(CorModels)]
  dp = ecp$ssnr@predpoints@SSNPoints[[
		which(ecp$ssnr@predpoints@ID == predsID)]]@point.data
	np = dim(dp)[1]
  coordsp = ecp$ssnr@predpoints@SSNPoints[[
		which(ecp$ssnr@predpoints@ID == predsID)]]@point.coords
	nearxy = knn(data = coords, query = coordsp, k = nNN)
	storePreds = matrix(NA, nrow = np, ncol = 2)
		for(i in 1:np) {
		DF1 = d1[nearxy$nn.idx[i,],]
		xy1 = coords[nearxy$nn.idx[i,],]
		distord <- order(as.integer(as.character(DF1[,"netID"])), DF1[,"pid"])
		DF1 = DF1[distord,]
		xy1 = xy1[distord,]
		DF2 = dp[i,]
		xy2 = coordsp[i,, drop = FALSE]
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
		qrV = qr(V)
		covb <- efe$covB
		z <- DF1[, response.col]
		n <- length(z)
		p <- dim(Xp)[2]
    Viz = solve(qrV, z)
    Vic = solve(qrV, Vpred)
    
		pred.out <- matrix(NA, nrow = 1, ncol = 2)
		pred.out[1,1] = sum(as.vector((Viz)) * Vpred) +
			Xp %*% efe$betaHat - t(Vic) %*% Xobs %*% efe$betaHat
			
		pred.out[1,2] = sqrt(sumparsil - sum(Vic * Vpred) +
			sum((covb %*% t(Xp))*t(Xp)) - 2*sum((covb %*% t(Xp))*(t(Xobs) %*% Vic)) +
			sum((covb %*% t(Xobs) %*% Vic)*(t(Xobs) %*% Vic)))		
		storePreds[i,] = pred.out	
	}
	storePreds
}

