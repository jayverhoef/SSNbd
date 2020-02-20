#-------------------------------------------------------------------------------
#
#           poba
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
poba <- function(ecp, efe, pul)
{
	d1 <- ecp$ssnr@obspoints@SSNPoints[[1]]@point.data
	coords <- ecp$ssnr@obspoints@SSNPoints[[1]]@point.coords
	theta <- ecp$estCovPar
	nobs <- dim(d1)[1]
	CorModels <- as.character(ecp$mfcall$CorModels)
	CorModels = CorModels[2:length(CorModels)]
	predsID = pul$predsID
  dp = ecp$ssnr@predpoints@SSNPoints[[
		which(ecp$ssnr@predpoints@ID == predsID)]]@point.data
	np = dim(dp)[1]
  coordsp = ecp$ssnr@predpoints@SSNPoints[[
		which(ecp$ssnr@predpoints@ID == predsID)]]@point.coords
	respCol = as.formula(ecp$mfcall[match("formula", names(ecp$mfcall), 0L)])[2]
	y = d1[,as.character(respCol)]

	v_star = pul$v_star[match(as.character(d1$pid),colnames(pul$yall))]
	Sig_pp_v = pul$Sig_pp_v
	
	Sig_op_v = rep(0, times = nobs)
	Sig_oo_vstar = rep(0, times = nobs)
	for(i in 1:nobs) {
		cat("\r", "Working on row", i, "out of", nobs)
			dmts = SSNbd:::dMatsEtc(ecp$ssn, CorModels, 'obs', d1[i,], coords[i,],
        ecp$mfcall$addfunccol, predsID, dp, coordsp)	
			Vopvec <- SSN:::makeCovMat(theta = theta, dist.hydro = dmts$dist.hydro,
					a.mat = dmts$a.mat, b.mat = dmts$b.mat, w.matrix = dmts$w.matrix,
					net.zero = dmts$net.zero, x.row = coords[i,1], y.row = coords[i,2],
					x.col = coordsp[,1], y.col = coordsp[,2],
					CorModels, useTailDownWeight = FALSE,
					use.nugget = FALSE,
					use.anisotropy = FALSE, dmts$REPs)	
			Sig_op_v[i] = mean(Vopvec) 

			dmts = SSNbd:::dMatsEtc(ecp$ssn, CorModels, 'obs', d1[i,], coords[i,],    
        ecp$mfcall$addfunccol, 'obs', d1, coords)	
			Voovec <- SSN:::makeCovMat(theta = theta, dist.hydro = dmts$dist.hydro,
					a.mat = dmts$a.mat, b.mat = dmts$b.mat, w.matrix = dmts$w.matrix,
					net.zero = dmts$net.zero, x.row = coords[i,1], y.row = coords[i,2],
					x.col = coords[,1], y.col = coords[,2],
					CorModels, useTailDownWeight = FALSE,
					use.nugget = FALSE,
					use.anisotropy = FALSE, dmts$REPs)	
			Sig_oo_vstar[i] = sum(Voovec*v_star)
				#ith term in the vector will have an additional nugget variance
				#v_star[i]*theta[attr(theta,'terms') == 'Nugget']
	}
	
		Yp_bar_hat = v_star %*% y
		Yp_bar_hat_se = sqrt(sum(Sig_oo_vstar*v_star) - 
			2*sum(Sig_op_v*v_star) + 
			mean(Sig_pp_v))

	list(block_pred = Yp_bar_hat, block_pred_se = Yp_bar_hat_se)
		
}
