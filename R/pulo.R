#-------------------------------------------------------------------------------
#
#           pulo
#
#-------------------------------------------------------------------------------

#' prediction at unobserved locations with fast methods for big data with SSN
#'
#' Prediction atunobserved locations with fast methods for big data with SSN
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

pulo <- function(ecp, efe, predsID, nNN, poba_prep = FALSE)
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
	nearxy = nabor:::knn(data = coords, query = coordsp, k = nNN)
	
  Vlist = efe$Vlist
  Sxx = Reduce('+',lapply(Vlist,function(x){x[['XViX']]}))
  ViXlist = lapply(Vlist,function(x){t(x[['ViX']])})
  ylist = lapply(Vlist,function(x){t(x[['y']])})
  Sxxi = solve(Sxx)
  QQ = NULL
  yall = NULL
  for(i in 1:length(ViXlist)) {
    QQ = cbind(QQ, Sxxi %*% ViXlist[[i]])
    yall = cbind(yall, ylist[[i]])
  }

	v_star = rep(0, times = nobs)
	Sig_pp_v = rep(0, times = np)
	predDF = matrix(NA, nrow = np, ncol = 3)
	for(i in 1:np) {
		cat("\r", "Working on prediction", i, "out of", np)
		DF1 = d1[nearxy$nn.idx[i,],]
		xy1 = coords[nearxy$nn.idx[i,],]
		distord <- order(as.integer(as.character(DF1[,"netID"])), DF1[,"pid"])
		DF1 = DF1[distord,]
		xy1 = xy1[distord,]
		DF2 = dp[i,]
		xy2 = coordsp[i,, drop = FALSE]
		#rbind Obs and Preds data.frames, take care about factors!
		Fcols = names(which(unlist(lapply(DF1[1,], is.factor))))
		if(length(Fcols)>0)
			for(j in 1:length(Fcols)){
				DF1[,Fcols[j]] = as.character(DF1[,Fcols[j]])
				DF2[,Fcols[j]] = as.character(DF2[,Fcols[j]])
			}
		DF12 = rbind(DF1, rep(NA, times = length(DF1[1,])))
		DF12[length(DF12[,1]), match(colnames(DF2), colnames(DF1))] = DF2
		rownames(DF12) = c(rownames(DF1), rownames(DF2))
		if(length(Fcols)>0)
			for(j in 1:length(Fcols))
				DF12[,Fcols[j]] = as.factor(DF12[,Fcols[j]])
		X12 = model.matrix(as.formula(ecp$mfcall$formula[c(1,3)]), DF12)
		ind = rownames(X12) %in% rownames(DF12)
		ind1 = ind[1:(length(ind)-1)]
		ind2 = ind[length(ind)]
		DF1 = DF1[ind1,]
		DF2 = DF2[ind2,]
		Xobs = X12[1:(length(X12[,1]) - 1),]
		Xp = X12[length(X12[,1]),, drop = FALSE]
		xy1 = xy1[ind1,, drop = FALSE]
		xy2 = xy2[ind2,, drop = FALSE]
	
	  if(ind2) {
			mf = ecp$mfcall
			m = match("ssn.object", names(mf), 0L)
			names(mf)[m] = "data"
			mf[[m]] = quote(DF1)
			m = match(c("formula", "data"), names(mf), 0L)
			mf = mf[c(1L, m)]
			mf$drop.unused.levels = TRUE
			mf[[1L]] = quote(stats::model.frame)
			mf = eval(mf, DF1)
			mt = attr(mf, "terms")
			y = model.response(mf, "numeric")
			ind = rownames(DF1) %in% names(y)

			dmts = SSNbd:::dMatsEtc(ecp$ssn, CorModels, 'Obs', DF1, xy1, 
        ecp$mfcall$addfunccol)
      if(is.null(ecp$mfcall$use.nugget)) {use.nugget = TRUE} else {
				use.nugget = ecp$mfcall$use.nugget}
			V <- SSN:::makeCovMat(theta = theta, dist.hydro = dmts$dist.hydro,
					a.mat = dmts$a.mat, b.mat = dmts$b.mat, w.matrix = dmts$w.matrix,
					net.zero = dmts$net.zero, x.row = xy1[,1], y.row = xy1[,2],
					x.col = xy1[,1], y.col = xy1[,2],
					CorModels, useTailDownWeight = FALSE,
					use.nugget = use.nugget,
					use.anisotropy = FALSE, dmts$REs)
			dmts = SSNbd:::dMatsEtc(ecp$ssn, CorModels, 'Obs', DF1, xy1, 
        ecp$mfcall$addfunccol, predsID, DF2, xy2)	
			Vpred <- SSN:::makeCovMat(theta = theta, dist.hydro = dmts$dist.hydro,
					a.mat = dmts$a.mat, b.mat = dmts$b.mat, w.matrix = dmts$w.matrix,
					net.zero = dmts$net.zero, x.row = xy1[,1], y.row = xy1[,2],
					x.col = xy2[,1], y.col = xy2[,2],
					CorModels, useTailDownWeight = FALSE,
					use.nugget = FALSE,
					use.anisotropy = FALSE, dmts$REPs)	
			
			# get the sum of partial sills
			sumparsil <- sum(theta[attr(theta,"type") == "parsill"])
			Vi = solve(V)
			qrV = qr(V)
			covb <- efe$covB
			z <- y
			n <- length(z)
			p <- dim(Xp)[2]
			Viz = solve(qrV, z)
			Vic = solve(qrV, Vpred)    
      
			if(poba_prep == TRUE) {
        dmts = SSNbd:::dMatsEtc(ecp$ssn, CorModels, predsID, dp, coordsp, 
          ecp$mfcall$addfunccol, predsID, DF2, xy2)	
        Vpvec <- SSN:::makeCovMat(theta = theta, dist.hydro = dmts$dist.hydro,
            a.mat = dmts$a.mat, b.mat = dmts$b.mat, w.matrix = dmts$w.matrix,
            net.zero = dmts$net.zero, x.row = coordsp[,1], y.row = coordsp[,2],
            x.col = xy2[,1], y.col = xy2[,2],
            CorModels, useTailDownWeight = FALSE,
            use.nugget = FALSE,
            use.anisotropy = FALSE, dmts$REPs)	
        Sig_pp_v[i] = mean(Vpvec)  
          #one term in the vector will have a nugget variance
          #theta[attr(theta,'terms') == 'Nugget']/np		
    
        cViN = rep(0, times = length(yall))
        cViN[match(names(z),colnames(yall))] = t(Vic)
        lambdai = Xp %*% QQ + cViN - t(Vic) %*% Xobs %*% QQ
        v_star = v_star + lambdai
      }
      
			predDF[i,] = c(
				pid = dp$pid[i], 
				sum(as.vector((Viz)) * Vpred) +
					Xp %*% efe$betaHat - t(Vic) %*% Xobs %*% efe$betaHat,
				sqrt(sumparsil - sum(Vic * Vpred) +
					sum((covb %*% t(Xp))*t(Xp)) - 2*sum((covb %*% t(Xp))*(t(Xobs) %*% Vic)) +
					sum((covb %*% t(Xobs) %*% Vic)*(t(Xobs) %*% Vic))) 
				)	
		}	
	}
	colnames(predDF) = c('pid', 'pred', 'predse')
	if(poba_prep == FALSE) outpt = list(predDF = as.data.frame(predDF))
	if(poba_prep == TRUE)	outpt = list(predDF = as.data.frame(predDF), yall= yall,
		v_star = v_star/np, Sig_pp_v = Sig_pp_v, predsID = predsID)
		
	outpt
}

