#-------------------------------------------------------------------------------
#
#           estFixEffSSNbd
#
#-------------------------------------------------------------------------------

#' estimates fixed effects with fast methods for big data with SSN
#'
#' estimates fixed effects with fast methods for big data with SSN
#'
#' @param estCovParSSN an object of class "estCovParSSNbd" created
#'   from estimating the covariance parameters. 

#'
#' @return a list
#'
#' @author Jay Ver Hoef
#' @export

estFixEffSSNbd = function(estCovParSSN)
{    
			CorModels = as.list(estCovParSSN$mfcall[['CorModels']])
			CorModels = CorModels[2:length(CorModels)]
			CorModels = unlist(CorModels)
			nSubSamp = length(estCovParSSN$distanceList)
			Vlist = foreach(i=1:nSubSamp) %dopar% {
      V = SSN:::makeCovMat(theta = estCovParSSN$estCovPar, 
				dist.hydro = estCovParSSN$distanceList[[i]]$netD, 
        a.mat = estCovParSSN$distanceList[[i]]$A, 
        b.mat = estCovParSSN$distanceList[[i]]$B, 
        w.matrix = estCovParSSN$distanceList[[i]]$W, 
        net.zero = estCovParSSN$distanceList[[i]]$net0, 
        x.row = estCovParSSN$distanceList[[i]]$xc, 
        y.row = estCovParSSN$distanceList[[i]]$yc, 
        x.col = estCovParSSN$distanceList[[i]]$xc, 
        y.col= estCovParSSN$distanceList[[i]]$yc, 
        useTailDownWeight = FALSE,
        CorModels = CorModels, 
        use.nugget = estCovParSSN$mfcall[['use.nugget']], 
        use.anisotropy = FALSE, 
        REs = estCovParSSN$distanceList[[i]]$Zs)
        qrV = qr(V)
        list(V = V, qrV = qrV, ViX = solve(qrV,
          estCovParSSN$distanceList[[i]]$X), 
          Viy = solve(qrV,estCovParSSN$distanceList[[i]]$y),
          logdet = sum(log(abs(diag(qr.R(qrV))))),
          XViX = crossprod(estCovParSSN$distanceList[[i]]$X,
            solve(qrV,estCovParSSN$distanceList[[i]]$X)),
          XViy = crossprod(estCovParSSN$distanceList[[i]]$y,
            solve(qrV,estCovParSSN$distanceList[[i]]$X)),
          yViy = crossprod(estCovParSSN$distanceList[[i]]$y,
          solve(qrV,estCovParSSN$distanceList[[i]]$y)))
    }

    Sxx = Reduce('+',lapply(Vlist,function(x){x[['XViX']]}))
    sxy = t(Reduce('+',lapply(Vlist,function(x){x[['XViy']]})))
    syy = t(Reduce('+',lapply(Vlist,function(x){x[['yViy']]})))

    betaHat = solve(Sxx,sxy)
    
    CijList <- foreach(i = 1:(nSubSamp -1)) %:% 
		foreach(j = (i+1):nSubSamp) %dopar% {
			makeSigijMats(ssnr = estCovParSSN$ssnr, 
			DFr = estCovParSSN$DF, 
			xy = estCovParSSN$xy, 
			CorModels = CorModels, 
			theta = estCovParSSN$estCovPar, 
			addfunccol = estCovParSSN$mfcall[['addfunccol']],
			subSampIndxCol = estCovParSSN$mfcall[['subSampIndxCol']], i, j,
			distPath = eval(estCovParSSN$mfcall[['ssn.object']])@path, 
			useTailDownWeight = FALSE, Vlist = Vlist)
		}

		CijList = unlist(CijList, recursive = FALSE)
		Wxx = Reduce('+',lapply(CijList,function(x){x[['XViCViX']]}))

		Sxxi = solve(Sxx)
		covB = Sxxi + 2*Sxxi %*% Wxx %*% Sxxi

    list(betaHat = betaHat, covB = covB, Vlist = Vlist)
}
