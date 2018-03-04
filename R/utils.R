distUpdater = function(DFr, y, X, xy, CorModels, addfunccol, subSampIndxCol, i,
	distPath, useTailDownWeight = FALSE)
{    
    DFi = DFr[DFr[subSampIndxCol] == i,] 
    yi = y[DFr[subSampIndxCol] == i]
    Xi = X[DFr[subSampIndxCol] == i,]
    nIDs <- sort(as.integer(as.character(unique(DFi[,"netID"]))))
    netDJ <- matrix(0, nrow = length(DFi[,1]), ncol = length(DFi[,1]))
    net0 <-  matrix(0, nrow = length(DFi[,1]), ncol = length(DFi[,1]))
    rn = NULL
    nsofar <- 0
    distord <- order(as.integer(as.character(DFi[,"netID"])),
      DFi[,"pid"])
    DFi = DFi[distord,]
    names(distord) <- rownames(DFi)[distord]
    for(k in nIDs) {
      workspace.name <- paste("dist.net", k, ".RData", sep = "")
      path <- file.path(distPath, "distance", "obs",
		    workspace.name)
	    if(!file.exists(path)) {
		    stop("Unable to locate required distance matrix")
	    }
	    file_handle <- file(path, open="rb")
	    distmat <- unserialize(file_handle)
	    ordpi <- order(as.numeric(rownames(distmat)))
	    close(file_handle)
      distmatk = distmat[rownames(distmat) %in% DFi$pid, 
        rownames(distmat) %in% DFi$pid, drop = FALSE]
	    nk <- dim(distmatk)[1]
	    netDJ[(nsofar + 1):(nsofar + nk),(nsofar + 1):
		    (nsofar + nk)] <- distmatk
	    net0[(nsofar + 1):(nsofar + nk),(nsofar + 1):(nsofar + nk)] <- 1
      rn = c(rn,rownames(distmatk))
	    nsofar <- nsofar + nk
    }
    xc = xy[DFr[subSampIndxCol] == i,1]
    yc = xy[DFr[subSampIndxCol] == i,1]
    rownames(netDJ) = rn
    colnames(netDJ) = rn
    rownames(net0) = rn
    colnames(net0) = rn
    A = pmax(netDJ,t(netDJ))
    B = pmin(netDJ,t(netDJ))
    netD = as.matrix(netDJ + t(netDJ))
    W = NULL
	  if(length(grep("tailup",CorModels)) | useTailDownWeight == TRUE) {
		  if(missing(addfunccol) || is.null(addfunccol) ||
			  length(addfunccol) == 0 || !(addfunccol %in% colnames(DFr)))
	      	stop("The addfunccol argument is missing or mis-specified")
        FCmat <- 1 - (B > 0)*1
	      W <- sqrt(
          pmin(
            outer(DFi[,addfunccol],
	              rep(1, times = dim(A)[1])), 
		          t(outer(DFi[, addfunccol],
                rep(1, times = dim(A)[1]))) ) /
		      pmax(
            outer(DFi[, addfunccol],
                rep(1, times = dim(A)[1])),
			        t(outer(DFi[, addfunccol],
                rep(1, times = dim(A)[1])))))*
			        FCmat*net0
    }
    Zs <- NULL
      REind <- which(names(DFi) %in% CorModels)
      if(length(REind)) {
        Zs <- list()
        REnames <- sort(names(DF)[REind])
        ## model matrix for a RE factor
        for(ii in 1:length(REind)) {
          DFi[,REnames[ii]] = as.factor(as.character(DFi[,REnames[ii]]))
          Zs[[REnames[ii]]] = model.matrix(~DFi[,REnames[ii]] - 1)
          Zs[[REnames[ii]]] = tcrossprod(Zs[[REnames[ii]]])
          rownames(Zs[[REnames[ii]]]) = DFi$pid
          colnames(Zs[[REnames[ii]]]) = DFi$pid
        }
      }

    list(y = yi, X = Xi, netDJ = netDJ, net0 = net0, xc = xc, yc = yc, 
      A = A, B = B, netD = netD, W = W, FCmat = FCmat, Zs = Zs)
}

distList = function(DFr, y, X, xy, CorModels, addfunccol, subSampIndxCol,
	distPath)
{
  nResamp = max(DFr[,subSampIndxCol])
  Dlists = foreach(i=1:nResamp) %dopar% {
    distUpdater(DFr = DFr, y = y, X = X, 
    xy = xy, CorModels = CorModels,  addfunccol = addfunccol, 
    subSampIndxCol = subSampIndxCol, i = i, distPath = distPath)
  }
 Dlists
}

m2LLstrbd <- function(theta, distLi,
	CorModels, use.nugget, use.anisotropy, useTailDownWeight,
	EstMeth = "REML", n, p, scale, maxrang = NULL)
{
	  theta1 <- SSN:::untrans.theta(theta = theta, scale = scale)
    nResamp = length(distLi)
    Vlists = foreach(i=1:nResamp) %dopar% {
      V = SSN:::makeCovMat(theta = theta1, dist.hydro = distLi[[i]]$netD, 
        a.mat = distLi[[i]]$A, b.mat = distLi[[i]]$B, 
        w.matrix = distLi[[i]]$W, net.zero = distLi[[i]]$net0, 
        x.row = distLi[[i]]$xc, y.row = distLi[[i]]$yc, 
        x.col = distLi[[i]]$xc, y.col= distLi[[i]]$yc, 
        useTailDownWeight = FALSE,
        CorModels = CorModels, 
        use.nugget = TRUE, use.anisotropy = FALSE, 
        REs = distLi[[i]]$Zs)
        qrV = qr(V)
        list(V = V, qrV = qrV, ViX = solve(qrV,distLi[[i]]$X), 
          Viy = solve(qrV,distLi[[i]]$y),
          logdet = sum(log(abs(diag(qr.R(qrV))))),
          XViX = crossprod(distLi[[i]]$X,solve(qrV,distLi[[i]]$X)),
          XViy = crossprod(distLi[[i]]$y,solve(qrV,distLi[[i]]$X)),
          yViy = crossprod(distLi[[i]]$y,solve(qrV,distLi[[i]]$y)))
    }

    Sxx = Reduce('+',lapply(Vlists,function(x){x[['XViX']]}))
    sxy = Reduce('+',lapply(Vlists,function(x){x[['XViy']]}))
    syy = Reduce('+',lapply(Vlists,function(x){x[['yViy']]}))
    ell = Reduce('+',lapply(Vlists,function(x){x[['logdet']]}))
    f1 = syy - crossprod(solve(Sxx,as.vector(sxy)),as.vector(sxy)) +  ell 
	
    if(EstMeth == "REML") f1 <- f1 + determinant(Sxx, logarithm = TRUE)$modulus
	  nmult <- (n - p)*log(2*pi)
	  if(EstMeth == "ML") nmult <- n*log(2*pi)
	  result <- f1 + nmult

    as.numeric(result)
}

makeSigijMats = function(DFr, xy, CorModels, theta, addfunccol, subSampIndxCol, i, j,
	distPath, useTailDownWeight = FALSE, Vlist)
{    
#		DFr = DF
#		resampIndxCol = 'resampIndx'
#		y = junk$y
#		X = junk$X
#		xy = junk$xy
#		CorModels = c("Exponential.tailup", "Exponential.taildown", 
#			"Exponential.Euclid", "locID")
#		addfunccol = 'afvArea'
#		i = 1
#		j = 2
#		distPath = ssn@path
#		Vlist = ssnbdout$Vlist

    REind <- which(names(DFr) %in% CorModels)
    if(length(REind)) {
      REnames <- sort(names(DF)[REind])
      ## model matrix for a RE factor
      for(ii in 1:length(REind)) {
        DFr[,REnames[ii]] = as.factor(as.character(DFr[,REnames[ii]]))
			}
		}
    DFi = DFr[DFr[subSampIndxCol] == i,] 
		DFj = DFr[DFr[subSampIndxCol] == j,]
    nIDs = sort(as.integer(as.character(unique(rbind(DFi,DFj)[,"netID"]))))
    netDJ = matrix(0, nrow = length(DFi[,1]), ncol = length(DFj[,1]))
    netDJt = t(netDJ)
    net0 =  matrix(0, nrow = length(DFi[,1]), ncol = length(DFj[,1]))
		net0t = t(net0)
    rn = NULL
		cn = NULL
    nsofari = nsofari1 = 0
		nsofarj = nsofarj1 = 0
    distordi <- order(as.integer(as.character(DFi[,"netID"])),
      DFi[,"pid"])
    distordj <- order(as.integer(as.character(DFj[,"netID"])),
      DFj[,"pid"])

    DFi = DFi[distordi,]
    DFj = DFj[distordj,]
    names(distordi) <- rownames(DFi)[distordi]
    names(distordj) <- rownames(DFj)[distordj]
    for(k in nIDs) {
      workspace.name <- paste0("dist.net", k, ".RData")
      path <- file.path(distPath, "distance", "obs",
		    workspace.name)
	    if(!file.exists(path)) {
		    stop("Unable to locate required distance matrix")
	    }
	    file_handle <- file(path, open="rb")
	    distmat <- unserialize(file_handle)
	    ordpi <- order(as.numeric(rownames(distmat)))
	    close(file_handle)
      distmatk = distmat[rownames(distmat) %in% DFi$pid, 
        rownames(distmat) %in% DFj$pid, drop = FALSE]
      distmatkt = distmat[rownames(distmat) %in% DFj$pid, 
        rownames(distmat) %in% DFi$pid, drop = FALSE]
	    nki <- dim(distmatk)[1]
	    nkj <- dim(distmatk)[2]
			if(nki > 0) nsofari1 = nsofari + nki
			if(nkj > 0) nsofarj1 = nsofarj + nkj
			if(nki > 0 & nkj > 0) {
			  netDJ[(nsofari+1):nsofari1,(nsofarj+1):nsofarj1] <- distmatk
			  net0[(nsofari+1):nsofari1,(nsofarj+1):nsofarj1] <- 1
			  netDJt[(nsofarj+1):nsofarj1,(nsofari+1):nsofari1] <- distmatkt
			  net0t[(nsofarj+1):nsofarj1,(nsofari+1):nsofari1] <- 1

			}
      rn = c(rn,rownames(distmatk))
      cn = c(cn,colnames(distmatk))
	    nsofari = nsofari1
			nsofarj = nsofarj1
    }
    xci = xy[DFr[subSampIndxCol] == i,1]
    yci = xy[DFr[subSampIndxCol] == i,2]
    xcj = xy[DFr[subSampIndxCol] == j,1]
    ycj = xy[DFr[subSampIndxCol] == j,2]
    rownames(netDJ) = rn
    colnames(netDJ) = cn
    rownames(net0) = rn
    colnames(net0) = cn
    rownames(netDJt) = cn
    colnames(netDJt) = rn
    rownames(net0t) = cn
    colnames(net0t) = rn
    A = pmax(netDJ,t(netDJt))
    B = pmin(netDJ,t(netDJt))
    netD = as.matrix(netDJ + t(netDJt))
    W = NULL
	  if(length(grep("tailup",CorModels)) | useTailDownWeight == TRUE) {
		  if(missing(addfunccol) || is.null(addfunccol) ||
			  length(addfunccol) == 0 || !(addfunccol %in% colnames(DFr)))
	      	stop("The addfunccol argument is missing or mis-specified")
        FCmat <- 1 - (B > 0)*1
	      W <- sqrt(
          pmin(
            outer(DFi[,addfunccol],
	              rep(1, times = dim(A)[2])), 
		          t(outer(DFj[, addfunccol],
                rep(1, times = dim(A)[1]))) ) /
		      pmax(
            outer(DFi[, addfunccol],
                rep(1, times = dim(A)[2])),
			        t(outer(DFj[, addfunccol],
                rep(1, times = dim(A)[1])))))*
			        FCmat*net0
    }
    Zs <- NULL
    if(length(REind)) {
      Zs <- list()
        ## model matrix for a RE factor
        for(ii in 1:length(REind)) {
          Zs[[REnames[ii]]] = tcrossprod(model.matrix(~DFi[,REnames[ii]] - 1),
						model.matrix(~DFj[,REnames[ii]] - 1))
          rownames(Zs[[REnames[ii]]]) = DFi$pid
          colnames(Zs[[REnames[ii]]]) = DFj$pid
        }
      }

     Cij = SSN:::makeCovMat(theta = theta, dist.hydro = netD, 
        a.mat = A, b.mat = B, 
        w.matrix = W, net.zero = net0, 
        x.row = xci, y.row = yci, 
        x.col = xcj, y.col = ycj, 
        useTailDownWeight = useTailDownWeight,
        CorModels = CorModels, 
        use.nugget = FALSE, use.anisotropy = FALSE, 
        REs = Zs)

				return(list(XViCViX = t(Vlist[[i]]$ViX) %*% Cij %*% Vlist[[j]]$ViX))
}
