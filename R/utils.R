################################################################################
################################################################################
#                          distUpdater
################################################################################
################################################################################


distUpdater = function(ssnr, DFr, y, X, xy, CorModels, addfunccol, subSampIndxCol, i,
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
#			workspace.name <- "dist.net2.bmat"
#			workspace.name <- paste0("dist.net", k, ".bmat")
#      workspace.name <- paste("dist.net", k, ".RData", sep = "")
#      path <- file.path(distPath, "distance", "obs",
#		    workspace.name)
#	    distMatPoint = fm.open(path)
#	    file_handle <- file(path, open="rb")
#	    distmat <- unserialize(file_handle)
#	    ordpi <- order(as.numeric(rownames(distmat)))
#	    close(file_handle)
#      distmatk = distmati[rownames(distmati) %in% DFi$pid, 
#        rownames(distmati) %in% DFi$pid, drop = FALSE]
      distmatk = getStreamDistMatInt(ssnr,
				DFi[DFi$netID == k, 'pid'], 'Obs')[[1]]
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
        REnames <- sort(names(DFr)[REind])
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

################################################################################
################################################################################
#                          distList
################################################################################
################################################################################

distList = function(ssnr, DFr, y, X, xy, CorModels, addfunccol, subSampIndxCol,
	distPath)
{
  nResamp = max(DFr[,subSampIndxCol])
  Dlists = foreach(i=1:nResamp) %dopar% {
    distUpdater(ssnr = ssnr, DFr = DFr, y = y, X = X, 
    xy = xy, CorModels = CorModels,  addfunccol = addfunccol, 
    subSampIndxCol = subSampIndxCol, i = i, distPath = distPath)
  }
 Dlists
}

################################################################################
################################################################################
#                          m2LLstrbd
################################################################################
################################################################################

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

################################################################################
################################################################################
#                          makeSigijMats
################################################################################
################################################################################

makeSigijMats = function(ssnr, DFr, xy, CorModels, theta, addfunccol, subSampIndxCol, i, j,
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
      REnames <- sort(names(DFr)[REind])
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
			if(sum(DFi$netID == k) > 0 & sum(DFj$netID == k) > 0) {
				distmat = getStreamDistMatInt(ssnr,
					DFi[DFi$netID == k, 'pid'], 'Obs', 
					DFj[DFj$netID == k, 'pid'], 'Obs')
				distmatk = distmat[[1]]
				distmatkt = distmat[[2]]
			}
	    nki <- sum(DFi$netID == k)
	    nkj <- sum(DFj$netID == k)
			nsofari1 = nsofari + nki
			nsofarj1 = nsofarj + nkj
			if(nki > 0 & nkj > 0) {
			  netDJ[(nsofari+1):nsofari1,(nsofarj+1):nsofarj1] <- distmatk
			  net0[(nsofari+1):nsofari1,(nsofarj+1):nsofarj1] <- 1
			  netDJt[(nsofarj+1):nsofarj1,(nsofari+1):nsofari1] <- distmatkt
			  net0t[(nsofarj+1):nsofarj1,(nsofari+1):nsofari1] <- 1
			}
			if(sum(DFi$netID == k) > 0)
				rn = c(rn,rownames(DFi[DFi$netID == k,,drop = F]))
			if(sum(DFj$netID == k) > 0)
				cn = c(cn,rownames(DFj[DFj$netID == k,,drop = F]))
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


dMatsEtc = function(ssn, CorModels, dname1, DF1, xy1, addfunccol = NULL,
	dname2 = NULL, DF2 = NULL, xy2 = NULL)
{	
#	ssn = ecp$ssn
#	dname1 ='Obs'
#	dname2 ='pred1km'
	useTailDownWeight = FALSE
	a.mat <- NULL
	b.mat <- NULL
	net.zero <- NULL
	w.matrix <- NULL
	dist.hydro <- NULL
	flow.con.mat = NULL
	dist.junc = NULL
	dist.junc.a = NULL
	dist.junc.b = NULL
	REs = NULL
	REPs = NULL
	rnames <- NULL
	cnames = NULL
	# if any "tail models"
	if(length(grep("tail",CorModels)) > 0 | useTailDownWeight == TRUE){
			n1 = dim(DF1)[1]
			n2 = dim(DF2)[1]
			if(is.null(dname2)) {
				dist.junc <- matrix(0, nrow = n1, ncol = n1)
				net.zero <-  matrix(0, nrow = n1, ncol = n1)
			}
			if(!is.null(dname2)) {
				npred = dim(DF2)[1]
				dist.junc.a <- matrix(0, nrow = n1, ncol = n2)
				dist.junc.b <- matrix(0, nrow = n2, ncol = n1)
				net.zero <-  matrix(0, nrow = n1, ncol = npred)
			}		
			nsofar = 0
			nsofari = 0
			nsofarj = 0
			nIDs <- sort(as.integer(as.character(unique(c(DF1$netID,DF2$netID)))))
			for(k in nIDs){
				if(is.null(dname2)) {
					distmat = getStreamDistMatInt(ssn,
						DF1[DF1$netID == k, 'pid'], dname1)[[1]]
					ni <- length(distmat[1,])
				}
				if(!is.null(dname2)) {
					distmat = getStreamDistMatInt(ssn,
						DF1[DF1$netID == k, 'pid'], dname1,
						DF2[DF2$netID == k, 'pid'], dname2)
					ni = dim(distmat[[1]])[1]
					nj = dim(distmat[[1]])[2]
				}
				rnames <- c(rnames,rownames(distmat[[1]]))
				cnames = c(cnames,colnames(distmat[[1]]))
				if(is.null(dname2)) {
					dist.junc[(nsofar + 1):(nsofar + ni),
						(nsofar + 1):(nsofar + ni)] <- distmat
					net.zero[(nsofar + 1):(nsofar + ni),
						(nsofar + 1):(nsofar + ni)] <- 1
					nsofar <- nsofar + ni
				}
				if(!is.null(dname2)) {
					dist.junc.a[(nsofari + 1):(nsofari + ni),
						(nsofarj + 1):(nsofarj + nj)] <- distmat[[1]]
					dist.junc.b[(nsofarj + 1):(nsofarj + nj),
						(nsofari + 1):(nsofari + ni)] <- distmat[[2]]
					net.zero[(nsofari + 1):(nsofari + ni),
						(nsofarj + 1):(nsofarj + nj)] <- 1
				}
			}
			if(is.null(dname2)) {
				rownames(dist.junc) = rnames
				colnames(dist.junc) = cnames
			}
			if(!is.null(dname2)) {
				rownames(dist.junc.a) = rnames
				colnames(dist.junc.a) = cnames
				rownames(dist.junc.b) = cnames
				colnames(dist.junc.b) = rnames
			}

			rownames(net.zero) = rnames
			colnames(net.zero) = cnames
			if(is.null(dname2)) {
				# maximum distance to common junction between two sites
				a.mat <- pmax(dist.junc,t(dist.junc))
				# minimum distance to common junction between two sites
				b.mat <- pmin(dist.junc,t(dist.junc))
				# hydrological distance
				dist.hydro <- as.matrix(dist.junc + t(dist.junc))*net.zero
			}
			if(!is.null(dname2)) {
				# creat A matrix (longest distance to junction of two points)
    		a.mat <- pmax(dist.junc.a,t(dist.junc.b))
    		# creat B matrix (shorted distance to junction of two points)
    		b.mat <- pmin(dist.junc.a,t(dist.junc.b))
    		# get hydrologic distance
    		dist.hydro <- as.matrix(dist.junc.a + t(dist.junc.b))
			}
			# binary flow connection matrix
			flow.con.mat <- 1 - (b.mat > 0)*1

			if(is.null(dname2)) {
				DF2 = DF1
				n2 = n1
			}
			w.matrix <- sqrt(pmin(outer(DF1[,addfunccol], rep(1, times = n2)),
				t(outer(DF2[,addfunccol],rep(1, times = n1)))) /
				pmax(outer(DF1[,addfunccol],rep(1, times = n2)),
				t(outer(DF2[,addfunccol],rep(1, times = n1)))))*flow.con.mat*net.zero
		}
# end "tail model parts
# if any random effects
		REind <- which(names(DF1) %in% CorModels)
		if(length(REind)) {
				REs <- list()
				REnames <- sort(names(DF1)[REind])
				## model matrix for a RE factor
				for(ii in 1:length(REind)) REs[[REnames[ii]]] <-
					model.matrix(~DF1[,REnames[ii]] - 1)
					rownames(REs[[REnames[ii]]]) <- DF1[,"pid"]
				## corresponding block matrix
				for(ii in 1:length(REind)) REs[[ii]] <- REs[[ii]] %*% t(REs[[ii]])
			if(!is.null(dname2)) {
					for(ii in 1:length(REind)) if(any(is.na(DF2[,REnames[ii]])))
					 stop("Cannot having missing values when creating random effects")
					REOs <- list()
					REPs <- list()
					ObsSimDF <- DF1
					PredSimDF <- DF2
					## model matrix for a RE factor
					for(ii in 1:length(REnames)){
						#we'll add "o" to observed levels and "p" to prediction
						# levels so create all possible levels
						plevels <- unique(c(levels(PredSimDF[,REnames[[ii]]]),
							paste("o",levels(ObsSimDF[,REnames[[ii]]]),sep = ""),
							paste("p",levels(PredSimDF[,REnames[[ii]]]),sep = "")))
						# sites with prediction levels same as observation levels
						pino <- PredSimDF[,REnames[[ii]]] %in% ObsSimDF[,REnames[[ii]]]
						#add "o" to observed levels
						ObsSimDF[,REnames[[ii]]] <- paste("o",
							ObsSimDF[,REnames[[ii]]], sep = "")
						ObsSimDF[,REnames[[ii]]] <- as.factor(as.character(
							ObsSimDF[,REnames[[ii]]]))
						#add all possible levels to prediction data frame
						levels(PredSimDF[,REnames[[ii]]]) <- plevels
								# add "o" to prediction sites with observation levels
						if(any(pino)) PredSimDF[pino,REnames[[ii]]] <- paste("o",
							PredSimDF[pino,REnames[[ii]]], sep = "")
						# add "p" to all predicition sites without observation levels
						if(any(!pino)) PredSimDF[!pino,REnames[[ii]]] <- paste("p",
							PredSimDF[!pino,REnames[[ii]]], sep = "")
						PredSimDF[,REnames[[ii]]] <- as.factor(as.character(
							PredSimDF[,REnames[[ii]]]))
						# now get down to just levels with "o" & "p" added
						blevels <- unique(c(levels(ObsSimDF[,REnames[[ii]]]),
							levels(PredSimDF[,REnames[[ii]]])))
						ObsSimDF[,REnames[[ii]]] <- factor(ObsSimDF[,REnames[[ii]]],
							levels = blevels, ordered = FALSE)
						PredSimDF[,REnames[[ii]]] <- factor(PredSimDF[,REnames[[ii]]],
							levels = blevels, ordered = FALSE)
						# now ordering of factors in Z matrices should be compatible
						# with obs x obs Z matrices
						REOs[[ii]] <- model.matrix(~ObsSimDF[,
							REnames[[ii]]] - 1)
						REPs[[ii]] <- model.matrix(~PredSimDF[,
							REnames[[ii]]] - 1)
						rownames(REOs[[ii]]) <- DF1[,"pid"]
						rownames(REPs[[ii]]) <- DF2[,"pid"]
						if(any(rownames(REs[[ii]])!=rownames(a.mat)))
						stop("rownames RE for obs do not match rownames of a.mat")
						if(any(rownames(REPs[[ii]])!=colnames(a.mat)))
						stop("rownames RE for preds do not match colnames of a.mat")
					}
					## corresponding block matrix
					for(ii in 1:length(REnames)) REPs[[ii]] <-
						REOs[[ii]] %*% t(REPs[[ii]])
				}
			}
# end random effects

	return(list(dist.junc = dist.junc, dist.junc.a = dist.junc.a, 
		dist.junc.b = dist.junc.b, net.zero = net.zero, a.mat = a.mat,
		b.mat = b.mat, dist.hydro = dist.hydro, w.matrix = w.matrix,
		REs = REs, REPs = REPs))
}

