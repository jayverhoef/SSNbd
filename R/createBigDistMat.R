createBigDistMat <- function(ssn, predpts = NULL, o.write = FALSE, amongpreds = FALSE,
                             no.cores = 1) {

    ## ssn = SpatialStreamNetwork object
    ## predpts = name of prediction points in character format
    ## o.write = should existing distance directories be overwritten
    ## amongpreds = should the pred-pred distance matrix be generated
    ## no.cores = number of cores to use and number of splits in dataset


    ## require(filematrix)
    ## require(foreach)
    ## require(doParallel)
    ## require(itertools)

  #############################################################################
  ## Check that arguments are valid and set up directories
  if(amongpreds && (missing(predpts) || is.null(predpts)))
  {
	stop("A named collection of prediction points must be specified via the predpts option when amongpreds is TRUE")
  }
  ##Check to see whether distance folder exists...
  if (!file.exists(file.path(ssn@path, "distance"))) {
    dir.create(file.path(ssn@path, "distance"))
  }

  ##Check whether an observation folder exists
  obs.path <- file.path(ssn@path, "distance", "obs")
  if (!file.exists(obs.path)) {
    dir.create(obs.path)
  }

  ## And then whether prediction folder exists
  if (!is.null(predpts)) {
    preds.path <- file.path(ssn@path, "distance", predpts)
    if(!file.exists(preds.path)) {
      dir.create(preds.path)
    }
    ## Get ID number from predictions
    count <- 0
    if(length(ssn@predpoints@ID) > 0) {
	for (m in 1:length(ssn@predpoints@ID)) {
	    if (ssn@predpoints@ID[m] == predpts) {
	         pred.num <- m
                 count <- count + 1}
	}
    }

    if (count==0) {
      stop(predpts, " does not exist in SSN")}

    if (count > 1) {
      stop("SSN contains more than one copy of ", predpts)}

    ssn@predpoints@SSNPoints[[pred.num]]@point.data$netID<- as.factor(ssn@predpoints@SSNPoints[[pred.num]]@point.data$netID)
  }

  if (is.null(predpts)) {
      pred.num <- 0}


  ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID<-
      as.factor(ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID)
  net.count <- length(levels(ssn@network.line.coords$NetworkID))
  warned.overwrite <- FALSE

      if (file.exists(file.path(ssn@path,"binaryID.db")) == FALSE)
        stop("binaryID.db is missing from ssn object")

      driver <- RSQLite::SQLite()
      connect.name <- file.path(ssn@path,"binaryID.db")

  connect <- dbConnect(SQLite(), connect.name)


      cl <- makeCluster(no.cores)
      registerDoParallel(cl)

      ## close connections upon function exit
      on.exit({
          dbDisconnect(connect)
          closeAllConnections()
      })
  #######################################################################
  ### FOR EACH NETWORK---------------------------------------------------
  for (i in 1:net.count) {
      net.num <- levels(ssn@network.line.coords$NetworkID)[i]
      ind.obs <- ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID == as.numeric(net.num)

      ## figure out how many observed sites there are in the network
      site.no <- nrow(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,])

      if(pred.num > 0) {
          ind.preds <- ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[,"NetworkID"] %in% net.num
          pred.site.no <- nrow(ssn@obspoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])
      } else {
          pred.site.no <- 0
      }

      net.name <- paste("net", net.num, sep = "")



      ## If observed sites > 0 #########################################
      if (site.no > 0) {

          ## get sorted pids to use as row and column names
	  obs.pids<- sort(as.numeric(rownames(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,])))

          ## Extract binaryID table
          bin.table <- dbReadTable(connect, net.name)

	  workspace.name1 <- paste(obs.path, "/dist.net", net.num, sep = "")

	  ## if(!o.write) {
	  ##       exists <- file.exists(file.path(ssn@path, "distance", "obs", workspace.name1))
	  ##       if (!missing(predpts) && !is.null(predpts))
	  ##       {
          ##           exists <- c(exists, file.exists(file.path(ssn@path,
          ##                       "distance", predpts, workspace.name.a)),
          ##                       file.exists(file.path(ssn@path, "distance",
          ##                       predpts, workspace.name.b)))
	  ##       }
	  ##       ##If they're all already there, warn but continue
	  ##       if(all(exists))
	  ##       {
	  ##       	if(!warned.overwrite) {
	  ##       	    warned.overwrite <- TRUE
	  ##       	    cat("Distance matrices already existed while o.write was set to FALSE. Not overwriting existing matrices\n")}
	  ##       	next
	  ##       }
	  ##       ##if some exist but some don't that's an error
	  ##       else if(any(exists) && any(!exists)) {
	  ##       	stop("o.write was set to FALSE and some (but not all) distance matrices already existed")}
	  ## }

          ######################################################################
          ## Create empty obs distance matrix
          current_distance_matrix <-
                fm.create(filenamebase = workspace.name1,
                          nrow = site.no, ncol = site.no, type = "double")


           rownames(current_distance_matrix) <- as.character(obs.pids)
           colnames(current_distance_matrix) <- as.character(obs.pids)

          ##Calculate distances between observed sites
          if (is.null(getDoParName())) {
              registerDoSEQ() # A little hack to avoid the foreach warning 1st time.
          }

          ##cl <- makeCluster(no.cores)
          ##registerDoParallel(cl)
          ## registerDoParallel(cores=no.cores)

          ## Create vector iterator
          itCol <- isplitVector(obs.pids, chunks = no.cores)

          ans1<- foreach(obs.pids.vec=itCol,.packages = c("SSN", "filematrix","itertools","iterators"),
                             .errorhandling = "pass") %dopar% {

                      amongObsBigDistMat(ssn = ssn, net.num = net.num, pids = obs.pids.vec,
                            bin.table = bin.table, workspace.name = workspace.name1)
                      finished1 <- TRUE
                      return(finished1)
                  }

          ## stopCluster(cl)
          ## registerDoSEQ()

      } ## END IF SITE.NO > 0

##################################################################################
      ## Calculate Observed to prediction site distances
      if(site.no>0 & pred.site.no > 0) {
          ## figure out how many prediction sites there are in the network
         ## ind.preds <- ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[,"NetworkID"] %in% net.num
         ## pred.site.no <- nrow(ssn@obspoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])

         #if(pred.site.no > 0) {

             ## Get pred pids
             pred.pids<- sort(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])))

             ## create empty distance matrices for obs-preds
             workspace.name.a <- paste(preds.path, "/dist.net", net.num, ".a", sep = "")
             workspace.name.b <- paste(preds.path, "/dist.net", net.num, ".b", sep = "")

             current_distance_matrix_a <-
                 fm.create(filenamebase = paste0(workspace.name.a),
                       nrow = pred.site.no, ncol = site.no, type = "double")

             rownames(current_distance_matrix_a) <- as.character(pred.pids)
             colnames(current_distance_matrix_a) <- as.character(obs.pids)

             current_distance_matrix_b <-
                    fm.create(filenamebase = paste0(workspace.name.b),
                              nrow = pred.site.no, ncol = site.no, type = "double")
             rownames(current_distance_matrix_b) <- as.character(pred.pids)
             colnames(current_distance_matrix_b) <- as.character(obs.pids)


             ## Create vector iterator
             itCol <- isplitVector(obs.pids, chunks = no.cores)

             ## cl <- makeCluster(no.cores)
             ## registerDoParallel(cl)

             ans2<- foreach(obs.pids.vec=itCol,.packages = c("SSN", "filematrix","itertools","iterators"),
                             .errorhandling = "pass") %dopar% {

                      amongObsPredsBigDistMat(ssn = ssn, net.num = net.num, pred.num, obs.pids = obs.pids.vec,
                                                    pred.pids=pred.pids, bin.table = bin.table, workspace.name.a = workspace.name.a,
                                                    workspace.name.b = workspace.name.b)
                      finished2 <- TRUE
                      return(finished2)
             }

             ## stopCluster(cl)
             ## registerDoSEQ()


      } ## END OBS-PRED SITE DISTANCES

##------------------------------------------------------
## ###############################################################################################

          ## Calculate distances between prediction sites ------------------------------------
          if (amongpreds & pred.num > 0) {

               pred.pids<- sort(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])))
               net.name <- paste("net", net.num, sep = "")
               bin.table <- dbReadTable(connect, net.name)

               among.path <- paste0(preds.path, "/dist.net", net.num)
               among_distance_matrix <- fm.create(filenamebase =among.path,
                     nrow = pred.site.no, ncol = pred.site.no, type = "double")
               rownames(among_distance_matrix) <- as.character(pred.pids)
               colnames(among_distance_matrix) <- as.character(pred.pids)


               ##Calculate distances between observed sites
               ## cl <- makeCluster(no.cores)
               ## registerDoParallel(cl)
          ##registerDoParallel(cores=no.cores)

               ## Create vector iterator
               itCol <- isplitVector(pred.pids, chunks = no.cores)

               ans3<- foreach(pred.pids.vec=itCol,.packages = c("SSN", "filematrix","itertools","iterators"),
                              .errorhandling = "pass") %dopar% {

                      amongPredsBigDistMat(ssn = ssn, net.num = net.num, pids = pred.pids.vec,
                            pred.num = pred.num, bin.table = bin.table,
                                           workspace.name = among.path)
                      finished3 <- TRUE
                      return(finished3)
               }
               ## stopCluster(cl)
               ## registerDoSEQ()

      }
  }
  stopCluster(cl)
  registerDoSEQ()

}

