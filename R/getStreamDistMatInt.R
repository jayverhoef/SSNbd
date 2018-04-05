#-------------------------------------------------------------------------------
#
#           getStreamDistMatInt
#
#-------------------------------------------------------------------------------

#' Extracts downstream-only hydrologic distances for two vectors of pids
#'
#' Extracts downstream-only hydrologic distances for two vectors of pids and
#' returns a list of distance matrices.
#'
#' @param ssn a SpatialStreamNetwork object created using the SSN package
#'
#' @param pidset1 a vector of pids found within a SpatialStreamNetwork object
#'
#' @param dataID1 a character string representing the name of the dataset found
#'   in the SpatialStreamNetwork object slot ID
#'
#' @param pidset2 a vector of pids found within a SpatialStreamNetwork object.
#'   Default = NULL.
#'
#' @param dataID2 a character string representing the name of the dataset found
#'   in the SpatialStreamNetwork object slot ID. Default = NULL.
#'
#'
#' @return a list matrices containing the of downstream-only hydrologic distance for
#'   pidset1 and pidset2. If pidset2 is omitted, the function returns the downstream-only
#'   hydrologic distances between pidset1 and pidset1. See documentation for
#'   createBigDistMat for more information about the format of the distance matrices.
#'
#' @author Erin Peterson
#' @export
#'
#' @examples
#' ## ssn<- importSSN("raw5.ssn", predpts = "preds")
#' ## createBigDistMat(ssn, predpts = "preds", o.write = TRUE, amongpreds = TRUE,
#' ##                  no.cores = 2)
#' ## pidset1 <- c(1, 3, 9, 28, 34, 51)
#' ## pidset2 <- c(4, 7, 30, 35)

#' ## distmats <- getStreamDistMatInt(ssn, pidset1, dataID1 = "preds",
#' ##                                  pidset2 = pidset2, dataID2 = "preds")


getStreamDistMatInt <- function(ssn, pidset1, dataID1, pidset2 = NULL,
                                dataID2 = NULL)
{
    ##library(filematrix)
    ## Check arguments
    if(class(ssn) != "SpatialStreamNetwork")
        stop("Object not of class SpatialStreamNetwork")

    ## Set dataID2 and pidset2 if necessary
    if(is.null(dataID2)) dataID2 <- dataID1
    if(is.null(pidset2)) pidset2 <-pidset1

    if(sum(pidset1 %in% pidset2) == length(pidset1) & sum(pidset2 %in% pidset1) == length(pidset2)){
        same <- TRUE
    } else { same <- FALSE}


    ## Fix typos
    if(dataID1 == "Obs") dataID1 = "obs"
    if(dataID2 == "Obs") dataID2 = "obs"

    ## Make sure that obs always comes first
    if(dataID1 != "obs" & dataID2 == "obs") {
        ## dataID2 <- dataID1
        ## dataID2 <- "obs"

        ## tmp <- pidset1
        ## pidset1 <- pidset2
        ## pidset2 <- tmp
        ## rm(tmp)
        stop("The observed pids must be passed as pidset1.")
    }

    ## IDENTIFY DATASETS----------------------------
    ## set object x1
    if(dataID1 == "obs") {
        x1 <- ssn@obspoints
    } else {
        x1 <- ssn@predpoints
    }

    ## set object x2
    if(dataID2 == "obs") {
        x2 <- ssn@obspoints
    } else if(!is.null(dataID2)){
        x2 <- ssn@predpoints
    }

    ## Get ID index for pidset1
    count <- 0
    if(dataID1 == "obs") {
        id.num1 <- 1
    } else {
        if(length(x1@ID) > 0) {
            for (m in 1:length(x1@ID)) {
                if (x1@ID[m] == dataID1) {
                    id.num1 <- m
                    count <- count + 1
                } else {
                    count <- count + 1}}
                if(count == 0) {
                    stop(paste0(dataID1, " not found in SSN"))}
        } else {
            stop("No prediction sites in SSN")}
    }

    ## get ID index for pidset2
    count <- 0
    if(dataID2 == "obs") {
        id.num2 <- 1
    } else {
        if(length(x2@ID) > 0) {
            for (m in 1:length(x2@ID)) {
                if (x2@ID[m] == dataID2) {
                    id.num2 <- m
                    count <- count +1
                } else {
                    count <- count + 1}}
                if(count == 0) {
                    stop(paste0(dataID2, " not found in SSN"))}
        } else {
            stop("No prediction sites in SSN")}
    }


    ## Get netIDs for pidset1 and pidset2
    ind1 <- x1@SSNPoints[[id.num1]]@point.data$pid %in% pidset1
    data1 <- cbind(x1@SSNPoints[[id.num1]]@point.data$pid[ind1], as.numeric(as.character(x1@SSNPoints[[id.num1]]@point.data$netID[ind1])))
    colnames(data1) <- c("pid", "netID")
    data1 <- data1[order(data1[,"netID"], data1[,"pid"]),]
    nIDs1 <- unique(data1[,"netID"])

    if(same == TRUE) {
        data2 <- data1
        nIDs2 <- nIDs1
    } else {
        ind2 <- x2@SSNPoints[[id.num2]]@point.data$pid %in% pidset2
        data2 <- cbind(x2@SSNPoints[[id.num2]]@point.data$pid[ind2], as.numeric(as.character(x2@SSNPoints[[id.num2]]@point.data$netID[ind2])))
        colnames(data2) <- c("pid", "netID")
        data2 <- data2[order(data2[,"netID"], data2[,"pid"]),]
        nIDs2 <- unique(data2[,"netID"])
    }
    ## Check to make sure all pids came from 1 dataset. Take this
    ## out later.
    if(length(pidset1) != nrow(data1)) {
        stop(paste("Some pids in pidset1 were not found in", dataID1))}

    if(length(pidset2) != nrow(data2)) {
        stop(paste("Some pids in pidset2 were not found in", dataID2))}

##############################################################
    ## Fill obs x obs OR preds x preds matrix
    if(same == TRUE) {
        ## Create empty distance matrix
        pid1.pid1 <- matrix(-1, nrow = length(pidset1),
                          ncol = length(pidset1))
        rownames(pid1.pid1) <- as.character(data1[,"pid"])
        colnames(pid1.pid1) <- as.character(data1[,"pid"])

        for(i in 1:length(nIDs1)) {
            file.name <-file.path(ssn@path,
                   "distance", dataID1, paste0("dist.net", nIDs1[i]))

            ## Check to make sure the file exists
            if (!file.exists(paste0(file.name, ".bmat"))) {
                stop("missing distance matrix")
            }
            dist <- fm.open(filenamebase = file.name, readonly = FALSE)
            ind.sub <- data1[,"netID"] == nIDs1[i]
            ind.whole <- colnames(dist) %in% data1[ind.sub, "pid"]

            ## Extract whole column because its faster
            tmp <- dist[,ind.whole]
            colnames(tmp)<- colnames(dist)[ind.whole]
            tmp <- tmp[,order(as.numeric(colnames(tmp)))]

            ## If tmp is a single column, make sure its a matrix
            if(!is.matrix(tmp)) {tmp <- as.matrix(tmp)}

            pid1.pid1[ind.sub, ind.sub]<- tmp[ind.whole,]
            ##pid1.pid1[ind.sub,ind.sub] <- dist[ind.whole,ind.whole]
            rm(tmp)
            ##close(dist)
        }
        rm(dist)
        ind <- pid1.pid1 == -1
        pid1.pid1[ind]<- NA
        dist.list <- list(pid1.pid1)
        attr(dist.list, "Matrix.Name") <- "dist"
    }
################################################################
    ## Fill in obs x obs OR preds x preds with different pids
    ## A = pidset1 x pidset2, B = pidset2 x pidset1
    ## All Matrices are TO x FROM
    if(same == FALSE & dataID1==dataID2){

        ## Create TO x FROM matrices
        ## pid1.pid2.a = TO pidset1 x FROM pidset2
        pid1.pid2.a <- matrix(-1, nrow = length(pidset1),
                              ncol = length(pidset2))
        rownames(pid1.pid2.a) <- as.character(data1[,"pid"])
        colnames(pid1.pid2.a) <- as.character(data2[,"pid"])

        ## TO pidset2 FROM pidset1
        pid2.pid1.b <- matrix(-1, nrow = length(pidset2),
                             ncol = length(pidset1))
        rownames(pid2.pid1.b) <- as.character(data2[,"pid"])
        colnames(pid2.pid1.b) <- as.character(data1[,"pid"])

        nIDs.all <- intersect(nIDs1, nIDs2)

        for(i in 1:length(nIDs.all)) {
            file.name <-file.path(ssn@path,
                     "distance", dataID1, paste0("dist.net", nIDs.all[i]))

            ## Check to make sure the files exist
            if (!file.exists(paste0(file.name, ".bmat"))) {
                stop("missing distance matrix")
            }
            ## t(pid1.pid2.a) = obs1 x preds2, pid2.pid1.b = preds2 x obs1
            dist <- fm.open(filenamebase = file.name, readonly = FALSE)

            ## Populate dist.a--------------------------------------------
            ind.sub1 <- data1[,"netID"] == nIDs.all[i]
            ind.sub2 <- data2[,"netID"] == nIDs.all[i]
            ind.whole1 <- colnames(dist) %in% data2[ind.sub2, "pid"]

            ## Extract whole column for speed
            tmp.c <- dist[,ind.whole1]

            ## Double check ordering
            colnames(tmp.c)<- colnames(dist)[ind.whole1]
            rownames(tmp.c) <- rownames(dist)
            tmp.c <- tmp.c[,order(as.numeric(colnames(tmp.c)))]

            ## If tmp.c is a single column, make sure its a matrix
            if(!is.matrix(tmp.c)) {tmp.c <- as.matrix(tmp.c)}

            ## Subset by row and populate pid1.pid2.a
            ind.whole2 <- rownames(tmp.c) %in% data1[ind.sub1, "pid"]
            pid1.pid2.a[ind.sub1, ind.sub2]<- tmp.c[ind.whole2,]
            ## close(dist)
            rm(ind.whole1, ind.whole2, tmp.c)

            ## Populate dist.b--------------------------------------------
            ind.whole1 <- colnames(dist) %in% data1[ind.sub1, "pid"]

            ## Extract whole column for speed
            tmp.c <- dist[,ind.whole1]

            ## Double check ordering
            colnames(tmp.c)<- colnames(dist)[ind.whole1]
            rownames(tmp.c) <- rownames(dist)
            tmp.c <- tmp.c[,order(as.numeric(colnames(tmp.c)))]

            ## If tmp.c is a single column, make sure its a matrix
            if(!is.matrix(tmp.c)) {tmp.c <- as.matrix(tmp.c)}

            ## Subset by row and populate pid1.pid2.a
            ind.whole2 <- rownames(as.matrix(tmp.c)) %in% data2[ind.sub2, "pid"]
            pid2.pid1.b[ind.sub2, ind.sub1]<- tmp.c[ind.whole2,]
            ##close(dist)
            rm(ind.whole1, ind.whole2, ind.sub1, ind.sub2, tmp.c)
        }

        rm(dist)
        ind <- pid1.pid2.a == -1
        pid1.pid2.a[ind]<- NA

        ind <- pid2.pid1.b == -1
        pid2.pid1.b[ind]<- NA
        dist.list <- list(pid1.pid2.a, pid2.pid1.b)
        attr(dist.list, "Matrix.Name") <- "dist.a, dist.b"
    }


    ################################################################
    ## Fill A: obs x preds and B: preds x obs
    ## Corresponds to pidset1 x pidset2, pidset2 x pidset1
    ## All Matrices are TO x From

    if(same == FALSE & dataID1 != dataID2) {

        ## ## Get data folder
        ## if(dataID1 == dataID2){
        ##     dist.dir <- dataID1
        ## } else {
        ##     if(dataID1 == "obs"){
        ##         dist.dir <- dataID2}
        ##     if(dataID2 == "obs") {
        ##         dist.dir <- dataID1}
        ##     if(dataID1 != "obs" & dataID2 != "obs") {
        ##         stop("extracting distances between pids from two prediction sites datasets is not currently supported.")}
        ## }


        ## Create TO x FROM matrices
        ## pid1.pid2.a = TO pidset1 x FROM pidset2
        pid1.pid2.a <- matrix(-1, nrow = length(pidset1),
                              ncol = length(pidset2))
        rownames(pid1.pid2.a) <- as.character(data1[,"pid"])
        colnames(pid1.pid2.a) <- as.character(data2[,"pid"])

        ## TO pidset2 FROM pidset1
        pid2.pid1.b <- matrix(-1, nrow = length(pidset2),
                             ncol = length(pidset1))
        rownames(pid2.pid1.b) <- as.character(data2[,"pid"])
        colnames(pid2.pid1.b) <- as.character(data1[,"pid"])

        nIDs.all <- intersect(nIDs1, nIDs2)

        for(i in 1:length(nIDs.all)) {
            file.name.a <-file.path(ssn@path,
                                    "distance", dataID2, paste0("dist.net", nIDs.all[i], ".a"))
            file.name.b <-file.path(ssn@path,
                   "distance", dataID2, paste0("dist.net", nIDs.all[i], ".b"))

            ## Check to make sure the files exist
            if (!file.exists(paste0(file.name.a, ".bmat"))|| !file.exists(paste0(file.name.b, ".bmat")) ) {
                stop("missing distance matrix")
            }
            ## t(pid1.pid2.a) = obs1 x preds2, pid2.pid1.b = preds2 x obs1
             dist.a <- fm.open(filenamebase = file.name.a, readonly = FALSE)

            ## Populate dist.a--------------------------------------------
            ind.sub1 <- data1[,"netID"] == nIDs.all[i]
            ind.sub2 <- data2[,"netID"] == nIDs.all[i]
            ind.whole1 <- colnames(dist.a) %in% data1[ind.sub1, "pid"]

            ## Extract whole column for speed
            tmp.c <- dist.a[,ind.whole1]

            ## Double check ordering
            colnames(tmp.c)<- colnames(dist.a)[ind.whole1]
            rownames(tmp.c) <- rownames(dist.a)
            tmp.c <- tmp.c[,order(as.numeric(colnames(tmp.c)))]

            ## If tmp.c is a single column, make sure its a matrix
            if(!is.matrix(tmp.c)) {tmp.c <- as.matrix(tmp.c)}

            ## Subset by row and populate pid1.pid2.a
            ind.whole2 <- rownames(tmp.c) %in% data2[ind.sub2, "pid"]
            pid1.pid2.a[ind.sub1, ind.sub2]<- t(tmp.c[ind.whole2,])
            ##close(dist.a)
            rm(ind.whole1, ind.whole2)

            ## Populate dist.b--------------------------------------------
            dist.b <- fm.open(filenamebase = file.name.b, readonly = FALSE)
            ind.whole1 <- colnames(dist.b) %in% data1[ind.sub1, "pid"]

            ## Extract whole column for speed
            tmp.c <- dist.b[,ind.whole1]

            ## Double check ordering
            colnames(tmp.c)<- colnames(dist.b)[ind.whole1]
            rownames(tmp.c) <- rownames(dist.b)
            tmp.c <- tmp.c[,order(as.numeric(colnames(tmp.c)))]

            ## If tmp.c is a single column, make sure its a matrix
            if(!is.matrix(tmp.c)) {tmp.c <- as.matrix(tmp.c)}

            ## Subset by row and populate pid1.pid2.a
            ind.whole2 <- rownames(tmp.c) %in% data2[ind.sub2, "pid"]
            pid2.pid1.b[ind.sub2, ind.sub1]<- tmp.c[ind.whole2,]
            ##close(dist.b)
            rm(ind.whole1, ind.whole2, ind.sub1, ind.sub2)
        }

        rm(dist.a, dist.b)
        ind.a <- pid1.pid2.a == -1
        pid1.pid2.a[ind.a]<- NA
        ind.b <- pid2.pid1.b == -1
        pid2.pid1.b[ind.b]<- NA
        dist.list <- list(pid1.pid2.a, pid2.pid1.b)
        attr(dist.list, "Matrix.Name") <- c("dist.a", "dist.b")
    }

    ## return list
    return(dist.list)

}



