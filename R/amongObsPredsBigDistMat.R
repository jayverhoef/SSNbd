

amongObsPredsBigDistMat <- function(ssn, net.num, pred.num, obs.pids, pred.pids, bin.table,
                                    workspace.name.a, workspace.name.b){

    ##require(filematrix)
    site.no <- length(obs.pids)
    pred.site.no <- length(pred.pids)

    obs.pids.ind <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% obs.pids
    pred.pids.ind <- ssn@predpoints@SSNPoints[[pred.num]]@point.data$netID %in% net.num

    ## CREATE DISTANCE MATRIX OUTSIDE OF THIS FUNCTION AND THEN ACCESS IT
    current_distance_matrix_a <- fm.open(filenamebase = workspace.name.a, readonly = FALSE)
    current_distance_matrix_b <- fm.open(filenamebase = workspace.name.b, readonly = FALSE)

    ## Create a data.frame of prediction point data
    locID.pred.data <- attributes(ssn@predpoints@SSNPoints[[1]]@network.point.coords)$locID[pred.pids.ind]
    pred.data <- as.data.frame(cbind(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[pred.pids.ind,])),
        as.numeric(levels(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[pred.pids.ind]))[ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[pred.pids.ind]],
        locID.pred.data, ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$DistanceUpstream[pred.pids.ind]))

    colnames(pred.data)<- c("pid","rid", "locID", "upDist")
    pred.data <- pred.data[order(pred.data$pid),]
    pred.data$locID <- as.factor(pred.data$locID)

    pred.data$binaryID <- bin.table$binaryID[match(pred.data$rid, bin.table$rid)]
    pred.data <- pred.data[order(pred.data[,"pid"]),]
    rownames(pred.data) <- pred.data$pid

    ##
    pred_by_locID <- pred.data[order(pred.data[,"locID"]),]
    pred_by_locID$pid <- as.numeric(pred_by_locID$pid)
    pred_by_locID$locID <- as.numeric(pred_by_locID$locID)
    pred_reordering <- order(pred_by_locID$pid)

    dup.pred <- !duplicated(pred_by_locID$locID)
    ##-----------

     ## Create a data.frame of observed point data
    locID.obs.data <- attributes(ssn@obspoints@SSNPoints[[1]]@network.point.coords)$locID[obs.pids.ind]
    obs.data <- as.data.frame(cbind(as.numeric(rownames(ssn@obspoints@SSNPoints[[1]]@network.point.coords[obs.pids.ind,])),
        as.numeric(levels(ssn@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[obs.pids.ind]))[ssn@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[obs.pids.ind]],
        locID.obs.data, ssn@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream[obs.pids.ind]))

    colnames(obs.data)<- c("pid","rid", "locID", "upDist")
    obs.data <- obs.data[order(obs.data$pid),]
    obs.data$locID <- as.factor(obs.data$locID)

    obs.data$binaryID <- bin.table$binaryID[match(obs.data$rid, bin.table$rid)]
    obs.data <- obs.data[order(obs.data[,"pid"]),]
    rownames(obs.data) <- obs.data$pid

    ##
    obs_by_locID <- obs.data[order(obs.data[,"locID"]),]
    obs_by_locID$pid <- as.numeric(obs_by_locID$pid)
    obs_by_locID$locID <- as.numeric(obs_by_locID$locID)
    obs_reordering <- order(obs_by_locID$pid)

    ##-----------

    ##locID values can be repeated, in which case they have the same distance data.
    locID.old <- -1
    for(ob in 1:(nrow(obs_by_locID))){

       ind.pid <- which(obs_by_locID$pid == obs.pids[ob])
       pid.ob <- obs_by_locID[ind.pid, "pid"]
       locID.ob <- obs_by_locID[ind.pid, "locID"]
       upDist.ob <- obs_by_locID[ind.pid, "upDist"]


        if(locID.ob != locID.old) {

            junk <- SSN:::get.rid.fc(pred_by_locID[dup.pred,"binaryID"],
                                     obs_by_locID$binaryID[ind.pid])

            ob.j <- getObsPredsRelationshipsDF(ssn, junk, dup.pred, pred_by_locID,
                                     obs_by_locID[ind.pid,], bin.table, pred.num)

            ob.j <-ob.j[pred_reordering,]

            ##obs fills in by column because it's working between obs to obs.
            ind.fc<-ob.j$fc==1

            dist.a <- ifelse(ind.fc, ob.j$upDist.j-upDist.ob, ob.j$upDist.j - ob.j$juncDist)
            ## WE CHANGE THIS SO THAT WE'RE ADDING COLUMNS (TRANSPOSE IT LATER FOR USE)
            col.ind <- colnames(current_distance_matrix_a) == as.character(pid.ob)
            current_distance_matrix_a[,col.ind] <- ifelse(dist.a<0, 0, dist.a)

            dist.b <- ifelse(ind.fc, upDist.ob - ob.j$upDist.j,
                           upDist.ob - ob.j$juncDist)
            current_distance_matrix_b[,col.ind] <- ifelse(dist.b<0, 0, dist.b)

            locID.old <- locID.ob

        } else {
            col.ind <- colnames(current_distance_matrix_a) == as.character(pid.ob)
            current_distance_matrix_a[,col.ind]<- ifelse(dist.a<0, 0, dist.a)

            current_distance_matrix_b[,col.ind] <- ifelse(dist.b<0, 0, dist.b)
            ##locID.old <- locID.ob
        }

   }
    close(current_distance_matrix_a)
    close(current_distance_matrix_b)
}

