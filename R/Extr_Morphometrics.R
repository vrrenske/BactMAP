##8-11-2017

##Renske van Raaphorst

##Extracting Morphometrics Data


##########################################################


##Dependencies:
    #library(R.matlab)

#' @export
extr_Morphometrics_cellList <- function(morphpath){
  morphdata <- R.matlab::readMat(morphpath)
  framenum <- 2*(1:(length(morphdata$frame)/2))
  for(u in framenum){
    cellListf <- data.frame(t(data.frame(morphdata$frame[u])))
    if(nrow(cellListf)!=0){
      cellListf$frame <- u/2
      if(u==2){
        cellList <- cellListf
      }
      else{
        cellList <- rbind(cellList, cellListf)
      }
    }
  }
  return(cellList)
}


spotrExtractMorphMESH <- function(cellList){

  for(n in 1:nrow(cellList)){
    meshcell <- data.frame(cellList$Xcont[n], cellList$Ycont[n], cellList$cellID[n], cellList$area[n],
                           cellList$pole1[n], cellList$pole2[n], cellList$frame[n])
    colnames(meshcell) <- c("X", "Y", "cell", "area", "pole1", "pole2", "frame")
    if("length"%in%colnames(cellList)){
      meshcell$max.length <- cellList$length[n]
      meshcell$max.width <- cellList$width[n]
    }

    meshcell$numpoint <- 1:nrow(meshcell)
    if(n==1){
      MESH <- meshcell
    }
    else{
      MESH <- rbind(MESH, meshcell)
    }
  }
  if("length"%in%colnames(cellList)!=T){
    M1 <- MESH[MESH$numpoint==MESH$pole1,]
    M2 <- MESH[MESH$numpoint==MESH$pole2,]
    M1 <- M1[order(M1$cell,M1$frame),]
    M2 <- M2[order(M2$cell, M2$frame),]
    M1$max.length <- sqrt((M1$X-M2$X)^2+(M1$Y-M2$Y)^2)
    MESH <- merge(MESH, M1[,c("cell", "frame", "max.length")])
  }
  return(MESH)

}

#' @export
extr_Morphometrics <- function(morphpath, mag, turncells = TRUE){
  C <- extr_Morphometrics_cellList(morphpath)
  M <- spotrExtractMorphMESH(C)
  listM <- list()
  listM$cellList <- C
  if(turncells==TRUE){
    listM$mesh <- meshTurn(M)
  }
  if(turncells==FALSE){
    listM$mesh <- M
  }

  if(!missing(mag)){
    if(is.numeric(unlist(get(magnificationList,envir=magEnv)[mag]))==FALSE){
      stop("Magnification conversion factor not recognized. Please use addPixels2um('pixelName', pixelsize) to add your conversion factor")
    }
    if(turncells==TRUE){
      listM$mesh$Xrotum <- listM$mesh$X_rot * unlist(get(magnificationList, envir=magEnv)[mag])
      listM$mesh$Yrotum <- listM$mesh$Y_rot * unlist(get(magnificationList, envir=magEnv)[mag])
      listM$mesh$maxwum <- listM$mesh$max.width * unlist(get(magnificationList, envir=magEnv)[mag])
    }
    listM$mesh$max_um <- listM$mesh$max.length * unlist(get(magnificationList, envir=magEnv)[mag])
    listM$mesh$area_um <- listM$mesh$area *  unlist(get(magnificationList, envir=magEnv)[mag])^2
    listM$pixel2um <- unlist(get(magnificationList, envir=magEnv)[mag])
  }

  return(listM)
}

