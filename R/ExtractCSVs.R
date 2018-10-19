#########################translation_csvs############################

#8/16/2017

#Renske van Raaphorst

#####################################################################

##MicrobeJ

extr.MicrobeJMESH <- function(dataloc, sep=","){
  MESH <- read.csv(dataloc, header=T, sep=sep)
  meshL <- list()
  meshL$cellList <- MESH
  IDlist <- data.frame(NAME.id = unique(MESH$NAME.id), cell = c(1:length(unique(MESH$NAME.id))))
  MESH <- merge(MESH, IDlist)
  MESH$frame <- MESH$POSITION
  MESH$cellID <- MESH$NAME.id
  MESH <- MESH[,c("X", "Y", "cell", "frame", "cellID")]
  MESH$num <- c(1:nrow(MESH))
  meshL$meshList <- MESH
  return(meshL)
}


extr.MicrobeJSpots <- function(spotloc ,mag, sep=","){
  SPOTS <- read.csv(spotloc, header=F, sep=sep)
  colnamesspot <- SPOTS[1,]
  colnamesspot <- colnamesspot[!is.na(colnamesspot)]
  SPOTS <- SPOTS[-1,]

  #comma problem
  if(ncol(SPOTS)>length(colnamesspot)){
    a <- which(lapply(1:ncol(SPOTS), function(x) is.logical(SPOTS[,x]))==TRUE)
    SPOTS[,a] <- NULL
    u <- t(SPOTS[1,])
    y <- which(substr(u, 1, 1) == "(")
    for(z in y){
      SPOTS[,z] <- paste(SPOTS[,z], SPOTS[,(z+1)], sep=";")
    }
    SPOTS[,(y+1)] <- NULL
  }
  colnames(SPOTS) <- colnamesspot

  spotL <- list()
  spotL$cellList <- SPOTS
 #SPOTS$x <- t(data.frame(strsplit(stringr::str_sub(SPOTS$LOCATION,2,-2),";"))[1,])
  #SPOTS$y <- t(data.frame(strsplit(stringr::str_sub(SPOTS$LOCATION,2,-2),";"))[2,])
  SPOTS$x <- as.numeric(as.character(SPOTS$LOCATION.x))/unlist(get(magnificationList, envir=magEnv)[mag])
  SPOTS$y <- as.numeric(as.character(SPOTS$LOCATION.y))/unlist(get(magnificationList, envir=magEnv)[mag])
  SPOTS$cellID <- SPOTS$PARENT.id
  SPOTS$frame <- SPOTS$POSITION.slice
  SPOTS <- SPOTS[,c("x", "y", "cellID", "frame")]
  spotL$spotList <- SPOTS
  return(spotL)
}

#' @export
extr.MicrobeJ <- function(dataloc, spotloc, mag, sepspot=",", sepmesh=","){
  if(missing(spotloc)!=T&missing(dataloc)!=T&missing(mag)){
    stop("Magnification conversion needed for proper intercellular measurements! MicrobeJ already converted the spot localizations for you, but not the contours.")
  }
  if(missing(mag)!=T&is.numeric(unlist(get(magnificationList,envir=magEnv)[mag]))==FALSE){
    stop("Magnification conversion factor not recognized. Please use addPixels2um('pixelName', pixelsize) to add your conversion factor")
  }
  outlist <- list()
  if(missing(dataloc)!=T){
    MESHout <- extr.MicrobeJMESH(dataloc, sepmesh)
    cellList <- MESHout$cellList
    MESH <- MESHout$meshList
    if(missing(spotloc)==T){
      cellList$num <- cellList$INDEX
      outlist$cellList <- cellList
      outlist$mesh <- MESH
    }
  }
  if(missing(spotloc)!=T){
    if(missing(mag)){
      mag <- "No_PixelCorrection"
    }
    spotsout <- extr.MicrobeJSpots(spotloc ,mag, sep=sepspot)
    SPOTS <- spotsout$spotList
    cellList2 <- spotsout$cellList
    if(missing(dataloc)==T){
      outlist$cellList <- cellList2
      IDframe <- data.frame(cellID = unique(SPOTS$cellID), cell=c(1:length(unique(SPOTS$cellID))))
      outlist$spotframe <- SPOTS
    }
  }
  if(missing(spotloc)!=T&missing(dataloc)!=T){
    IDframe <- unique(MESH[,c("cellID", "cell")])
    SPOTS <- merge(SPOTS, IDframe)
    listbox <- spotsInBox(SPOTS, MESH)
    outlist$spotframe <- SPOTS
    if(missing(mag)){
      mag <- "No_PixelCorrection"
    }
    spot_mesh <- mergeframes(listbox$spots_relative, listbox$mesh, mag)
    outlist$spots_relative <- spot_mesh
    outlist$pixel2um <- unlist(get(magnificationList, envir=magEnv)[mag])
    outlist$mesh <- listbox$mesh
    cellList3 <- list()
    cellList3$Mesh <- cellList
    cellList3$Spots <- cellList2
    outlist$cellList <- cellList3
  }
  return(outlist)
}

##ISBatch
#' @export
extr.ISBatch <- function(dataloc){
  if(substr(dataloc, nchar(dataloc)-3, nchar(dataloc))==".txt"){
    SPOTS <- read.table(dataloc, header=T, sep="\t")
  }
  if(substr(dataloc, nchar(dataloc)-3, nchar(dataloc))==".csv"){
    SPOTS <- read.csv(dataloc, header=T)
  }
  SPOTS$frame  <- SPOTS$slice
  SPOTS$slice <- NULL
  if("trajectory"%in%colnames(SPOTS)){
    spotminimal <- SPOTS[,c("x", "y", "frame", "displacement_sq", "trajectory", "trajectory_length")]
  }else{
    spotminimal <- SPOTS[,c("x","y","frame")]
  }

  listout <- list()
  listout$cellList <- SPOTS
  listout$spotframe <- spotminimal
  return(listout)
}


#ObjectJ - mostly manual entry.
#' @export
extr.ObjectJ <- function(dataloc, mag="No_PixelCorrection"){
  #prepare list for output
  outlist <- list()
  #read out .txt file
  oj <- read.table(dataloc, header=T, sep="\t")

  #get a frame count from 1:x
  celldats <- data.frame(Frame=unique(oj$Frame))

  celldats <- data.frame(Frame = celldats[order(celldats$Frame),])

  celldats$f <- c(1:nrow(celldats))

  oj <- merge(celldats, oj)

  oj$frame <- oj$f

  oj$f <- NULL
  oj$Frame <- NULL
  #separate chain dataframe
  oj_chain <- oj[!is.na(oj$ChainAxis),]

  oj_chain$chainID <- oj_chain$n

  if(length(unique(oj_chain$GFPfluor)>1)){
    oj_chain <- oj_chain[,c("chainID", "ChainAxis", "ChainDia", "frame", "GFPfluor", "GFPMax", "GFPMean", "GFPMid")]
  }
  else{
    oj_chain <- oj_chain[,c("chainID", "ChainAxis", "ChainDia", "frame")]
  }

  #save full file as cellList
  outlist$cellList <- oj

  #save chainfile as chainframe
  outlist$chainframe <- oj_chain

  #take now only the cells
  oj <- oj[!is.na(oj$CellAxis),]

  #combine all x/y axes of the bounding boxes around the cells.
  oj2 <- reshape2::melt(oj, measure.vars=c("X1","X2","X3","X4","X5","X6","X7","X8", "X9","X10","X11"), value.name="Xum",
                        variable.name="index_X")

  oj2 <- reshape2::melt(oj2, measure.vars=c("Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8", "Y9","Y10","Y11"), value.name="Yum",
                        variable.name="index_Y")

  oj2$index_X <- lapply(oj2$index_X, function(x) strsplit(as.character(x), "X")[[1]][[2]])

  oj2$index_Y <- lapply(oj2$index_Y, function(x) strsplit(as.character(x), "Y")[[1]][[2]])

  oj2$index_X <- as.numeric(oj2$index_X)

  oj2$index_Y <- as.numeric(oj2$index_Y)

  oj2 <- oj2[oj2$index_X==oj2$index_Y,]

  oj2$num <- oj2$index_X
  oj2$max.length <- oj2$CellAxis
  oj2$cell <- oj2$n
  oj2$chainID <- oj2$ParentID

  oj3 <- unique(oj2[,c("cell", "max.length", "chainID", "GFPfluor", "GFPMean", "GFPMax", "GFPMid", "allSepta", "manSepta", "manTrace", "Thr", "startX", "startY", "barX", "barY", "barW", "barH", "SeptaPos", "frame")])

  outlist$GFPframe <- oj3

  oj2 <- oj2[,c("cell", "frame","max.length", "chainID", "Xum", "Yum", "num")]
  p2um <- as.numeric(get(magnificationList, envir=magEnv)[mag])

  oj2$X <- oj2$Xum/p2um
  oj2$Y <- oj2$Yum/p2um

  oj2 <- oj2[order(oj2$frame, oj2$cell, oj2$num),]
  outlist$mesh <- oj2

  return(outlist)
}


choosecolumns <- function(){
  u <- "i"
  columnamelist <- list()
  while(u!="#"){
    u <- readline()
    if(u!="#"){
      columnamelist <- append(columnamelist, u)
    }
  }
  return(unlist(columnamelist))
}

spotrExtractSpotsObjectJ <- function(SF){
  nmax <- (ncol(SF)-4)/3
  for(n in 1:nmax){
    SN <- SF[,c("n", "len", "seq", "mirr", paste("type",n, sep=""), paste("pp",n,sep=""), paste("off",n, sep=""))]
    colnames(SN) <- c("cell", "length", "spotorder", "mirr", "type", "l", "d")
    if(n==1){Sframe <- SN}
    if(n>1){
      SN <-SN[!is.na(SN$L),]
      Sframe <- rbind(Sframe,SN)}
  }
  Sframe <- Sframe[order(Sframe$cell, Sframe$L),]
  return(Sframe)
}
