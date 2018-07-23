#########################translation_csvs############################

#8/16/2017

#Renske van Raaphorst

#####################################################################

##MicrobeJ

extr.MicrobeJMESH <- function(dataloc){
  MESH <- read.table(dataloc, header=T, sep=",")
  meshL <- list()
  meshL$cellList <- MESH
  IDlist <- data.frame(NAME.id = unique(out$cellList$Mesh$NAME.id), cell = c(1:length(unique(out$cellList$Mesh$NAME.id))))
  MESH <- merge(MESH, IDlist)
  MESH$frame <- MESH$POSITION
  MESH$cellID <- MESH$NAME.id
  MESH <- MESH[,c("X", "Y", "cell", "frame", "cellID")]
  MESH$num <- c(1:nrow(MESH))
  meshL$meshList <- MESH
  return(meshL)
}


extr.MicrobeJSpots <- function(spotloc ,mag){
  SPOTS <- read.csv(spotloc, header=F)
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
extr.MicrobeJ <- function(dataloc, spotloc, mag){
  if(missing(mag)){
    stop("Magnification conversion needed for proper intercellular measurements! MicrobeJ already converted the spot localizations for you, but not the contours.")
  }
  outlist <- list()
  if(missing(dataloc)!=T){
    MESH <- extr.MicrobeJMESH(dataloc)$meshList
    cellList <- extr.MicrobeJMESH(dataloc)$cellList
    if(missing(spotloc)==T){
      outlist$cellList <- cellList
    }
  }
  if(missing(spotloc)!=T){
    SPOTS <- extr.MicrobeJSpots(spotloc ,mag)$spotList
    cellList2 <- extr.MicrobeJSpots(spotloc , mag)$cellList
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
    spot_mesh <- mergeframes(listbox$REP, listbox$Mfull, mag)
    outlist$spots_relative <- spot_mesh
    outlist$pixel2um <- unlist(get(magnificationList, envir=magEnv)[mag])
    outlist$mesh <- listbox$Mfull
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
  spotminimal <- SPOTS[,c("x","y","frame")]
  listout <- list()
  listout$cellList <- SPOTS
  listout$spotframe <- spotminimal
  return(listout)
}


#ObjectJ - mostly manual entry.
#' @export
extr.ObjectJ <- function(dataloc, spots="X", constr, xcolumns, ycolumns, xconstr, yconstr, spotcol,spotloc){
  if(substr(dataloc, nchar(dataloc)-3, nchar(dataloc))==".txt"){
    OJ <- read.table(dataloc, header=T, sep="\t")
  }
  if(substr(dataloc, nchar(dataloc)-3, nchar(dataloc))==".csv"){
    OJ <- read.csv(dataloc, header=T)
  }

  OJ$cell <- OJ$n
  OJ$max.length <- round(OJ$Axis,2)
  OJ$max.width <- OJ$Dia
  OJ$n <- NULL
  OJ$Dia <- NULL


  if(spots=="X"&missing(spotloc)){
    spots <- readline("Does the file contain spot coordinates? Y/N: ")
  }

  if(spots=="Y"|spots=="y"){
    if(missing(xcolumns)|missing(ycolumns)){
    print("Give the name(s) of the column(s) containing the x localization of the spots.\nPress <enter> between names, end with #<enter>")
    xcolumns <- choosecolumns()
    print("Give the name(s) of the column(s) containing the y localization of the spots.\nPress <enter> between names, end with #<enter>" )
    ycolumns <- choosecolumns()
    }
    n <- length(xcolumns)
    if(n!=length(ycolumns)){
      stop("Number of x and y columns does not match.")
    }

    OJ[,xcolumns[1]][is.na(OJ[,xcolumns[1]])] <- 0
    OJ <- tidyr::gather(OJ, "xspot", "l", xcolumns[1:n], na.rm=T)
    OJ$l[OJ$l==0] <- NA
    OJ[,ycolumns[1]][is.na(OJ[,ycolumns[1]])] <- 0
    OJ <- tidyr::gather(OJ, "yspot", "d", ycolumns[1:n], na.rm=T)
    OJ$d[OJ$d==0] <- NA
    if(missing(spotcol)){
    spotcol <- readline("If included, give the name of the column containing the number of spots per cell.\nOtherwise, press #<enter>")
    }
    if(spotcol!="#"){
      colnames(OJ)[colnames(OJ)==spotcol] <- "spotn"

    }
  }
  if(missing(spotloc)!=T){
    if(substr(spotloc, nchar(spotloc)-3, nchar(spotloc))==".csv"){
      SF <- read.csv(spotloc, header=T)
    }
    if(substr(spotloc, nchar(spotloc)-3, nchar(spotloc))==".txt"){
      SF <- read.table(spotloc, header=T, sep="\t")
    }
    SF <- spotrExtractSpotsObjectJ(SF)
    OJ <- merge(OJ, SF, all=T)
  }


  if(missing(constr)){
  constr <- readline("Does the file contain constriction point coordinates? Y/N: ")
  }
  if(constr=="Y"|constr=="y"){
    if(missing(xconstr)|missing(yconstr)){
    print("Give the name(s) of the column(s) containing the distance of the constriction point to the Axis.\nPress <enter> between names, end with #<enter>")
    xconstr <- choosecolumns()
    print("Give the name(s) of the column(s) containing the distance of the constriction point to the Axis.\nPress <enter> between names, end with #<enter>")
    yconstr <- choosecolumns()
    }
    i <- length(xconstr)
    if(i!=length(yconstr)){
      stop("Number of x and y columns does not match.")
    }

    OJ[,xconstr[1]][is.na(OJ[,xconstr[1]])] <- 0
    OJ <- tidyr::gather(OJ, "xconstr", "LC", xconstr[1:i], na.rm=T)
    OJ$LC[OJ$LC==0] <- NA
    OJ[,yconstr[1]][is.na(OJ[,yconstr[1]])] <- 0
    OJ <- tidyr::gather(OJ, "yconstr", "DC", yconstr[1:i], na.rm=T)
    OJ$DC[OJ$DC==0] <- NA
  }

  OJ$max.length <- OJ$Axis
  OJ$Axis <- NULL
  return(OJ)

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
