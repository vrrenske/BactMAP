####transformation of OUFTI into shape which works for MT#################

##################source file cellList####################################

sprOuftipickfile <- function(cellListfile){
#pick a file
  if(missing(cellListfile)){
  cellListfile = file.choose()
  }
  cellList <- utils::read.delim(cellListfile, header=FALSE)
  #listname=basename(cellListfile)
#find line where the relevant meshdata starts
  startnum <- which((substring(as.character(cellList$V1), 1, 5) == "frame"))
#get column names
  listcols <- as.character(cellList[startnum,])
#to get rid of the ";" in the end of the line and just in case of all of them..
  listcols <- gsub(";", "", listcols)
  listcols <- strsplit(listcols, ",", fixed=T)
#remove annoying bug of missing ", " in end when there is an object.
#if(U=="n"){
  #listcols <- listcols[[1]][-13]
#}else(
  listcols <- listcols[[1]]
#)

  cellList <- cellList[(startnum+1):nrow(cellList),]
  cellList <- data.frame(cellList)
  cellList <- suppressWarnings(data.frame(do.call("rbind", strsplit(as.character(cellList$cellList), ",", fixed=T))))
  colnames(cellList) <- listcols
#remove annoying bug of missing ", " in end when there is an object.
#needed with old version of OUFTI
  if("#1"%in%cellList$cellId){
    cellList$cellId[(as.character(cellList$cellId))==as.character(cellList$frameNumber)] <- NA
    cellList$cellId <- as.numeric(cellList$cellId)
    cellList$cellId[is.na(cellList$cellId)] <- cellList$objects[is.na(cellList$cellId)]
    cellList$objects[grepl(" ", cellList$objects)==FALSE] <- NA
  }
#get rid of hashtag.
  cellList$frameNumber <- cellList$frameNumber
  cellList$frameNumber <- gsub("#", "", cellList$frameNumber)
  cellList$frameNumber <- as.numeric(cellList$frameNumber)
  cellList <- cellList[cellList$cellId>0,]
  return(cellList)
  }

#######function to extract and reshape meshes out of .csv file; name .csv file = cellList

sprOuftimesh <- function(cellList){

  mesh <- cellList[c("frameNumber", "cellId", "length", "mesh")]
  mesh$frame <- mesh$frameNumber
  mesh$cell <- as.numeric(mesh$cellId)
  mesh$max.length <- as.numeric(mesh$length)
  mesh$length <- NULL
  mesh$frameNumber <- NULL
  mesh$cellId <- NULL
  x <- 1
  for(n in 1:nrow(mesh)){
    if(!is.na(mesh$max.length[n])){
      dat <- data.frame(t(do.call('rbind', strsplit(as.character(mesh$mesh[n]), ';', fixed=TRUE))))
      colnames(dat) <- "cord"
      dat$cord <- as.character(dat$cord)
      for(z in 1:nrow(dat)){
        if(substring(dat$cord[z], 1, 1) == " "){
          dat$cord[z] <- substring(dat$cord[z], 2)
        }
      }
      dat <- data.frame(t(do.call('rbind', strsplit(as.character(dat$cord), ' ', fixed=TRUE))))
      colnames(dat) <- c("x0", "y0", "x1", "y1")
      dat$num <- 0:(nrow(dat)-1)
      dat$x0 <- as.numeric(as.character(dat$x0))
      dat$y0 <- as.numeric(as.character(dat$y0))
      dat$x1 <- as.numeric(as.character(dat$x1))
      dat$y1 <- as.numeric(as.character(dat$y1))
      xdist <- dat$x0 - dat$x1
      ydist <- dat$y0 - dat$y1
      widths <- sqrt(xdist^2 + ydist^2)
      max_width <- max(widths)
      dat$max.width <- max_width
      dat$frame <- mesh$frame[n]
      dat$cell <- mesh$cell[n]
      dat$max.length <- mesh$max.length[n]
      dat$length <- dat$max.length/max(dat$num)*dat$num
      #print(mesh$cell[n])
      if(x==1){
        MESH <- dat
      }
      if(x>1){MESH <- rbind(dat,MESH)}
      x <- x + 1
    }
  }
  minPointTwo <- stats::quantile(MESH$max.length, probs=seq(0, 0.02, 0.02))[[2]]
  return(MESH[MESH$max.length>minPointTwo,])
}


######################getting the meshfile "M"########################################################
sprOuftiM <- function(cellList, MESH){
  M <- cellList[c("frameNumber", "cellId", "ancestors", "descendants", "length", "area")]
  M$frame <- as.numeric(as.character(M$frameNumber))
  M$cell <- as.numeric(as.character(M$cellId))
  M$frameNumber <- NULL
  M$cellId <- NULL
  M$length <- as.numeric(as.character(M$length))
  M$area <- as.numeric(as.character(M$area))
  M <- merge(M, MESH[c("frame", "cell", "max.width")], all=T)
  M <- M[M$max.width>1&!is.na(M$frame),]
}



######################spot file  ######################################################
sprOuftispot <- function(cellList){
  REP <- cellList[c("frameNumber", "cellId", "spots")]
  REP$frame <- as.numeric(REP$frameNumber)
  REP$cell <- as.numeric(REP$cellId)
  REP <- REP[!is.na(REP$cell),]
  REP$frameNumber <- NULL
  REP$cellId <- NULL
  u <- 0
  for(n in 1:nrow(REP)){
    if(REP$spots[n] != " "){
      dat3 <- data.frame(t(do.call('rbind', strsplit((REP$spots[n]), ';', fixed=TRUE))))
      colnames(dat3) <- "sp"
      dat3$cnt <- 1:nrow(dat3)
      dat3 <- dat3[1:5,]
      dat3$sp <- as.character(dat3$sp)
      dat3$sp <- gsub("-", " -", dat3$sp)
      for(z in 1:nrow(dat3)){
        if(substring((dat3$sp[z]), 1, 1)== " "){
          dat3$sp[z] <- substring((dat3$sp[z]), 2)
        }
      }
      dat3 <- data.frame(t(do.call('rbind', strsplit((dat3$sp), " ", fixed=TRUE))))
      colnames(dat3) <- c("l", "d", "x", "y", "positions")
      dat3$frame <- REP$frame[n]
      dat3$cell <- REP$cell[n]
      if(nrow(dat3)>0){
        u <- u+1
      }
      if(u==1){
        repN <- dat3
      } else { repN <- rbind(repN,dat3)}
    }
  }
  repN$l <- as.numeric(as.character(repN$l))
  repN$d <- as.numeric(as.character(repN$d))
  repN$x <- as.numeric(as.character(repN$x))
  repN$y <- as.numeric(as.character(repN$y))
  return(repN)
}


  sprOuftiobject <- function(cellList){
    OBJ <- cellList[c("frameNumber", "cellId", "objects")]
    OBJ$frame <- as.numeric(as.character(OBJ$frameNumber))
    OBJ$cell <- as.numeric(as.character(OBJ$cellId))
    OBJ <- OBJ[!is.na(OBJ$cell),]
    OBJ$frameNumber <- NULL
    OBJ$cellId <- NULL
    OBJ$objects <- as.character(OBJ$objects)
    for(n in 1:nrow(OBJ)){
      if(!is.na(OBJ$objects[n])){
         if(OBJ$objects[n] != " "){
             dat3 <- data.frame(t(do.call('rbind', strsplit((OBJ$objects[n]), ';', fixed=TRUE))))
             colnames(dat3) <- "ob"
             dat3$cnt <- 1:nrow(dat3)
             dat3 <- dat3[dat3$ob != " ",]
             dat3$ob <- as.character(dat3$ob)
             dat3 <- data.frame(t(do.call('rbind', strsplit((dat3$ob), " ", fixed=TRUE))))
             numobjects <- c(1:(ncol(dat3)/5))
             selectobj <- (numobjects*5)-4
             selectobj <- append(selectobj, selectobj+1)
             selectobj <- selectobj[order(selectobj)]
             dat3 <- dat3[,selectobj]
             dat3 <- dat3[as.character(dat3$X1)!="",]
             colnames(dat3) <- rep(c("ob_x", "ob_y"), (ncol(dat3)/2))


             for(i in numobjects){
               datpart <- dat3[,(2*i-1):(2*i)]
               datpart <- unique(datpart)
               datpart$ob_x <- as.numeric(datpart$ob_x)
               datpart$ob_y <- as.numeric(datpart$ob_y)
               datpart <- datpart[!is.na(datpart$ob_x)&!is.na(datpart$ob_y),]
               datpart$obnum <- i
               datpart$obpath <- c(1:nrow(datpart))
               bbox <- shotGroups::getMinBBox(data.frame(x= datpart$ob_x, y=datpart$ob_y))
               datpart$obwidth <- min(c(bbox$height, bbox$width))
               datpart$oblength <- max(c(bbox$height, bbox$width))
               datpart$obarea <- sp::Polygon(data.frame(x=datpart$ob_x, y=datpart$ob_y))@area
               if(i==1){
                 datfull <- datpart
               }
               if(i>1){
                 datfull <- rbind(datfull, datpart)
               }
             }
             datfull$frame <- OBJ$frame[n]
             datfull$cell <- OBJ$cell[n]
             datfull$obID <- paste(datfull$cell, datfull$frame, datfull$obnum, sep="_")
             if(exists("OBJn")==FALSE){
                OBJn <- datfull
             } else { OBJn <- rbind(OBJn,datfull)}
             }
      }
    }




    return(OBJn)
  }

#
extr_OuftiCSV <- function(dataloc){
  if(missing(dataloc)){
    dataloc <- file.choose()
  }
  C <- sprOuftipickfile(dataloc)
  MESH <- sprOuftimesh(C)
  outlist <- list()
  outlist$cellList <- C
  outlist$mesh <- MESH
  if(length(unique(C$spots))>1){
    SPOTS <- sprOuftispot(C)
    outlist$spotframe <- SPOTS
  }
  if(length(unique(C$objects))>1){
    OBJ <- suppressWarnings(sprOuftiobject(C))
    outlist$objectframe <- OBJ
  }
  return(outlist)
}

