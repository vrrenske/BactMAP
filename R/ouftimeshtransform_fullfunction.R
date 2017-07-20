####transformation of OUFTI into shape which works for MT#################

##################source file cellList####################################

sprOuftipickfile <- function(){
#pick a file
  cellListfile = file.choose()
  cellList <-read.delim(cellListfile, header=FALSE)
  listname=basename(cellListfile)
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
  cellList <- data.frame(do.call("rbind", strsplit(as.character(cellList$cellList), ",", fixed=T)))
  colnames(cellList) <- listcols
#remove annoying bug of missing ", " in end when there is an object.
#needed with old version of OUFTI
  #if(U=="n"){
    #cellList$cellId[(as.character(cellList$cellId))==as.character(cellList$frameNumber)] <- NA
   # cellList$cellId <- as.numeric(as.character(cellList$cellId))
   # cellList$objects <- as.numeric(as.character(cellList$objects))
   # cellList$cellId[is.na(cellList$cellId)] <- cellList$objects[is.na(cellList$cellId)]
   # cellList$objects[cellList$objects==cellList$cellId] <- NA
  #}
#get rid of hashtag.
  cellList$frameNumber <- as.character(cellList$frameNumber)
  cellList$frameNumber <- gsub("#", "", cellList$frameNumber)
  cellList$frameNumber <- as.numeric(cellList$frameNumber)
  return(cellList)
  }

#######function to extract and reshape meshes out of .csv file; name .csv file = cellList

sprOuftimesh <- function(cellList){

  mesh <- cellList[c("frameNumber", "cellId", "length", "mesh")]
  mesh$frame <- as.numeric(as.character(mesh$frameNumber))
  mesh$cell <- as.numeric(as.character(mesh$cellId))
  mesh$max_length <- as.numeric(as.character(mesh$length))
  mesh$length <- NULL
  mesh$frameNumber <- NULL
  mesh$cellId <- NULL
  x <- 1
  for(n in 1:nrow(mesh)){
    if(!is.na(mesh$max_length[n])){
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
      dat$max_length <- mesh$max_length[n]
      dat$length <- dat$max_length/max(dat$num)*dat$num
      print(mesh$cell[n])
      if(x==1){
        MESH <- dat
      }
      if(x>1){MESH <- rbind(dat,MESH)}
      x <- x + 1
    }
  }
  return(MESH[MESH$max_length>1,])
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
  REP$frame <- as.numeric(as.character(REP$frameNumber))
  REP$cell <- as.numeric(as.character(REP$cellId))
  REP <- REP[!is.na(REP$cell),]
  REP$frameNumber <- NULL
  REP$cellId <- NULL
  REP$spots <- as.character(REP$spots)
  u <- 0
  for(n in 1:nrow(REP)){
    if(REP$spots[n] != " "){
      dat3 <- data.frame(t(do.call('rbind', strsplit((REP$spots[n]), ';', fixed=TRUE))))
      colnames(dat3) <- "sp"
      dat3$cnt <- 1:nrow(dat3)
      dat3 <- dat3[dat3$sp != " ",]
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
               datpart$obnum <- i
               if(i==1){
                 datfull <- datpart
               }
               if(i>1){
                 datfull <- rbind(datfull, datpart)
               }
             }
             datfull$frame <- OBJ$frame[n]
             datfull$cell <- OBJ$cell[n]
             if(exists("OBJn")==FALSE){
                OBJn <- datfull
             } else { OBJn <- rbind(OBJn,datfull)}
             }
    }
    OBJn$ob_x <- as.numeric(as.character(OBJn$ob_x))
    OBJn$ob_y <- as.numeric(as.character(OBJn$ob_y))
    return(OBJn)
  }

spotrTransform <- function(sourceprogram = "Oufti", outputtype = "M"){
  if(sourceprogram=="Oufti"){
    C <- sprOuftipickfile()
    MESH <- sprOuftimesh(C)
    if("S"%in%outputtype){
      REP <- sprOuftispot(C)
      return(REP)
    }
    if("O"%in%outputtype){
      OBJ <- sprOuftiobject(C)
      return(OBJ)
    }
    else(return(MESH))
  }
}

