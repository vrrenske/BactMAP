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
    message("found spot data. saving as 'spotframe'.")
    SPOTS <- sprOuftispot(C)
    outlist$spotframe <- SPOTS
  }
  if(length(unique(C$objects))>1){
    message("found object data. saving as 'objectframe'.")
    OBJ <- suppressWarnings(sprOuftiobject2(C))
    outlist$objectframe <- OBJ
  }
  return(outlist)
}


sprOuftiobject2 <- function(cellList){
  OBJ <- cellList[c("frameNumber", "cellId", "objects")] %>%
    dplyr::transmute(frame = as.numeric(as.character(.data$frameNumber)),
                     cell = as.numeric(as.character(.data$cellId)),
                     objects = as.character(.data$objects)) %>%
    dplyr::filter(!is.na(.data$cell),
                  !is.na(.data$objects),
                  .data$objects != " ")

  OBJ <- lapply(c(1:nrow(OBJ)), function(x) perRowObject(OBJ[x,])) %>%
    dplyr::bind_rows()
  return(OBJ)
}


perRowObject <- function(obRow){
  dat3 <- data.frame(t(do.call('rbind', strsplit((obRow$objects), ';', fixed=TRUE))))
  colnames(dat3) <- "ob"
  dat3$cnt <- 1:nrow(dat3)
  dat3 <- dat3 %>%
    dplyr::filter(.data$ob != " ") %>%
    dplyr::mutate(ob = as.character(.data$ob))
  dat3 <- data.frame(t(do.call('rbind', strsplit((dat3$ob), " ", fixed=TRUE))))
  numobjects <- c(1:(ncol(dat3)/5))
  selectobj <- (numobjects*5)-4
  selectobj <- append(selectobj, selectobj+1)
  selectobj <- selectobj[order(selectobj)]
  dat3 <- dat3[,selectobj]
  colnames(dat3) <- rep(c("ob_x", "ob_y"), (ncol(dat3)/2))
  datfull <- lapply(numobjects, function(x) perObject(x, dat3)) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(frame = obRow$frame,
                  cell = obRow$cell,
                  obID = paste(.data$frame, .data$cell, .data$obnum, sep="_"))
  return(datfull)

}


perObject <- function(i, dat3){
  datpart <- dat3[,(2*i-1):(2*i)] %>%
    dplyr::distinct() %>%
    dplyr::mutate(ob_x = as.numeric(.data$ob_x),
                  ob_y = as.numeric(.data$ob_y)) %>%
    dplyr::filter(!is.na(.data$ob_x)) %>%
    dplyr::filter(!is.na(.data$ob_y)) %>%
    dplyr::mutate(obnum = i,
                  obpath = dplyr::row_number())
  bbox <- shotGroups::getMinBBox(data.frame(point.x= datpart$ob_x, point.y=datpart$ob_y))
  datpart <- datpart %>%
    mutate(obwidth = min(c(bbox$height, bbox$width)),
           oblength = max(c(bbox$height, bbox$width)),
           obarea = sp::Polygon(data.frame(x=datpart$ob_x, y=datpart$ob_y))@area)
  return(datpart)
}



