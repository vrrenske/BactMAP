##Renske van Raaphorst
##7/14/2017

##Helperfunctions: functions necessary to make other functions working properly.
##other package dependencies:

#merge data functions
#library(ggplot2)
#library(ggthemes)
#library(gridExtra)

#' @export
#'
#' @title Function to get a dataframe per cell out of your mesh, spots, object or tiff dataset.
#'
#' \code{onePerCell()} returns a dataframe with one table per cell, simplified especially for people who wish
#' to make their own plots (e.g. histograms, etc) based on the summarized cell data. Automatically, it groups
#' the input dataset by unique cell (cell & frame), but it is possible to use this function with other column names.
#'
#' @param inputframe the dataset (in form of a dataframe) to be summarized
#' @param by the grouping value. default = c("cell", "frame")
#' @return the function will return a dataframe (tibble) with the data summarized by the grouping value \code{by}.
#' it will do this by returning the common values per cell (mean intensity, etc), removing values that are useless
#' in one-column context (x/y coordinates of outlines), and calculating the mean, median and standard deviation of
#' numerical values it doesn't recognize. if there are character values with different values per grouped variable,
#' the function will take ONLY THE FIRST VARIABLE.
#'
#' @examples
#' \dontrun{
#' df <- data.frame("cell"= c(rep(1, 5), rep(2,5)), "frame"=rep(1, 10), "channel"="GFP", "intensity"=c(1,3,4,5,6,1,2,7,3,10))
#'
#' onePerCell(df)
#' }
onePerCell <- function(inputframe, by=c("cell","frame")){

  #remove columns with x/y coordinates
  cols <- colnames(inputframe)
  cols <- cols[!cols%in%c("X",
                          "Y",
                          "ob_out_x",
                          "ob_out_y",
                          "num",
                          "obpath",
                          "X_rot",
                          "Y_rot",
                          "Xrot_micron",
                          "Yrot_micron"
  )]

  inputframe <- inputframe %>% select(cols)

  for(n in by){
    inputframe <- inputframe %>% dplyr::group_by(.data[[n]], add=TRUE)
  }



  #summarize
  inputframe_num <- inputframe %>%
    dplyr::summarize_if(is.double, list(mean = mean, sd=sd)) %>%
    ungroup()
  inputframe_num <- inputframe_num[,colSums(inputframe_num!=0)>0]

  inputframe_char <- inputframe %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::summarize_if(is.character, funs(dplyr::first))
  inputframe_out <- dplyr::full_join(inputframe_num, inputframe_char) %>% ungroup()

  #check for ambiguous stuff
  chcar <- inputframe_char %>% ungroup() %>% select_if(is.character) %>% colnames()
  message(paste("Character columns : '", chcar, "' are not used to subset this dataset and instead only the first occuring
                variable of each column per group was kept. To subset by these columns, add them to the input parameter 'by'.",
                collapse = " "))
  return(inputframe_out)
}


#' @export
#'
#'
checkVersionCompatible <- function(oldDataFrame, returnDataFrame=TRUE){
  if("Xrotum"%in%colnames(oldDataFrame)){
    colnames(oldDataFrame)[colnames(oldDataFrame)=="Xrotum"] <- "Xrot_micron"
    colnames(oldDataFrame)[colnames(oldDataFrame)=="Yrotum"] <- "Yrot_micron"
    message("Found old variable names 'Xrotum' and 'Yrotum'.")
    if(returnDataFrame==TRUE){
      message("Changed old variables for new variables 'Xrot_micron' and 'Yrot_micron'")
      return(oldDataFrame)
    }
  }
  if("Xrotum"%in%colnames(oldDataFrame)==F){
    message("No compatibility problems found.")
    if(returnDataFrame==TRUE){
      return(oldDataFrame)
    }
  }
}


##Set the pixel to um conversion
#' @export
addPixels2um <- function(pixelName, pixels2um){

  if(missing(pixels2um)){
    pixels2um <- readline("Give the conversion factor from pixel to um:  ")
  }
  if(missing(pixelName)){
    pixelName <- readline("Give the name of your conversion factor um:  ")
  }

  newp <- c(as.numeric(pixels2um))
  names(newp) <- pixelName
  newML <- append(get(magnificationList, envir=magEnv), newp)
  assign(magnificationList, newML, envir=magEnv)
  print("Currently loaded magnification converters:")
  print(get(magnificationList, envir=magEnv))

}

#' @export
getPixels2um <- function(){
  print("Currently loaded magnification converters:")
  print(get(magnificationList, envir=magEnv))
}



##merge spotfiles with only raw coordinates with mesh file with only raw data. add mesh length/width while on it.
#' @export
spotsInBox <- function(spotdata, meshdata, Xs = "x", Ys = "y", Xm = "X", Ym = "Y", meshInOutput=FALSE){
  if (!requireNamespace("shotGroups", quietly = TRUE)) {
    inp <- readline("Package 'shotGroups' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("shotGroups")}else{stop("Canceled")}
  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    inp <- readline("Package 'sp' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("sp")}else{stop("Canceled")}
  }

   #rewrite colnames if not the same as suggested
  pb <- utils::txtProgressBar(min=0, max=100, title="Total Progress SpotsInBox:")
  if(Xs!="x"){
    colnames(spotdata)[colnames(spotdata)==Xs] <- "x"
  }
  if(Ys!="y"){
    colnames(spotdata)[colnames(spotdata)==Ys] <- "y"
  }
  if(Xm!="X"){
    colnames(meshdata)[colnames(meshdata)==Xm] <- "X"
  }
  if(Ym!="Y"){
    colnames(meshdata)[colnames(meshdata)==Ym] <- "Y"
  }

  if("max.width"%in%colnames(meshdata)==T){u <- 1}
  if("max.width"%in%colnames(meshdata)==F){u<-2}

  if("length"%in%colnames(meshdata)==T){a <- 1}
  if("length"%in%colnames(meshdata)==F){a<-2} #if length and max width are already defined, don't touch them.
  utils::setTxtProgressBar(pb, 5)

#for-loop to re-write
  #cellframe_meshdata <- paste(meshdata$cell, meshdata$frame, sep="_")
 # o <- meshdata %>%
  #  dplyr::group_by(.data$frame, .data$cell) %>%
 #   getSpotsInBox(spotdatap = spotdata[spotdata$frame==meshdata$frame,],
       #           u=u,
     #             a=a,
       #           returnMESH = FALSE)


  if (!requireNamespace("pbapply", quietly = TRUE)) {
    o <- lapply(unique(meshdata$frame), function(x) lapply(unique(meshdata[meshdata$frame==x,]$cell), function(y) getSpotsInBox(meshp=meshdata[meshdata$frame==x&meshdata$cell==y,],
                                                                                                                                spotdatap=spotdata[spotdata$frame==x,],
                                                                                                                                u,
                                                                                                                                a,
                                                                                                                                returnMESH=meshInOutput)))
  }
  if (requireNamespace("pbapply", quietly = TRUE)) {
    o <- pbapply::pblapply(unique(meshdata$frame), function(x) lapply(unique(meshdata[meshdata$frame==x,]$cell), function(y) getSpotsInBox(meshp=meshdata[meshdata$frame==x&meshdata$cell==y,],
                                                                                                                     spotdatap=spotdata[spotdata$frame==x,],
                                                                                                                     u,
                                                                                                                   a,
                                                                                                                   returnMESH=meshInOutput)))
  }

  utils::setTxtProgressBar(pb, 45)

  outs <- list()
  outs$spots_relative <- do.call('rbind', lapply(o, function(x) do.call('rbind', lapply(x, function(y) y$REP))))
  utils::setTxtProgressBar(pb,65)
  names(outs) <- c("spots_relative")

  if(meshInOutput==TRUE){
    outs$mesh <-  do.call('rbind', lapply(o, function(x) do.call('rbind', lapply(x, function(y) y$MESH))))
    utils::setTxtProgressBar(pb, 85)
    names(outs) <- c("spots_relative", "mesh")
  }
  utils::setTxtProgressBar(pb, 100)
  return(outs) #return datasets as list of dataframes

}

getSpotsInBox <- function(meshp, spotdatap, u, a, returnMESH=FALSE){

  box <- suppressWarnings(shotGroups::getMinBBox(data.frame(x= meshp$X, y=meshp$Y))) #bounding box of cell
  lengthwidth <- c(box$width, box$height)

  if(u==2){
    meshp$max.width <- min(lengthwidth)
  }
  if(a==2){
    meshp$max.length <- max(lengthwidth) #take length/width if not already defined
  }

  pinps <- suppressWarnings(sp::point.in.polygon(spotdatap$x, spotdatap$y, meshp$X,meshp$Y)) #find spot/object coordinates inside cell
  pinps <- data.frame("x"=spotdatap$x, "y"=spotdatap$y, "pip"=pinps)
  if(nrow(pinps)>0){
    pinps <- pinps[pinps$pip==1,]
  }

  pts <- data.frame(box$pts) #get midpoint of the bounding box + define median lines
  shadowpts <- rbind(pts[2:4,], pts[1,])
  distances <- (pts+shadowpts)/2 #coordinates of median lines

  d1 <- data.frame(x = abs(distances$x-distances$x[1]), y =abs(distances$y-distances$y[1]))
  d1$pointn <- 1:4
  d1 <- d1[-1,]
  d1$dist <- sqrt(d1$x^2+d1$y^2)

  un <- data.frame(table(round(d1$dist, 5)))

  partdist <- un$Var1[un$Freq==1]
  partner <- d1$pointn[round(d1$dist,5)==partdist]
  rest <- d1$pointn[round(d1$dist,5)!=partdist]

  if(round(d1$dist[d1$pointn==partner],5)==round(min(lengthwidth),5)){
    widthline <- distances[c(1,partner),]
    lengthline <- distances[rest,]
  }
  if(round(d1$dist[d1$pointn==partner],5)==round(max(lengthwidth),5)){
    widthline <- distances[rest,]
    lengthline <- distances[c(1,partner),] #pick which line is length/width
  }

  mp <- c(mean(lengthline$x), mean(lengthline$y)) #midpoint
  X_cor <- meshp$X-mp[1]
  Y_cor <- meshp$Y-mp[2]
  angle <- (-box$angle)*pi/180 #angle to lay cell flat on x axis

  meshp$X_rot <- X_cor * cos(angle) - Y_cor * sin(angle)
  meshp$Y_rot <- X_cor * sin(angle) + Y_cor * cos(angle) #rotate cell

  if(nrow(pinps)>0){ #rotate spot/object points
    Lc <- pinps$x-mp[1]
    Dc <- pinps$y-mp[2]
    pinps$l <- -(Lc*cos(angle)-Dc*sin(angle))
    pinps$d <- -(Lc*sin(angle) +Dc*cos(angle))
    pinps$max.width <- unique(meshp$max.width)
    if("max_length"%in%colnames(meshp)){pinps$max.length <- unique(meshp$max_length)
    mesh$max.length <- mesh$max_length
    mesh$max_length <- NULL}
    else{pinps$max.length <- unique(meshp$max.length)}
    pinps$cell <- unique(meshp$cell)
    pinps$frame <- unique(meshp$frame)
  }
  #if(i==min.i&&n==min.n){ #first data frame in dataset
  #  if(nrow(pinps)>0){
  #    REP <- pinps
  #    }
  #  Mfull <- meshp
  #  }
  if("trajectory"%in%colnames(spotdatap)==T){
    pinps <- merge(spotdatap[,c("x", "y", "trajectory", "displacement_sq", "trajectory_length")], pinps)
  }

  if(returnMESH==TRUE){
    return(list("REP" = pinps, "MESH" = meshp))
  }
  if(returnMESH==FALSE){
    return(list("REP"=pinps))
  }
}


getpointsaround <- function(datsaround, angle){
  xlist <- c()
  ylist <- c()
  xlist[1] <- (datsaround$X_corRM)*cos(angle) - (datsaround$Y_corRM)*sin(angle)
  xlist[2]<- (datsaround$X_corRM)*cos(angle) - (datsaround$Y_corRMa)*sin(angle)
  xlist[3] <- (datsaround$X_corRMa)*cos(angle) - (datsaround$Y_corRMa)*sin(angle)
  xlist[4] <- (datsaround$X_corRMa)*cos(angle) - (datsaround$Y_corRM)*sin(angle)

  ylist[1] <- (datsaround$X_corRM)*sin(angle) + (datsaround$Y_corRM)*cos(angle)
  ylist[2] <- (datsaround$X_corRM)*sin(angle) + (datsaround$Y_corRMa)*cos(angle)
  ylist[3] <- (datsaround$X_corRMa)*sin(angle) + (datsaround$Y_corRMa)*cos(angle)
  ylist[4] <- (datsaround$X_corRMa)*sin(angle) + (datsaround$Y_corRM)*cos(angle)

  datforpoint <- data.frame(xt = xlist, yt=ylist, pointN=datsaround$pointN)
  return(datforpoint)
}

turnraws <- function(rawdatafile, i, n, mp, angle){
  rawr <- rawdatafile[rawdatafile$frame==i&rawdatafile$cell==n,]
  X_corR <- rawr$x - mp[1]
  Y_corR <- rawr$y - mp[2]
  rawr$X_rot <- X_corR * cos(angle) - Y_corR * sin(angle)
  rawr$Y_rot <- X_corR * sin(angle) + Y_corR * cos(angle)
  rawr$pointN <- c(1:nrow(rawr))
  datsaround <- data.frame(X_corRM = X_corR - 0.5, X_corRMa = X_corR + 0.5, Y_corRM = Y_corR - 0.5, Y_corRMa = Y_corR + 0.5, pointN = rawr$pointN)
  datsaround <- lapply(1:nrow(datsaround), function(x) getpointsaround(datsaround[x,], angle))
  datsaround <- do.call(rbind, datsaround)
  rawr <- merge(rawr, datsaround)
  return(rawr)
}

turncell <- function(MESHp, u, rawdatafile, a, n, i, ars){
  if (!requireNamespace("shotGroups", quietly = TRUE)) {

    inp <- readline("Package 'shotGroups' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("shotGroups")}else{stop("Canceled")}
  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    inp <- readline("Package 'sp' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("sp")}else{stop("Canceled")}
  }
  box <- suppressWarnings(shotGroups::getMinBBox(data.frame(x= MESHp$X, y=MESHp$Y))) #bounding box of cell
  lengthwidth <- c(box$width, box$height)
  if(ars==2){
    MESHp$area <- sp::Polygon(cbind(x=MESHp$X, y=MESHp$Y))@area
  }
  if(u==2){
    MESHp$max.width <- min(lengthwidth)
  }
  if(a==1){
    MESHp$max.length <- max(lengthwidth) #take length/width if not already defined
  }
  pts <- data.frame(box$pts) #get midpoint of the bounding box + define median lines
  shadowpts <- rbind(pts[2:4,], pts[1,])
  distances <- (pts+shadowpts)/2 #coordinates of median lines

  d1 <- data.frame(x = abs(distances$x-distances$x[1]), y =abs(distances$y-distances$y[1]))
  d1$pointn <- 1:4
  d1 <- d1[-1,]
  d1$dist <- sqrt(d1$x^2+d1$y^2)

  un <- data.frame(table(round(d1$dist, 5)))

  partdist <- un$Var1[un$Freq==1]
  partner <- d1$pointn[round(d1$dist,5)==partdist]
  rest <- d1$pointn[round(d1$dist,5)!=partdist]

  if(round(d1$dist[d1$pointn==partner],5)==round(min(lengthwidth),5)){
    widthline <- distances[c(1,partner),]
    lengthline <- distances[rest,]
  }
  if(round(d1$dist[d1$pointn==partner],5)==round(max(lengthwidth),5)){
    widthline <- distances[rest,]
    lengthline <- distances[c(1,partner),] #pick which line is length/width
  }

  mp <- c(mean(lengthline$x), mean(lengthline$y)) #midpoint
  X_cor <- MESHp$X-mp[1]
  Y_cor <- MESHp$Y-mp[2]
  angle <- (180-box$angle)*pi/180 #angle to lay cell flat on x axis

  MESHp$angle <- angle
  MESHp$Xmid <- mp[1]
  MESHp$Ymid <- mp[2]
  MESHp$X_rot <- X_cor * cos(angle) - Y_cor * sin(angle)
  MESHp$Y_rot <- X_cor * sin(angle) + Y_cor * cos(angle) #rotate cell
  if(MESHp$Y_rot[[1]]>0){
    MESHp$Y_rot <-  MESHp$X_rot * sin(pi) + MESHp$Y_rot * cos(pi)
    MESHp$X_rot <- MESHp$X_rot * cos(pi) - MESHp$Y_rot * sin(pi)
    MESHp$angle <- MESHp$angle - pi
  }


  if(!missing(rawdatafile)){
    message(paste("Turning raw data for cell", n))

    rawr <- turnraws(rawdatafile, i, n, mp, angle)
    return(list(mesh=MESHp, rawdat=rawr))
  }
  if(missing(rawdatafile)){return(MESHp)}
}


meshTurn <- function(MESH, Xm="X", Ym="Y", rawdatafile){
  if(Xm!="X"){colnames(MESH)[colnames(MESH)==Xm] <- "X"}
  if(Ym!="Y"){colnames(MESH)[colnames(MESH)==Ym] <- "Y"}

  if("max.width"%in%colnames(MESH)){u <- 1}else{u<-2}
  if("max.length"%in%colnames(MESH)){a <- 2}else{a<-1}
  if("area"%in%colnames(MESH)){ars <- 1}else{ars<-2}
  if("x0"%in%colnames(MESH)){
    MESH <- spotrXYMESH(MESH)
  }
  if(!missing(rawdatafile)){
    Rlist <- list()
  }
  Mlist <- list()

  #if length and max width are already defined, don't touch them.
  for(i in unique(MESH$frame)){ #per frame
    #print(paste("Turning meshes for frame", i))
    M <- MESH[MESH$frame==i,]
    if(!missing(rawdatafile)){
      Mlistboth <- lapply(unique(M$cell), function(x) turncell(M[M$cell==x,], u, rawdatafile,a, x, i, ars=ars))
      MlistF <- lapply(Mlistboth, function(x) x$mesh)
      RlistF <- lapply(Mlistboth, function(x) x$rawdat)
      Rlist[[i]] <- do.call(rbind, RlistF)
    }
    if(missing(rawdatafile)){MlistF <- lapply(unique(M$cell), function(x) turncell(M[M$cell==x,], u, a=a, n=x, i=i, ars=ars))}
    Mlist[[i]] <- do.call(rbind,MlistF)

  }

  Mfull <- do.call(rbind, Mlist)
  if(missing(rawdatafile)){
    return(Mfull) #return datasets as list of dataframes
  }else{
    rawdata_turned <- do.call(rbind, Rlist)
    rawdata_turned <- merge(rawdata_turned, unique(Mfull[,c("cell", "frame", "max.length", "max.width", "area")]))
    return(list(mesh = Mfull, rawdata_turned = rawdata_turned))}
}

spotrXYMESH <- function(MESH, x_1="x1", y_1="y1",x_0="x0", y_0="y0" ){
  u <- colnames(MESH)
  MESH <- MESH[!is.na(MESH$cell),]
  MESH0 <- MESH[,u[u!=x_1&u!=y_1]]
  MESH1 <- MESH[,u[u!=x_0&u!=y_0]]
  MESH1$xy <- 1
  MESH0$xy <- 0
  colnames(MESH1) <- gsub("1", "", colnames(MESH1))
  colnames(MESH0) <- gsub("0", "", colnames(MESH0))
  ##need to fix lapply function - check dapply? mapply? otherwise first make list of each row.
  MESH1$cf <- paste(MESH1$cell, MESH1$frame, sep="_")
  MESH1 <- do.call('rbind', lapply(unique(MESH1$cf), function(x) mesh1Fun(MESH1[MESH1$cf==x,])))
  MESH1$cf <- NULL
  MESH <- merge(MESH0, MESH1, all=T)
  colnames(MESH)[colnames(MESH)=="x"] <- "X"
  colnames(MESH)[colnames(MESH)=="y"] <- "Y"
  MESH <- MESH[order(MESH$frame, MESH$cell, MESH$num),]
  return(MESH)
}

mesh1Fun <- function(MESH1){
  MESH1$num <- 1 + max(MESH1$num, na.rm=T)*2 - MESH1$num
  return(MESH1)
}

mergeframes <- function(REP, MESH, mag="No_PixelCorrection", cutoff=T, maxfactor=2, minfactor=0.5, remOut=T, ouf=F){

  #REP<- REP[(0<REP$l),]
  if("rel.l" %in% colnames(REP)){
    REP <-REP[(REP$rel.l<1),]
  }
  if("max.length"%in%colnames(MESH)!=T&"max_length"%in%colnames(MESH)!=T){
    MESH$max.length <- MESH$length
  }
  if("max_length"%in%colnames(MESH)){
    MESH$max.length <- MESH$max_length
  }
  if("area"%in%colnames(MESH)){
    M <- unique(MESH[,c("cell", "frame", "max.length", "max.width", "area")])
  }else{M <- unique(MESH[,c("cell", "frame", "max.length", "max.width")])}
  M <- M[order(M$max.length),]
  M$cellnum <- c(1:nrow(M))

  #merging
  MR <- merge(M,REP, all=T)

  #remove MR's cells which have NA's in the cell area
  MR <- MR[!is.na(MR$max.length),]
  #remove duplicated rows
  MR <- MR[!duplicated(MR$l)|is.na(MR$l)&!duplicated(MR$d)|is.na(MR$d),]

  MR <- spotMR(MR)
  MR$totalspot[is.na(MR$totalspot)] <- 0

  #if needed: remove smallest and largest ones (cutoff: smaller than 1/2 mean and larger than 2x mean)
  if(cutoff==T){
    MR <- MR[MR$max.length<(maxfactor*mean(MR$max.length)),]
    MR <- MR[MR$max.length>(minfactor*mean(MR$max.length)),]
  }

  MR <- MR[order(MR$max.length),]

  #make column with row numbers per cell length.
  MR$num <- c(1:nrow(MR))

  pix2um <- unlist(get(magnificationList, envir=magEnv)[mag])
  MR <- LimDum(MR, pix2um, ouf=ouf)
  return(MR)
}

#add spotnumbers!
spotMR <- function(dat){
  if("spot" %in% colnames(dat)){
    #NA in spots replaced by "0"
    dat$spot[is.na(dat$spot)] <- 0
  }else{
    dat <-dat[order(dat$frame, dat$cell,dat$max.length),]
    dat$spot <- 1
    for(n in 1:(nrow(dat)-1)){
      if(dat$max.length[n+1]==dat$max.length[n]){
        dat$spot[n+1] <- dat$spot[n] + 1
      } }
    dat$spot[is.na(dat$spot)] <- 0 }
  if("totalspot"%in%colnames(dat)){
    dat$totalspot[is.na(dat$totalspot)] <- 0
  }else{dat$cellframe <- paste(dat$cell, dat$frame, sep=".")
  spotn <- data.frame(cellframe = unique(dat$cellframe), totalspot = unlist(lapply(unique(dat$cellframe), function(x)max(dat$spot[dat$cellframe==x]))))
  dat <- merge(dat, spotn)
  }
  return(dat)
}

centrefun <- function(dat, xie="ob_x", yie="ob_y"){
  if (!requireNamespace("shotGroups", quietly = TRUE)) {
    inp <- readline("Package 'shotGroups' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("shotGroups")}else{stop("Canceled")}
  }
  dat <- dat[!is.na(dat$ob_x),]
  dat$centre_x <- NA
  dat$centre_y <- NA
  if (requireNamespace("pbapply", quietly = TRUE)) {
  dat <- do.call('rbind', pbapply::pblapply(unique(dat$obID), function(x) takeObjectCentre(dat[dat$obID==x,], xie, yie)))
  }
  if (!requireNamespace("pbapply", quietly = TRUE)) {
    dat <- do.call('rbind', lapply(unique(dat$obID), function(x) takeObjectCentre(dat[dat$obID==x,], xie, yie)))
  }
  return(dat)
}

takeObjectCentre <- function(dat, xie, yie){
  woei <- as.matrix(dat[c(xie, yie)])
  cen <- shotGroups::getMinCircle(woei)$ctr
  dat$centre_x <- cen[1]
  dat$centre_y <- cen[2]
  return(dat)
}

#add object centre to mesh file and turn accordingly
#' @importFrom dplyr %>%
midobject <- function(MESH, OBJ, p2um){
  OBJ$angle <- NULL
  MESH <- MESH %>%
    dplyr::left_join(OBJ) %>%
    dplyr::mutate(xccor = .data$centre_x - .data$Xmid,
                  yccor = .data$centre_y - .data$Ymid,
                  xcorcor = .data$ob_x - .data$Xmid,
                  ycorcor = .data$ob_y - .data$Ymid,
                  Lmid = .data$xccor * cos(.data$angle)-.data$yccor*sin(.data$angle),
                  Dum = .data$xccor * sin(.data$angle)+.data$yccor*cos(.data$angle),
                  ob_out_x = .data$xcorcor * cos(.data$angle) - .data$ycorcor * sin(.data$angle),
                  ob_out_y = .data$xcorcor * sin(.data$angle) + .data$ycorcor * cos(.data$angle)) %>%
    dplyr::select(.data$frame, .data$cell, .data$obpath, .data$obnum, .data$obID, .data$max.length, .data$max.width, .data$Dum, .data$Lmid, .data$ob_out_x, .data$ob_out_y) %>%
    dplyr::distinct() %>%
    dplyr::mutate(num = dplyr::dense_rank(.data$max.length)) %>%
    LimDum(p2um) %>%
    dplyr::mutate(ob_out_x = .data$ob_out_x * p2um, ob_out_y = .data$ob_out_y * p2um)
  return(MESH)

  #MESH$xccor <- MESH$centre_x - MESH$Xmid
  #MESH$yccor <- MESH$centre_y - MESH$Ymid
  #MESH$xcorcor <- MESH$ob_x - MESH$Xmid
 # MESH$ycorcor <- MESH$ob_y - MESH$Ymid

 # MESH$Lmid <- MESH$xccor * cos(MESH$angle) - MESH$yccor * sin(MESH$angle)
 # MESH$Dum <- MESH$xccor * sin(MESH$angle) + MESH$yccor * cos(MESH$angle)
 # MESH$ob_out_x <- MESH$xcorcor * cos(MESH$angle) - MESH$ycorcor * sin(MESH$angle)
 # MESH$ob_out_y <- MESH$xcorcor * sin(MESH$angle) + MESH$ycorcor * cos(MESH$angle)
  #MO <- MESH[,c("frame", "cell", "obpath", "obnum", "obID", "max.length", "max.width", "Dum", "Lmid", "ob_out_x", "ob_out_y")]
#  MO <- unique(MO)

  #MOnum <- unique(MO[,c("frame", "cell", "max.length", "obnum")])
  #MOnum <- MOnum[order(MOnum$max.length),]
  #MOnum$num <- 1:nrow(MOnum)
  #MO <- merge(MOnum, MO)
  #MO <- LimDum(MO, p2um)
  #MO$ob_out_x <- MO$ob_out_x*p2um
  #MO$ob_out_y <- MO$ob_out_y*p2um
  #return(MO)
}

################################################################################################
#plot preparation
#quartiles, maxima, etc.
LimDum <- function(MR, pix2um, remOut=T, ouf=F){
  if("l"%in%colnames(MR)==TRUE&ouf==F){
    MR$Lmid<-MR$l*pix2um
  }
  if("l"%in%colnames(MR)==TRUE&ouf==T){
    MR$Lmid<-(MR$l-(MR$max.length/2))*pix2um
  }
  if("l"%in%colnames(MR)!=TRUE){
    MR$Lmid <- MR$Lmid*pix2um
    }
  MR$pole1<- -MR$max.length*0.5*pix2um
  MR$pole2<- -MR$pole1
  if("d"%in%colnames(MR)){MR$Dum <- MR$d*pix2um}
  if("d"%in%colnames(MR)==FALSE){MR$Dum <- MR$Dum*pix2um}
  MR$max_um <- MR$max.length*pix2um
  MR$maxwum <- MR$max.width*pix2um
  if(remOut==T){
    MR <- MR[abs(MR$Lmid)<MR$pole2|is.na(MR$Lmid),]
    MR <- MR[abs(MR$Dum)<(MR$max.width/2)|is.na(MR$Lmid),]
  }
  return(MR)
}

mL <- list("100x_TIRF" = 0.05204891, "100x_DVMolgen" = 0.0645500, "No_PixelCorrection" = 1, "100x_FRAP" = 0.0499548)
magnificationList <- "magnificationList"
magEnv <- new.env()
assign(magnificationList, mL, envir=magEnv)



#'@export
micron <- function(){
  return("\u00b5m")
}


#'@export
#' @title Orient your cells by the side in which your spots or objects (predominantly) are located
#'
#' \code{orientCells()} takes a \code{spots_relative} or \code{object_relative} dataset and uses the relative localization of each spot or mid-point of each object
#' to flip all related cells & content to one side based on this localization.
#'
#' @param dataTurn has to be a \code{spots_relative} or \code{object_relative} dataset
#' @param cellData a list of corresponding datasets that need to be turned as well. can be 1 dataset (for instance \code{mesh}) or a list of datasets (like \code{list(mesh, object_relative)})
#' @return the function will return a list of the datasets, where the relative X/Y coordinates of the cells and internal features are all oriented to one side.
#'
#' @examples
#' \dontrun{
#'
#' output <- orientCells(dataTurn = myData$spots_relative, cellData = list(myData$mesh, myData$object_relative))
#'}
#'
orientCells <- function(dataTurn, cellData, removeNonOriented = FALSE){
  #check if the dataset is relative
  if(!"Lmid"%in%colnames(dataTurn)){
    if(!"l"%in%colnames(dataTurn)){
      stop("No relative localizations found in the 'dataTurn' input. Please use a dataframe with relative spot or object locatizations.")
    }
    if("l"%>%colnames(dataTurn)){
      #rename l & d
      dataTurn <- dataTurn %>%
        dplyr::rename(Lmid=.data$l, Dum = .data$d)
    }
  }
  #in case of objects
  if("obID"%in%colnames(dataTurn)){
    dataTurn <- dataTurn %>%
      dplyr::group_by(.data$frame, .data$cell) %>%
      dplyr::mutate(totalspot = dplyr::n_distinct(.data$obID),
                    spot = obID) %>%
      dplyr::ungroup()
  }

  #classify location of spot/object middle with variable 'pol' (handy if you have more than 1 spot)
  dataTurn <- dataTurn %>%
    dplyr::mutate(pol = dplyr::if_else(.data$Lmid<0, -1, 1, 0))

  if(max(dataTurn$totalspot)>1){
    message("More than one spot per cell detected. Function takes the side most spots are on (in case of a tie, cell is classified as 'non-polarized' by setting variable 'pol' to 0).")
    dataSmall <- dataTurn %>%
      dplyr::distinct(.data$frame, .data$cell, .data$spot, .data$pol) %>%
      dplyr::group_by(.data$frame, .data$cell) %>%
      dplyr::mutate(pol = mean(.data$pol)) %>%
      dplyr::mutate(pol = if_else(pol==-1|pol==1, .data$pol, 0, 0))
    dataTurn <- dataTurn %>%
      dplyr::select(-.data$pol) %>%
      left_join(dataSmall)
  }

  if(removeNonOriented == TRUE){
    dataTurn <- dataTurn[dataTurn$pol!=0,]
  }

  listPol <- dataTurn %>%
    dplyr::select(.data$frame, .data$cell, .data$pol)

  dataTurn <- turnEach(dataTurn)

  if(is.data.frame(cellData)==FALSE){
    if(is.list(cellData)==FALSE){
      stop("Does not recognize the cellData input. Either provide one dataframe or a list() of dataframes.")
    }
    cellData <- lapply(cellData, function(x) turnEach(x, listPol))
  }else{
    cellData <- turnEach(cellData, listPol)
  }

  return(list(dataTurn = dataTurn, cellData = cellData))

}


turnEach <- function(partCD, listPol){
  if("spot"%in%colnames(partCD)){
    dType <- "spot"
    if(!"Lmid"%in%colnames(partCD)){
      if(!"l"%in%colnames(partCD)){
        stop("The spot dataframe does not contain relative localizations. Please add a spots_relative dataframe.")
      }
      partCD <- partCD %>%
        dplyr::rename(Lmid = .data$l, Dum = .data$d)
    }
  }
  if("obID"%in%colnames(partCD)){
    dType <- "object"
  }
  if("Xmid"%in%colnames(partCD)){
    dType <- "mesh"
  }
  if(missing(dType)){
    stop("Did not recognize the input dataframe in cellData.")
  }

  if(!missing(listPol)){
    partCD <- partCD %>%
      dplyr::right_join(listPol)
  }


  if(dType == "spot"){
    partCD <- partCD %>%
      dplyr::mutate(
        Lmid = dplyr::if_else(.data$pol!=0, .data$Lmid*.data$pol, .data$Lmid, .data$Lmid),
      )
  }
  if(dType == "object"){
    partCD <- partCD %>%
      dplyr::mutate(
        Lmid = dplyr::if_else(.data$pol!=0, .data$Lmid*.data$pol, .data$Lmid, .data$Lmid),
        ob_out_x = dplyr::if_else(.data$pol!=0, .data$ob_out_x*.data$pol, .data$ob_out_x, .data$ob_out_x)
      )
  }
  if(dType == "mesh"){
    partCD <- partCD %>%
      dplyr::mutate(
        X_rot = dplyr::if_else(.data$pol!=0, .data$X_rot*.data$pol, .data$X_rot, .data$X_rot),
        Xrot_micron = dplyr::if_else(.data$pol!=0, .data$Xrot_micron*.data$pol, .data$Xrot_micron, .data$Xrot_micron)
      )
    if("xt"%in%colnames(partCD)){
      partCD <- partCD %>%
        dplyr::mutate(
          xt = dplyr::if_else(.data$pol!=0, .data$xt*.data$pol, .data$xt, .data$xt)
        )
    }
  }

  return(partCD)

}

