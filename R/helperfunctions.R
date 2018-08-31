##Renske van Raaphorst
##7/14/2017

##Helperfunctions: functions necessary to make other functions working properly.
##other package dependencies:

#merge data functions
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(scales)

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

##error message + solution when pixel2um is not set
#spotrPixelerror <- function(){
  #errormessage <- readline(caption="The conversion factor from pixel to um is not indicated. Please use 'spotrSetpixel2um' by pressing 'a'.
  #        if you don't want to convert from pixel to um, press 'b'./n")
 # if(errormessage=="a"|errormessage=="A"){
 #   conv == spotrSetpixel2um()
 # }
 # if(errormessage=="b"|errormessage=="B"){
 #   conv == 1
 # }
  #else(print("Did not receive 'a' nor 'b'. If you want to convert from pixel to um, use the function 'spotrSetpixel2um()' manually./n"))
  #return(conv)
#}

##merge spotfiles with only raw coordinates with mesh file with only raw data. add mesh length/width while on it.
#' @export
spotsInBox <- function(spotfile, MESH, Xs = "x", Ys = "y", Xm = "X", Ym = "Y"){
  q <- 0
  b <- 0
   #rewrite colnames if not the same as suggested
  if(Xs!="x"){
    colnames(spotfile)[colnames(spotfile)==Xs] <- "x"
  }
  if(Ys!="y"){
    colnames(spotfile)[colnames(spotfile)==Ys] <- "y"
  }
  if(Xm!="X"){
    colnames(MESH)[colnames(MESH)==Xm] <- "X"
  }
  if(Ym!="Y"){
    colnames(MESH)[colnames(MESH)==Ym] <- "Y"
  }

  if("max.width"%in%colnames(MESH)==T){u <- 1}
  if("max.width"%in%colnames(MESH)==F){u<-2}

  if("length"%in%colnames(MESH)==T){a <- 1}
  if("length"%in%colnames(MESH)==F){a<-2} #if length and max width are already defined, don't touch them.

  min.i <- min(MESH$frame)
  for(i in unique(MESH$frame)){ #per frame
    min.n <- min(MESH$cell[MESH$frame==i])
    spotfilep <- spotfile[spotfile$frame==i,]
    for(n in unique(MESH$cell[MESH$frame==i])){ #per cell
      MESHp <- MESH[MESH$cell==n&MESH$frame==i,] #define part of the frame to run faster

      box <- suppressWarnings(shotGroups::getMinBBox(data.frame(x= MESHp$X, y=MESHp$Y))) #bounding box of cell
      lengthwidth <- c(box$width, box$height)

      if(u==2){
        MESHp$max.width <- min(lengthwidth)
      }
      if(a==2){
        MESHp$max.length <- max(lengthwidth) #take length/width if not already defined
      }

      pinps <- suppressWarnings(SDMTools::pnt.in.poly(spotfilep[,c("x","y")], data.frame(MESHp$X,MESHp$Y))) #find spot/object coordinates inside cell
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
      X_cor <- MESHp$X-mp[1]
      Y_cor <- MESHp$Y-mp[2]
      angle <- (-box$angle)*pi/180 #angle to lay cell flat on x axis

      MESHp$X_rot <- X_cor * cos(angle) - Y_cor * sin(angle)
      MESHp$Y_rot <- X_cor * sin(angle) + Y_cor * cos(angle) #rotate cell

      if(nrow(pinps)>0){ #rotate spot/object points
        Lc <- pinps$x-mp[1]
        Dc <- pinps$y-mp[2]
        pinps$l <- -(Lc*cos(angle)-Dc*sin(angle))
        pinps$d <- -(Lc*sin(angle) +Dc*cos(angle))
        pinps$max.width <- unique(MESHp$max.width)
        if("max_length"%in%colnames(MESHp)){pinps$max.length <- unique(MESHp$max_length)
                                            MESH$max.length <- MESH$max_length
                                            MESH$max_length <- NULL}
        else{pinps$max.length <- unique(MESHp$max.length)}
        pinps$cell <- n
        pinps$frame <- i
      }
      #if(i==min.i&&n==min.n){ #first data frame in dataset
      #  if(nrow(pinps)>0){
      #    REP <- pinps
    #    }
      #  Mfull <- MESHp
    #  }

          if(nrow(pinps)>0){
          if(q==0){
            REP <- pinps
            q <- 1
          }
          if(q!=0){
            REP <- rbind(REP, pinps)
          }
        }
        if(b==0){
          Mfull <- MESHp
          b <- 1
        }
        if(b!=0){
         #bind the rest to it
          Mfull <- rbind(Mfull, MESHp)
        }


    }


  }
  outs <- list()
  REP <- spotMR(REP)
  outs$spots_relative <- REP
  outs$mesh <- Mfull
  return(outs) #return datasets as list of dataframes

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

turncell <- function(MESHp, u, rawdatafile, a, n, i){
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

  if(missing(rawdatafile)!=T){
    print(paste("Turning raw data for cell", n))

    rawr <- turnraws(rawdatafile, i, n, mp, angle)
    return(list(mesh=MESHp, rawdat=rawr))
  }
  if(missing(rawdatafile)==T){return(MESHp)}
}


meshTurn <- function(MESH, Xm="X", Ym="Y", rawdatafile){
  if(Xm!="X"){colnames(MESH)[colnames(MESH)==Xm] <- "X"}
  if(Ym!="Y"){colnames(MESH)[colnames(MESH)==Ym] <- "Y"}

  if("max.width"%in%colnames(MESH)==T){u <- 1}
  else{u<-2}
  if("max.length"%in%colnames(MESH)==T){a <- 2}
  else{a<-1}
  if("area"%in%colnames(MESH)==T){ars <- 1}
  else{ars<-2}
  if("x0"%in%colnames(MESH)){
    MESH <- spotrXYMESH(MESH)
  }
  if(missing(rawdatafile)!=T){
    Rlist <- list()
  }
  Mlist <- list()

  #if length and max width are already defined, don't touch them.
  for(i in unique(MESH$frame)){ #per frame
    #print(paste("Turning meshes for frame", i))
    M <- MESH[MESH$frame==i,]
    if(missing(rawdatafile)!=T){
      Mlistboth <- lapply(unique(M$cell), function(x) turncell(M[M$cell==x,], u, rawdatafile,a, x, i))
      MlistF <- lapply(Mlistboth, function(x) x$mesh)
      RlistF <- lapply(Mlistboth, function(x) x$rawdat)
      Rlist[[i]] <- do.call(rbind, RlistF)
    }
    if(missing(rawdatafile)==T){MlistF <- lapply(unique(M$cell), function(x) turncell(M[M$cell==x,], u, a=a, n=x, i=i))}
    Mlist[[i]] <- do.call(rbind,MlistF)

  }

  Mfull <- do.call(rbind, Mlist)
  if(missing(rawdatafile)==T){
    return(Mfull) #return datasets as list of dataframes
  }
  else{
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
  for(n in unique(MESH1$frame)){
    #print(n)
    for(i in unique(MESH1$cell[MESH1$frame==n])){
      MESH1$num[MESH1$cell==i&MESH1$frame==n] <- 1 + max(MESH1$num[MESH1$cell==i&MESH1$frame==n], na.rm=T)*2 - MESH1$num[MESH1$cell==i&MESH1$frame==n]
    }
  }
  MESH <- merge(MESH0, MESH1, all=T)
  colnames(MESH)[colnames(MESH)=="x"] <- "X"
  colnames(MESH)[colnames(MESH)=="y"] <- "Y"
  MESH <- MESH[order(MESH$frame, MESH$cell, MESH$num),]
  return(MESH)
}


mergeframes <- function(REP, MESH, mag="100x_LeicaVeening", cutoff=T, maxfactor=2, minfactor=0.5, remOut=T, ouf=F){

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
  }
  else{M <- unique(MESH[,c("cell", "frame", "max.length", "max.width")])}
  M <- M[order(M$max.length),]
  M$cellnum <- c(1:nrow(M))

  #merging
  MR <- merge(M,REP, all=T)

  #remove MR's cells which have NA's in the cell area
  MR <- MR[!is.na(MR$max.length),]
  #remove duplicated rows
  MR <- MR[!duplicated(MR$l)&!duplicated(MR$d),]

  MR <- spotMR(MR)

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
  } else {
    dat <-dat[order(dat$frame, dat$cell,dat$max.length),]
    dat$spot <- 1
    for(n in 1:(nrow(dat)-1)){
      if(dat$max.length[n+1]==dat$max.length[n]){
        dat$spot[n+1] <- dat$spot[n] + 1
      } }
    dat$spot[is.na(dat$spot)] <- 0
  }
  dat$cellframe <- paste(dat$cell, dat$frame, sep=".")
  spotn <- data.frame(cellframe = unique(dat$cellframe), totalspot = unlist(lapply(unique(dat$cellframe), function(x)max(dat$spot[dat$cellframe==x]))))
  dat <- merge(dat, spotn)
  return(dat)
}

centrefun <- function(dat, xie="ob_x", yie="ob_y"){
  dat <- dat[!is.na(dat$ob_x),]
  dat$centre_x <- NA
  dat$centre_y <- NA
  for(n in unique(dat$frame)){
    for(u in unique(dat$cell)){
      if(nrow(dat[dat$cell==u&dat$frame==n,])>0){
        for(i in unique(dat$obnum[dat$cell==u&dat$frame==n])){
          woei <- as.matrix(dat[c(xie, yie)][dat$cell==u&dat$frame==n&dat$obnum==i,])
          cen <- shotGroups::getMinCircle(woei)$ctr
          dat$centre_x[dat$cell==u&dat$frame==n&dat$obnum==i] <- cen[1]
          dat$centre_y[dat$cell==u&dat$frame==n&dat$obnum==i] <- cen[2]
        }
      }
    }
  }
  return(dat)
}

#add object centre to mesh file and turn accordingly
midobject <- function(MESH, OBJ, p2um){
  MESH <- merge(MESH, OBJ, all=T)
  MESH$xccor <- MESH$centre_x - MESH$Xmid
  MESH$yccor <- MESH$centre_y - MESH$Ymid
  MESH$xcorcor <- MESH$ob_x - MESH$Xmid
  MESH$ycorcor <- MESH$ob_y - MESH$Ymid
  #MESHp$X_rot <- X_cor * cos(angle) - Y_cor * sin(angle)
  MESH$Lmid <- MESH$xccor * cos(MESH$angle) - MESH$yccor * sin(MESH$angle)
  MESH$Dum <- MESH$xccor * sin(MESH$angle) + MESH$yccor * cos(MESH$angle)
  MESH$ob_out_x <- MESH$xcorcor * cos(MESH$angle) - MESH$ycorcor * sin(MESH$angle)
  MESH$ob_out_y <- MESH$xcorcor * sin(MESH$angle) + MESH$ycorcor * cos(MESH$angle)
  MO <- MESH[,c("frame", "cell", "obpath", "obnum", "obID", "max.length", "max.width", "Dum", "Lmid", "ob_out_x", "ob_out_y")]
  MO <- unique(MO)

  MOnum <- unique(MO[,c("frame", "cell", "max.length", "obnum")])
  MOnum <- MOnum[order(MOnum$max.length),]
  MOnum$num <- 1:nrow(MOnum)
  MO <- merge(MOnum, MO)
  MO <- LimDum(MO, p2um)
  MO$ob_out_x <- MO$ob_out_x*p2um
  MO$ob_out_y <- MO$ob_out_y*p2um
  return(MO)
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
  MR$max.length <- MR$max.length*pix2um
  MR$max.width <- MR$max.width*pix2um
  if(remOut==T){
    MR <- MR[abs(MR$Lmid)<MR$pole2,]
    MR <- MR[abs(MR$Dum)<(MR$max.width/2),]
  }
  return(MR)
}

mL <- list("100x_TIRF" = 0.0499538, "100x_DVMolgen" = 0.0645500, "No_PixelCorrection" = 1, "100x_FRAP" = 0.0499548)
magnificationList <- "magnificationList"
magEnv <- new.env()
assign(magnificationList, mL, envir=magEnv)




