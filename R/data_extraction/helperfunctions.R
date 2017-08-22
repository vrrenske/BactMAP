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
spotrSetpixel2um <- function(){
  pixel2um <- readline(caption="Give the conversion factor from pixel to um./n")
  return(pixel2um)
}

##error message + solution when pixel2um is not set
spotrPixelerror <- function(){
  errormessage <- readline(caption="The conversion factor from pixel to um is not indicated. Please use 'spotrSetpixel2um' by pressing 'a'.
          if you don't want to convert from pixel to um, press 'b'./n")
  if(errormessage=="a"|errormessage=="A"){
    conv == spotrSetpixel2um()
  }
  if(errormessage=="b"|errormessage=="B"){
    conv == 1
  }
  else(print("Did not receive 'a' nor 'b'. If you want to convert from pixel to um, use the function 'spotrSetpixel2um()' manually./n"))
  return(conv)
}

##merge spotfiles with only raw coordinates with mesh file with only raw data. add mesh length/width while on it.

spotrSpotsinBox <- function(spotfile, MESH, Xs = "x", Ys = "y", Xm = "X", Ym = "Y"){
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
  else{u<-2}
  if("length"%in%colnames(MESH)==T){a <- 1}
  else{a<-2} #if length and max width are already defined, don't touch them.

  min.i <- min(MESH$frame)
  for(i in unique(MESH$frame)){ #per frame
    min.n <- min(MESH$cell[MESH$frame==i])

    for(n in unique(MESH$cell[MESH$frame==i])){ #per cell
      MESHp <- MESH[MESH$cell==n&MESH$frame==i,] #define part of the frame to run faster

      box <- shotGroups::getMinBBox(data.frame(x= MESHp$X, y=MESHp$Y)) #bounding box of cell
      lengthwidth <- c(box$width, box$height)

      if(u==2){
        MESHp$max.width <- min(lengthwidth)
      }
      if(a==1){
        MESHp$length <- max(lengthwidth) #take length/width if not already defined
      }

      pinps <- SDMTools::pnt.in.poly(spotfile[,c("x","y")], data.frame(MESHp$X,MESHp$Y)) #find spot/object coordinates inside cell
      pinps <- pinps[pinps$pip==1,]

      pts <- data.frame(box$pts) #get midpoint of the bounding box + define median lines
      shadowpts <- rbind(pts[2:4,], pts[1,])
      distances <- (pts+shadowpts)/2 #coordinates of median lines

      d1 <- data.frame(x = abs(distances$x-distances$x[1]), y =abs(distances$y-distances$y[1]))
      d1$pointn <- 1:4
      d1 <- d1[-1,]
      d1$dist <- sqrt(d1$x^2+d1$y^2)

      un <- data.frame(table(round(d1$dist, 4)))

      partdist <- un$Var1[un$Freq==1]
      partner <- d1$pointn[round(d1$dist,4)==partdist]
      rest <- d1$pointn[round(d1$dist,4)!=partdist]

      if(round(d1$dist[d1$pointn==partner],4)==round(min(lengthwidth),4)){
        widthline <- distances[c(1,partner),]
        lengthline <- distances[rest,]
      }
      if(round(d1$dist[d1$pointn==partner],4)==round(max(lengthwidth),4)){
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
        pinps$L <- Lc*cos(angle)-Dc*sin(angle)
        pinps$D <- Lc*sin(angle) +Dc*cos(angle)
        pinps$max.width <- unique(MESHp$max.width)
        if("max_length"%in%colnames(MESHp)){pinps$length <- unique(MESHp$max_length)}
        else{pinps$length <- unique(MESHp$length)}
        pinps$cell <- n
        pinps$frame <- i
      }
      if(i==min.i&&n==min.n){ #first data frame in dataset
        REP <- pinps
        Mfull <- MESHp
      }
      else{
        REP <- rbind(REP, pinps) #bind the rest to it
        Mfull <- rbind(Mfull, MESHp)
      }
    }


  }
  return(list(REP, Mfull)) #return datasets as list of dataframes

}



spotrMeshTurn <- function(MESH, Xm="X", Ym="Y"){
  if(Xm!="X"){colnames(MESH)[colnames(MESH)==Xm] <- "X"}
  if(Ym!="Y"){colnames(MESH)[colnames(MESH)==Ym] <- "Y"}

  if("max.width"%in%colnames(MESH)==T){u <- 1}
  else{u<-2}
  if("length"%in%colnames(MESH)==T){a <- 1}
  else{a<-2}
  if("x0"%in%colnames(MESH)){
    spotrXYMESH(MESH)
  }

  #if length and max width are already defined, don't touch them.
  min.i <- min(MESH$frame)
  for(i in unique(MESH$frame)){ #per frame
    min.n <- min(MESH$cell[MESH$frame==i])

    for(n in unique(MESH$cell[MESH$frame==i])){ #per cell
      MESHp <- MESH[MESH$cell==n&MESH$frame==i,] #define part of the frame to run faster

      box <- shotGroups::getMinBBox(data.frame(x= MESHp$X, y=MESHp$Y)) #bounding box of cell
      lengthwidth <- c(box$width, box$height)

      if(u==2){
        MESHp$max.width <- min(lengthwidth)
      }
      if(a==1){
        MESHp$length <- max(lengthwidth) #take length/width if not already defined
      }

      pts <- data.frame(box$pts) #get midpoint of the bounding box + define median lines
      shadowpts <- rbind(pts[2:4,], pts[1,])
      distances <- (pts+shadowpts)/2 #coordinates of median lines

      d1 <- data.frame(x = abs(distances$x-distances$x[1]), y =abs(distances$y-distances$y[1]))
      d1$pointn <- 1:4
      d1 <- d1[-1,]
      d1$dist <- sqrt(d1$x^2+d1$y^2)

      un <- data.frame(table(round(d1$dist, 4)))

      partdist <- un$Var1[un$Freq==1]
      partner <- d1$pointn[round(d1$dist,4)==partdist]
      rest <- d1$pointn[round(d1$dist,4)!=partdist]

      if(round(d1$dist[d1$pointn==partner],4)==round(min(lengthwidth),4)){
        widthline <- distances[c(1,partner),]
        lengthline <- distances[rest,]
      }
      if(round(d1$dist[d1$pointn==partner],4)==round(max(lengthwidth),4)){
        widthline <- distances[rest,]
        lengthline <- distances[c(1,partner),] #pick which line is length/width
      }

      mp <- c(mean(lengthline$x), mean(lengthline$y)) #midpoint
      X_cor <- MESHp$X-mp[1]
      Y_cor <- MESHp$Y-mp[2]
      angle <- (180-box$angle)*pi/180 #angle to lay cell flat on x axis

      MESHp$X_rot <- X_cor * cos(angle) - Y_cor * sin(angle)
      MESHp$Y_rot <- X_cor * sin(angle) + Y_cor * cos(angle) #rotate cell

      if(i==min.i&&n==min.n){ #first data frame in dataset
        Mfull <- MESHp
      }
      else{#bind the rest to it
        Mfull <- rbind(Mfull, MESHp)
      }
    }


  }
  return(Mfull) #return datasets as list of dataframes

}

spotrXYMESH <- function(MESH, x_1="x1", y_1="y1",x_0="x0", y_0="y0" ){
  u <- colnames(MESH)
  MESH0 <- MESH[,u[u!=x_1&u!=y_1]]
  MESH1 <- MESH[,u[u!=x_0&u!=y_0]]
  MESH1$xy <- 1
  MESH0$xy <- 0
  colnames(MESH1) <- gsub("1", "", colnames(MESH1))
  colnames(MESH0) <- gsub("0", "", colnames(MESH0))
  for(n in unique(MESH1$frame)){
    for(i in unique(MESH1$cell)){
      MESH1$num[MESH1$cell==i&MESH1$frame==n] <- 1 + max(MESH1$num[MESH1$cell==i&MESH1$frame==n], na.rm=T)*2 - MESH1$num[MESH1$cell==i&MESH1$frame==n]
    }
  }
  MESH <- merge(MESH0, MESH1, all=T)
  MESH <- MESH[order(MESH$frame, MESH$cell, MESH$num),]
  return(MESH)
}


