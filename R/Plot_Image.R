##collection of raw image plotting functions for BactMAP

#upload your image stack (one color)
#' @export
extr.OriginalStack <- function(picloc){
  im <- tiff::readTIFF(picloc, all=T) #if you want the best resolution, it needs to be a .tiff file
  im <- lapply(im, function(x) raster::raster(x))
  imdatframe <- lapply(im, function(x) as.data.frame(as(x, "SpatialPixelsDataFrame"))) #get values
  imdatframe <- lapply(imdatframe, function(x) changecols(x))
  nx <- length(unique(imdatframe[[1]]$x))
  ny <- length(unique(imdatframe[[1]]$y))
  imdatframe <- lapply(imdatframe, function(x) changeres(x, nx, ny))
  return(imdatframe)
}


#' @export
extr.OriginalCells <- function(imdatframe, mesh){
    allcellslist <- lapply(unique(mesh$frame), function(x) pipperframe(imdatframe, mesh, x))
    allcellsframe <- do.call(rbind, allcellslist)
    allcellsframe$cell <- allcellsframe$pip
    allcellsframe$pip <- NULL
    outlist <- meshTurn(mesh, rawdatafile=allcellsframe)
  return(outlist)
}


#' @export
plotCellsTime <- function(celdat,
                           updown = T,
                           movie = F,
                           viridisoption = "magma",
                           cellN,
                           minf,
                           maxf){
  if(missing(minf)){
    minf <- min(celdat$frame)
  }
  if(missing(maxf)){
    maxf <- max(celdat$frame)
  }
  #when no cell number is indicated:return a list of plots/movie objects
  if(missing(cellN)){
    plotout <- lapply( unique(celdat$cell), function(x) plotcellsframelist(celdat[celdat$cell==x,], maxframes=maxf, minframes=minf, updown, movie, viridisoption) + ggplot2::ggtitle(x))
  }
  #and just plot one plot if cellN exists
  if(missing(cellN)!=T){
    if(length(cellN)==1){
      plotout <- plotcellsframelist(celdat[celdat$cell==cellN,], maxf, minf, updown, movie, viridisoption) + ggplot2::ggtitle(cellN)
    }
    if(length(cellN)>1){
      plotout <- lapply(cellN, function(x) plotcellsframelist(celdat[celdat$cell==x,], maxf, minf,updown,movie,viridisoption) + ggplot2::ggtitle(x))
    }
  }
  return(plotout)
}

###################################|(*-*)/#############################################################
changecols <- function(x){
  colnames(x) <- c("values" , "x", "y")
  return(x)
}
changeres <- function(dat, nx, ny){
  dat$x <- dat$x*nx
  dat$y <- max(dat$y)-dat$y
  dat$y <- dat$y*ny
  return(dat)
}



pinping <- function(dat, mesh, x){
  print(paste("Cell", x))
  mesh <- mesh[mesh$cell==x,]
  minmeshx <- min(mesh$X)-2
  minmeshy <- min(mesh$Y)-2
  maxmeshx <- max(mesh$X)+2
  maxmeshy <- max(mesh$Y)+2
  dat <- dat[dat$x>minmeshx&dat$y>minmeshy&dat$x<maxmeshx&dat$y<maxmeshy,]
  p <- SDMTools::pnt.in.poly(dat[,c("x","y")], mesh[mesh$cell==x,][,c("X","Y")])
  p$pip[p$pip!=0] <- x
  datje <- merge(dat, p[p$pip!=0,])
  return(datje)
}

pipperframe <- function(dat, mesh, y){
  mesh <- mesh[mesh$frame==y,]
  dat <- dat[[y]]
  print(paste("Finding & saving the raw data per cell for frame", y))
  datjeslist <- lapply(unique(mesh$cell), function(x) pinping(dat, mesh, x))
  datjesframe <- do.call(rbind, datjeslist)
  datjesframe$frame <- y
  return(datjesframe)
}

##################Plotting cells in a tower/row/movie per cell, per frame

plotcellsframelist <- function(TRframe, maxframes, minframes, updown=F, movie=F, viridisoption="magma"){
  nframes <- length(TRframe$frame)
  if(missing(minframes)){
    minframes <- min(TRframe$frame)
  }
  if(missing(maxframes)){
    maxframes <- max(TRframe$frame)
  }
  if(nframes>0){

    minTRframe <- min(TRframe$frame)
    maxTRframe <- max(TRframe$frame)
    if(minTRframe>minframes){
      minframes<-minTRframe
    }
    if(maxframes>maxTRframe){
      maxframes <- maxTRframe
    }


    p <- ggplot2::ggplot(TRframe, aes(frame=frame)) +
      ggplot2::geom_polygon(ggplot2::aes(x=xt,y=yt,fill=values,group=pointN),color=NA) +
      ggimage::theme_transparent() + ggplot2::coord_fixed() +
      viridis::scale_fill_viridis(option=viridisoption) +
      ggplot2::xlim(c(-19,19)) +
      ggplot2::ylim(c(-12, 12)) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL) +
      ggplot2::theme(strip.background=ggplot2::element_rect(fill="transparent", colour=NA),
            strip.text=ggplot2::element_blank(),
            legend.position="none",
            axis.text=ggplot2::element_blank(),
            axis.ticks=ggplot2::element_blank(),
            plot.margin=ggplot2::unit(c(0,0,0,0),"mm"),
            panel.spacing = ggplot2::unit(0, "mm"),
            panel.background=ggplot2::element_rect(fill="transparent", colour=NA),
            plot.background=ggplot2::element_rect(fill="transparent", colour=NA),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank())

    if(updown==T&movie==F){

      p <- p + ggplot2::facet_grid(frame~.)
    }

    if(updown==F&movie==F){
      p <- p + ggplot2::facet_grid(~frame)
    }
    return(p)
  }

  if(nframes==0){
    return(NA)
  }

}

###plot raw image

#' @export
plotRaw <- function(tiffdata, meshdata, pointsdata, frameN=1, xrange, yrange, viridisoption="inferno", meshcolor="white", spotcolor="yellow"){

  plotcells <- ggplot2::ggplot(tiffdata[[frameN]]) + #plot raw image
    ggplot2::geom_raster(ggplot2::aes(x=x,y=y,fill=values)) + #use geom_raster to remake image out of dataframe
    ggplot2::theme_classic() + #simple theme, no backgrounds
    viridis::scale_fill_viridis(option=viridisoption) + #well-working color scheme for gradient values
    ggplot2::theme(legend.position="none") # remove legend for easy viewing

  #add x and/or y range + fixed coordinates when indicated:
  if(missing(xrange)!=T&missing(yrange)!=T){
    plotcells <- plotcells + ggplot2::coord_fixed(xlim=xrange, ylim=xrange)  #sub-set of the image frame to zoom in
  }
  if(missing(xrange)!=T&missing(yrange)==T){
    plotcells <- plotcells + ggplot2::coord_fixed(xlim=xrange)
  }
  if(missing(xrange)==T&missing(yrange)!=T){
    plotcells <- plotcells + ggplot2::coord_fixed(ylim=yrange)
  }
  if(missing(xrange)==T&missing(yrange)==T){
    plotcells <- plotcells + ggplot2::coord_fixed()
  }

  #add mesh data when given:
  if(missing(meshdata)!=T){
    plotcells <- plotcells + #plot made above
      geom_path(data=meshdata[meshdata$frame==frameN,], aes(x=X,y=Y, group=cell), color=meshcolor) #add outline of cells, only frame one, white color
  }
  if(missing(pointsdata)!=T){
    plotcells <- plotcells +
      geom_point(data=pointsdata[pointsdata$frame==frameN,], aes(x=x,y=y), shape=1, color=spotcolor) # add yellow empty dots of our spot localizations on top
  }
  return(plotcells)
}


########plot time series (or by cell length) raw data kymographs
#' @export

prepForKymo <- function(turnedCells, dimension="length", bins = 25, sizeAV=FALSE){

  turnedCells <- turnedCells[,c("values", "frame", "cell", "X_rot","Y_rot", "max.length", "max.width")]
  groupL <- list()
  for(n in unique(turnedCells$frame)){

    U <- turnedCells[turnedCells$frame==n,]
    if(dimension=="length"){
      U$group <- unlist(lapply(unique(U$cell), function(x) cut(U$X_rot[U$cell==x], breaks=bins, labels= c(1:bins))))
    }
    if(dimension=="width"){
      U$group <- unlist(lapply(unique(U$cell), function(x) cut(U$Y_rot[U$cell==x], breaks=bins, labels= c(1:bins))))
    }

    Umeans <- lapply(unique(U$cell), function(x) aggregate(U[U$cell==x,], by=list(U$group[U$cell==x]), FUN=mean, na.rm=T))
    Umeansall <- Umeans[[1]]
    for(x in 2:length(Umeans)){
      Umeansall <- rbind(Umeansall, Umeans[[x]])
    }


    if(sizeAV==TRUE){

      Umeansall$a <- Umeansall$frame + 0.5
      Umeansall$b <- Umeansall$frame - 0.5
      if(dimension=="length"){
        Umeansall$d <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$X_rot[Umeansall$cell==x] + (max(Umeansall$X_rot[Umeansall$cell==x])/bins)*1.15))
        Umeansall$c <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$X_rot[Umeansall$cell==x] - (max(Umeansall$X_rot[Umeansall$cell==x])/bins)*1.15))
      }
      if(dimension=="width"){
        Umeansall$d <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$Y_rot[Umeansall$cell==x] + (max(Umeansall$Y_rot[Umeansall$cell==x])/bins)*1.15))
        Umeansall$c <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$Y_rot[Umeansall$cell==x] - (max(Umeansall$Y_rot[Umeansall$cell==x])/bins)*1.15))
      }


      Umeansall <- reshape2::melt(Umeansall, measure.vars=c("a", "b"), value.name="x_coords")
      Umeansall$frameh <- Umeansall$variable
      Umeansall$variable <- NULL
      Umeansall <- reshape2::melt(Umeansall, measure.vars=c("d", "c"), value.name="y_coords")
      Umeansall$frameh <- as.character(Umeansall$frameh)
      Umeansall$frameh[Umeansall$frameh=="b"&Umeansall$variable == "c"] <- "c"
      Umeansall$frameh[Umeansall$frameh=="a"&Umeansall$variable == "c"] <- "d"
      Umeansall$variable <- NULL
      Umeansall <-  Umeansall[order( Umeansall$frame,  Umeansall$cell,  Umeansall$Group.1,  Umeansall$frameh),]

    }

    Umeansall$group <- Umeansall$Group.1
    Umeansall$Group.1 <- NULL

    groupL[[n]] <- Umeansall

  }

  groupL <- do.call("rbind", groupL)

  groupL <- groupL[order(groupL$max.length),]
  cellnumdatframe<- data.frame(cellnum.length = 1:length(unique(groupL$max.length)), max.length = unique(groupL$max.length))
  groupL <- merge(groupL, cellnumdatframe)

  groupL <- groupL[order(groupL$max.width),]
  cellnumdatframe<- data.frame(cellnum.width = 1:length(unique(groupL$max.width)), max.width=unique(groupL$max.width))
  groupL <- merge(groupL, cellnumdatframe)

  return(groupL)
}


##############plot#########
#' @export

bactKymo <- function(groupL, timeD = FALSE, dimension = "length", bins=25, sizeAV=FALSE, cells="all", prep=TRUE){

  if(prep==TRUE){
    groupL <- prepForKymo(groupL, dimension=dimension, bins=bins, sizeAV=sizeAV)
  }

  if(timeD==FALSE){

    if(dimension=="length"){
      plot1 <- ggplot2::ggplot(groupL, ggplot2::aes(x=cellnum.length,y=group, fill=values)) +
        ggplot2::geom_raster() +
        ggplot2::coord_fixed(ratio=20) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("n(th) cell ordered by cell length") +
        ggplot2::ylab("bin (by cell length)") +
        viridis::scale_fill_viridis(name="Fluorescence\nIntensity")
    }

    if(dimension=="width"){
      plot1 <- ggplot2::ggplot(groupL, ggplot2::aes(x=cellnum.width,y=group, fill=values)) +
        ggplot2::geom_raster() +
        ggplot2::coord_fixed(ratio=20) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("n(th) cell ordered by cell width") +
        ggplot2::ylab("bin (by cell width)") +
        viridis::scale_fill_viridis(name="Fluorescence\nIntensity")
    }
  }

  if(timeD==TRUE&sizeAV==FALSE){
    if(cells=="all"){
      plot1 <- lapply(unique(groupL$cell), function(x) ggplot2::ggplot(groupL[groupL$cell==x,]) +
                        ggplot2::geom_raster(ggplot2::aes(x=frame, y=group, fill=values)) +
                        ggplot2::theme_minimal() +
                        ggplot2::xlab("Time (frames)") +
                        ggplot2::ylab(paste("bin (by cell ", dimension, ")", sep="")) +
                        viridis::scale_fill_viridis()
      )
    }
    if(is.numeric(cells)==TRUE&length(cells)>1){
      plot1 <- lapply(cells, function(x) ggplot2::ggplot(groupL[groupL$cell==x,]) +
                        ggplot2::geom_raster(ggplot2::aes(x=frame, y=group, fill=values)) +
                        ggplot2::theme_minimal() +
                        ggplot2::xlab("Time (frames)") +
                        ggplot2::ylab(paste("bin (by cell ", dimension, ")", sep="")) +
                        viridis::scale_fill_viridis()
      )
    }
    if(is.numeric(cells)==TRUE&length(cells)==1){
      plot1 <- ggplot2::ggplot(groupL[groupL$cell==cells,]) +
        ggplot2::geom_raster(ggplot2::aes(x=frame, y=group, fill=values)) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("Time (frames)") +
        ggplot2::ylab(paste("bin (by cell ", dimension, ")", sep="")) +
        viridis::scale_fill_viridis()
    }
    if(cells!="all" & is.numeric(cells)==FALSE){
      stop("'cells' must be either be a character string 'all', a number or a vector of numbers c(x, x,)")
    }
  }

  if(timeD==TRUE&sizeAV==TRUE){
    if(cells=="all"){
        plot1 <- lapply(unique(ggplot2::groupL$cell), function(x) ggplot2::ggplot(groupL[groupL$cell==x,]) +
                          ggplot2::geom_polygon(ggplot2::aes(x=x_coords, y=y_coords, fill=values, group=X_rot)) +
                          ggplot2::theme_minimal() +
                          ggplot2::xlab("Time (frames)") +
                          ggplot2::ylab(paste("bin (by cell ", dimension, ")", sep="")) +
                          viridis::scale_fill_viridis() +
                          ggtitle(x)
        )
        plot1 <- !is.na(plot1)
    }
    if(is.numeric(cells)==TRUE&length(cells)>1){
      plot1 <- lapply(cells, function(x) ggplot2::ggplot(groupL[groupL$cell==x,]) +
                        ggplot2::geom_polygon(ggplot2::aes(x=x_coords, y=y_coords, fill=values, group=X_rot)) +
                        ggplot2::theme_minimal() +
                        ggplot2::xlab("Time (frames)") +
                        ggplot2::ylab(paste("bin (by cell ", dimension, ")", sep="")) +
                        viridis::scale_fill_viridis() +
                        ggtitle(x)
      )
    }
    if(is.numeric(cells)==TRUE&length(cells)==1){
      plot1 <- ggplot2::ggplot(groupL[groupL$cell==cells,]) +
        ggplot2::geom_polygon(ggplot2::aes(x=x_coords, y=y_coords, fill=values, group=X_rot)) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("Time (frames)") +
        ggplot2::ylab(paste("bin (by cell ", dimension, ")", sep="")) +
        viridis::scale_fill_viridis()
    }
    if(cells!="all" & is.numeric(cells)==FALSE){
      stop("'cells' must be either be a character string 'all', a number or a vector of numbers c(x, x,)")
    }
  }

  return(plot1)

}

