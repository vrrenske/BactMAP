##9-1-2015- last update 13-6-2016
##Renske van Raaphorst

#preparation of dataset containing spots for plotting.

#might be good to clean up first:

#before running the script, open datasets:
#"REP": output of peakfitter.
#"M": the excel output file of the meshes used for peakfitter
#aand start



#!! --> if using .txt extension: make sure the "length" cells are notated as numeric, with the same amount of decimal points
#       for both data frames before importing.
#!! --> also check your decimal seperator for the tab delimited .txt files.

###############################################################################################
#' @export
addPalette <- function(palList, palName){
  if(missing(palName)){
    palName <- readline("Give the name of your new Palette:\n")
  }
  if(missing(palList)){
    low <- readline("Give the color (hex) for the lowest value:\n")
    med <- readline("Give the color (hex) for the midpoint value:\n")
    hig <- readline("Give the color (hex) for the highest value:\n")
    palList <- list(low, med, hig)
  }
  palList <- list(palList)
  names(palList) <- palName
  colopts2 <- append(get(colopts, envir=colEnv), palList)
  assign(colopts, colopts2, envir=colEnv)
  showCurrentPalettes(get(colopts, envir=colEnv))

}

#' @export
showCurrentPalettes <- function(colchoice = get(colopts, envir=colEnv)){
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      inp <- readline("Package 'gridExtra' and 'mvtnorm' needed for this function to work. Press 'y' to install them, or any other key to cancel.")
      if(inp=="y"|inp=="Y"){install.packages(c("gridExtra", "mvtnorm"))}else{stop("Canceled")}
    }
    else{inp <- readline("Package 'gridExtra' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){install.packages("gridExtra")}else{stop("Canceled")}
    }
  }
  if(!requireNamespace("mvtnorm", quietly = TRUE)) {
    inp <- readline("Package 'mvtnorm' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){install.packages("mvtnorm")}else{stop("Canceled")}
  }
  plotlist <- list()
  namelist <- names(colchoice)
  for(n in 1:length(colchoice)){
    p <- colorchoiceplot(colchoice[[n]], n, namelist[n])
    plotlist <- append(plotlist,list(p))
  }
  return(do.call(gridExtra::grid.arrange,c(plotlist, ncol=2)))
}

colorchoiceplot <- function(colchoice, nums, pname){

  z <- mvtnorm::rmvnorm(100, mean=c(3,5), sigma=matrix(c(1,0.5,0.5,2), nrow=2))
  z <- data.frame(z)
  return(ggplot2::ggplot(z, ggplot2::aes(x=X1, y=X2)) + ggplot2::stat_density2d(ggplot2::aes(fill=..density..), geom="raster", contour=FALSE) + ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint=0.06) + ggplot2::theme_minimal() + ggplot2::theme(legend.position="none") + ggplot2::ggtitle(pname) + ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::xlim(0,5) + ggplot2::ylim(0,10) + ggplot2::theme(axis.text=ggplot2::element_blank()))

}

#' @export
pal2Default <- function(){
  assign(colopts, cols, envir=colEnv)
  message("Colorpalettes back to default")
}

#wat basisplotfuncties
densityplot <- function(plot){
  return(plot + ggplot2::stat_density2d(ggplot2::aes(fill=..density..), geom="raster", contour = FALSE))
}

LWplot <- function(plot, u="black", maxn){
  return(plot + #ggplot2::geom_rect(xmin=0, xmax=maxn, ymin=-1.5, ymax=1.5, fill=u)
           ggplot2::theme_minimal()  + ggplot2::scale_x_continuous(limits=c(0,maxn)))
}

heatmap <- function(pdens, mp, colchoice, viridis = F){
  if(viridis==F){
    return(pdens + ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint =mp, space = "Lab") + ggplot2::theme_minimal() + ggplot2::theme(legend.position="none", panel.background=ggplot2::element_rect(fill=colchoice[1])))
  }
  if(viridis==T){
    return(pdens + ggplot2::scale_fill_viridis_c(option=colchoice) + ggplot2::theme_minimal() + ggplot2::theme(legend.position="none", panel.background=ggplot2::element_rect(fill="black")))
  }
  }

#function for goodlooking x/y coordinate plot:
#makes a plot sized as the max cell (width/length: xmax) of the quartile inside a plot
#which has a grey background as large as the largest cell in the dataset
#so all quartile plots will have the same dimensions.
#the title, y axis and x axis will also be drawn.
coplot <- function(pheat, xmax, ymax, u="black"){
 return(pheat + ggplot2::xlab("Length (\u00B5m)") + ggplot2::ylab("Width (\u00B5m)") + ggplot2::coord_fixed(xlim = c(-xmax,xmax), ylim=c(-ymax,ymax)) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = u), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()))
}

plotlocation_histograms <- function(){}



######################### mid-points ################################################################################

superfun <- function(dat, bins,mag){
  dat <- dat[order(dat$frame, dat$cell, dat$num),]
  dat <- unique(dat)
  dat$av <- 0
  if("max.length"%in%colnames(dat)!=T){
    dat$max.length <- dat$length
  }
  dat$av <- dat$X_rot/dat$max.length*100
  dat <- dat[!is.na(dat$av),]

  cutpoints<-quantile(dat$av,(0:bins)/bins)

  dat$binned1 <-cut(dat$av,cutpoints, include.lowest=TRUE, labels = 1:bins)

  xb <- split(dat$X_rot, dat$binned)
  ybp <- split(dat$Y_rot[dat$Y_rot>0], dat$binned[dat$Y_rot>0])
  ybm <- split(dat$Y_rot[dat$Y_rot<0], dat$binned[dat$Y_rot<0])
  xmeans <- lapply(xb, function(x)mean(x))
  ypmeans <- lapply(ybp, function(x)mean(x))
  ymmeans <- lapply(ybm, function(x)mean(x))
  Lb <- split(dat$max.length, dat$binned1)
  maxLmeans <- lapply(Lb, function(x)mean(x))
  meanframe <- data.frame(x=unlist(xmeans), yp=unlist(ypmeans), ym=unlist(ymmeans),max.length=unlist(maxLmeans))
  meanframe <- reshape2::melt(meanframe, c("x","max.length"), value.name="y")
  meanframe$variable <- NULL
  meanframe <- meanframe * mag
  meanframe <- meanframe[abs(meanframe$x)<=0.5*meanframe$max.length,]
  #add num as an order identifyer for drawing a path or a polygon as cell
  meanframe$num <- c(1:(nrow(meanframe)/2), nrow(meanframe):(nrow(meanframe)/2+1))
  meanframe <- rbind(meanframe, meanframe[1,])
  meanframe[nrow(meanframe),]$num <- nrow(meanframe)
  meanframe <- meanframe[order(meanframe$num),]
  return(meanframe)
}

#or two, by quartiles of the number of cells:
#' @export
createPlotlist <- function(spotdata,  meshdata, groups =4 , colorpalette="GreenYellow", mag="No_PixelCorrection", AllPlot=F, Xm="X", Ym="Y", viridis=FALSE){
  if (!requireNamespace("MASS", quietly = TRUE)) {
  inp_P <- readline("Package 'MASS' needed for this function to work. Press 'y' to install, or any other key to cancel.")
  if(inp_P=="y"|inp_P=="Y"){install.packages("MASS")}else{stop("Canceled")}
  }
  if(AllPlot==T){
    if (!requireNamespace("grid", quietly = TRUE)) {
      inp_P <- readline("Package 'grid' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp_P=="y"|inp_P=="Y"){install.packages("grid")}else{stop("Canceled")}
    }
  }
  if(missing(mag)!=T&is.numeric(unlist(get(magnificationList,envir=magEnv)[mag]))==FALSE){
    stop("Magnification conversion factor not recognized. Please use addPixels2um('pixelName', pixelsize) to add your conversion factor")
  }
  if("Lmid"%in%colnames(spotdata)==T){
    MR <- spotdata
    MR$cellnum <- MR$num
  }
  if("l"%in%colnames(spotdata)==T&"Lmid"%in%colnames(spotdata)==F){
    MR <- LimDum(spotdata, pix2um = unlist(get(magnificationList, envir=magEnv)[mag]))
    MR <- spotMR(MR)
    #if(length(which(unique(meshdata$cell)%in%spotdata$cell!=T))==0){MR <- spotdata}else{
      MR <- merge(MR, unique(meshdata[,c("cell",
                                              "frame",
                                              "max.width",
                                              "max.length")]),
                  all=TRUE)
      #}
    MR$totalspot[is.na(MR$totalspot)] <- 0
    MR <- MR[order(MR$max.length, MR$max.width),]
    MR <- MR[!is.na(MR$cell),]
    MR$num <- c(1:nrow(MR))
    MR$cellnum <- MR$num
  }
  if("Lmid"%in%colnames(spotdata)==F&"l"%in%colnames(spotdata)==F){
    message("Did not find cell data in spot dataset. Running spotsInBox to connect spots to meshes...")
    combineframes <- spotsInBox(spotdata, meshdata, Xm=Xm,Ym=Ym)
    meshdata <- combineframes$mesh
    spotdata <- combineframes$spots_relative

    MR <- LimDum(spotdata, pix2um = unlist(get(magnificationList, envir=magEnv)[mag]))
    MR <- spotMR(MR)
    #if(length(which(unique(meshdata$cell)%in%spotdata$cell!=T))==0){MR <- spotdata}else{
    MR <- merge(MR, unique(meshdata[,c("cell",
                                       "frame",
                                       "max.width",
                                       "max.length")]),
                all=TRUE)
    #}

    MR$totalspot[is.na(MR$totalspot)] <- 0
    MR <- MR[order(MR$max.length, MR$max.width),]
    MR <- MR[!is.na(MR$cell),]
    MR$num <- c(1:nrow(MR))
    MR$cellnum <- MR$num
  }
  if("obID"%in%colnames(spotdata)){
    MR <- unique(MR[,c("num", "frame", "cell", "max.length", "max.width", "Dum", "Lmid", "pole1", "pole2", "cellnum")])
    MR <- MR[!is.na(MR$cell),]
  }

  if(viridis==FALSE){
    colc <- get(colopts, envir = colEnv)[colorpalette][[1]]
  }

  if(viridis==TRUE){
    colc <- colorpalette
  }

  MR$q1 <- cut(MR$cellnum, breaks=groups, labels = 1:groups)


  xmax <- 0.5*max(MR$max.length, na.rm=TRUE)
  ymax <- 0.5*max(MR$max.width, na.rm=TRUE)

  #Length and width corrected for average length per quartile for plotting the coordinate plots
  mframe <- suppressWarnings(aggregate(MR, by=list(MR$q1), FUN=mean, na.rm=TRUE, na.action=na.omit))
  meansL <- mframe$max.length
  meansW <- mframe$max.width

  #length correction
  MR$Lcor <- (MR$Lmid/MR$max.length*meansL[as.numeric(MR$q1)])
  #width correction
  MR$Dcor <- (MR$Dum/MR$max.width*meansW[as.numeric(MR$q1)])

###############################################################################################################################################
  ##plotting! -> coordinate plots
  allMRs <- split(MR, MR$q1)
  plist <- list()
  for(n in 1:groups){
    p <- ggplot2::ggplot(allMRs[[n]], ggplot2::aes(x=Lcor, y=Dcor))
    p <- densityplot(p)
    plist <- append(plist, list(p))
  }

#plotting! -> L and D ordered by cell length
  pL <- ggplot2::ggplot(MR, ggplot2::aes(x=num, y=Lmid))
  if(viridis==TRUE){
    pL <- LWplot(pL, "black", max(MR$num, na.rm=T))
  }
  if(viridis==FALSE){
    pL <- LWplot(pL, colc[[1]], max(MR$num, na.rm=T))
  }
  #pLpoint <- pL + ggplot2::geom_point() + ggplot2::ggtitle("Spot location on length axis ordered by cell length") + ggplot2::xlab("Spot - ordered by cell length") + ggplot2::ylab("Y-position (\u03BCm)") + ggplot2::theme_bw()
  pLD <- densityplot(pL) + ggplot2::ggtitle("Spot location on length axis ordered by cell length") + ggplot2::xlab("Spot - ordered by cell length") + ggplot2::ylab("Y-position (\u03BCm)") + ggplot2::geom_line(data=MR, ggplot2::aes(x=num,y=pole1),colour="white") + ggplot2::geom_line(data=MR, ggplot2::aes(x=num,y=pole2),colour="white")

  pW <- ggplot2::ggplot(MR, ggplot2::aes(x=num, y=Dum))
  if(viridis == TRUE){
    pW <- LWplot(pW, "black", max(MR$num,na.rm=T))
  }
  if(viridis==FALSE){
    pW <- LWplot(pW, colc[[1]], max(MR$num,na.rm=T))
  }
  #pWpoint <- pW + ggplot2::geom_point() + ggplot2::ggtitle("Spot location on width axis ordered by cell length") + ggplot2::xlab("Spot - ordered by cell length") + ggplot2::ylab("X-position (\u03BCm)") + ggplot2::theme_bw()
  pWD <- densityplot(pW) + ggplot2::ggtitle("Spot location on width axis ordered by cell length") + ggplot2::xlab("Spot - ordered by cell length") + ggplot2::ylab("X-position (\u03BCm)")

#make heatmap using the half max densities:

  ##for all quartiles (or n-tiles)
  allMRLmids <- split(MR[,c("Lmid","Dum")][!is.na(MR$Lmid),], MR$q1[!is.na(MR$Lmid)])
  #get maximum half max density of all groups and apply for each plot to have same scaling.
  mp<- max(sapply(allMRLmids, function(x) median(range(MASS::kde2d(x$Lmid,x$Dum)$z))))
  #create n-tile plots saved in a list.
  phlist <- list()
  for(n in 1:groups){
    p <- heatmap(plist[[n]], mp, colc, viridis)
    p <- coplot(p, xmax, ymax)
    phlist <- append(phlist, list(p))
  }

  u <- list()
  #do the same procedure for length/width heatmap plots
  mpL <- MASS::kde2d(MR$num[!is.na(MR$Lmid)], MR$Lmid[!is.na(MR$Lmid)])
  mpL1 <- median(range(mpL$z))
  pLD <- heatmap(pLD, mpL1, colc,viridis)
  if(viridis==FALSE){
    pLD <- pLD + ggplot2::theme_minimal() + ggplot2::theme(panel.background=ggplot2::element_rect(fill=colc[[1]]), panel.grid = ggplot2::element_blank())
  }
  if(viridis==TRUE){
    pLD <- pLD + ggplot2::theme_minimal() + ggplot2::theme(panel.background=ggplot2::element_rect(fill="black"), panel.grid = ggplot2::element_blank())
  }
  u$lengthplot <- pLD

  mpW <- MASS::kde2d(MR$num[!is.na(MR$Dum)], MR$Dum[!is.na(MR$Dum)])
  mpW1 <- median(range(mpW$z))
  pWD <- heatmap(pWD, mpW1, colc,viridis)
  if(viridis==TRUE){
    pWD <- pWD + ggplot2::theme_minimal() + ggplot2::theme(panel.background=ggplot2::element_rect(fill="black"), panel.grid = ggplot2::element_blank())
  }

  if(viridis==FALSE){
    pWD <- pWD + ggplot2::theme_minimal() + ggplot2::theme(panel.background=ggplot2::element_rect(fill=colc[[1]]), panel.grid = ggplot2::element_blank())
  }
  u$widthplot <- pWD

  if(missing(meshdata)){
    u$qplots <- phlist
  }

  #################add meshdata data###############################################################################################
  if(!missing(meshdata)){

    if("X_rot"%in%colnames(meshdata)!=T){
    message("Turning the cells... this may take a while.")

    meshdata <- meshTurn(meshdata, Xm, Ym)
    message("Finished turning the cells")
    }
    p2um <- as.numeric(get(magnificationList, envir=magEnv)[mag])
    if("max_um"%in%colnames(meshdata)!=T){
    meshdata$max_um <- meshdata$max.length*p2um
    meshdata$maxwum <- meshdata$max.width*p2um
    }


    meshdata <- merge(meshdata, MR[,c("cell", "frame", "q1")], all=T)
    MESHlist <- split(meshdata, meshdata$q1)
    message("Calculating mean cell outlines..")
    if(nrow(unique(meshdata[meshdata$max.length==max(meshdata$max.length),]))< 12){
      means <- lapply(MESHlist, function(x) superfun(x, 12, p2um))
    }
    if(nrow(unique(meshdata[meshdata$max.length==max(meshdata$max.length),]))>11){
     means <- lapply(MESHlist, function(x)superfun(x, round(min(x$max.length), digits=0), p2um))
    }
    u$mean_outlines <- means
    message("Finished calculating mean cell outlines")

    phmlist <- list()
    for(n in 1:groups){
      p <- phlist[[n]] + ggplot2::geom_path(data=means[[n]], ggplot2::aes(x=x,y=y), colour="white") + ggplot2::coord_fixed(xlim=c(min(means[[groups]]$x)*1.2, max(means[[groups]]$x*1.2)), ylim=c(min(means[[groups]]$y)*1.2, max(means[[groups]]$y*1.2)))
      phmlist <- append(phmlist, list(p))
    }
    if(AllPlot==F){
      u$qplots <- phmlist
    }
  }
    #and create the plots:
    if(AllPlot==T){
      phmalist <- list()
      for(n in 1:groups){
        if(missing(meshdata)!=TRUE){
          p1_all <- suppressMessages(allplot(phmlist[[n]], allMRs[[n]], max(means[[groups]]$x*1.2), max(means[[groups]]$y*1.2), empty))
        }
        if(missing(meshdata)==TRUE){
          p1_all <- suppressMessages(allplot(phlist[[n]], allMRs[[n]], mean(allMRs[[groups]]$max.length)*1.2, mean(allMRs[[groups]]$max.width)*1.2, empty))
        }
        phmalist <- append(phmalist,list(p1_all))
        u$qplots <- phmalist
      }
    }

    #in case you want it for the whole thing instead of only quartiles:

    pall <- ggplot2::ggplot(MR, ggplot2::aes(x=Lcor, y=Dcor))
    pall <- densityplot(pall)
    mppall <- mean(range(MASS::kde2d(MR$Lcor[!is.na(MR$Lcor)&!is.na(MR$Dcor)],MR$Dcor[!is.na(MR$Dcor)&!is.na(MR$Lcor)])$z))
    pall <- heatmap(pall, mppall, colc, viridis)


    if(nrow(unique(meshdata[meshdata$max.length==max(meshdata$max.length),]))< 12){
      meantotal <- superfun(meshdata, 12, p2um)
    }
    if(nrow(unique(meshdata[meshdata$max.length==max(meshdata$max.length),]))> 11){
      meantotal <- superfun(meshdata, 30, p2um)
    }
    pall <- coplot(pall,max(meantotal$x)*1.2, max(meantotal$y)*1.2)
    pall <- pall + ggplot2::geom_path(data=meantotal, ggplot2::aes(x=x,y=y), colour="white")
    if(AllPlot==F){
      u$plottotal <- pall
    }
    if(AllPlot==T){
    pall_all <- suppressMessages(allplot(pall, MR, max(meantotal$x)*1.2, max(meantotal$y)*1.2, empty))
    u$plottotal <- pall_all

    }




  hislist <- list()
  #save all histograms (L coordinates) of the quartiles too:
  for(n in 1:groups){
    if(viridis==FALSE){
      p1his <- ggplot2::ggplot(allMRs[[n]], ggplot2::aes(x=Lcor)) + ggplot2::geom_density(fill=colc[[2]], color=colc[[2]]) + ggplot2::theme_minimal() + ggplot2::labs(x="Length(\u03BCm)") + ggplot2::coord_cartesian(xlim=c(min(MR$Lcor, na.rm=T),max(MR$Lcor, na.rm=T)))
    }
    if(viridis==TRUE){
      singopt <- viridissinglecols[colc]
      p1his <- ggplot2::ggplot(allMRs[[n]], ggplot2::aes(x=Lcor)) + ggplot2::geom_density(fill=singopt, color=singopt) + ggplot2::theme_minimal() + ggplot2::labs(x="Length(\u03BCm)") + ggplot2::coord_cartesian(xlim=c(min(MR$Lcor,na.rm=T), max(MR$Lcor,na.rm=T)))

    }
    hislist <- append(hislist, list(p1his))
  }

  u$histograms <- hislist
  u$spotdata <- MR
  if(!missing(meshdata)){
  	u$meshdata <- meshdata
  }
  u$pixel2um <- p2um
  u$data_summary <- makePlotListSummary_2(MR, groups=groups)
  message("Done plotting.")
  return(u)
}

###########################################double histograms!! yay!!!######################################3333
##allplot function combines histograms and density plot.

allplot <- function(plot, data, xmax, ymax, empty){

  #prepare seperate plots: histograms (hL, hD) and modified coordinate plots(remove legend )
  p1D <- plot + ggplot2::theme(legend.position = "none")
  p1hL <- ggplot2::ggplot(data, ggplot2::aes(x=Lcor)) + ggplot2::geom_histogram(fill="black") + ggplot2::coord_cartesian(xlim = c(-xmax, xmax)) + ggplot2::theme_minimal() + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  p1hD <- ggplot2::ggplot(data, ggplot2::aes(x=Dcor)) + ggplot2::geom_histogram(fill="black") + ggplot2::coord_flip(xlim = c(-ymax, ymax)) + ggplot2::theme_minimal() + ggplot2::theme(axis.title.y = ggplot2::element_blank())

  #align the plots properly before putting them together
  p1Dg <- ggplot2::ggplotGrob(p1D)
  p1hLg <- ggplot2::ggplotGrob(p1hL)
  p1hDg <- ggplot2::ggplotGrob(p1hD)

  maxWidth = grid::unit.pmax(p1Dg$widths[2:5], p1hLg$widths[2:5])
  p1Dg$widths[2:5] <- as.list(maxWidth)
  p1hLg$widths[2:5] <- as.list(maxWidth)

  #put the grids together using gridarrange
  return(gridExtra::arrangeGrob(p1hLg, empty, p1Dg, p1hD, ncol=2, nrow=2, widths=c(10*xmax, 2.5*xmax), heights=c(2.5*xmax, 10*ymax)))
}

##before using the function:
#create mockup plot to make space
empty <- ggplot2::ggplot(data.frame(u=1), ggplot2::aes(u,u)) +
  ggplot2::theme(
    plot.background = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )

###################################################

# Function to change the color and/or "groups" of createPlotList function.

#' @export
changePlotlist <- function(plotlist, changecolor = TRUE, changegrouping = FALSE, colorpalette="ColdBlue", viridis=FALSE, groups = 4){
  if(changecolor == TRUE){
    if(viridis==TRUE){
      colchoice <- colorpalette
    }
    #for the plots which just need a new palette; find mid point for manual gradient option & add new palette.
    if(viridis==FALSE){
      colchoice <- get(colopts, envir = colEnv)[colorpalette][[1]]
      #lengthplot
      mpL <- MASS::kde2d(plotlist$spotdata$num[!is.na(plotlist$spotdata$Lmid)], plotlist$spotdata$Lmid[!is.na(plotlist$spotdata$Lmid)])
      mpL1 <- median(range(mpL$z))
      plotlist$lengthplot <- suppressMessages(plotlist$lengthplot + ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint =mpL1, space = "Lab"))
      #widthplot
      mpW <- MASS::kde2d(plotlist$spotdata$num[!is.na(plotlist$spotdata$Dum)], plotlist$spotdata$Dum[!is.na(plotlist$spotdata$Dum)])
      mpW1 <- median(range(mpW$z))
      plotlist$widthplot <- suppressMessages(plotlist$widthplot + ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint =mpW1, space = "Lab"))
      #totalplot
      if(ggplot2::is.ggplot(plotlist$plottotal)==TRUE){
        mpA <- mean(range(MASS::kde2d(plotlist$spotdata$Lcor[!is.na(plotlist$spotdata$Lcor)&!is.na(plotlist$spotdata$Dcor)],plotlist$spotdata$Dcor[!is.na(plotlist$spotdata$Dcor)&!is.na(plotlist$spotdata$Lcor)])$z))
        plotlist$plottotal <- suppressMessages(plotlist$plottotal + ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint =mpA, space = "Lab"))
      }
      #listofplots
      if(changegrouping ==FALSE & ggplot2::is.ggplot(plotlist$qplots[[1]])==TRUE){
        plotlist$qplots <- suppressMessages(lapply(plotlist$qplots, function(x) x + ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint =mpA, space = "Lab")))
      }

    }
  }

    #remake totalplot when it's a Grob. heatmap function takes viridis so no need to separate on this.
    if(ggplot2::is.ggplot(plotlist$plottotal)!=TRUE){
      pall <- ggplot2::ggplot(plotlist$spotdata, ggplot2::aes(x=Lcor, y=Dcor))
      pall <- densityplot(pall)
      mppall <- mean(range(MASS::kde2d(plotlist$spotdata$Lcor[!is.na(plotlist$spotdata$Lcor)&!is.na(plotlist$spotdata$Dcor)],plotlist$spotdata$Dcor[!is.na(plotlist$spotdata$Dcor)&!is.na(plotlist$spotdata$Lcor)])$z))
      pall <- heatmap(pall, mppall, colchoice, viridis)
      pall <- coplot(pall,max(plotlist$spotdata$max.length)*0.5, max(plotlist$spotdata$max.width)*0.5)
      if(length(plotlist$meshdata$X)!=0){
        if(nrow(unique(plotlist$meshdata[plotlist$meshdata$max.length==max(plotlist$meshdata$max.length),]))< 12){
          meantotal <- superfun(plotlist$meshdata, 12, plotlist$pixel2um)
        }
        if(nrow(unique(plotlist$meshdata[plotlist$meshdata$max.length==max(plotlist$meshdata$max.length),]))> 11){
          meantotal <- superfun(plotlist$meshdata, 30, plotlist$pixel2um)
        }

        pall <- pall + ggplot2::geom_path(data=meantotal, ggplot2::aes(x=x,y=y), colour="white")
      }

      pall_all <- suppressMessages(allplot(pall, plotlist$spotdata, max(plotlist$spotdata$max.length)*0.5, max(plotlist$spotdata$max.width)*0.5, empty))
      plotlist$plottotal <- pall_all
    }

  if(changecolor==TRUE){
    #viridis option; bit easier, just add the palette with the option.
    if(viridis==TRUE){
      #lengthplot
      plotlist$lengthplot <- suppressMessages(plotlist$lengthplot + ggplot2::scale_fill_viridis_c(option=colorpalette))
      #widthplot
      plotlist$widthplot <- suppressMessages(plotlist$widthplot + ggplot2::scale_fill_viridis_c(option=colorpalette))
      #totalplot
      if(ggplot2::is.ggplot(plotlist$plottotal)==TRUE){
        plotlist$plottotal <- suppressMessages(plotlist$plottotal + ggplot2::scale_fill_viridis_c(option=colorpalette))
      }
      #listofplots
      if(changegrouping ==FALSE & ggplot2::is.ggplot(plotlist$qplots[[1]])==TRUE){
        plotlist$qplots <- suppressMessages(lapply(plotlist$qplots, function(x) x + ggplot2::scale_fill_viridis_c(option=colorpalette)))
      }
    }
  }

  if(changegrouping==TRUE){
    plotlist$spotdata$q1 <-  cut(plotlist$spotdata$cellnum, breaks=groups, labels = 1:groups)
    if(length(plotlist$meshdata$X)!=0){
      plotlist$meshdata$q1 <- NULL
      plotlist$meshdata <- merge(plotlist$meshdata, plotlist$spotdata[,c("cell", "frame", "q1")])
    }

  }


  #list of qplots. make in same way as allplot. only run if either the grouping has been changed or when the plot has to be remade becaus it was an AllPlot
  if(changegrouping==TRUE|ggplot2::is.ggplot(plotlist$qplots[[1]])!=TRUE){
    spotsplit <- split(plotlist$spotdata, plotlist$spotdata$q1)
    spotsplitLmids <- split(plotlist$spotdata[,c("Lmid","Dum")][!is.na(plotlist$spotdata$Lmid),], plotlist$spotdata$q1[!is.na(plotlist$spotdata$Lmid)])
    #get maximum half max density of all groups and apply for each plot to have same scaling.
    mp <- max(sapply(spotsplitLmids, function(x) median(range(MASS::kde2d(x$Lmid,x$Dum)$z))))

    #if mesh is there; create new mean outlines:
    if(length(plotlist$meshdata$X)!=0){
      MESHlist <- split(plotlist$meshdata, plotlist$meshdata$q1)
      if(nrow(unique(plotlist$meshdata[plotlist$meshdata$max.length==max(plotlist$meshdata$max.length),]))< 12){
        means <- lapply(MESHlist, function(x) superfun(x, 12, plotlist$pixel2um))
      }
      if(nrow(unique(plotlist$meshdata[plotlist$meshdata$max.length==max(plotlist$meshdata$max.length),]))>11){
        means <- lapply(MESHlist, function(x)superfun(x, round(min(x$max.length), digits=0), plotlist$pixel2um))
      }
      plotlist$mean_outlines <- means
    }
    #create 1:groups plots saved in a list.
    phlist <- list()

    xmax <- 0.5*max(plotlist$spotdata$max.length, na.rm=TRUE)
    ymax <- 0.5*max(plotlist$spotdata$max.width, na.rm=TRUE)

    for(n in 1:groups){
      p <- ggplot2::ggplot(spotsplit[[n]], ggplot2::aes(x=Lcor, y=Dcor))
      p <- densityplot(p)
      p <- heatmap(p, mp, colchoice, viridis)
      p <- coplot(p, mean(plotlist$spotdata$max.length[plotlist$spotdata$q1==groups])*0.6, mean(plotlist$spotdata$max.width[plotlist$spotdata$q1==groups])*0.6)
      if(length(plotlist$meshdata$X)!=0){
        p <- p + ggplot2::geom_path(data=plotlist$mean_outlines[[n]], ggplot2::aes(x=x,y=y), colour="white")
      }
      phlist <- append(phlist, list(p))
    }

    if(ggplot2::is.ggplot(plotlist$qplots[[1]])==TRUE){
      plotlist$qplots <- phlist
    }
    if(ggplot2::is.ggplot(plotlist$qplots[[1]])!=TRUE){
      phmalist <- list()
      for(n in 1:groups){
          p1_all <- suppressMessages(allplot(phlist[[n]], spotsplit[[n]], mean(spotsplit[[groups]]$max.length)*0.6, mean(spotsplit[[groups]]$max.width)*0.6, empty))
        phmalist <- append(phmalist,list(p1_all))
      }
      plotlist$qplots <- phmalist
    }
  }


  if(viridis==FALSE){
      #listofhistograms
      plotlist$histograms <- ggplot2::ggplot(plotlist$spotdata, ggplot2::aes(x=Lcor)) + ggplot2::geom_density(fill=colchoice[2], color=colchoice[2]) + ggplot2::theme_minimal() + ggplot2::facet_grid(q1~.)
    }

  if(viridis==TRUE){
      #listofhistograms
      plotlist$histograms <- suppressWarnings(ggplot2::ggplot(plotlist$spotdata, ggplot2::aes(x=Lcor)) + ggplot2::geom_density(fill=as.character(viridissinglecols[colchoice]), color=as.character(viridissinglecols[colchoice])) + ggplot2::theme_minimal() + ggplot2::facet_grid(q1~.))
    }

  return(plotlist)
}

makePlotListSummary_1 <- function(spotdata, timelapse=FALSE){

  totalcells <- nrow(unique(spotdata[,c("cell", "frame")]))
  localization_x <- summary(spotdata$Lmid)
  localization_y <- summary(spotdata$Dum)
  max_length <- summary(unique(spotdata[,c("cell", "frame","max_um")])$max_um)
  max_width <- summary(unique(spotdata[,c("cell", "frame","maxwum")])$maxwum)
  spotscell <- summary(as.factor(unique(spotdata[,c("cell", "frame","totalspot")])$totalspot))
  return(list("Total_amount_of_cells" = totalcells,
                   "Spot_distance_from_midcell_length_axis" = localization_x,
                   "Spot_distance_from_midcell_width_axis" = localization_y,
                   "Cell_length" = max_length,
                   "Cell_width" = max_width,
                   "Distribution_spots_per_cell" = spotscell))
}

makePlotListSummary_2 <- function(spotdata, groups){
  if(groups==1){
    return(makePlotListSummary_1(spotdata))
  }
  if(groups>1){
    Full_dataset <- makePlotListSummary_1(spotdata)
    grouplist <- lapply(c(1:groups), function(x) makePlotListSummary_1(spotdata[spotdata$q1==x,]))
    names(grouplist) <- paste("Group", c(1:groups), sep="_")
    grouplist$Full_dataset <- Full_dataset
    return(grouplist)
  }
}

###################################################

# color palette lists & options. possible to add one yourself manually

cols <-  list(OrangeHot=list("#000000", "#D55E00", "#F0E442"), GreenYellow = list("#000000", "#009E73", "#F0E442"), ColdBlue = list("#000000", "#0072B2", "#FFFFFF"), YellowHot = list("#000000", "#e69f00", "#F0E442"), RedHot = list("#000000", "#FF0000", "#FFFF00"), WhiteOrange = list("#FFFFFF", "#F0E442", "#D55E00"))
colopts <- "colopts"
colEnv <- new.env()
assign(colopts, cols, envir=colEnv)


singcollist <- list("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#D55E00", "#CC79A7")

#'@export
colorsCUD <- function(colornames){
  if(!missing(colornames)){
    outlist <- singcollist[1:length(colornames)]
    names(outlist) <- colornames
    return(unlist(outlist))
  }
  if(missing(colornames)){
   return(unlist(singcollist))
  }
}

xaxislist <- c("Cell Length (\u03BCm)", "Cell Width (\u03BCm)", "Cell Area (\u03BCm\u00B2)",
               "Spot location in length axis (\u03BCm)", "Spot location on width axis (\u03BCm)")


viridissinglecols <- list(magma = "#B63679FF", A = "#B63679FF", inferno = "#BB3754FF", B ="#BB3754FF", plasma = "#CC4678FF", C = "#CC4678FF", cividis = "#7C7B78FF", E = "#7C7B78FF", viridis = "#21908CFF", D = "#21908CFF")

