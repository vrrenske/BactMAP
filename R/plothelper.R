##9-1-2015- last update 13-6-2016
##Renske van Raaphorst

#preparation of dataset containing spots for plotting.

#might be good to clean up first:

#before running the script, open datasets:
#"REP": output of peakfitter.
#"M": the excel output file of the meshes used for peakfitter
#aand start
library(ggplot2)
library(gridExtra)
library(scales)
library(ggthemes)
library(mvtnorm)
library(MASS)



#!! --> if using .txt extension: make sure the "length" cells are notated as numeric, with the same amount of decimal points
#       for both data frames before importing.
#!! --> also check your decimal seperator for the tab delimited .txt files.

###############################################################################################
#' @export
addPalette <- function(palList, palName, env = parent.frame()){
  if(missing(palName)){
    palName <- readline("Give the name of your new Palette:\n")
  }
  if(missing(palList)){
    low <- readline("Give the color (hex) for the lowest value:\n")
    med <- readline("Give the color (hex) for the midpoint value:\n")
    hig <- readline("Give the color (hex) for the highest value:\n")
    palList <- list(low, med, hig)
  }
  palList <-
  colopts2 <- append(colopts, palList)
  names(colopts2)[length(colopts2)] <- palName
  env$colopts <- colopts2
  showCurrentPalettes()

}

#' @export
showCurrentPalettes <- function(colchoice = colopts){
  plotlist <- list()
  for(n in 1:length(colchoice)){
    p <- colorchoiceplot(colchoice[[n]], n)
    plotlist <- append(plotlist,list(p))
  }
  return(do.call(gridExtra::grid.arrange,c(plotlist, ncol=2)))
}

colorchoiceplot <- function(colchoice, nums){
  namelist <- names(colopts)
  z <- rmvnorm(100, mean=c(3,5), sigma=matrix(c(1,0.5,0.5,2), nrow=2))
  z <- data.frame(z)
  return(ggplot2::ggplot(z, ggplot2::aes(x=X1, y=X2)) + ggplot2::stat_density2d(ggplot2::aes(fill=..density..), geom="raster", contour=FALSE) + ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint=0.06) + ggplot2::theme_minimal() + ggplot2::theme(legend.position="none") + ggplot2::ggtitle(namelist[nums]) + ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::xlim(0,5) + ggplot2::ylim(0,10) + ggplot2::theme(axis.text=ggplot2::element_blank()))

}

#wat basisplotfuncties
densityplot <- function(plot){
  return(plot + ggplot2::stat_density2d(ggplot2::aes(fill=..density..), geom="raster", contour = FALSE))
}

LWplot <- function(plot, u="black", maxn){
  return(plot + #ggplot2::geom_rect(xmin=0, xmax=maxn, ymin=-1.5, ymax=1.5, fill=u)
           ggplot2::theme_minimal()  + ggplot2::scale_x_continuous(limits=c(0,maxn)))
}

heatmap <- function(pdens, mp, colchoice, u = "black", viridis = F){
  if(viridis==F){
  return(pdens + ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint =mp, space = "Lab") + theme_minimal() + theme(legend.position="none", panel.background=ggplot2::element_rect(fill=colchoice[1])))
  }
  if(viridis==T){
    return(pdens + viridis::scale_fill_viridis(option=colchoice,  space = "Lab") + theme_minimal() + theme(legend.position="none", panel.background=ggplot2::element_rect(fill=u)))
  }
  }

#function for goodlooking x/y coordinate plot:
#makes a plot sized as the max cell (width/length: xmax) of the quartile inside a plot
#which has a grey background as large as the largest cell in the dataset
#so all quartile plots will have the same dimensions.
#the title, y axis and x axis will also be drawn.
coplot <- function(pheat, xmax, ymax, xqmax, u="black"){
 return(pheat + ggplot2::xlab("Length (\u00B5m)") + ylab("Width (\u00B5m)") + coord_cartesian(xlim = c(-xmax,xmax), ylim=c(-ymax,ymax)) + geom_vline(xintercept = xqmax) + geom_vline(xintercept=-xqmax) + geom_hline(yintercept = ymax) + geom_hline(yintercept = -ymax) + theme(panel.background = element_rect(fill = u), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
}

plotlocation_histograms <- function(){}



######################### mid-points ################################################################################

superfun <- function(dat, bins,mag){
  dat <- dat[order(dat$frame, dat$cell, dat$num, dat$max.length, dat$xy),]
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
  return(meanframe)
}

#or two, by quartiles of the number of cells:
#' @export
createPlotlist <- function(REP, inp, MESH, colorpalette="GreenYellow", mag="100x_LeicaVeening", AllPlot=T, Xm="X", Ym="Y", viridis=FALSE){
  if("Lmid"%in%colnames(REP)!=T){
  MR <- mergeframes(REP, MESH, mag)
  }
  if("Lmid"%in%colnames(REP)==T){
    MR <- REP
    rm(REP)
  }

  colc <- colopts[colorpalette][[1]]
  MR$q1 <- cut(MR$cellnum, breaks=inp, labels = 1:inp)


  xmax <- 0.5*max(MR$max.length, na.rm=TRUE)
  ymax <- 0.5*max(MR$max.width, na.rm=TRUE)

  #Length and width corrected for average length per quartile for plotting the coordinate plots
  mframe <- aggregate(MR, by=list(MR$q1), FUN=mean, na.rm=T)
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
  for(n in 1:inp){
    p <- ggplot2::ggplot(allMRs[[n]], aes(x=Lcor, y=Dcor))
    p <- densityplot(p)
    plist <- append(plist, list(p))
  }

#plotting! -> L and D ordered by cell length
  pL <- ggplot2::ggplot(MR, aes(x=num, y=Lmid))
  pL <- LWplot(pL, colc[[1]], max(MR$num, na.rm=T))
  pLpoint <- pL + ggplot2::geom_point() + ggplot2::ggtitle("Spot location on length axis ordered by cell length") + ggplot2::xlab("n\u209C\u2095 (ordered by cell length)") + ggplot2::ylab("Y-position (\u03BCm)") + ggplot2::theme_bw()
  pLD <- densityplot(pL) + ggplot2::ggtitle("Spot location on length axis ordered by cell length") + ggplot2::xlab("n\u209C\u2095 (ordered by cell length)") + ggplot2::ylab("Y-position (\u03BCm)") + ggplot2::geom_line(data=MR, ggplot2::aes(x=num,y=pole1),colour="white") + ggplot2::geom_line(data=MR, ggplot2::aes(x=num,y=pole2),colour="white")

  pW <- ggplot2::ggplot(MR, ggplot2::aes(x=num, y=Dum))
  pW <- LWplot(pW, colc[[1]], max(MR$num,na.rm=T))
  pWpoint <- pW + geom_point() + ggplot2::ggtitle("Spot location on width axis ordered by cell length") + ggplot2::xlab("(n\u209C\u2095 cell (ordered by cell length)") + ggplot2::ylab("X-position (\u03BCm)") + ggplot2::theme_bw()
  pWD <- densityplot(pW) + ggplot2::ggtitle("Spot location on width axis ordered by cell length") + ggplot2::xlab("n\u209C\u2095 (ordered by cell length)") + ggplot2::ylab("X-position (\u03BCm)") + ggplot2::geom_hline(yintercept=ymax) + ggplot2::geom_hline(yintercept=-ymax) + ggplot2::coord_cartesian(ylim=c(-ymax,ymax))

#make heatmap using the half max densities:

  ##for all quartiles (or n-tiles)
  allMRLmids <- split(MR[,c("Lmid","Dum")][!is.na(MR$Lmid),], MR$q1[!is.na(MR$Lmid)])
  #get maximum half max density of all groups and apply for each plot to have same scaling.
  mp<- max(sapply(allMRLmids, function(x) median(range(MASS::kde2d(x$Lmid, x$Dum)$z))))
  #create n-tile plots saved in a list.
  phlist <- list()
  for(n in 1:inp){
    p <- heatmap(plist[[n]], mp, colc, viridis)
    p <- coplot(p, xmax, ymax, max(allMRs[[n]]$max.length,na.rm=T)*0.5)
    phlist <- append(phlist, list(p))
  }

  u <- list()
  #do the same procedure for length/width heatmap plots
  mpL <- MASS::kde2d(MR$num[!is.na(MR$Lmid)], MR$Lmid[!is.na(MR$Lmid)])
  mpL1 <- median(range(mpL$z))
  pLD <- heatmap(pLD, mpL1, colc,viridis)
  pLD <- pLD + theme_minimal() + theme(panel.background=element_rect(fill=colc[[1]]), panel.grid = element_blank())
  u$lengthplot <- pLD

  mpW <- MASS::kde2d(MR$num[!is.na(MR$Dum)], MR$Dum[!is.na(MR$Dum)])
  mpW1 <- median(range(mpW$z))
  pWD <- heatmap(pWD, mpW1, colc,viridis)
  pWD <- pWD + theme_minimal() + theme(panel.background=element_rect(fill=colc[[1]]), panel.grid = element_blank())
  u$widthplot <- pWD

  if(missing(MESH)){
    u$qplots <- phlist
  }

  #################add mesh data###############################################################################################
  if(!missing(MESH)){
    if("X_rot"%in%colnames(MESH)!=T){
    print("Turning the cells... this may take a while.")
    MESH <- meshTurn(MESH, Xm, Ym)
    print("Finished turning the cells")
    }
    p2um <- as.numeric(magnificationList[mag])
    if("max_um"%in%colnames(MESH)!=T){
    MESH$max_um <- MESH$max.length*p2um
    MESH$maxwum <- MESH$max.width*p2um
    }


    MESH <- merge(MESH, MR[,c("cell", "frame", "q1")], all=T)
    MESHlist <- split(MESH, MESH$q1)
    print("Calculating mean cell outlines..")
    means <- lapply(MESHlist, function(x)superfun(x, 30, p2um))
    print("Finished calculating mean cell outlines")

    phmlist <- list()
    for(n in 1:inp){
      p <- phlist[[n]] + ggplot2::geom_point(data=means[[n]], ggplot2::aes(x=x,y=y), colour="white") + ggplot2::coord_fixed()
      phmlist <- append(phmlist, list(p))
    }
    if(AllPlot==F){
      u$qplots <- phmlist
    }
    #and create the plots:
    if(AllPlot==T){
      phmalist <- list()
      for(n in 1:inp){
      p1_all <- allplot(phmlist[[n]], allMRs[[n]], xmax, ymax, empty)
      phmalist <- append(phmalist,list(p1_all))
      u$qplots <- phmalist
      }
    }

    #in case you want it for the whole thing instead of only quartiles:

    pall <- ggplot2::ggplot(MR, aes(x=Lcor, y=Dcor))
    pall <- densityplot(pall)
    mppall <- mean(range(MASS::kde2d(MR$Lcor[!is.na(MR$Lcor)&!is.na(MR$Dcor)],MR$Dcor[!is.na(MR$Dcor)&!is.na(MR$Lcor)])$z))
    pall <- heatmap(pall, mppall, colc)
    pall <- coplot(pall,xmax, ymax, max(MR$max.length)*0.5)


    meantotal <- superfun(MESH, 30, p2um)
    pall <- pall + geom_point(data=meantotal, aes(x=x,y=y), colour="white")
    if(AllPlot==F){
      u$plottotal <- pall
    }
    if(AllPlot==T){
    pall_all <- allplot(pall, MR, xmax, ymax, empty)
    u$plottotal <- pall_all

    }


  }

  hislist <- list()
  #save all histograms (L coordinates) of the quartiles too:
  for(n in 1:inp){
    p1his <- ggplot2::ggplot(allMRs[[n]], aes(x=Lcor)) + geom_density(fill=colc[[1]][2], color=colc[[1]][2]) + theme_bw() + labs(x="Length(\u03BCm)") + coord_cartesian(xlim=c(-1.5,1.5))
    hislist <- append(hislist, list(p1his))
  }

  u$histograms <- hislist
  u$spotdata <- MR
  if(!missing(MESH)){
  	u$meshdata <- MESH
  }
  print("Done plotting.")
  return(u)
}

###########################################double histograms!! yay!!!######################################3333
##allplot function combines histograms and density plot.

allplot <- function(plot, data, xmax, ymax, empty){

  #prepare seperate plots: histograms (hL, hD) and modified coordinate plots(remove legend )
  p1D <- plot + theme(legend.position = "none")
  p1hL <- ggplot2::ggplot(data, aes(x=Lcor)) + geom_histogram(fill="black") + coord_cartesian(xlim = c(-xmax, xmax)) + theme_minimal() +theme(axis.title.x = element_blank())
  p1hD <- ggplot2::ggplot(data, aes(x=Dcor)) + geom_histogram(fill="black") + coord_flip(xlim = c(-ymax, ymax)) + theme_minimal() + theme(axis.title.y = element_blank())

  #align the plots properly before putting them together
  p1Dg <- ggplotGrob(p1D)
  p1hLg <- ggplotGrob(p1hL)
  p1hDg <- ggplotGrob(p1hD)

  maxWidth = grid::unit.pmax(p1Dg$widths[2:5], p1hLg$widths[2:5])
  p1Dg$widths[2:5] <- as.list(maxWidth)
  p1hLg$widths[2:5] <- as.list(maxWidth)

  #put the grids together using gridarrange
  return(arrangeGrob(p1hLg, empty, p1Dg, p1hD, ncol=2, nrow=2, widths=c(10*xmax, 2.5), heights=c(2, 10*ymax)))
}

##before using the function:
#create mockup plot to make space
empty <- ggplot2::ggplot(data.frame(u=1), ggplot2::aes(u,u)) +
  ggplot2::theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

colopts <- list(OrangeHot=list("#000000", "#D55E00", "#F0E442"), GreenYellow = list("#000000", "#009E73", "#F0E442"), ColdBlue = list("#000000", "#0072B2", "#FFFFFF"), YellowHot = list("#000000", "#e69f00", "#F0E442"), RedHot = list("#000000", "#FF0000", "#FFFF00"), WhiteOrange = list("#FFFFFF", "#F0E442", "#D55E00"))

singcollist <- list("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#D55E00", "#CC79A7")

xaxislist <- c("Cell Length (\u03BCm)", "Cell Width (\u03BCm)", "Cell Area (\u03BCm\u00B2)",
               "Spot location in length axis (\u03BCm)", "Spot location on width axis (\u03BCm)")



