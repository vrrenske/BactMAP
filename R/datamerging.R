##15-05-13##
##Renske van Raaphorst##

#########################merging GFP and RFP data############################################

#' @export
spotrMerge <- function(datasetlist, samecells=TRUE, namesetlist=list("GFP","RFP", "CFP", "YFP", "DAPI", "CY5"), groups = 4){
##we need to get both datasets in one. only, we need to keep them apart. that's why we need to add an extra column
#indicating which spot is gfp (or any name you want) and which is RFP (idem).
  numsets <- length(datasetlist)
  datasetlist <- lapply(c(1:numsets), function(x) nameSetList(datasetlist[[x]], namesetlist[[x]]))
  GR <- do.call(merge, c(datasetlist, all=T))

###################################################################################################################


##now you might want to go back to the prep and plotting script to make the quartile partitions.
##however, we have more cells now because of the GFP & RFP. we can circumvent that by cutting
##the four partitions by length:
  if(samecells==TRUE){
  GR <- GR[order(GR$max.length),]
  GR$num2 <- c(1:nrow(GR))
# by quartiles of the number of cells:
  GR$q1 <- cut(GR$num2, breaks=groups, labels = 1:groups)
  }
#when 2 different datasets are chosen, groups are grouped by their original quartiles (by length)
  if(samecells==FALSE){
    GR$q1 <- cut(GR$num, breaks=groups, labels = 1:groups)
  }


#prep for the code: define maxima.
  xmax <- 0.5*max(GR$max.length, na.rm=TRUE)
  ymax <- 0.5*max(GR$max.width, na.rm=TRUE)

####################################################################################################################

#Length and width corrected for average length per quartile for plotting the coordinate plots

  for(n in 1:groups){

    meansL <- mean(GR$max.length[GR$q1==n], na.rm=TRUE)
    meansW <- mean(GR$max.width[GR$q1==n], na.rm=TRUE)
    QR <- GR[GR$q1==n,]
    QR$Lcor <- QR$Lmid/QR$max.length*meansL
    QR$Dcor <- QR$Dum/QR$max.width*meansW
    if(n==1){
      QRall <- QR
    }
    if(n>1){
      QRall <- rbind(QR, QRall)
    }
}
QRall$q1 <- as.numeric(as.character(QRall$q1))
return(QRall)
}
#########################################################################################################################

#plotting
#' @export
plotMultiChannel <- function(dataset, cpalette=bm_Colors, colorpicks){
  if(missing(colorpicks)==TRUE){
    cdots <- lapply(1:length(unique(dataset$q1)), function(x) plotallsides(dataset[dataset$q1==x,], x, cpalette))
    histos <- hisfun(dataset, cpalette=cpalette)
  }
  if(missing(colorpicks)==FALSE){
    cdots <- lapply(1:length(unique(dataset$q1)), function(x) plotallsides(dataset[dataset$q1==x,], x, cpalette=cpalette, colorpicks=colorpicks))
    histos <- hisfun(dataset, cpalette=cpalette, colorpicks=colorpicks)
  }
  plotlist <- list(dotplots = cdots, histograms=histos)
  return(plotlist)
}

#plotfunction for a 2-color dotplot
cdot2 <- function(dataset, quartile, cpalette){
  plot <- ggplot2::ggplot(dataset[dataset$q1==quartile,], ggplot2::aes(x=Lcor, y= Dcor))
  return(plot + ggplot2::geom_point(aes(colour=color), size=4, alpha=0.3) + ggplot2::theme_minimal() + ggplot2::xlab("Length(\u00B5m)") + ggplot2::ylab("Width(\u00B5m)") + ggplot2::scale_color_manual(values=cpalette))
}


###############################################################################################################################
#adding. here's a copy of the allplot function. I modified the histograms to become 2 color density plots instead, having the same color
#codes.

allplotmulticolor <- function(plot, data, xmax, ymax, empty=empty2, xqmax, cpalette){

  #prepare seperate plots: histograms (hL, hD) and modified coordinate plots(remove legend )
  p1D <- plot + ggplot2::theme(legend.position = "none") + ggplot2::coord_cartesian(xlim = c(-xmax, xmax), ylim = c(-ymax,ymax)) + ggplot2::geom_vline(xintercept=xqmax, alpha=0.4) + ggplot2::geom_vline(xintercept=-xqmax, alpha=0.4)
  p1hL <- ggplot2::ggplot(data, aes(x=Lcor, colour=color), alpha=0.6) + ggplot2::geom_density(aes(fill=color), alpha = 0.3) + ggplot2::coord_cartesian(xlim = c(-xmax, xmax)) + ggplot2::theme_minimal() +ggplot2::theme(axis.title.x = element_blank(), legend.position="none") + ggplot2::scale_color_manual(values=cpalette) + ggplot2::scale_fill_manual(values=cpalette)
  p1hD <- ggplot2::ggplot(data, aes(x=Dcor, colour=color), alpha=0.6) + ggplot2::geom_density(aes(fill=color), alpha = 0.3) + ggplot2::coord_flip(xlim = c(-ymax, ymax)) + ggplot2::theme_minimal() + ggplot2::theme(axis.title.y = element_blank(), legend.position="none") + ggplot2::scale_color_manual(values=cpalette) + ggplot2::scale_fill_manual(values=cpalette)

  #align the plots properly before putting them together
  p1Dg <- ggplotGrob(p1D)
  p1hLg <- ggplotGrob(p1hL)
  p1hDg <- ggplotGrob(p1hD)

  maxWidth = grid::unit.pmax(p1Dg$widths[2:5], p1hLg$widths[2:5])
  p1Dg$widths[2:5] <- as.list(maxWidth)
  p1hLg$widths[2:5] <- as.list(maxWidth)

  #put the grids together using gridarrange
  return(gridExtra::grid.arrange(p1hLg, empty2, p1Dg,  p1hDg, ncol=2, nrow=2, widths=c(10*xmax, 2.5), heights=c(2, 10*ymax)))
}

#again you need the mockup plot in the corner:
empty2 <- ggplot2::ggplot(data.frame(u=1), ggplot2::aes(u,u)) +
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

#plotting

plotallsides <- function(dataset, quartile, cpalette = bm_Colors, colorpicks, empty = empty){
  xmax <- 0.5*max(dataset$max.length, na.rm=TRUE)
  ymax <- 0.5*max(dataset$max.width, na.rm=TRUE)

  if(missing(colorpicks)==TRUE){
    for(n in 1:length(unique(dataset$color))){
      names(cpalette)[n] = unique(dataset$color)[n]
    }
    cpalette <- cpalette[1:length(unique(dataset$color))]
  }
  if(missing(colorpicks)!=TRUE){
    for(n in 1:length(colorpicks)){
      if(colorpicks[n]%in%names(cpalette)==TRUE){
        colorpicks[n] <- as.character(cpalette[colorpicks[n]])
      }
    }
  }
  p <- cdot2(dataset, quartile, cpalette)
  pG <- allplotmulticolor(p, dataset[dataset$q1 == quartile,], xmax, ymax, empty, max(dataset$max.length[dataset$q1==quartile], na.rm=T)*0.5, cpalette)
  return(pG)
}




##################################################################################################################
#only histograms
hisfun <- function(dataset, colorpicks, cpalette = bm_Colors){
  xmax <- 0.5*max(dataset$max.length, na.rm=TRUE)
  if(missing(colorpicks)==TRUE){
    for(n in 1:length(unique(dataset$color))){
      names(cpalette)[n] = unique(dataset$color)[n]
    }
    cpalette <- cpalette[1:length(unique(dataset$color))]
  }
  if(missing(colorpicks)!=TRUE){
    for(n in 1:length(colorpicks)){
      if(colorpicks[n]%in%names(cpalette)==TRUE){
        colorpicks[n] <- as.character(cpalette[colorpicks[n]])
      }
    }
  }
  return(ggplot2::ggplot(dataset, ggplot2::aes(x=Lcor, colour=color), alpha=0.6) + ggplot2::geom_density(aes(fill=color), alpha = 0.3) + ggplot2::coord_cartesian(xlim = c(-xmax, xmax)) + ggplot2::theme_minimal() + ggplot2::theme(legend.position="none") + ggplot2::scale_color_manual(values=cpalette) + ggplot2::scale_fill_manual(values=cpalette)  + ggplot2::xlab("Location length-axis (\u00B5m from mid-point)") + facet_grid(q1~.))
}

##############################################################################################################

##Single Colors Standardly Loaded Into BactMap - from  http://jfly.iam.u-tokyo.ac.jp/color/:

bm_Colors <-  c("bm_BlueGreen" = "#009E73", "bm_Orange" = "#E69F00", "bm_SkyBlue" = "#56B4E9",
                "bm_Yellow" = "#F0E442", "bm_Blue" = "#0072B2", "bm_Vermillion" = "#D55E00",
                "bm_RedPurple" = "#CC79A7", "bm_Grey" = "#999999")


