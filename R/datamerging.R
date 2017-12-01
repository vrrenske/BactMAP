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

#plotfunction for a 2-color dotplot
cdot2 <- function(dataset, quartile, colorpalette = c("Green", "Red", "Blue", "Yellow", "Dark Blue"), namesetlist){
  colorpalette <- colorpalette[1:length(namesetlist)]
  names(colorpalette) <- unlist(namesetlist)
  plot <- ggplot2::ggplot(dataset[dataset$q1==quartile,], ggplot2::aes(x=Lcor, y= Dcor))
  return(plot + ggplot2::geom_point(aes(colour=color), size=4, alpha=0.3) + ggplot2::theme_minimal() + ggplot2::xlab("Length(\u00B5m)") + ggplot2::ylab("Width(\u00B5m)") + ggplot2::scale_color_manual(values=colorpalette))
}


###############################################################################################################################
#adding. here's a copy of the allplot function. I modified the histograms to become 2 color density plots instead, having the same color
#codes.

allplot <- function(plot, data, xmax, ymax, empty, xqmax){

  #prepare seperate plots: histograms (hL, hD) and modified coordinate plots(remove legend )
  p1D <- plot + ggplot2::theme(legend.position = "none") + ggplot2::coord_cartesian(xlim = c(-xmax, xmax), ylim = c(-ymax,ymax)) + ggplot2::geom_vline(xintercept=xqmax, alpha=0.4) + ggplot2::geom_vline(xintercept=-xqmax, alpha=0.4)
  p1hL <- ggplot2::ggplot(data, aes(x=Lcor, colour=color), alpha=0.6) + ggplot2::geom_density(aes(fill=color), alpha = 0.3) + ggplot2::coord_cartesian(xlim = c(-xmax, xmax)) + ggplot2::theme_minimal() +ggplot2::theme(axis.title.x = element_blank(), legend.position="none") + ggplot2::scale_color_manual(values=c(nameset1=colorpalette[1], nameset2=colorpalette[2])) + ggplot2::scale_fill_manual(values=c(nameset1=colorpalette[1], nameset2=colorpalette[2]))
  p1hD <- ggplot2::ggplot(data, aes(x=Dcor, colour=color), alpha=0.6) + ggplot2::geom_density(aes(fill=color), alpha = 0.3) + ggplot2::coord_flip(xlim = c(-ymax, ymax)) + ggplot2::theme_minimal() + ggplot2::theme(axis.title.y = element_blank(), legend.position="none") + ggplot2::scale_color_manual(values=c(nameset1=colorpalette[1], nameset2=colorpalette[2])) + ggplot2::scale_fill_manual(values=c(nameset1=colorpalette[1], nameset2=colorpalette[2]))

  #align the plots properly before putting them together
  p1Dg <- ggplotGrob(p1D)
  p1hLg <- ggplotGrob(p1hL)
  p1hDg <- ggplotGrob(p1hD)

  maxWidth = grid::unit.pmax(p1Dg$widths[2:5], p1hLg$widths[2:5])
  p1Dg$widths[2:5] <- as.list(maxWidth)
  p1hLg$widths[2:5] <- as.list(maxWidth)

  #put the grids together using gridarrange
  return(grid::arrangeGrob(p1hLg, empty, p1Dg, p1hD, ncol=2, nrow=2, widths=c(3, 1), heights=c(1, 2)))
}

#again you need the mockup plot in the corner:
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

#plotting

plotallsides <- function(dataset, quartile, colorpalette){
  p <- cdot2(dataset, quartile, colorpalette)
  pG <- allplot(p, dataset[dataset$q1 == quartile,], xmax, ymax, empty, max(dataset$length[dataset$q1==quartile], na.rm=T)*0.5)
  return(pG)
}



##################################################################################################################
#only histograms
hisfun <- function(dataset, quartile, colorpalette){
  return(ggplot2::ggplot(dataset[dataset$q1==quartile,], ggplot2::aes(x=Lcor, colour=color), alpha=0.6) + ggplot2::geom_density(aes(fill=color), alpha = 0.3) + ggplot2::coord_cartesian(xlim = c(-xmax, xmax)) + ggplot2::theme_minimal() +ggplot2::theme(legend.position="none") + ggplot2::scale_color_manual(values=c("GFP"=colorpalette[1], "RFP"=colorpalette[2])) + ggplot2::scale_fill_manual(values=c("GFP"=colorpalette[1], "RFP"=colorpalette[2]))  + ggplot2::xlab("Location length-axis (\u00B5m from mid-point)"))
}

##############################################################################################################


