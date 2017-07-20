##15-05-13##
##Renske van Raaphorst##

#########################merging GFP and RFP data############################################

spotrMerge <- function(dataset1, dataset2, samecells=TRUE, nameset1="GFP", nameset2="RFP", groups = 4){
##we need to get both datasets in one. only, we need to keep them apart. that's why we need to add an extra column
#indicating which spot is gfp (or any name you want) and which is RFP (idem).
  dataset1$color <- nameset1
  dataset2$color <- nameset2
  GR <- merge(dataset1, dataset2,all=T)

###################################################################################################################


##now you might want to go back to the prep and plotting script to make the quartile partitions.
##however, we have more cells now because of the GFP & RFP. we can circumvent that by cutting
##the four partitions by length:
  if(samecells==TRUE){
  GR <- GR[order(GR$length),]
  GR$num2 <- c(1:nrow(GR))
# by quartiles of the number of cells:
  GR$q1 <- cut(GR$num2, breaks=groups, labels = 1:groups)
  }
#when 2 different datasets are chosen, groups are grouped by their original quartiles (by length)
  if(samecells==FALSE){
    GR$q1 <- cut(GR$num, breaks=groups, labels = 1:groups)
  }


#prep for the code: define maxima.
  xmax <- 0.5*max(GR$length, na.rm=TRUE)
  ymax <- 0.5*max(GR$max.width, na.rm=TRUE)

####################################################################################################################

#Length and width corrected for average length per quartile for plotting the coordinate plots

  for(n in 1:groups){

    meansL <- mean(GR$length[GR$q1==n], na.rm=TRUE)
    meansW <- mean(GR$max.width[GR$q1==n], na.rm=TRUE)
    QR <- GR[GR$q1==n,]
    QR$Lcor <- QR$Lmid/QR$length*meansL
    QR$Dcor <- QR$Dum/QR$max.width*meansW
    if(n==1){
      QRall <- QR
    }
    if(n>1){
      QRall <- rbind(QR, QRall)
    }
}

return(QRall)
}
#########################################################################################################################

#plotting

#plotfunction for a 2-color dotplot
cdot2 <- function(dataset, quartile, colorpalette){
  plot <- ggplot(dataset[dataset$q1==quartile,], aes(x=Lcor, y= Dcor))
  return(plot + geom_point(aes(colour=color), size=4, alpha=0.3) + theme_minimal() + xlab("Length(\u00B5m)") + ylab("Width(\u00B5m)") + scale_color_manual(values=c(nameset1=colorpalette[1], nameset2=colorpalette[2])))
}


###############################################################################################################################
#adding. here's a copy of the allplot function. I modified the histograms to become 2 color density plots instead, having the same color
#codes.

allplot <- function(plot, data, xmax, ymax, empty, xqmax){

  #prepare seperate plots: histograms (hL, hD) and modified coordinate plots(remove legend )
  p1D <- plot + theme(legend.position = "none") + coord_cartesian(xlim = c(-xmax, xmax), ylim = c(-ymax,ymax)) + geom_vline(xintercept=xqmax, alpha=0.4) + geom_vline(xintercept=-xqmax, alpha=0.4)
  p1hL <- ggplot(data, aes(x=Lcor, colour=color), alpha=0.6) + geom_density(aes(fill=color), alpha = 0.3) + coord_cartesian(xlim = c(-xmax, xmax)) + theme_minimal() +theme(axis.title.x = element_blank(), legend.position="none") + scale_color_manual(values=c(nameset1=colorpalette[1], nameset2=colorpalette[2])) + scale_fill_manual(values=c(nameset1=colorpalette[1], nameset2=colorpalette[2]))
  p1hD <- ggplot(data, aes(x=Dcor, colour=color), alpha=0.6) + geom_density(aes(fill=color), alpha = 0.3) + coord_flip(xlim = c(-ymax, ymax)) + theme_minimal() + theme(axis.title.y = element_blank(), legend.position="none") + scale_color_manual(values=c(nameset1=colorpalette[1], nameset2=colorpalette[2])) + scale_fill_manual(values=c(nameset1=colorpalette[1], nameset2=colorpalette[2]))

  #align the plots properly before putting them together
  p1Dg <- ggplotGrob(p1D)
  p1hLg <- ggplotGrob(p1hL)
  p1hDg <- ggplotGrob(p1hD)

  maxWidth = grid::unit.pmax(p1Dg$widths[2:5], p1hLg$widths[2:5])
  p1Dg$widths[2:5] <- as.list(maxWidth)
  p1hLg$widths[2:5] <- as.list(maxWidth)

  #put the grids together using gridarrange
  return(arrangeGrob(p1hLg, empty, p1Dg, p1hD, ncol=2, nrow=2, widths=c(3, 1), heights=c(1, 2)))
}

#again you need the mockup plot in the corner:
empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(
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

#plotting

plotallsides <- function(dataset, quartile, colorpalette){
  p <- cdot2(dataset, quartile, colorpalette)
  pG <- allplot(p, dataset[dataset$q1 == quartile,], xmax, ymax, empty, max(dataset$length[dataset$q1==quartile], na.rm=T)*0.5)
  return(pG)
}



##################################################################################################################
#only histograms
hisfun <- function(dataset, quartile, colorpalette){
  return(ggplot(dataset[dataset$q1==quartile,], aes(x=Lcor, colour=color), alpha=0.6) + geom_density(aes(fill=color), alpha = 0.3) + coord_cartesian(xlim = c(-xmax, xmax)) + theme_minimal() +theme(legend.position="none") + scale_color_manual(values=c("GFP"=colorpalette[1], "RFP"=colorpalette[2])) + scale_fill_manual(values=c("GFP"=colorpalette[1], "RFP"=colorpalette[2]))  + xlab("Location length-axis (\u00B5m from mid-point)"))
}

##############################################################################################################


