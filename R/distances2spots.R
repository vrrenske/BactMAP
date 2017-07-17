#COLOCALIZATION

#23-7-2015

#Renske van Raaphorst

#pixel2um = conversion factor for pixel to um. (Deltavision 100x = 0.064)
#Plot single cell with detected spots:
#data = MESH (mesh.Rda) and GR (merged dataframe resulting from data merging RFP/GFP)

spotrSinglecellplot <- function(outlinedata, spotdata,framenum, cellnum, val1 = "GFP", val2 = "RFP"){
  oo <- outlinedata[outlinedata$slice==framenum&outlinedata$cell==cellnum,]
  so <- spotdata[spotdata$frame==framenum&spotdata$cell==cellnum,]
  totalframe <- data.frame(x=0, y=0, Lmid=0, Dum=0, color = "0", num =0)
  for(n in 1:nrow(so)){
    dat <- circleFun(c(so$Lmid[n],so$Dum[n]),so$fwhm_x[n]*pixel2um, so$fwhm_y[n]*pixel2um,npoints = 100, so$color[n])
    totalframe <- merge(totalframe, dat, all=T)
  }
  totalframe <- totalframe[totalframe$color!="0",]
  totalframe <- totalframe[order(totalframe$Lmid, totalframe$num),]
  return(ggplot(oo, aes(x=y0_rot*-pixel2um, y=x0_rot*-pixel2um)) +geom_point() + geom_point(aes(x=y1_rot*-pixel2um, y=x1_rot*-0.064)) + geom_point(data=so, aes(x=Lmid, y=Dum, color=color)) + geom_polygon(data=totalframe, aes(x=x, y=y, fill= color, linetype=as.factor(Lmid)), alpha=0.4) + theme_minimal() + scale_colour_manual(
    values = c(val1 = "green",val2 = "red")) + scale_fill_manual(values=c(val1="green", val2="red")) + xlab("length (um)") + ylab("width (um)") + ggtitle(paste("frame", framenum, ", cell", cellnum)))
}

spotrCircleFun <- function(center = c(0,0),diameterx = 1, diametery = 1, npoints = 100, color1){
  rx = diameterx / 2
  ry = diametery /2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + rx * cos(tt)
  yy <- center[2] + ry * sin(tt)
  return(data.frame(x = xx, y = yy, Lmid = diameterx, Dum = diametery, color = color1, num = 1:100))
}

##################measuring distances between spots############################

#distinguishing between GFP & RFP spots in 1 frame

rncols <- function(dat, name){
  dat <- dat[,c("frame", "cell", "q1", "fwhm_x","fwhm_y", "Lmid", "Dum")]
  colnames(dat)[4:7] <- paste0(colnames(dat)[4:7], "_",name)
  return(dat)
}

#merge RFP/GFP

mergeandmeasure <- function(dat1, dat2){
  comdat <- merge(dat1, dat2, all = T)
  comdat <- comdat[order(comdat$frame, comdat$cell),]
  comdat$xdist <- comdat$Lmid_GFP - comdat$Lmid_RFP
  comdat$ydist <- comdat$Dum_GFP - comdat$Dum_RFP
  comdat$sdist <- sqrt(comdat$xdist^2 + comdat$ydist^2)
  comdat <- comdat[!is.na(comdat$sdist),]
  comdat <- test(comdat)
  comdat <- groupdists(comdat)
  return(comdat)
}

test <- function(comdat){
  comdat$x <- 0
    for(n in 2:(nrow(comdat))){
      if(comdat$Lmid_GFP[n]==comdat$Lmid_GFP[(n-1)]|comdat$Lmid_RFP[n]==comdat$Lmid_RFP[(n-1)]){
        print(comdat$cell[n])
        if (comdat$sdist[n] > comdat$sdist[(n-1)]){
          print("BOOH")
          comdat$x[n] <- "nope"}
        else
          comdat$x[n-1] <- "nope"
    }
  }
  comdat <- comdat[comdat$x!="nope",]
  comdat$x <- NULL
  return(comdat)
}

groupdists <- function(comdat){
  comdat$group <- "no overlap"
  comdat$group[comdat$sdist < (comdat$fwhm_x_GFP + comdat$fwhm_x_RFP)*pixel2um] <- "< fwhm1+fwhm2"
  comdat$group[comdat$sdist < comdat$fwhm_x_GFP*pixel2um|comdat$sdist < comdat$fwhm_x_RFP*pixel2um] <- "< fwhm"
  return(comdat)
}

spotrPlotdists <- function(distancedata, dataname){
  comdat$u <- dataname
  return(ggplot(comdat, aes(x=u, fill=group)) + geom_bar(aes(y=(..count..)/sum(..count..)*100), width=0.3) + ylab("percentage of spots") + theme_minimal() + scale_fill_colorblind() + coord_flip())
}
  
spotrMeasuredistances <- function(dat1, dat2, val1="GFP", val2="RFP"){
  G <- rncols(dat1, val1)
  R <- rncols(dat2, val2)
  combi <- mergeandmeasure(G, R)
  return(combi)
}
  
 
spotrGetdistancefreqs <- function(distancedata){
  freqs <- data.frame(table(distancedata$group))
  freqs$sum <- sum(freqs$Freq)
  freqs$percentages <- freqs$Freq/freqs$sum*100
  return(freqs)
}



#standard dotplot function - should also work for other variables.
spotrDotplot <- function(data1, xaxis=q1, yaxis = sdist, bins = 0.012, xlabel = "", ylabel =""){
dotplot1 <- ggplot(data1, aes(x=xaxis, y=yaxis)) + geom_dotplot(aes(fill=xaxis, color=xaxis), binaxis="y", stackdir="center", binwidth=bins, xlab(xlabel), ylab(ylabel))
return(dotplot1)
}


#plot spots in one cell over time using colour code:

timecellplot <- function(dat, num, Mdat){
  return(ggplot(dat[dat$cell==num,], aes(x=Lmid, y=Dum, color=frame)) + geom_point() + scale_colour_gradient2(low = "#CC79A7", mid= "#F0E442", high = "#009E73", midpoint = 30, space = "Lab", guide = "colourbar") + theme_minimal() 
         + geom_point(data=Mdat[Mdat$cell==num&Mdat$slice==1,], aes(x=y0_rot*-pixel2um, y=x0_rot*-pixel2um), colour="black") + geom_point(data=Mdat[Mdat$cell==num&Mdat$slice==1,], aes(x=y1_rot*-pixel2um, y=x1_rot*-pixel2um),colour="black"))
}

timecelllineplot <- function(dat, num, Mdat){
  return(ggplot(dat[dat$cell==num,], aes(x=Lmid, y=Dum, color=frame)) + geom_path(aes(line=as.factor(spot))) + scale_colour_gradient2(low = "#CC79A7", mid= "#F0E442", high = "#009E73", midpoint = 30, space = "Lab", guide = "colourbar") + theme_minimal() 
         + geom_point(data=Mdat[Mdat$cell==num&Mdat$slice==1,], aes(x=y0_rot*-pixel2um, y=x0_rot*-pixel2um),colour="black") +geom_point() + geom_point(data=Mdat[Mdat$cell==num&Mdat$slice==1,], aes(x=y1_rot*-pixel2um, y=x1_rot*-pixel2um),colour="black"))
}



