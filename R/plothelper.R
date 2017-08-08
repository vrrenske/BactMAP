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

colorchoiceplot <- function(colchoice, nums){
  namelist <- c("Orange Hot", "Green Hot", "Cold Blue", "Yellow Hot", "Red Hot", "White Orange")
  z <- rmvnorm(100, mean=c(3,5), sigma=matrix(c(1,0.5,0.5,2), nrow=2))
  z <- data.frame(z)
  ggplot2::ggplot(z, aes(x=X1, y=X2)) + stat_density2d(aes(fill=..density..), geom="raster", contour=FALSE) + scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint=0.06) + theme_minimal() + theme(legend.position="none") + ggtitle(namelist[nums]) + xlab("") + ylab("")
}

#wat basisplotfuncties
densityplot <- function(plot){
  return(plot + ggplot2::stat_density2d(aes(fill=..density..), geom="raster", contour = FALSE))
}

LWplot <- function(plot, u="black", maxn){
  return(plot + ggplot2::geom_rect(xmin=0, xmax=maxn, ymin=-1.5, ymax=1.5, fill=u) + theme_minimal()  + scale_x_continuous(limits=c(0,maxn)))
}

heatmap <- function(pdens, mp, colchoice){
  return(pdens + ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint =mp, space = "Lab") + theme_minimal() + theme(legend.position="none"))
}

#function for goodlooking x/y coordinate plot:
#makes a plot sized as the max cell (width/length: xmax) of the quartile inside a plot
#which has a grey background as large as the largest cell in the dataset
#so all quartile plots will have the same dimensions.
#the title, y axis and x axis will also be drawn.
coplot <- function(pheat, xmax, ymax, xqmax){
 return(pheat + ggplot2::xlab("Length (\u00B5m)") + ylab("Width (\u00B5m)") + coord_cartesian(xlim = c(-xmax,xmax), ylim=c(-ymax,ymax)) + geom_vline(xintercept = xqmax) + geom_vline(xintercept=-xqmax) + geom_hline(yintercept = ymax) + geom_hline(yintercept = -ymax) + theme(panel.background = element_rect(fill = "dark grey"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
}


mergeframes <- function(REP, M){

  REP<- REP[(0<REP$l),]
  if("rel.l" %in% colnames(REP)){
    REP <-REP[(REP$rel.l<1),]
  }

  M <- M[order(M$length),]
  M$cellnum <- c(1:nrow(M))

  #merging
  MR <- merge(M,REP, all=T)

  #remove MR's cells which have NA's in the cell area.
  MR <- MR[!is.na(MR$length),]
  #remove duplicated rows - something with the van Oijen plugin :S
  MR <- MR[!duplicated(MR$l)&!duplicated(MR$d),]

  MR <- spotMR(MR)

  #if needed: remove smallest and largest ones (cutoff: smaller than 1/2 mean and larger than 2x mean)
  MR <- MR[MR$length<(2*mean(MR$length)),]
  MR <- MR[MR$length>(0.5*mean(MR$length)),]

  MR <- MR[order(MR$length),]

  #make column with row numbers per cell length.
  MR$num <- c(1:nrow(MR))

  MR <- LimDum(MR)
  }

#add spotnumbers!
spotMR <- function(dat){
if("spot" %in% colnames(dat)){
  #NA in spots replaced by "0"
  dat$spot[is.na(dat$spot)] <- 0
} else {
  dat <-dat[order(dat$frame, dat$cell,dat$length),]
  dat$spot <- 1
  for(n in 1:(nrow(dat)-1)){
    if(dat$length[n+1]==dat$length[n]){
      dat$spot[n+1] <- dat$spot[n] + 1
    } }
  dat$spot[is.na(dat$spot)] <- 0
}
return(dat)
}



################################################################################################
#plot preparation
#quartiles, maxima, etc.
LimDum <- function(MR){
  MR$Lmid<-(MR$l-0.5*MR$length)
  MR$pole1<- -MR$length*0.5
  MR$pole2<- -MR$pole1
  MR$Dum <- MR$d
  MR <- MR[abs(MR$Lmid)<MR$pole2,]
  MR <- MR[abs(MR$Dum)<MR$max.width/2,]
  return(MR)
}
##make quartile partitions
#measure one, by length:
#MR$q1 <- cut(MR$length, breaks=4, labels = 1:4)
############TURNING THE MESH TO THE RIGHT ANGLE#########################################################
meshturn <- function(MESH){
  MESH <- MESH[!is.na(MESH$x0),]
  if("frame" %in% colnames(MESH)){ MESH$slice <- MESH$frame}
  MESH <- MESH[order(MESH$slice, MESH$cell, MESH$length),]
  MESH$Xmid <- 0
  MESH$Ymid <- 0
  MESH <- MESH[MESH$x0 >1,]

  for(n in 1:max(MESH$slice)){
    for(z in as.integer(levels(as.factor(MESH$cell[MESH$slice==n])))){
      X_1 <- MESH$x0[MESH$slice==n&MESH$cell==z][1]
      Y_1 <- MESH$y0[MESH$slice==n&MESH$cell==z][1]
      X <- MESH$x0[MESH$slice==n&MESH$cell==z]
      X_2 <- tail(X[!is.na(X)], n=1)
      Y <- MESH$y0[MESH$slice==n&MESH$cell==z]
      Y_2 <- tail(Y[!is.na(Y)], n=1)
      Xmid <- X_1 + 0.5*(X_2-X_1)
      Ymid <- Y_1 + 0.5*(Y_2-Y_1)
      MESH$Xmid[MESH$slice==n&MESH$cell==z] <- Xmid
      MESH$Ymid[MESH$slice==n&MESH$cell==z] <- Ymid

    }
  }

  ######################### mid-point correction ######################################################################


  MESH$x0_cor <- MESH$x0 - MESH$Xmid
  MESH$y0_cor <- MESH$y0 - MESH$Ymid
  MESH$x1_cor <- MESH$x1 - MESH$Xmid
  MESH$y1_cor <- MESH$y1 - MESH$Ymid

  ########################## rotation #################################################################################

  #adding a column with the rotation angle for each cell
  MESH$angle <- 0
  for(n in 1:max(MESH$slice)){
    for(z in as.integer(levels(as.factor(MESH$cell[MESH$slice==n])))){
      Xmin_c <- MESH$x0_cor[MESH$slice==n&MESH$cell==z][1]
      Ymin_c <- MESH$y0_cor[MESH$slice==n&MESH$cell==z][1]
      MESH$angle[MESH$cell==z&MESH$slice==n] <- atan(Xmin_c/Ymin_c)
    }
  }

  #rotating the coordinates
  MESH$x0_rot <- MESH$x0_cor * cos(MESH$angle) - MESH$y0_cor * sin(MESH$angle)
  MESH$x1_rot <- MESH$x1_cor * cos(MESH$angle) - MESH$y1_cor * sin(MESH$angle)
  MESH$y0_rot <- MESH$x0_cor * sin(MESH$angle) + MESH$y0_cor * cos(MESH$angle)
  MESH$y1_rot <- MESH$x1_cor * sin(MESH$angle) + MESH$y1_cor * cos(MESH$angle)
  MESH <- MESH[!is.na(MESH$y1_rot),]
  #average of all cells
  MESH$x0_rot <- -1*abs(MESH$x0_rot)
  MESH$x1_rot <- abs(MESH$x1_rot)


  return(MESH)
}

######################### mid-points ################################################################################
mfun <- function(points, b, binlist){
  means <- c()
  for(q in 1:b){
    mq <- mean(points[binlist==q])
    means[q] <- mq
  }
  return(means)
}

superfun <- function(dat, bins,mag){
  dat$av <- 0
  dat$av <- dat$y0_rot/dat$max_length*100

  cutpoints<-quantile(dat$av,(0:bins)/bins)
  dat$binned <-cut(dat$av,cutpoints, include.lowest=TRUE, labels = 1:bins)

  x0means <- mfun(dat$x0_rot, bins, dat$binned)
  x1means <- mfun(dat$x1_rot, bins, dat$binned)
  y0means <- mfun(dat$y0_rot, bins, dat$binned)
  y1means <- mfun(dat$y1_rot, bins, dat$binned)
  maxLmeans <- mfun(dat$max_length, bins, dat$binned)

  meanframe <- data.frame(x0means, x1means, y0means, y1means, maxLmeans)
  meanframe <- meanframe * mag
  colnames(meanframe) <- c("x0", "x1", "y0", "y1", "max_length")
  meanframe <- meanframe[abs(meanframe$y0)<=0.5*meanframe$max_length&abs(meanframe$y1)<=meanframe$max_length,]
  return(meanframe)
}

#or two, by quartiles of the number of cells:
plotsout <- function(MR, inp, MESH, colc, mag){

  MR$q1 <- cut(MR$cellnum, breaks=inp, labels = 1:inp)

  xmax <- 0.5*max(MR$length, na.rm=TRUE)
  ymax <- 0.5*max(MR$max.width, na.rm=TRUE)


  #Length and width corrected for average length per quartile for plotting the coordinate plots
  meansL <- c(mean(MR$length[MR$q1==1], na.rm=TRUE),mean(MR$length[MR$q1==2], na.rm=TRUE), mean(MR$length[MR$q1==3], na.rm=TRUE))
  if(inp>3) meansL[4] <- mean(MR$length[MR$q1==4], na.rm=TRUE)
  if(inp>4) meansL[5] <- mean(MR$length[MR$q1==5], na.rm=TRUE)
  if(inp>5) meansL[6] <- mean(MR$length[MR$q1==6], na.rm=TRUE)
  meansW <- c(mean(MR$max.width[MR$q1==1], na.rm=TRUE), mean(MR$max.width[MR$q1==2], na.rm=TRUE), mean(MR$max.width[MR$q1==3], na.rm=TRUE))
  if(inp>3) meansW[4] <- mean(MR$max.width[MR$q1==4], na.rm=TRUE)
  if(inp>4) meansW[5] <- mean(MR$max.width[MR$q1==5], na.rm=TRUE)
  if(inp>5) meansW[6] <- mean(MR$max.width[MR$q1==6], na.rm=TRUE)

  #length
  MR$Lcor <- (MR$Lmid/MR$length*meansL[1])
  #width
  MR$Dcor <- (MR$Dum/MR$max.width*meansW[1])

  #seperate frames:
  Q1 <- MR[MR$q1==1,]
  Q2 <- MR[MR$q1==2,]
  Q3 <- MR[MR$q1==3,]
  if(inp>3){
  Q4 <- MR[MR$q1==4,]
  }
  if(inp>4){
    Q5 <- MR[MR$q1==5,]
  }
  if(inp>5){
    Q6 <- MR[MR$q1==6,]
  }

  #length
  Q2$Lcor <- (Q2$Lmid/Q2$length*meansL[2])
  Q3$Lcor <- (Q3$Lmid/Q3$length*meansL[3])
  if(inp>3) Q4$Lcor <- (Q4$Lmid/Q4$length*meansL[4])
  if(inp>4) Q5$Lcor <- (Q5$Lmid/Q5$length*meansL[5])
  if(inp>5) Q6$Lcor <- (Q6$Lmid/Q6$length*meansL[6])

  #width
  Q2$Dcor <- (Q2$Dum/Q2$max.width*meansW[2])
  Q3$Dcor <- (Q3$Dum/Q3$max.width*meansW[3])
  if(inp>3) Q4$Dcor <- (Q4$Dum/Q4$max.width*meansW[4])
  if(inp>4) Q5$Dcor <- (Q5$Dum/Q5$max.width*meansW[5])
  if(inp>5) Q6$Dcor <- (Q6$Dum/Q6$max.width*meansW[6])

###############################################################################################################################################
  ##plotting! -> coordinate plots
  p1 <- ggplot2::ggplot(Q1, aes(x=Lcor, y=Dcor))
  p1 <- densityplot(p1)
  p1 <- coplot(p1, xmax, ymax, max(Q1$length,na.rm=T)*0.5)


  p2 <- ggplot2::ggplot(Q2, aes(x=Lcor, y=Dcor))
  p2 <- densityplot(p2)
  p2 <- coplot(p2, xmax, ymax, max(Q2$length,na.rm=T)*0.5)

  p3 <- ggplot2::ggplot(Q3, aes(x=Lcor, y=Dcor))
  p3 <- densityplot(p3)
  p3 <- coplot(p3, xmax, ymax, max(Q3$length,na.rm=T)*0.5)

  if(inp>3){
  p4 <- ggplot2::ggplot(Q4, aes(x=Lcor, y=Dcor))
  p4 <- densityplot(p4)
  p4 <- coplot(p4, xmax, ymax, max(Q4$length, na.rm=TRUE)*0.5)
  }
  if(inp>4){
    p5 <- ggplot2::ggplot(Q5, aes(x=Lcor, y=Dcor))
    p5 <- densityplot(p5)
    p5 <- coplot(p4, xmax, ymax, max(Q5$length, na.rm=TRUE)*0.5)
  }
  if(inp>5){
    p6 <- ggplot2::ggplot(Q6, aes(x=Lcor, y=Dcor))
    p6 <- densityplot(p6)
    p6 <- coplot(p6, xmax, ymax, max(Q6$length, na.rm=TRUE)*0.5)
  }

#plotting! -> L and D ordered by cell length
  pL <- ggplot2::ggplot(MR, aes(x=num, y=Lmid))
  pL <- LWplot(pL, "black", max(MR$num, na.rm=T))
  pLpoint <- pL + geom_point() + ggtitle("Spot location on length axis ordered by cell length") + xlab("nth cell (ordered by cell length)") + ylab("Y-position (\u03BCm)") + theme_bw()
  pLD <- densityplot(pL) + ggtitle("Spot location on length axis ordered by cell length") + xlab("nth cell (ordered by cell length)") + ylab("Y-position (\u03BCm)") + geom_line(data=MR, aes(x=num,y=pole1),colour="white") + geom_line(data=MR, aes(x=num,y=pole2),colour="white")

  pW <- ggplot2::ggplot(MR, aes(x=num, y=Dum))
  pW <- LWplot(pW, "black", max(MR$num,na.rm=T))
  pWpoint <- pW + geom_point() + ggtitle("Spot location on width axis ordered by cell length") + xlab("nth cell (ordered by cell length)") + ylab("X-position (\u03BCm)") + theme_bw()
  pWD <- densityplot(pW) + ggtitle("Spot location on width axis ordered by cell length") + xlab("nth cell (ordered by cell length)") + ylab("X-position (\u03BCm)") + geom_hline(yintercept=ymax) + geom_hline(yintercept=-ymax) + coord_cartesian(ylim=c(-ymax,ymax))

#make heatmap using the half max densities:

  mp1 <- kde2d(Q1$Lmid[!is.na(Q1$Lmid)], Q1$Dum[!is.na(Q1$Dum)])
  mp2 <- kde2d(Q2$Lmid[!is.na(Q2$Lmid)], Q2$Dum[!is.na(Q2$Dum)])
  mp3 <- kde2d(Q3$Lmid[!is.na(Q3$Lmid)], Q3$Dum[!is.na(Q3$Dum)])
  if(inp>3) mp4 <- kde2d(Q4$Lmid[!is.na(Q4$Lmid)], Q4$Dum[!is.na(Q4$Dum)])
  if(inp>4) mp5 <- kde2d(Q5$Lmid[!is.na(Q5$Lmid)], Q5$Dum[!is.na(Q5$Dum)])
  if(inp>5) mp6 <- kde2d(Q6$Lmid[!is.na(Q6$Lmid)], Q6$Dum[!is.na(Q6$Dum)])

  mplist <- c(median(range(mp1$z)), median(range(mp2$z)), median(range(mp3$z)))
  if(inp>3)mplist[4] <- median(range(mp4$z))
  if(inp>4)mplist[5] <- median(range(mp5$z))
  if(inp>5)mplist[6] <- median(range(mp6$z))

  mp <- max(mplist)
  p1 <- heatmap(p1, mp, colc)
  p2 <- heatmap(p2, mp, colc)
  p3 <- heatmap(p3, mp, colc)
  if(inp>3)p4 <- heatmap(p4, mp, colc)
  if(inp>4)p5 <- heatmap(p5, mp, colc)
  if(inp>5)p6 <- heatmap(p6, mp, colc)


  mpL <- kde2d(MR$num[!is.na(MR$Lmid)], MR$Lmid[!is.na(MR$Lmid)])
  mpL1 <- median(range(mpL$z))

  pLD <- heatmap(pLD, mpL1, colc)

  mpW <- kde2d(MR$num[!is.na(MR$Dum)], MR$Dum[!is.na(MR$Dum)])
  mpW1 <- median(range(mpW$z))

  pWD <- heatmap(pWD, mpW1, colc)

  MESH <- meshturn(MESH)

  MESH$max_um <- MESH$max_length*mag
  MESH$maxwum <- MESH$max.width*mag



  MESH$q <- "0"
  MESH$q[MESH$max_um<=max(MR$length[MR$q1==1], na.rm=T)] <- "Q1"
  MESH$q[MESH$max_um<=max(MR$length[MR$q1==2], na.rm=T)&MESH$max_um>max(MR$length[(MR$q1)==1],na.rm=T)] <- "Q2"
  MESH$q[MESH$max_um<=max(MR$length[(MR$q1)==3],na.rm=T)&MESH$max_um>max(MR$length[(MR$q1)==2],na.rm=T)] <- "Q3"
  if(inp>3)MESH$q[MESH$max_um<=max(MR$length[(MR$q1)==4],na.rm=T)&MESH$max_um>max(MR$length[(MR$q1)==3],na.rm=T)] <- "Q4"
  if(inp>4)MESH$q[MESH$max_um<=max(MR$length[(MR$q1)==5],na.rm=T)&MESH$max_um>max(MR$length[(MR$q1)==4],na.rm=T)] <- "Q5"
  if(inp>5)MESH$q[MESH$max_um<=max(MR$length[(MR$q1)==6],na.rm=T)&MESH$max_um>max(MR$length[(MR$q1)==5],na.rm=T)] <- "Q6"

  meanq1 <- superfun(MESH[MESH$q=="Q1",], 30,mag)
  meanq2 <- superfun(MESH[MESH$q=="Q2",], 30, mag)
  meanq3 <- superfun(MESH[MESH$q=="Q3",], 30, mag)
  if(inp>3)meanq4 <- superfun(MESH[MESH$q=="Q4",], 30, mag)
  if(inp>4)meanq5 <- superfun(MESH[MESH$q=="Q5",], 30, mag)
  if(inp>5)meanq6 <- superfun(MESH[MESH$q=="Q6",], 30, mag)
  meantotal <- superfun(MESH, 30, mag)

  #made meanq's by running script above where MESH is replaced by MESH[MESH$q=="Q1",] etc and meanq <- meanframe
  p1 <- p1 + ggplot2::geom_point(data=meanq1, aes(x=y0,y=x0), colour="white") + geom_point(data=meanq1, aes(x=y1,y=x1), colour="white") + coord_fixed()
  p2 <- p2 + ggplot2::geom_point(data=meanq2, aes(x=y0,y=x0), colour="white") + geom_point(data=meanq2, aes(x=y1,y=x1), colour="white") + coord_fixed()
  p3 <- p3 + ggplot2::geom_point(data=meanq3, aes(x=y0,y=x0), colour="white") + geom_point(data=meanq3, aes(x=y1,y=x1), colour="white") + coord_fixed()
  if(inp>3)p4 <- p4 + ggplot2::geom_point(data=meanq4, aes(x=y0,y=x0), colour="white") + geom_point(data=meanq4, aes(x=y1,y=x1), colour="white") + coord_fixed()
  if(inp>4)p5 <- p5 + ggplot2::geom_point(data=meanq5, aes(x=y0,y=x0), colour="white") + geom_point(data=meanq5, aes(x=y1,y=x1), colour="white") + coord_fixed()
  if(inp>5)p6 <- p6 + ggplot2::geom_point(data=meanq6, aes(x=y0,y=x0), colour="white") + geom_point(data=meanq6, aes(x=y1,y=x1), colour="white") + coord_fixed()
  #and create the plots:
  #p1_all <- allplot(p1, Q1, xmax, ymax, empty)
  #p2_all <- allplot(p2, Q2, xmax, ymax, empty)
  #p3_all <- allplot(p3, Q3, xmax, ymax, empty)
  #if(inp>3)p4_all <- allplot(p4, Q4, xmax, ymax, empty)
  #if(inp>4)p5_all <- allplot(p5, Q5, xmax, ymax, empty)
  #if(inp>5)p6_all <- allplot(p6, Q6, xmax, ymax, empty)

  #in case you want it for the whole thing instead of only quartiles:
  #pall <- ggplot2::ggplot(MR, aes(x=Lcor, y=Dcor))
  #pall <- densityplot(pall)
  #pall <- coplot(pall,xmax, ymax, max(MR$length)*0.5)
  #mppall <- kde2d(MR$Lcor[!is.na(MR$Lcor)&!is.na(MR$Dcor)],MR$Dcor[!is.na(MR$Dcor)&!is.na(MR$Lcor)])

  #mpp <- mean(range(mppall$z))
  #pall <- heatmap(pall, mpp, colc)
  #pall <- pall + geom_point(data=meantotal, aes(x=y0,y=x0), colour="white") + geom_point(data=meantotal, aes(x=y1,y=x1), colour="white")
  #pall_all <- allplot(pall, MR, xmax, ymax, empty)


  u <- list()

  u$p1 <- p1#_all
  u$p2 <- p2#_all
  u$p3 <- p3#_all
  if(inp>3)u$p4 <- p4#_all
  if(inp>4)u$p5 <- p5#_all
  if(inp>5)u$p6 <- p6#_all


  #save all histograms (L coordinates) of the quartiles too:
  p1his <- ggplot2::ggplot(Q1, aes(x=Lcor)) + geom_density(fill=colc[2], color=colc[2]) + theme_bw() + labs(x="Length(\u03BCm)") + coord_cartesian(xlim=c(-1.5,1.5))
  p2his <- ggplot2::ggplot(Q2, aes(x=Lcor)) + geom_density(fill=colc[2], color=colc[2]) + theme_bw() + labs(x="Length(\u03BCm)") + coord_cartesian(xlim=c(-1.5,1.5))
  p3his <- ggplot2::ggplot(Q3, aes(x=Lcor)) + geom_density(fill=colc[2], color=colc[2]) + theme_bw() + labs(x="Length(\u03BCm)") + coord_cartesian(xlim=c(-1.5,1.5))
  if(inp>3)p4his <- ggplot2::ggplot(Q4, aes(x=Lcor)) + geom_density(fill=colc[2], color=colc[2]) + theme_bw() + labs(x="Length(\u03BCm)") + coord_cartesian(xlim=c(-1.5,1.5))
  if(inp>4)p5his <- ggplot2::ggplot(Q5, aes(x=Lcor)) + geom_density(fill=colc[2], color=colc[2]) + theme_bw() + labs(x="Length(\u03BCm)") + coord_cartesian(xlim=c(-1.5,1.5))
  if(inp>5)p6his <- ggplot2::ggplot(Q6, aes(x=Lcor)) + geom_density(fill=colc[2], color=colc[2]) + theme_bw() + labs(x="Length(\u03BCm)") + coord_cartesian(xlim=c(-1.5,1.5))

  u$h1 <- p1his
  u$h2 <- p2his
  u$h3 <- p3his
  if(inp>3)u$h4 <- p4his
  if(inp>4)u$h5 <- p5his
  if(inp>5)u$h6 <- p6his

  #u$alls <- pall_all
  u$L <- pLD
  u$D <- pWD

  return(u)
}

###########################################double histograms!! yay!!!######################################3333
##allplot function combines histograms and density plot.

allplot <- function(plot, data, xmax, ymax, empty){

  #prepare seperate plots: histograms (hL, hD) and modified coordinate plots(remove legend )
  p1D <- plot + theme_bw() + theme(legend.position = "none")
  p1hL <- ggplot2::ggplot(data, aes(x=Lcor)) + geom_histogram() + coord_cartesian(xlim = c(-xmax, xmax)) + theme_bw() +theme(axis.title.x = element_blank())
  p1hD <- ggplot2::ggplot(data, aes(x=Dcor)) + geom_histogram() + coord_flip(xlim = c(-ymax, ymax)) + theme_bw() + theme(axis.title.y = element_blank())

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


colopts <- list(list("#000000", "#D55E00", "#F0E442"), list("#000000", "#009E73", "#F0E442"), list("#000000", "#0072B2", "#FFFFFF"), list("#000000", "#e69f00", "#F0E442"), list("#000000", "#FF0000", "#FFFF00"), list("#FFFFFF", "#F0E442", "#D55E00"))

singcollist <- list("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#D55E00", "#CC79A7")

xaxislist <- c("Cell Length (\u03BCm)", "Cell Width (\u03BCm)", "Cell Area (\u03BCm\u00B2)",
               "Spot location in length axis (\u03BCm)", "Spot location on width axis (\u03BCm)")




