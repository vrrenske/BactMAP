#Timelapse averaging

#To get the division percentage/cell and plot the growth per division group

#18.10.31 Renske van Raaphorst


#division function:


mark_Division <- function(timelapse){
  timelapse$division <- 1
  timelapse <- timelapse[order(timelapse$cell, timelapse$frame),]

  for(n in 1:(nrow(timelapse)-1)){
    if(timelapse$cell[n] == timelapse$cell[n+1]){
      if(timelapse$max.length[n+1]<3/4*timelapse$max.length[n]){
        timelapse$division[n+1] <- timelapse$division[n]+1
      }
      else{
        timelapse$division[n+1]<- timelapse$division[n]
      }
    }
    else{
      timelapse$division[n+1] <- 1
    }
  }
  timelapse$division <- as.numeric(timelapse$division)
  return(timelapse)
}

#percentage function:

mark_Percentage <- function(timelapse, av=TRUE){

  division_time <- aggregate(timelapse$frame,
                             by=list(timelapse$cell,
                                     timelapse$division),
                             FUN=min)
  colnames(division_time) <- c("cell", "division", "min_frame")
  division_time$max_frame <- aggregate(timelapse$frame,
                                       by=list(timelapse$cell, timelapse$division),
                                       FUN=max)$x
  division_time$av_length <- aggregate(timelapse$max.length,
                                       by=list(timelapse$cell, timelapse$division),
                                       FUN=mean)$x
  division_time$length_var <- aggregate(timelapse$max.length,
                                        by=list(timelapse$cell,
                                                timelapse$division),
                                        FUN=sd)$x
  division_time$division_time <- division_time$max_frame-division_time$min_frame +1
  division_time <- division_time[!is.na(division_time$length_var),]

  ##cutoff based on standard deviation of the mean cell length
  cutoff <- sd(division_time$length_var, na.rm=T)
  division_time <- division_time[division_time$length_var>cutoff,]


  ##cutoff based on the total division time.take out based on SD?
  meantime <- mean(division_time$division_time)
  sdtime <- sd(division_time$division_time)

  division_time <- division_time[division_time$division_time>(meantime-sdtime)&division_time$division_time<(sdtime+meantime),]

  division_time$fulldivision <- TRUE

  #put the information back into the "test" dataframe
  timelapse <- merge(timelapse, division_time, all=T)
  timelapse$fulldivision[timelapse$fulldivision!=TRUE] <- FALSE

  #then calculate the percentage
 # min_frame <- aggregate(timelapse$frame,
 #                        by=list(timelapse$cell, timelapse$division),
    #                     FUN=min)
  #colnames(min_frame) <- c("cell", "division", "min_frame")

  #timelapse <- merge(timelapse, min_frame)

  timelapse$percentage <- (timelapse$frame-timelapse$min_frame)/(timelapse$division_time-1)*100

  #there are a few cells which are "found back" after a few frames. in general those segmentations are not very good + they mess with the percentage so I remove them.
  timelapse <- timelapse[timelapse$percentage<100,]

  #binning: use "cut"
  timelapse$percentage_binned <- cut(timelapse$percentage,
                                     breaks=10,
                                     labels=paste(c(0:9)*10,c(1:10)*10, sep="-"))

  if(av==TRUE){
    mean_test <- suppressWarnings(aggregate(timelapse[colnames(timelapse)!="percentage_binned"],
                           by=list("percentage_binned" = timelapse$percentage_binned),
                           FUN=mean,na.rm=T))
    return(list("timelapse"=timelapse, "mean_by_percentage"=mean_test))
  }

  if(av==FALSE){
    return(timelapse)
  }

}

#final function for use
#' @export
perc_Division <- function(timelapse, av=TRUE, plotgrowth =TRUE){
  timelapseS <- unique(timelapse[,c("cell", "frame", "max.length")])
  timelapseS <- mark_Division(timelapseS)
  if(av==TRUE){
    out <- mark_Percentage(timelapseS, av=TRUE)
    out$timelapse <- merge(timelapse, out$timelapse)
  }
  if(av!=TRUE){
    out <- list()
    out$timelapse  <- mark_Percentage(timelapseS, av=FALSE)
    out$timelapse <- merge(timelapse, out$timelapse)
  }

  if(plotgrowth!=TRUE){return(out)}
  if(plotgrowth==TRUE){
    out$plot_growth <- plot_GrowthTime(out$timelapseS)
    if(av==TRUE){
      out$plot_avgrowth <- plot_avGrowth(out$timelapseS, out$mean_by_percentage)
    }
    return(out)
  }

}

#' @export
plot_GrowthTime <- function(timelapse, divisionmarked=TRUE, facet_division=TRUE, variable="max.length"){
  #get division percentage. of course without plotting or we would get caught in a loop!!
  if(divisionmarked==FALSE){
    timelapse <- perc_Division(timelapse, av=F, plotgrowth=F)$timelapse
  }
  #make column same as column name indicated.
  timelapse$u <- timelapse[,variable]
  p <- ggplot2::ggplot(timelapse[!is.na(timelapse$division),],
                       ggplot2::aes(x=percentage, y=u)) +
                        ggplot2::geom_line(ggplot2::aes(group=cell), alpha=0.5) +
                        ggplot2::theme_classic() +
                        ggplot2::ylab(variable)
  if(facet_division==TRUE){
    p <- p + ggplot2::facet_wrap(~division)
  }

  return(p)
}

#' @export
plot_avGrowth <- function(timelapse, mean_by_percentage, variable="max.length"){
  timelapse$u <- timelapse[,variable]
  mean_by_percentage$u <- mean_by_percentage[,variable]
  p <- ggplot2::ggplot(timelapse, ggplot2::aes(x=percentage_binned, y=u)) + ggplot2::geom_jitter(alpha=0.8) +
                                                                      ggplot2::ylab(variable) +
                                                                      ggplot2::geom_line(data=mean_by_percentage,
                                                                                         ggplot2::aes(x=as.numeric(percentage_binned), y=u), color="orange", size=2) +
                                                                      ggplot2::theme_classic()
  return(p)
}
