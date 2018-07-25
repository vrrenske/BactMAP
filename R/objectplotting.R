##24.7.2018
##Renske van Raaphorst

##function to plot object projections (now only possible with oufti output)

#' @export
plotobjects <- function(obdat, meshdat, groups=1, cellcolor="black", objectcolor="green", transparancy=0.1, getdata=FALSE){
  if(getdata==TRUE){outlist <- list()}

  obdat <- obdat[order(obdat$max.length, obdat$obnum, obdat$obpath),]
  if(groups>1){
    obdat$lnum <- c(1:nrow(obdat))
    obdat$q1 <- cut(obdat$lnum, breaks=groups, labels=c(1:groups))
  }
  if(groups==1){
    obdat$q1 <- 1
  }

  if(getdata==TRUE){outlist$object_data <- obdat}
  if(missing(meshdat)==TRUE){
    outplot <- ggplot2::ggplot(obdat, ggplot2::aes(x=ob_out_x, y=ob_out_y, group=obID)) +
      ggplot2::geom_polygon(alpha=transparancy, fill=objectcolor) +
      ggplot2::facet_grid(q1~.) +
      ggplot2::coord_fixed() +
      ggplot2::theme_classic() +
      ggplot2::xlab("X coordinates (micron)") +
      ggplot2::ylab("Y coordinates (micron)")

    if(getdata==TRUE){outlist$objectplot <- outplot}
  }
  if(missing(meshdat)!=TRUE){
    meshdat <- merge(meshdat, obdat[,c("frame", "cell", "q1")])
    meshdat <- meshdat[order(meshdat$frame, meshdat$cell, meshdat$num),]
    p2um <- unlist(get(magnificationList, envir=magEnv)[mag])
    outplot <- ggplot2::ggplot(meshdat) +
      ggplot2::geom_polygon(ggplot2::aes(x=X_rotum, y=Y_rotum, group=cell), alpha=transparancy, fill=cellcolor) +
      ggplot2::geom_polygon(data=obdat, ggplot2::aes(x=ob_out_x, y=ob_out_y, group=obID), alpha=transparancy, fill=objectcolor) +
      ggplot2::facet_grid(q1~.) +
      ggplot2::coord_fixed() +
      ggplot2::theme_classic() +
      ggplot2::xlab("X coordinates (micron)") +
      ggplot2::ylab("Y coordinats (micron)")

    if(getdata==TRUE){
      outlist$objectplot <- outplot
      outlist$mesh_data <- meshdat
    }

  }
  if(getdata==TRUE){return(outlist)}
  if(getdata==FALSE){return(outplot)}
}
