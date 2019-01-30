#plot tracks
#12.18.18
#R van Raaphorst
#bactMAP


plotTracks <- function(meshdata,
                       spotdata,
                       objectdata,
                       tracks=TRUE,
                       ignore_singles=FALSE,
                       movie=FALSE,
                       timepalette_lines="viridis",
                       timepalette_fill="magma",
                       dimension = c("length", "width"),
                       turn_cells = TRUE,
                       mag,
                       cell = "all"
){
  outplot <- ggplot2::ggplot() + theme_minimal()
  if(cell!="all"&is.numeric(cell)==FALSE){stop("Cell number not recognized. Give cell number or keep default (cell='all') to return all cells")}

      if(missing(meshdata)!=TRUE){
        if("length"%in%dimension&"width"%in%dimension){
        if(cell!="all"){
          meshdata <- meshdata[meshdata$cell%in%cell,]
        }
        if(turn_cells==TRUE){
          if("X_rot"%in%colnames(meshdata)){
            if("Xrotum"%in%colnames(meshdata)!=T){
              if(missing(mag)){stop("Please specify the pixel-micron conversion value 'mag'.")}
                meshdata$Xrotum <- meshdata$X_rot*unlist(get(magnificationList,envir=magEnv)[mag])
                meshdata$Yrotum <- meshdata$Y_rot*unlist(get(magnificationList,envir=magEnv)[mag])
              }
            }
            meshpart <- ggplot2::geom_polygon(data=meshdata, ggplot2::aes(x=Xrotum, y=Yrotum, color=frame), fill="black", alpha=0.01)
        }
        meshpart <- ggplot2::geom_polygon(data=meshdata, ggplot2::aes(x=X, y=Y, color=frame), fill="black", alpha=0.01)
        }
      }

      if(missing(objectdata)!=TRUE){
        if(cell!="all"){
          objectdata <- objectdata[objectdata$cell%in%cell,]
        }
        if(turn_cells==TRUE){
          if("ob_out_x"%in%colnames(objectdata)!=TRUE){
            if("ob_x"%in%colnames(objectdata)!=TRUE){
              stop("Cannot find object coordinate data (ob_out_x/ob_out_y for turned cells, ob_x/ob_y for raw coordinates.")
              if(missing(meshdata)){
                stop("Cannot find turned object data or mesh information to connect the objects to the cell localizations.
                     Replace object dataframe for a dataframe (object_relative) containing this information or add mesh data to convert dataframe.")
              }
              if(missing(mag)){
                stop("Need pixel to micron conversion factor 'mag' to correctly connect mesh to object data.")
              }

              objectdata <- suppressWarnings(centrefun(objectframe))
              objectdata <- suppressWarnings(midobject(meshdata, objectframe, get(magnificationList,envir=magEnv)[mag]))
            }
          }
          objectpart <- ggplot2::geom_polygon(data=objectdata, ggplot2::aes(x=ob_out_x,y=ob_out_y, group=obID, fill=frame), color="black", alpha=0.2)
        }
        if(turn_cells!=TRUE){
          if("ob_x"%in%colnames(objectdata)!=TRUE){
            stop("Cannot find object coordinates 'ob_x'/'ob_y'.")
          }
          objectpart <- ggplot2::geom_polygon(data=objectdata, ggplot2::aes(x=ob_x, y=ob_y, group=obID, fill=frame), color="black", alpha=0.2)
        }

        objectpart <- objectpart + scale_fill_viridis_c(option=timepalette_fill)
      }

      if(missing(spotdata)!=TRUE){
        if(cell!="all"){
          spotdata <- spotdata[spotdata$cell%in%cell,]
        }
        if(tracks==TRUE){if(ignore_singles==TRUE){spotdata <- spotdata[spotdata$trajectory!=-1,]}}
        if(turn_cells==TRUE){
          if("Lmid"%in%colnames(spotdata)!=T){
            if("l"%in%colnames(spotdata)!=T){stop("Data doesn't include relative spot positions. Please use spotsInBox() to relate spot positions to cells.")}
            if(missing(mag)){stop("Please specify the pixel-micron conversion value 'mag'.")}
            spotdata$Lmid <- spotdata$l * unlist(get(magnificationList,envir=magEnv)[mag])
            spotdata$Dum <- spotdata$d * unlist(get(magnificationList,envir=magEnv)[mag])
          }
          spotpart <- ggplot2::geom_point(data=spotdata, ggplot2::aes(x=Lmid, y=Dum, color=frame), size=2)
          if(tracks==TRUE){
            spotpart <- spotpart + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes(x=Lmid, y=Dum, group=trajectory, color=frame), size=1)
          }
        }
        if(turn_cells==FALSE){
          spotpart <- ggplot2::geom_point(data=spotdata, ggplot::aes(x=x,y=y, color=frame), size=2)
          if(tracks==TRUE){
            spotpart <- spotpart + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes(x=x,y=y,group=trajectory, color=frame),size=1)
          }

        }
      }

      if(missing(meshdata)!=T|missing(spotdata)!=T){
        plotout <- plotout + ggplot2::scale_color_viridis_c(option = timepalette_lines)
      }
      plotout <- plotout + ggplot2::coord_fixed()
    }


