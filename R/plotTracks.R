#plot tracks
#12.18.18
#R van Raaphorst
#bactMAP

#'@export
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
                       cell = "all",
                       transparency = 0.2
){
  outplot <- ggplot2::ggplot() + ggplot2::theme_minimal()
  if(cell[1]!="all"&is.numeric(cell)==FALSE){stop("Cell number not recognized. Give cell number or keep default (cell='all') to return all cells")}

  if(missing(meshdata)!=TRUE){
    if(cell[1]!="all"){
      meshdata <- meshdata[meshdata$cell%in%cell,]
    }
    if("length"%in%dimension&"width"%in%dimension){

      if(turn_cells==TRUE){
        if("X_rot"%in%colnames(meshdata)){
          if("Xrotum"%in%colnames(meshdata)!=T){
            if(missing(mag)){stop("Please specify the pixel-micron conversion value 'mag'.")}
            meshdata$Xrotum <- meshdata$X_rot*unlist(get(magnificationList,envir=magEnv)[mag])
            meshdata$Yrotum <- meshdata$Y_rot*unlist(get(magnificationList,envir=magEnv)[mag])
          }
        }
        outplot <- outplot +  ggplot2::geom_polygon(data=meshdata, ggplot2::aes(x=Xrotum, y=Yrotum, color=frame, group=frame), fill="black", alpha=transparency/5)
      }
      if(turn_cells==FALSE){
        outplot <- outplot + ggplot2::geom_polygon(data=meshdata, ggplot2::aes(x=X, y=Y, color=frame, group=frame), fill="black", alpha=transparency/5)
      }
    }
  }

  if(missing(objectdata)!=TRUE){
    if(cell[1]!="all"){
      objectdata <- objectdata[objectdata$cell%in%cell,]
    }
    if(turn_cells==TRUE|("length"%in%dimension&!"width"%in%dimension)|(!"length"%in%dimension&"width"%in%dimension)){
      if("ob_out_x"%in%colnames(objectdata)!=TRUE){
        if("ob_x"%in%colnames(objectdata)!=TRUE){
          stop("Cannot find object cobjectdatardinate data (ob_out_x/ob_out_y for turned cells, ob_x/ob_y for raw cobjectdatardinates.")
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
      if("length"%in%dimension&"width"%in%dimension){
        outplot <- outplot + ggplot2::geom_polygon(data=objectdata, ggplot2::aes(x=ob_out_x,y=ob_out_y, group=(paste(obID,frame)), fill=frame), color="black", alpha=transparency)
      }
      if("length"%in%dimension&!"width"%in%dimension){
        Om <- aggregate(objectdata[,c("ob_out_x", "pole1", "pole2")], by=list("cell"=objectdata$cell, "frame"=objectdata$frame, "obID"=objectdata$obID), FUN=max)
        Omin <- aggregate(objectdata[,c("ob_out_x", "pole1", "pole2")], by=list("cell"=objectdata$cell, "frame"=objectdata$frame, "obID"=objectdata$obID), FUN=min)
        Om$obxmax <- Om$ob_out_x
        Om$ob_out_x <- NULL
        Omin$obxmin <- Omin$ob_out_x
        Omin$ob_out_x <- NULL
        objectdata <- merge(Om, Omin)

        objectdata$framemin <- objectdata$frame - 0.1
        objectdata$framemax <- objectdata$frame + 0.1
        outplot <- outplot + ggplot2::geom_polygon(data=objectdata, ggplot2::aes(xmin = obxmin, xmax=obxmax, ymin=framemin, ymax = framemax, group=paste(obID, frame), fill=frame))
        if(missing(spotdata)){
          outplot <- outplot + ggplot2::geom_path(data=objectdata, ggplot2::aes(x=pole1, y=frame)) + ggplot2::geom_path(data=objectdata, ggplot2::aes(x=pole2, y=frame))
          outplot <- outplot + ggplot2::xlab("location on cell length axis (\u00b5m)") + ggplot2::ylab("Time (Frames)")
        }
      }
      if(!"length"%in%dimension&"width"%in%dimension){
        Om <- aggregate(objectdata[,c("ob_out_y", "maxwum")], by=list("cell"=objectdata$cell, "frame"=objectdata$frame, "obID"=objectdata$obID), FUN=max)
        Omin <- aggregate(objectdata[,c("ob_out_y", "maxwum")], by=list("cell"=objectdata$cell, "frame"=objectdata$frame, "obID"=objectdata$obID), FUN=min)
        Om$obymax <- Om$ob_out_y
        Om$ob_out_y <- NULL
        Omin$obymin <- Omin$ob_out_y
        Omin$ob_out_y <- NULL
        objectdata <- merge(Om, Omin)

        objectdata$framemin <- objectdata$frame - 0.1
        objectdata$framemax <- objectdata$frame + 0.1
        outplot <- outplot + ggplot2::geom_polygon(data=objectdata, ggplot2::aes(xmin = obymin, xmax=obymax, ymin=framemin, ymax = framemax, group=paste(obID, frame), fill=frame))
        if(missing(spotdata)){
          objectdata$pole1 <- objectdata$maxwum/2
          objectdata$pole2 <- objectdata$pole1*-1
          outplot <- outplot + geom_path(data=objectdata, aes(x=pole1, y=frame)) + geom_path(data=objectdata, aes(x=pole2, y=frame))
          outplot <- outplot + ggplot2::xlab("location on cell length axis (\u00b5m)") + ggplot2::ylab("Time (Frames)")}
      }
    }
    if(turn_cells!=TRUE&"length"%in%dimension&"width"%in%dimension){
      if("ob_x"%in%colnames(objectdata)!=TRUE){
        stop("Cannot find object cobjectdatardinates 'ob_x'/'ob_y'.")
      }
      outplot <- outplot + ggplot2::geom_polygon(data=objectdata, ggplot2::aes(x=ob_x, y=ob_y, group=paste(obID,frame), fill=frame), color="black", alpha=transparency)
    }

    outplot <- outplot + scale_fill_viridis_c(option=timepalette_fill)
  }

  if(missing(spotdata)!=TRUE){
    if(cell[1]!="all"){
      spotdata <- spotdata[spotdata$cell%in%cell,]
    }
    if(tracks==TRUE){if(ignore_singles==TRUE){spotdata <- spotdata[spotdata$trajectory!=-1,]}}
    if(turn_cells==TRUE|("length"%in%dimension&!"width"%in%dimension)|(!"length"%in%dimension&"width"%in%dimension)){
      if("Lmid"%in%colnames(spotdata)!=T){
        if("l"%in%colnames(spotdata)!=T){
          if(missing(meshdata)!=T){
            A <- readline("Data doesn't include relative spot positions. Press 'Y'+enter to start spotsInBox() to relate spot positions to cells or any other key to stop the function.")
            if(A=="Y"|A=="y"){
              spotdata <- spotsInBox(spotdata, meshdata)$spots_relative
            }else{stop("Function stopped.")}
          }
        }
        if(missing(mag)){stop("Please specify the pixel-micron conversion value 'mag'.")}
        spotdata$Lmid <- spotdata$l * unlist(get(magnificationList,envir=magEnv)[mag])
        spotdata$Dum <- spotdata$d * unlist(get(magnificationList,envir=magEnv)[mag])
      }
      if("length"%in%dimension&"width"%in%dimension){
        outplot <- outplot + ggplot2::geom_point(data=spotdata, ggplot2::aes(x=Lmid, y=Dum, color=frame), size=2)
        if(tracks==TRUE){
          outplot <- outplot + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes(x=Lmid, y=Dum, group=trajectory, color=frame), size=1)
        }
      }
      if("length"%in%dimension&!"width"%in%dimension){
        outplot <- outplot + ggplot2::geom_point(data=spotdata, ggplot2::aes(x=Lmid, y=frame, color=frame), size=2)
        if(tracks==TRUE){
          outplot <- outplot + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes(x=Lmid, y=frame, group=trajectory, color=frame), size=1)
        }
        if(!"pole1"%in%colnames(spotdata)){
          if(!"maxum"%in%colnames(spotdata)){
            spotdata$maxum <- spotdata$max.length*unlist(get(magnificationList,envir=magEnv)[mag])
          }
          spotdata$pole1 <- 0.5*spotdata$maxum
          spotdata$pole2 <- -spotdata$pole1
          outplot <- outplot + ggplot2::geom_path(data=spotdata, ggplot2::aes(x=pole1, y=frame)) + ggplot2::geom_path(data=spotdata, ggplot2::aes(x=pole2, y=frame))
        }

      }
    }
    if(turn_cells==FALSE&"length"%in%dimension&"width"%in%dimension){
      outplot <- outplot + ggplot2::geom_point(data=spotdata, ggplot2::aes(x=x,y=y, color=frame), size=2)
      if(tracks==TRUE){
        outplot <- outplot + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes(x=x,y=y,group=trajectory, color=frame),size=1)
      }

    }
  }

  if(missing(meshdata)!=T|missing(spotdata)!=T){
    outplot <- outplot + ggplot2::scale_color_viridis_c(option = timepalette_lines)
  }
  if(turn_cells==TRUE&"length"%in%dimension&"width"%in%dimension){
    outplot <- outplot + ggplot2::cobjectdatard_fixed() + facet_wrap(~cell)
  }
  if(turn_cells==FALSE|("length"%in%dimension&!"width"%in%dimension)|(!"length"%in%dimension&"width"%in%dimension)){
    outplot <- outplot + facet_wrap(~cell, scales="free")
  }
  return(outplot)
}


