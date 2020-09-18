##8.25.2020
##Ilastik import of masks & information


#'@export
#'@title Function to take Ilastik cell information and cell masks
#'@description To be able to use this function, you need one directory with the Object Predictions (saved as TIFF) and Object Tables (saved as CSV).
#'Save them with the names the software gives them ("*_Object Predictions.TIFF" and "_table.csv").
#'
getFilesIlastik <- function(directory,name_classes = c("attack", "inside"), smoothing = TRUE, mag = "No_PixelCorrection"){

  #get the directory if missing.
  if(missing(directory)){
    directory <- getwd()
  }

  #read out mask files & tables. have to be TIFF (mask) & csv(Tables)
  masks <- list.files(directory, pattern = "Object Predictions")

  tables <- list.files(directory, pattern = "table") %>%
    lapply(function(x) read.csv(paste(directory, "/", x, sep="")))
  tables <- lapply(c(1:length(tables)), function(x)perTable(x, tables))
  #change column names to standard bactmap names
  tables <- tables %>% dplyr::bind_rows() %>%
    rename(cell = .data$object_id,
           Xmid = .data$Center.of.the.object_0,
           Ymid = .data$Center.of.the.object_1,
           area = .data$Size.in.pixels,
           class = .data$Predicted.Class)
  #continue with masks; read out
  masks <- lapply(masks, function(x) raster::raster(tiff::readTIFF(paste(directory, "/", x, sep=""))))
  #check amount of pixels per image (NB function implies SQUARE images)
  pixels <- masks[[1]]@ncols
  #convert masks to polygons & smooths them if set to TRUE
  pols <- lapply(masks, function(x) raster::rasterToPolygons(x, fun=function(x){x>0}, dissolve=TRUE))
  if(smoothing==TRUE){
    pols <- lapply(pols, function(x) smoothr::smooth(x, method="ksmooth", smoothness=4, n=2))
  }

  # per class input (standard "attack" and "inside"): turn separate polygon data points to one single data frame containing X, Y coordinates, frame, cell identifier & class.
  if(length(name_classes)>1){
    eachpol <- lapply(c(1:length(name_classes)), function(n)
      lapply(c(1:length(pols)), function(y) lapply(c(1:length(pols[[y]]@polygons[[n]]@Polygons)), function(x) data.frame(X= pols[[y]]@polygons[[n]]@Polygons[[x]]@coords[,1],
                                                                                                                                Y= pols[[y]]@polygons[[n]]@Polygons[[x]]@coords[,2],
                                                                                                                                frame = y,
                                                                                                                                cellIDpol = x,
                                                                                                                                class = name_classes[n]))) %>% bind_rows()) %>%
      bind_rows()
  }else{
    eachpol <- lapply(c(1:length(pols)), function(y) lapply(c(1:length(pols[[y]]@polygons[[1]]@Polygons)), function(x) data.frame(X= pols[[y]]@polygons[[1]]@Polygons[[x]]@coords[,1],
                                                                                                                                                Y= pols[[y]]@polygons[[1]]@Polygons[[x]]@coords[,2],
                                                                                                                                                frame = y,
                                                                                                                                                cellIDpol = x,
                                                                                                                                                class = name_classes))) %>%
      bind_rows()
  }

  #polygon conversion puts data from 0-->1 so correcting for the amount of pixels. + the image is flipped, so also flipping the Y axis.
  eachpol <- eachpol %>%
    dplyr::mutate(X = X*pixels,
                    Y = ((Y*pixels)*-1)+pixels)

  #average the x/y coordinates of the outlines to get one single point referring to it.
  av <- eachpol %>%
    dplyr::group_by(frame, cellIDpol, class) %>%
    dplyr::summarize(X=mean(X), Y=mean(Y)) %>% dplyr::ungroup()


  #find corresponding spot (some subgrouping to make sure each calculation is not done too often)
  av <- lapply(name_classes, function(n)
    lapply(unique(tables$frame[tables$class==n]), function(f)
      lapply(unique(tables$cell[tables$frame==f&tables$class==n]), function(c)
        polyPerFrameClass(av[av$frame==f&av$class==n,], tables[tables$frame==f&tables$class==n&tables$cell==c,])
        ) %>% dplyr::bind_rows()
      ) %>% dplyr::bind_rows()
    ) %>% dplyr::bind_rows()

  #add cell identifier to original X/Y coordinate set
  eachpol <- eachpol %>%
    left_join(av) %>%
    right_join(tables) %>%
    mutate(cellID = paste(.data$frame, .data$cell, .data$class, sep="_")) %>%
    select(-.data$cellIDpol) %>%
    filter(!is.na(.data$X))


  #add length and width of the minimal bounding box. note: bbox-dimensions, to properly discern length and width a medial axis / real mesh calculation would be much better!!
  bblist <- lapply(unique(eachpol$cellID), function(x) as.numeric(suppressWarnings(shotGroups::getMinBBox(data.frame(x= eachpol[eachpol$cellID==x,]$X, y=eachpol[eachpol$cellID==x,]$Y))[c("width","height")])))
  lengthlist <- lapply(c(1:length(bblist)), function(x) max(bblist[[x]]))
  widthlist <- lapply(c(1:length(bblist)), function(x) min(bblist[[x]]))
  MESH <- data.frame("cellID"=unique(eachpol$cellID), "max.length" = unlist(lengthlist), "max.width"=unlist(widthlist)) %>%
    right_join(eachpol) %>%
    arrange(.data$max.length) %>%
    mutate(num = row_number()) %>% meshTurn()

  midLine <- lapply(unique(MESH$cellID), function(x) make_Axis(MESH[MESH$cellID==x,])) %>% dplyr::bind_rows()

  MESH <- MESH %>% select(-.data$max.length) %>%
    dplyr::left_join(distinct(midLine[,c("cell","frame", "class", "max.length")])) %>%
    filter(.data$max.length<3*median(.data$max.length))
  outlist <- list()
  outlist$mesh <- MESH
  outlist$midLines <- midLine

  outlist$mesh$Xrot_micron <- outlist$mesh$X_rot * unlist(get(magnificationList, envir=magEnv)[mag])
  outlist$mesh$Yrot_micron <- outlist$mesh$Y_rot * unlist(get(magnificationList, envir=magEnv)[mag])
  outlist$mesh$max_um <- outlist$mesh$max.length * unlist(get(magnificationList, envir=magEnv)[mag])
  outlist$mesh$maxwum <- outlist$mesh$max.width * unlist(get(magnificationList, envir=magEnv)[mag])
  outlist$mesh$area_um <- outlist$mesh$area * unlist(get(magnificationList, envir=magEnv)[mag])

  return(outlist)

}

#make data frame of each table including the frame number.
perTable <- function(x, tables){
  tables <- tables[[x]]
  tables$frame <- x
  tables$User.Label <- NULL
  return(tables)
}


polyPerFrameClass <- function(av, tables){
  pplist <- sp::point.in.polygon(point.x = av$X, point.y=av$Y, pol.x = c(rep(tables$Bounding.Box.Maximum_0, 2), rep(tables$Bounding.Box.Minimum_0, 2)), pol.y= c(tables$Bounding.Box.Maximum_1, tables$Bounding.Box.Minimum_1, tables$Bounding.Box.Minimum_1, tables$Bounding.Box.Maximum_1))
  av$pp <- pplist
  av <- av[av$pp==1,]
  av$cell <- tables$cell
  av <- av %>% dplyr::select(-.data$X, -.data$Y)
  return(av[av$pp==1,])
}







##make medial-axis line based on voronoi diagram (based on described method used in Microbetracker)
##This will work the best when the cells are already turned based on their minimal bounding box! So use this only later in the process.

#needs package 'deldir'


make_Axis <- function(oneCell){

  ##get coordinates of voronoi using the deldir function
  vD <- deldir::deldir(oneCell$X_rot, oneCell$Y_rot)$dirsgs
  ##take only the vertices that are not on the edge of the diagram (in practice, the ones inside the cell)
  vD <- vD[vD$bp2==TRUE,] %>%
    select(.data$x1, .data$y1) %>%
    distinct() %>%
    filter(.data$x1>quantile(.data$x1, c(0.25))) %>%
    filter(.data$x1<quantile(.data$x1, c(0.75)))

  vDp <- data.frame(predict(smooth.spline(vD$x1, vD$y1), seq(min(oneCell$X_rot), max(oneCell$X_rot), by=0.4))) %>%
    dplyr::mutate(cell=unique(oneCell$cell),
                  frame= unique(oneCell$frame),
                  class=unique(oneCell$class),
                  dist = sqrt((.data$x-dplyr::lead(.data$x))^2+(.data$y-dplyr::lead(.data$y))^2),
                  max.length=sum(.data$dist, na.rm=T))
  return(vDp)

}

