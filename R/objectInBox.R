##ObjectInBox


#'@export

objectInBox <- function(objectdata, meshdata, mag = "No_PixelCorrection"){
  pixel2um <- unlist(get(magnificationList, envir=magEnv)[mag])
  message("Matching Object with Cell..")
  object_rel <- spotsInBox(objectdata,meshdata, Xs="ob_x", Ys="ob_y")$spots_relative
  colnames(object_rel)[colnames(object_rel)=="x"] <- "ob_x"
  colnames(object_rel)[colnames(object_rel)=="y"] <- "ob_y"
  objectdata <- merge(objectdata, object_rel[,c("ob_x", "ob_y", "cell", "frame")])
  objectdata <- objectdata[order(objectdata$obID, objectdata$cell),]
  message("Counting Objects per Cell..")
  ObN <- do.call('rbind',
                 lapply(unique(objectdata$frame),
                        function(x)
                          do.call('rbind',
                                  lapply(unique(objectdata$cell[objectdata$frame==x]),
                                         function(y)
                                           data.frame("frame"=x,
                                                      "cell"=y,
                                                      "obID"=unique(objectdata$obID[objectdata$frame==x&objectdata$cell==y]),
                                                      "obnum"=c(1:length(unique(objectdata$obID[objectdata$frame==x&objectdata$cell==y])))
                                           )
                                  )
                          )

                 )
  )
  objectdata <- merge(objectdata, ObN)
  message("Marking the Object centrepoints..")
  OM <- suppressWarnings(centrefun(objectdata))
  message("Putting the objects in the correct orientation..")
  OM <- suppressWarnings(midobject(meshdata, OM, pixel2um))
  message("Done.")
  return(OM)
}
