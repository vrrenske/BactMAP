##ObjectInBox


#'@export

objectInBox <- function(objectdata, meshdata, mag = "No_PixelCorrection"){
  pixel2um <- unlist(get(magnificationList, envir=magEnv)[mag])
  message("Matching Object with Cell..")
  object_rel <- spotsInBox(objectdata, meshdata, Xs="ob_x", Ys="ob_y")$spots_relative

  object_rel <- object_rel %>%
    dplyr::rename(ob_x = x,
                  ob_y = y) %>%
    dplyr::right_join(objectdata) %>%
    dplyr::group_by(obID) %>%
    dplyr::mutate(pip=sum(pip, na.rm=T)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(pip)) %>%
    dplyr::filter(pip>0) %>%
    dplyr::group_by(obID) %>%
    dplyr::mutate(max.length = mean(max.length, na.rm=T),
                  max.width = mean(max.width, na.rm=T),
                  cell = mean(cell, na.rm=T)
                  ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-l, -d) %>%
    dplyr::group_by(cell) %>%
    dplyr::mutate(obnum = dplyr::dense_rank(obID)) %>%
    dplyr::ungroup()

  message("Marking the Object centrepoints..")
  OM <- suppressWarnings(centrefun(object_rel))
  message("Putting the objects in the correct orientation..")
  OM <- suppressWarnings(midobject(meshdata, OM, pixel2um))
  message("Done.")
  return(OM)
}
