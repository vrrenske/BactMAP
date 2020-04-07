##ObjectInBox

#'@importFrom rlang .data
#'@export

objectInBox <- function(objectdata, meshdata, mag = "No_PixelCorrection"){
  pixel2um <- unlist(get(magnificationList, envir=magEnv)[mag])
  message("Matching Object with Cell..")
  object_rel <- spotsInBox(objectdata, meshdata, Xs="ob_x", Ys="ob_y")$spots_relative

  object_rel <- object_rel %>%
    dplyr::rename(ob_x = .data$x,
                  ob_y = .data$y) %>%
    dplyr::right_join(objectdata) %>%
    dplyr::group_by(.data$obID) %>%
    dplyr::mutate(pip=sum(.data$pip, na.rm=T)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(.data$pip)) %>%
    dplyr::filter(.data$pip>0) %>%
    dplyr::group_by(.data$obID) %>%
    dplyr::mutate(max.length = mean(.data$max.length, na.rm=T),
                  max.width = mean(.data$max.width, na.rm=T),
                  cell = mean(.data$cell, na.rm=T)
                  ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$l, -.data$d) %>%
    dplyr::group_by(.data$cell) %>%
    dplyr::mutate(obnum = dplyr::dense_rank(.data$obID)) %>%
    dplyr::ungroup()

  message("Marking the Object centrepoints..")
  OM <- suppressWarnings(centrefun(object_rel))
  message("Putting the objects in the correct orientation..")
  OM <- suppressWarnings(midobject(meshdata, OM, pixel2um))
  message("Done.")
  return(OM)
}
