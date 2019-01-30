###########FTSZ ANGLE MEASUREMENT#########################################################################

#3-2-2016
#Renske van Raaphorst

#edited 5/22/2018 for bactMAP integration

##Goal: to obtain the angle of the longest axis of a polygon compared to the polygon it is sitting in.
##Practical goal: to measure the angle of FtsZ perpundicular to the long cell axis to see if the localization
##is skewed.


anglefun <- function(dat, xie, yie){
  if (!requireNamespace("shotGroups", quietly = TRUE)) {
    inp <- readline("Package 'shotGroups' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){install.packages("shotGroups")}else{stop("Canceled")}
  }
  dat$angle <- NA
  for(n in unique(dat$frame)){
    for(u in unique(dat$cell[dat$frame==n])){
      if(nrow(dat[dat$cell==u&dat$frame==n,])>0){
        if("obnum"%in%colnames(dat)==TRUE){
          for(z in dat$obnum[dat$cell==u&dat$frame==n]){
            woei <- as.matrix(dat[c(xie, yie)][dat$cell==u&dat$frame==n&dat$obnum==z,])
            an <- shotGroups::getMinBBox(woei)$angle
            dat$angle[dat$cell==u&dat$frame==n&dat$obnum==z] <- an
          }
        }
        if("obnum"%in%colnames(dat)!=TRUE){
        woei <- as.matrix(dat[c(xie, yie)][dat$cell==u&dat$frame==n,])
        an <- shotGroups::getMinBBox(woei)$angle
        dat$angle[dat$cell==u&dat$frame==n] <- an
        }
        }
    }
  }
  return(dat)
}

compareangle <- function(M, O){
  M <- unique(M[c("cell", "frame")])
  M$angle <- 0
  #not the original angle from the dataframe, but take the "turned" angle = 0 because the cells are laying on their backs now.
  O <- unique(O[c("cell", "frame", "angle", "obnum")])
  colnames(O)[3] <- "obj_an"
  MO <- merge(M, O, all=T)
  MO$dif <- MO$angle - MO$obj_an
  MO$angle <- NULL #column with 0 is a little useless
  return(MO)
}

#########CODE###########################################################################
#' @export
compare_ObjectAngle <- function (MESH, OBJ){
  MESH <- MESH[!is.na(MESH$Xrotum),]
  #MESH <- anglefun(MESH, "Xrotum", "Yrotum")
  OBJ <- OBJ[!is.na(OBJ$ob_out_x),]
  OBJ <- anglefun(OBJ, "ob_out_x", "ob_out_y")
  fincomp <- compareangle(MESH, OBJ)
  fincomp$dif <- abs(fincomp$dif)
  fincomp$dif90 <- abs(90-fincomp$dif)
  return(fincomp)
}

