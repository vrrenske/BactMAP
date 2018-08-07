##Supersegger

##3.08.2018
##Renske van Raaphorst


##Combination of functions to: * get mesh outlines from cell data * make it compatible with bactmap format


#function to get all cell files into R as a list per frame number
upperdir <- function(frameN, loc){
  cellnameList <- list.files(path=paste(loc, "xy", frameN, "\\cell", sep=""), pattern="cell")
  matlist <- lapply(cellnameList, function(x) R.matlab::readMat(paste(loc, "xy", frameN, "\\cell\\", x, sep="")))
  return(matlist)
}

#getting a list of cell files from frame 0:frames-1
readallsegcells <- function(frames, loc){
  allcells <- lapply(c(0:(frames-1)), function(x) upperdir(x, loc))
  return(allcells)
}

#getting the cell mask out of the cell lists
celmask <- function(matf){
  ras <- matf$CellA[[1]][[1]][[4]]
  ras[ras==0] <- NA
  ras <- raster::raster(ras,
                        xmn=min(matf$CellA[[1]][[1]][[2]]),
                        xmx = max(matf$CellA[[1]][[1]][[2]]),
                        ymn = min(matf$CellA[[1]][[1]][[3]]),
                        ymx = max(matf$CellA[[1]][[1]][[3]])
  )
  rasp <- raster::rasterToPolygons(ras, dissolve=TRUE)
  rasp <- as.data.frame(rasp@polygons[[1]]@Polygons[[1]]@coords)
  return(rasp)

}

#getting a list of length x (amount of cells) containing the cell mask for each cell (y)
getallmasks <- function(matlist){
  matlength <- length(matlist)
  celmasklist <- lapply(1:matlength, function(x) lapply(matlist[[x]], function(y) celmask(y)))
  return(celmasklist)
}


#turn the dataframes
flipallcells <- function(celllists){
  cellflip <- lapply(1:length(celllists), function(x) lapply(celllists[[x]], function(y) as.data.frame(t(as.data.frame(y$CellA[[1]][[1]])))))
  return(cellflip)
}

#add frame number
addcellnumframe <- function(cellflip, x){
  cellflip$frame <- x
  return(cellflip)
}

#cell number
addcellnumcell <- function(cellflip, y){
  cellflip$cell <- y
  return(cellflip)
}

#final combination function

bindallcellsandmeshes <- function(cellflip, cellmask){
  cellflipout <- list()
  for(n in 1:length(cellflip)){
    for(u in 1:length(cellflip[[n]])){
      cellflip[[n]][[u]]$cell <- u
      cellmask[[n]][[u]]$cell <- u
      cellflip[[n]][[u]]$cellmask <- list(cellmask[[n]][[u]])
      cellmask[[n]][[u]]$num <- c(1:nrow(cellmask[[n]][[u]]))
    }
    cellflipframe <- do.call('rbind', cellflip[[n]])
    cellflipframe$frame <- n
    cellmaskframe <- do.call('rbind', cellmask[[n]])
    cellmaskframe$frame <- n
    if(n==1){
      cellflipout$cellList <- cellflipframe
      cellflipout$MESH <- cellmaskframe
    }
    else{
      cellflipout$cellList <- rbind(cellflipout$cellList, cellflipframe)
      cellflipout$MESH <- rbind(cellflipout$MESH, cellmaskframe)
    }
  }
  return(cellflipout)
}


#####################################################
##extract function for exporting.

#' @export
extr.SuperSeggerCells <- function(loc, frames){
  segcells <- readallsegcells(frames=frames, loc=loc)
  segmasks <- getallmasks(segcells)
  finalframe <- bindallcellsandmeshes(flipallcells(segcells), segmasks)
  finalframe$MESH <- meshTurn(finalframe$MESH, "x", "y")
  return(finalframe)
}
