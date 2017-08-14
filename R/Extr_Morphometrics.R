##8-11-2017

##Renske van Raaphorst

##Extracting Morphometrics Data


##########################################################


##Dependencies:
    #library(R.matlab)


spotrExtractMorph <- function(morphpath){
  morphdata <- R.matlab::readMat(morphpath)
  framenum <- 2*(1:(length(morphdata$frame)/2))
  for(u in framenum){
    cellListf <- data.frame(t(data.frame(morphdata$frame[u])))
    cellListf$frame <- u/2
    if(u==2){
      cellList <- cellListf
    }
    else{
      cellList <- rbind(cellList, cellListf)
    }
  }
  return(cellList)
}


spotrExtractMESH <- function(cellList){

  for(n in 1:nrow(cellList)){
    meshcell <- data.frame(cellList$Xcont[n], cellList$Ycont[n], cellList$cellID[n], cellList$area[n],
                           cellList$pole1[n], cellList$pole2[n], cellList$frame[n])
    colnames(meshcell) <- c("Xcont", "Ycont", "cellID", "area", "pole1", "pole2", "frame")
    if("length"%in%colnames(cellList)){
      meshcell$length <- cellList$length[n]
      meshcell$width <- cellList$width[n]
    }
    if(n==1){
      MESH <- meshcell
    }
    else{
      MESH <- rbind(MESH, meshcell)
    }
  }
  return(MESH)

}


