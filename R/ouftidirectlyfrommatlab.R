##7/19/2017

##import functions when importing right from Matlab file.

##Renske van Raaphorst

##part of shinyspots package

############################################################

#dependencies
library(R.matlab)

############################################################

##extract cellList as a data frame with lists of data for the mesh, spots, etc.
spotrInMatOufti <- function(matfile){
  matlist <- R.matlab::readMat(matfile)
  matcellList <- matlist$cellList[[1]]
  matcelnums <- matlist$cellList[[2]]
  for(n in 1:length(matcellList)){
    for(y in 1:length(matcellList[[n]][[1]])){
      cell <- matcelnums[[n]][[1]][[y]]
      celldatas <- as.data.frame(matcellList[[n]][[1]][[y]][[1]])
      celldatas <- as.data.frame(t(celldatas))
      celldatas$cell <- cell
      celldatas$frame <- n
      if(y==1&&n==1){
        cellList <- celldatas
      }
      else{
        cellList <- rbind(cellList, celldatas)
      }
    }
  }
  return(cellList)
}

#make seperate data frame of the mesh - combined with parts of the meshturn function
spotrInMatGetMesh <- function(cellList){
  for(n in 1:nrow(cellList)){
    meshline <- as.data.frame(cellList$mesh[n])
    if(nrow(meshline)>1){
      colnames(meshline) <- c("x0", "y0", "x1", "y1")
      cell <- cellList$cell[n]
      meshline$cell <- cell
      frame <- cellList$frame[n]
      meshline$frame <- frame
      meshline$num <- 0:(nrow(meshline)-1)
      xdist <- meshline$x0 - meshline$x1
      ydist <- meshline$y0 - meshline$y1
      widths <- sqrt(xdist^2 + ydist^2)
      max_width <- max(widths)
      meshline$max.width <- max_width
      for(z in 1:nrow(meshline)){
        if(z==1){
          meshline$length[z] <- 0
          meshline$steplength[z] <- 0
        }
        if(z > 1){
          xdistL0 <- meshline$x0[z]-meshline$x0[z-1]
          ydistL0 <- meshline$y0[z]-meshline$y0[z-1]
          distL0 <- sqrt(xdistL0^2 + ydistL0^2)
          angL0 <- atan(ydistL0/xdistL0)
          angd <- atan(ydist[z]/xdist[z])
          anglength <- pi-abs(angL0)-abs(angd)
          steplength <- sin(anglength)*distL0
          meshline$steplength[z] <- steplength
          meshline$length[z] <- sum(meshline$steplength[1:z])
        }
        if(z==max(meshline$num+1)){
          angK <- atan(xdist[z-1]/ydist[z-1])
          anglengthEnd <- (pi/2)-abs(angL0-abs(angK))
          steplengthK <- sin(anglengthEnd)*distL0
          meshline$steplength[z] <- steplengthK
          meshline$length[z] <- sum(meshline$steplength[1:z])

        }
      }
      meshline$max_length <- sum(meshline$steplength)
      if(n==1){MESH <- meshline}
      if(n>1){MESH <- rbind(MESH, meshline)}
    }
  }
  return(MESH)
}

spotrInMatGetSpots <- function(cellList){
  spots <- cellList$spots
  u <- 1
  for(n in 1:length(spots)){
      l <- spots[[n]][1]
      d <- spots[[n]][2]
      x <- spots[[n]][3]
      y <- spots[[n]][4]
      position <- spots[[n]][5]
      adj_Rsquared <- spots[[n]][6]
      CI_xy <- spots[[n]][7]
      spotframe <- data.frame("l"=t(l[[1]]), "d"=t(d[[1]]), "x"=t(x[[1]]), "y"=t(y[[1]]), "position" = t(position[[1]]), "adj_Rsquared" = t(adj_Rsquared[[1]]), "CI_xy" =paste(CI_xy[[1]]))
      if(nrow(spotframe)>0){
      spotframe$frame <- cellList$frame[n]
      spotframe$cell <- cellList$cell[n]
        if(u!=1){
          spottotal <- rbind(spottotal, spotframe)
        }
        else{spottotal <- spotframe
             u <- u+1}
        }
      }
  return(spottotal)
}



matfile <- file.choose()
cellList <- spotrInMatOufti(matfile)
REP <- spotrInMatGetSpots(cellList)
MESH <- spotrInMatGetMesh(cellList)

