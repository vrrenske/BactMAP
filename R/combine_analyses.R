####28.08.2018

###Renske van Raaphorst

##BactMAP


#Merge datasets, both manually after already analyzing the dataset & already when uploading.
#Plot a multicolor or multi-condition plot. Either next to each other or overlay depending on the type.


combineAnalyses <- function(listofdataframes, listofconditions, listofchannels){

}

#1. combine.analyses: function(listofdataframes, listofconditions, listofcolors)
combineDataframes <- function(listofdataframes, listofconditions, listofchannels){
  if(!missing(listofchannels)){
    if(length(listofdataframes)!=length(listofchannels)){
      stop("The amount of input dataframes does not match the amount of channels.")
    }
  }
  if(!missing(listofconditions)){
    if(length(listofdataframes)!=length(listofconditions)){
      stop("The amount of input dataframes does not match the amount of conditions.")
    }
  }
  listofdataframes <- lapply(c(1:length(listofdataframes)), function(x) addcolumn(listofdataframes[[x]], listofconditions[x], "condition") )
  listofdataframes <- lapply(c(1:length(listofdataframes)), function(x) addcolumn(listofdataframes[[x]], listofchannels[x], "channel"))

  finalframe <- do.call('rbind', listofdataframes)
  return(finalframe)

  }

addColumn <- function(datframe, value, name){
  datframe[,name] <- value
  return(datframe)
}
#2. plotoverlay(combinedanalyses, by=c("condition", "color" or "both)), type= length/width/projection, quantiles=n)

#3. batch upload (mesh=c(filepath, "type"), spots="same" or c(filepath, "type"), object etc)
##for this: ask users to order their files according to condition/color/analysis type.
