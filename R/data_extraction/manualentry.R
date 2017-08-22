####manual upload script

##8-11-17

#Renske


#########################################################################################

datafile <- file.choose(caption="Choose .txt, .csv or .xls(x) file containing your data")

if(datafile[-3]=="txt"){
  delimeter <- readline("Seperator: ")
  titlesyesno <- readline("Column titles in dataset? Y/N: ")
  if(titlesyesno=="Y"|titlesyestno=="y"){
    datlist <- read.table(datafile, sep = delimeter, title=T)

    print(paste("Current column names are: ", colnames(datlist)))

    if("cell"%in%colnames(datlist)!=T&"cellID"%in%colnames(datlist)!=T){
      cellcol <- readline("Name column containing cell ID's: ")
      datlist$cell <- datlist[,cellcol]
      datlist[,cellcol] <- NULL
    }
    if("cellID"%in%colnames(datlist)){
        datlist$cell <- datlist$cellID
        datlist$cellID <- NULL
    }
    if("frame"%in%colnames(datlist))
  }
  if(titlesyesno=="N"|titlesyesno=="n"){
    datlist <- read.table(datafile, sep = delimeter, title=F)

  }


}
