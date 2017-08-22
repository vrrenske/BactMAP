###DATA SUMMARY CREATION ####


##Renske van Raaphorst

##08-14-2017

#################################################################################################


#upload (or not)
if(exists("cellList")){
  U <- readline("Use cellList currently open in environment? Y/N: ")
  }
if(U=="N"|U=="n"|exists("cellList")==F){
  datafile <- file.choose()
  NUM <- readline("Data origin:\n\t 1. Oufti\n\t 2. Microbetracker\n\t 3. MicrobeJ\n\t 4. Morphometrics\n\t 5. ObjectJ\n\t 6. ISBatch\u\t 7. Utrack\u\t 8. Manual Entry")
  if(NUM=="1"){
    cellList <- spotrInMatOufti(datafile)
    s <- readline("Extract spot data? Y/N: ")
    if(s=="Y"|s=="y"){REP <- spotrInMatGetSpots(cellList)}
    MESH <- spotrInMatGetMesh(cellList)
  }
  if(NUM=="2"){
    print("Function will be added soon!")
  }
  if(NUM=="3"){
    print("Function will be added soon!")
  }
  if(NUM=="4"){
    cellList <- spotrExtractMorph(datafile)
    MESH <- spotrExtractMES(cellList)
    s <- readline("Also upload spot data? Y/N: ")
    if(s=="Y"|s=="y"){
      spotdatafile <- file.choose()
      spotNUM <- readline("Data origin:\n\t 1. Oufti\n\t 2. MicrobeJ\n\t 3. ObjectJ\n\t 4. ISBatch\u\t 5. Utrack\u\t 6. Manual Entry")
      if(spotNUM=="1"){
        REP <- spotrInMatGetSpots(spotrInMatOufti(spotdatafile))
        ##need function putspotsinboxes here
      }
      if(spotNUM=="2"){
        print("Function will be added soon!")
      }
      if(spotNUM=="3"){
        print("Function will be added soon!")
      }
      if(spotNUM=="4"){

      }

    }
  }
  if(NUM=="5"){
    print("Function will be added soon!")
  }
  if(NUM=="6"){
    print("Function will be added soon!")
  }
  if(NUM=="7"){
    source("Utrackprocessing.R")
  }
  if(NUM=="8"){
    source("manualentry.R")
  }
}

#




