#
# Renske van Raaphorst 23 september 2016
# First test to load a file, plot a custom dot/violin plot. 
# Server File
#

library(shiny)

source("functionstransform.R")
source("plothelper.R")

function(input, output) {

  #### Input for tab 1  
  dataInput <- reactive({
    inFile <- input$file1
    if(is.null(inFile))
      return(NULL)
    C <- read.delim(inFile$datapath, header=FALSE)
    C <- prepareouf(C, input$U)
    return(C)
  }
  )
  
  meshInput <- reactive({
    M <- ouftimesh(dataInput())
    return(M)
  })
  
  spotInput <- reactive({
    if(is.null(input$Y)) 
      return(NULL)
    else{S <- ouftispot(dataInput(), as.numeric(input$mag))
    return(S)}
  })
  
  mInput <- reactive({
    msmall <- ouftiM(dataInput(), meshInput(), as.numeric(input$mag))
    return(msmall)
  })
  
  oInput <- reactive({
    if(is.null(input$U))
      return(NULL)
    else{O <- ouftiobject(dataInput())
    return(O)}
  })
  
  plotOutsstuff <- reactive({
    MR <- mergeframes(spotInput(), mInput())
    MR <- spotMR(MR)
    MR <- LimDum(MR)
    p1 <- plotsout(MR, (as.numeric(input$radio)+2), meshInput(), colopts[[as.numeric(input$colorchoice)]], as.numeric(input$mag))
    return(p1)
  })
  
  Mmid <- reactive({
    MR <- mergeframes(spotInput(), mInput())
    MR <- LimDum(MR)
    mpL <- kde2d(MR$num[!is.na(MR$Lmid)], MR$Lmid[!is.na(MR$Lmid)])
    mpL1 <- median(range(mpL$z))
    return(mpL1)
  })
  
  plothisinp <- reactive({
    M <- mergeframes(spotInput(),mInput())
    if(input$Vars[1]==1)M$u <- M$length
    if(input$Vars[1]==2)M$u <- M$max.width
    if(input$Vars[1]==3)M$u <- M$area*as.numeric(input$mag)^2
    if(input$Vars[1]==4)M$u <- M$Lmid
    if(input$Vars[1]==5)M$u <- M$Dum
    if(length(input$Vars)>1){
      if(input$Vars[2]==1)M$t <- M$length
      if(input$Vars[2]==2)M$t <- M$max.width
      if(input$Vars[2]==3)M$t <- M$area*as.numeric(input$mag)^2
      if(input$Vars[2]==4)M$t <- M$Lmid
      if(input$Vars[2]==5)M$t <- M$Dum
    }
    return(M)
  })

  plothis1 <- reactive({
    M <- plothisinp()
    M$b <- "a"
    z <- ggplot(M, aes(x=u)) + theme_minimal()
    if(input$Typep==1){
      z <- z + geom_histogram(binwidth=input$bins/200, color=singcollist[as.numeric(input$Cols3)], fill=singcollist[as.numeric(input$Cols3)]) + xlab(xaxislist[as.numeric(input$Vars[1])])
    }
    if(input$Typep==2){
      z <- z + geom_density(adjust=input$bins/50, color=singcollist[as.numeric(input$Cols3)], fill=singcollist[as.numeric(input$Cols3)]) + xlab(xaxislist[as.numeric(input$Vars[1])])
    }
    if(input$Typep==3){
      z <- z+ geom_dotplot(binwidth = input$bins/900,  stackdir = "center", color=singcollist[as.numeric(input$Cols3)], fill=singcollist[as.numeric(input$Cols3)]) + coord_flip() + ylab(xaxislist[as.numeric(input$Vars[1])]) +  xlab("Count") + theme(axis.text.x=element_blank())
        }
    if(input$Typep==4){
      z <- ggplot(M, aes(factor(b), u)) + geom_violin(scale="area",adjust=input$bins/50, color=singcollist[as.numeric(input$Cols3)], fill=singcollist[as.numeric(input$Cols3)]) + ylab(xaxislist[as.numeric(input$Vars[1])]) + xlab("Count") +theme_minimal() + theme(legend.position="none", axis.text.x=element_blank())
    }
    return(z)
  })
  
  plothis2 <- reactive({
    if(length(input$Vars)==1){return(NULL)}
    else{
    M <- plothisinp()
    M$b <- "count"
    z <- ggplot(M, aes(x=t)) + theme_minimal()
    if(input$Typep==1){
      z <- z + geom_histogram(binwidth=input$bins2/200, color=singcollist[as.numeric(input$Cols4)], fill=singcollist[as.numeric(input$Cols4)]) + xlab(xaxislist[as.numeric(input$Vars[2])])
    }
    if(input$Typep==2){
      z <- z + geom_density(adjust=input$bins2/50, color=singcollist[as.numeric(input$Cols4)], fill=singcollist[as.numeric(input$Cols4)]) + xlab(xaxislist[as.numeric(input$Vars[2])])
    }
    if(input$Typep==3){
      z <- z+ geom_dotplot(binwidth = input$bins2/500,  stackdir = "center", color=singcollist[as.numeric(input$Cols4)], fill=singcollist[as.numeric(input$Cols4)]) + coord_flip()  + ylab(xaxislist[as.numeric(input$Vars[2])]) + xlab("Count") + theme(axis.text.x=element_blank())
    }
    if(input$Typep==4){
      z <- ggplot(M, aes(factor(b), t)) + geom_violin(scale="area",adjust=input$bins2/50, color=singcollist[as.numeric(input$Cols4)], fill=singcollist[as.numeric(input$Cols4)]) + ylab(xaxislist[as.numeric(input$Vars[2])])+ xlab("Count") + theme_minimal() + theme(legend.position="none", axis.text.x=element_blank())
    }
    return(z)
    }
  })
  
  plotcombo <- reactive({
    if(length(input$Vars)>1){
    M <- plothisinp()
    z <- ggplot(M, aes(x=u, y=t)) + theme_minimal() + geom_point(color=singcollist[as.numeric(input$Cols3)], fill=singcollist[as.numeric(input$Cols3)])  + xlab(xaxislist[as.numeric(input$Vars[1])]) + ylab(xaxislist[as.numeric(input$Vars[2])])
    return(z) 
    }
    if(length(input$Vars)==1) return(NULL)
    })
  
  #### display & download transformed table tab 2, subtab 1
  
  output$contents <- renderTable(meshInput())

  
  output$downloadData <- downloadHandler(filename = function(){paste(input$file1, "meshdata.csv", sep="_")},
                                         content = function(file){write.csv(meshInput(), file)
                                           }
                                         )
  ####display & download spot data in tab 2, subtab 2
  
  output$spots <- renderTable(spotInput())
  
  output$downloadSpots <- downloadHandler(filename = function(){paste(input$file1, "spotdata.csv", sep="_")},
                                         content = function(file){write.csv(spotInput(), file)
                                         }
  )
  
  ####display & download spot data in tab 2, subtab 3
  
  output$lengthwidth <- renderTable(mInput())
  
  output$downloadM <- downloadHandler(filename = function(){paste(input$file1, "lengthwidthdata.csv", sep="_")},
                                          content = function(file){write.csv(mInput(), file)
                                          }
  )
  
  ####display & download spot data in tab 2, subtab 4
  
  output$objects <- renderTable(oInput())
  
  output$downloadObs <- downloadHandler(filename = function(){paste(input$file1, "Objectdata.csv", sep="_")},
                                          content = function(file){write.csv(oInput(), file)
                                          }
  )
  
  ####display the choices for the color scheme
  output$colchoiceplot1 <- renderPlot(colorchoiceplot(colopts[[1]], 1))
  
  output$colchoiceplot2 <- renderPlot(colorchoiceplot(colopts[[2]],2))
  
  output$colchoiceplot3 <- renderPlot(colorchoiceplot(colopts[[3]],3))
  
  output$colchoiceplot4 <- renderPlot(colorchoiceplot(colopts[[4]],4))
  
  output$colchoiceplot5 <- renderPlot(colorchoiceplot(colopts[[5]],5))
  
  output$colchoiceplot6 <- renderPlot(colorchoiceplot(colopts[[6]],6))
  
  ####display the quartile plots + histograms
  
  output$Q1 <- renderPlot(plotOutsstuff()$p1)
  output$Q2 <- renderPlot(plotOutsstuff()$p2)
  output$Q3 <- renderPlot(plotOutsstuff()$p3)
  output$Q4 <- renderPlot({if(is.null(plotOutsstuff()$p4))return(NULL) else plotOutsstuff()$p4})
  output$Q5 <- renderPlot({if(is.null(plotOutsstuff()$p5))return(NULL) else plotOutsstuff()$p5})
  output$Q6 <- renderPlot({if(is.null(plotOutsstuff()$p6))return(NULL) else plotOutsstuff()$p6})
  
  output$h1 <- renderPlot(plotOutsstuff()$h1)
  output$h2 <- renderPlot(plotOutsstuff()$h2)
  output$h3 <- renderPlot(plotOutsstuff()$h3)
  output$h4 <- renderPlot({if(is.null(plotOutsstuff()$h4))return(NULL) else plotOutsstuff()$h4})
  output$h5 <- renderPlot({if(is.null(plotOutsstuff()$h5))return(NULL) else plotOutsstuff()$h5})
  output$h6 <- renderPlot({if(is.null(plotOutsstuff()$h6))return(NULL) else plotOutsstuff()$h6})
  
  output$downloadplots <- downloadHandler(filename = function(){paste(input$file1, "groupplots.pdf", sep="_")}, 
                                          content=function(file){if(input$radio==1)ggsave(arrangeGrob(plotOutsstuff()$p1, 
                                                                                                                   plotOutsstuff()$p2, 
                                                                                                                   plotOutsstuff()$p3, 
                                                                                                                   ncol=1), 
                                                                                                                   filename=file, width=10, height=12)
                                                                 if(input$radio==2)ggsave(arrangeGrob(plotOutsstuff()$p1, 
                                                                                                                   plotOutsstuff()$p2, 
                                                                                                                   plotOutsstuff()$p3, 
                                                                                                                   plotOutsstuff()$p4,
                                                                                                                   ncol=1), 
                                                                                                                   filename=file, width=10, height=16)
                                            
                                                                 if(input$radio==3)ggsave(arrangeGrob(plotOutsstuff()$p1, 
                                                                                                                       plotOutsstuff()$p2, 
                                                                                                                       plotOutsstuff()$p3, 
                                                                                                                       plotOutsstuff()$p4,
                                                                                                                       plotOutsstuff()$p5,
                                                                                                                       ncol=1), 
                                                                                                                       filename=file, width=10, height=20)
                                                                 if(input$radio==4)ggsave(arrangeGrob(plotOutsstuff()$p1, 
                                                                                                                         plotOutsstuff()$p2, 
                                                                                                                         plotOutsstuff()$p3, 
                                                                                                                         plotOutsstuff()$p4,
                                                                                                                         plotOutsstuff()$p5,
                                                                                                                         plotOutsstuff()$p6,
                                                                                                                         ncol=1), 
                                                                                                             filename=file, width=10, height=24)
                                            }
                                          
  )

  output$downloadhist <- downloadHandler(filename = function(){paste(input$file1, "histplots.pdf", sep="_")}, 
                                          content=function(file){if(input$radio==1)ggsave(arrangeGrob(plotOutsstuff()$h1, 
                                                                                                      plotOutsstuff()$h2, 
                                                                                                      plotOutsstuff()$h3, 
                                                                                                      ncol=1), 
                                                                                          filename=file, width=10, height=12)
                                            if(input$radio==2)ggsave(arrangeGrob(plotOutsstuff()$h1, 
                                                                                 plotOutsstuff()$h2, 
                                                                                 plotOutsstuff()$h3, 
                                                                                 plotOutsstuff()$h4,
                                                                                 ncol=1), 
                                                                     filename=file, width=10, height=16)
                                            
                                            if(input$radio==3)ggsave(arrangeGrob(plotOutsstuff()$h1, 
                                                                                 plotOutsstuff()$h2, 
                                                                                 plotOutsstuff()$h3, 
                                                                                 plotOutsstuff()$h4,
                                                                                 plotOutsstuff()$h5,
                                                                                 ncol=1), 
                                                                     filename=file, width=10, height=20)
                                            if(input$radio==4)ggsave(arrangeGrob(plotOutsstuff()$h1, 
                                                                                 plotOutsstuff()$h2, 
                                                                                 plotOutsstuff()$h3, 
                                                                                 plotOutsstuff()$h4,
                                                                                 plotOutsstuff()$h5,
                                                                                 plotOutsstuff()$h6,
                                                                                 ncol=1), 
                                                                     filename=file, width=10, height=24)
                                          }
                                          
  )
  
  output$Lplot <- renderPlot(plotOutsstuff()$L)
  output$Wplot <- renderPlot(plotOutsstuff()$D)
  
  output$Plotvar1 <- renderPlot(plothis1())
  output$Plotvar2 <- renderPlot(plothis2())
  output$Combpplot <- renderPlot(plotcombo())
  
  
  output$Downloaddimensions <- downloadHandler(filename = function(){paste(input$file1, "dimensions.pdf", sep="_")},
                                      content = function(file){
                                        if(length(input$Vars)==1){
                                          ggsave(plothis1(), filename=file, height=4, width=3)
                                        }
                                        else ggsave(arrangeGrob(plothis1(), plothis2(), plotcombo(), ncol=1), filename=file, height = 16, width = 6)
                                      }
  )
}