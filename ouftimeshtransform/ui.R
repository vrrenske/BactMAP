#
# Renske van Raaphorst 23 sept 2016
# UI interface for first violin/dotplot app
#

library(shiny)

# Define UI for application that draws a histogram
navbarPage("SpotProcessR",
           tabPanel("Upload file",
           h3("Analysis of single Oufti File"),
           sidebarLayout(
            sidebarPanel(
              fileInput('file1', 'Choose .csv or .txt File',
                        accept=c('.csv', '.txt')),
              checkboxInput("U", "Object analysis", value = FALSE, width = NULL),
              checkboxInput("Y", "Spot analysis", value = FALSE, width = NULL),
              textInput("mag", "pixel to um (use . and no ,)")
            ),
            mainPanel(
            h1("SpotProcessR app in Construction"),
            p("Please upload your segmentation and spot/object detection data.Afterwards, SpotProcessR will order 
              & visualize your data. You can download the ordered data sets, plots and compare your datasets to 
              each other."),
            strong("More features to be added soon!",style = "color:#D55E00")
            )
           )
           ),
           
           tabPanel("Download data",
                      navlistPanel(
                        tabPanel("Mesh Data", 
                                 downloadButton('downloadData', "Download the Meshes here (.csv)"),
                                 tableOutput('contents')
                        ),
                        tabPanel("Spot Data",
                          downloadButton('downloadSpots', "Download the Spots here (.csv)"),
                          tableOutput('spots')
                        ),
                        tabPanel("Length/width Data",
                          downloadButton('downloadM', "Download only length/width data here (.csv)"),
                          tableOutput('lengthwidth')
                          
                        ),
                        tabPanel("Object Data",
                          downloadButton('downloadObs', "Download the Objects here (.csv)"),
                          tableOutput('objects')
                          
                        )
                      )
           ),
           
           tabPanel("Plot Spot Coordinates",
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("radio",
                                     "Amount of Groups",
                                     choices = c("3"=1, "4"= 2, "5"=3, "6"=4), selected = 2),
                        selectInput("colorchoice",
                                    "Pick color to download plot",
                                    choices = c("Orange Hot" = 1, "Green Hot" = 2, "Cold Blue" = 3, "Yellow Hot" = 4, "Red Hot" = 5, "White Orange" = 6),
                                    selected = 2),
                        downloadButton("downloadplots", "Download heat maps here"),
                        downloadButton("downloadhist", "Dowload histograms here"),
                        downloadButton("downloadL", "Download the length heat map here"),
                        downloadButton("downloadW", "Dowload the width heat map here")
                        
                      ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Pick your color",
                                   fluidRow(
                                     verticalLayout(
                                     splitLayout(
                                       plotOutput("colchoiceplot1"),
                                       plotOutput("colchoiceplot2"),
                                       plotOutput("colchoiceplot3")),
                                     splitLayout(
                                       plotOutput("colchoiceplot4"),
                                       plotOutput("colchoiceplot5"),
                                       plotOutput("colchoiceplot6")
                                     )
                                     )
                                     )
                                   ),
                          tabPanel("Group Plots",
                                   fluidRow(
                                   splitLayout(
                                     verticalLayout(
                                       plotOutput("Q1"),
                                       plotOutput("Q2"),
                                       plotOutput("Q3"),
                                       plotOutput("Q4"),
                                       plotOutput("Q5"),
                                       plotOutput("Q6")
                                     ),
                                     verticalLayout(
                                       plotOutput("h1"),
                                       plotOutput("h2"),
                                       plotOutput("h3"),
                                       plotOutput("h4"),
                                       plotOutput("h5"),
                                       plotOutput("h6")
                                     )
                                    )
                                   )
                                   ),
                          tabPanel("Plot Spots over Length/Width Axis",
                                   fluidRow(
                                     splitLayout(
                                       verticalLayout(
                                        plotOutput("Lplot"),
                                        plotOutput("Wplot")
                                       )
                                   )
                                   )
                                   
                          )
                              
                                  
                        )
                      )
                    )
                    ),
           
           tabPanel("Plot Dimensions",
                    sidebarLayout(
                      sidebarPanel(
                        checkboxGroupInput("Vars",
                          "Pick 1 or 2 variables you want to plot",
                                    choices = c("Cell Length"=1, "Cell Width" = 2, "Cell Area"= 3,
                                    "Spot location in length axis" = 4, "Spot location on width axis" = 5),
                                    selected = 1),
                        selectInput("Cols3",
                                    "Pick color for 1th variable",
                                    choices = c("Grey"=1, "Ocra"=2, "Light Blue"=3, "Green"=4, "Yellow"=5,
                                                "Black" = 6, "Dark Blue"=7, "Orange"=8, "Pink"=9),
                                    selected = 1),
                        selectInput("Cols4",
                                    "Pick color for 2nd variable",
                                    choices = c("Grey"=1, "Ocra"=2, "Light Blue"=3, "Green"=4, "Yellow"=5,
                                                "Black" = 6, "Dark Blue"=7, "Orange"=8, "Pink"=9),
                                    selected = 1),
                        sliderInput("bins", "Select bin size for variable 1",
                                    min=0, max=50, value=30),
                        sliderInput("bins2", "Select bin size for variable 2",
                                    min=0, max=50, value=30),
                        radioButtons("Typep", "Pick plot type (single plot)",
                                     choices = c("Histogram"=1, "Density"=2, "Dotplot"=3, "Violin"=4),
                                     selected = 2),
                        downloadButton("Downloaddimensions", "Download the plot(s) here")
                        
                      ),
                      mainPanel(
                        fluidPage(
                          verticalLayout(
                            splitLayout(
                            plotOutput("Plotvar1"),
                            plotOutput("Plotvar2")
                            ),
                            plotOutput("Combpplot")
                          )
                          
                        )
                        
                      )
                    )
           )
           
           
)







