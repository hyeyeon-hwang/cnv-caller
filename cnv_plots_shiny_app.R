# install.packages("shiny")
library(shiny)
library(scatterD3)

library(magrittr) # pipe
library(tibble) # as_tibble(), add_column()
library(dplyr)
library(tidyr) # unite()

# Shiny tutorial pdf
# https://ibiostat.be/seminar/uploads/introdcution-r-shiny-package-20160330.pdf

# Tutorial app 9: users can upload their own data file 

# scatterD3 shiny app source code
# https://github.com/juba/scatterD3_shiny_app/blob/master/server.R

# Data
readData <- function(inputData, sampleNames) {
  data <- read.delim(inputData, header = TRUE, col.names = sampleNames) %>%
    #tibble::as_tibble() %>%
    tibble::add_column(chrm = .$chr, .before = 1) %>%
    tidyr::unite("pos", chr:start, sep = ":") %>%
    tidyr::unite("pos", pos:end, sep = "-")
  return(data)
}
data <- readData(
  inputData = "Dec06_cnv_output.txt",
  sampleNames = c("chr", "start", "end",
                  "cJLKD", "c6978", "c6980", "c7015", "c7016",
                  "e7005", "e7006", "e7007", "e7008"))


ui <- fluidPage(
  titlePanel("CNV caller visualization"),
  sidebarLayout(
    sidebarPanel(
      selectInput("scatterD3_chrm", 
                  label = "Chromosome position (x variable)", 
                  choices = unique(data$chrm), 
                  selected = "chr21"),
      helpText("Format of the x variable in the plot is position chr:start-end"),
      selectInput("scatterD3_y", 
                  label = "CNV of sample (y variable)", 
                  choices = names(data), 
                  selected = "cJLKD")
    ),
    mainPanel(scatterD3Output("scatterPlot"))
  )
)


# Shiny app server
server <- function(input, output) {
  # dataD3 <- reactive({
  #   data[which(data$chrm == input$scatterD3_chrm), ]
  # })

  data_x <- reactive({
    data$pos[which(data$chrm == input$scatterD3_chrm)]
  })
  data_y <- reactive({
    data[which(data$chrm == input$scatterD3_chrm), input$scatterD3_y]
  })
  
  output$scatterPlot <- renderScatterD3({
    scatterD3(x = data_x(),
              y = data_y(),
              xlab = "Positions of chromosome",
              ylab = "Copy number estimation of sample",
              point_size = 18, point_opacity = 0.6,
              hover_size = 20, hover_opacity = 1,
              lines = data.frame(slope = c(0, 0), 
                                 intercept = c(2, 3), 
                                 stroke = c("green", "red"),
                                 stroke_width = c(2, 2)),
              ylim = c(0, 5))
  })
  # default ylim = c(0,5), user can zoom out to see full graph with outliers
  
  #scatterD3(x = data[which(data$chrm == input$scatterD3_chrm), input$scatterD3_x],#data_x(), #data()[,input$scatterD3_x],#input$xcol, 
  #          y = data[which(data$chrm == input$scatterD3_chrm), input$scatterD3_y])#data_y())#data()[,input$scatterD3_y])#input$ycol)
}

shinyApp(ui = ui, server = server)

