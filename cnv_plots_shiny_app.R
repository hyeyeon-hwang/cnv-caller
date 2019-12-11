library(shiny)
library(scatterD3)

library(magrittr) # pipe
library(tibble) # as_tibble(), add_column()
library(dplyr)
library(tidyr) # unite()
library(glue) # glue() 

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
  inputData = "Dec11_cnv_bins_full.txt",
  sampleNames = c("chr", "start", "end",
                  "cJLKD", "c6978", "c6980", "c7015", "c7016",
                  "e7005", "e7006", "e7007", "e7008"))


ui <- fluidPage(
  titlePanel("Copy number variation (CNV) caller visualization"),
  sidebarLayout(
    sidebarPanel(
      h4(strong("x variable:")),
      h4("Chromosome bin"),
      selectInput("scatterD3_chrm", 
                  label = "Select chromosome.",
                  choices = unique(data$chrm), 
                  selected = "chr21"),
      helpText("Format of the x variable in the plot is position chr:start-end"),
      br(),
      h4(strong("y variable:")),
      h4("Copy number of chromosome bin for sample"),
      selectInput("scatterD3_y", 
                  label = "Select sample.", 
                  choices = names(data)[-c(1:2)], 
                  selected = names(data)[3]),
      checkboxInput("checkbox_yrange", "Set y variable range to [0, 6]", value = FALSE)
    ),
    mainPanel(
      span(textOutput("plotTitle"), style = "font-weight: bold"),
      scatterD3Output("scatterPlot"))
  )
)


# Shiny app server
server <- function(input, output) {
  data_x <- reactive({
    data$pos[which(data$chrm == input$scatterD3_chrm)]
  })
  data_y <- reactive({
    data[which(data$chrm == input$scatterD3_chrm), input$scatterD3_y]
  })
  
  output$plotTitle <- renderText({
    paste("Copy number plot of bins in chromosome ", input$scatterD3_chrm, "and sample ", input$scatterD3_y)
  })
  
  output$scatterPlot <- renderScatterD3({
    if (input$checkbox_yrange == TRUE) {
      ylim_var <- c(0,6)
    } else {
      ylim_var <- NULL
    }
    
    scatterD3(x = data_x(),
              y = data_y(),
              xlab = glue::glue("Bin in chromosome {input$scatterD3_chrm}"),
              ylab = glue::glue("Copy number of chromosome bin in sample {input$scatterD3_y}"),
              point_size = 18, point_opacity = 0.6,
              hover_size = 20, hover_opacity = 1,
              lines = data.frame(slope = c(0, 0), 
                                 intercept = c(2, 3), 
                                 stroke = c("black", "red"),
                                 stroke_width = c(2, 2)),
              caption = "This is the caption for scatterD3",
              ylim = ylim_var
              )
  })

}

shinyApp(ui = ui, server = server)

