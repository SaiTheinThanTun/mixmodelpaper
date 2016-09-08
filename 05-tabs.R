# 05-tabs.R

library(shiny)

ui <- fluidPage(
  h2("Identify  artemisinin resistance from parasite clearance half-life data"),
  #title = "Random generator",
  tabsetPanel(              
    tabPanel(title = "Introduction",
             p("Assuming that the parasite clearance half-lives are in log-normal distribution, and that the values for 
sensitive and resistant populations each assume unimodal distribution. When the sensitive population has 
a geomatric half-life mean of : ", 
strong(textOutput("senmuO", inline=TRUE)), 
"and the standard deviation of ",
strong(textOutput("sensdO",inline=T)),
"and the resistant population has a geomatric half-life mean of : ",
strong(textOutput("resmuO",inline=T)),
"and the standard deviation of ",
strong(textOutput("ressdO",inline=T)),
"then the distribution will look like the following."
               ),
      # plotOutput("norm"),
      # actionButton("renorm", "Resample"),
      fluidRow(
        column(5,
               plotOutput(outputId = "densityplot")
        ),
        column(5,
               plotOutput(outputId = "ROC")
        )
        
      ),
      fluidRow(
        column(4,
               h3("Sensitive Distribution"),
               sliderInput(inputId = "senmu",
                           label = "Mean half-life",
                           value = 3, min = 1, max = 6.5, step = .5
               ),
               sliderInput(inputId = "sensd",
                           label = "SD",
                           value = 1.45, min = 1, max = 2.1
               )
        ),
        column(4,
               h3("Resistant Distribution"),
               sliderInput(inputId = "resmu",
                           label = "Mean half-life",
                           value = 6.5, min = 5, max = 10, step = .5
               ),
               sliderInput(inputId = "ressd",
                           label = "SD",
                           value = 1.22, min = 1, max = 2.1
               ),
               sliderInput(inputId = "prop_resist",
                           label = "Proportion resistant",
                           value = .1, min = 0, max = 1
               )
        ),
        column(3,
               numericInput(inputId = "nn",
                            label = "Sample Size:",
                            value = 200
               ),
               
               
               sliderInput(inputId = "cutoff",
                           label = "Cut-off half-life value",
                           value = 5, min = 0, max = 10, step=.5
               ),
               h4("Histogram options"),
               checkboxInput(inputId = "bStacked",
                             label = "Stacked histogram",
                             value = TRUE
               ),
               checkboxInput(inputId = "bDensity",
                             label = "Density",
                             value = TRUE
               )
        )
        
      )
    ),
    tabPanel(title = "Uniform data",
      plotOutput("unif"),
      actionButton("reunif", "Resample")
    ),
    tabPanel(title = "Chi Squared data",
      plotOutput("chisq"),
      actionButton("rechisq", "Resample")
    )
  )
)

server <- function(input, output) {
  
  rv <- reactiveValues(
    norm = rnorm(500), 
    unif = runif(500),
    chisq = rchisq(500, 2))
  
  observeEvent(input$renorm, { rv$norm <- rnorm(input$bb) })
  observeEvent(input$reunif, { rv$unif <- runif(500) })
  observeEvent(input$rechisq, { rv$chisq <- rchisq(500, 2) })
  
  output$norm <- renderPlot({
    hist(rv$norm, breaks = 30, col = "grey", border = "white",
      main = "500 random draws from a standard normal distribution")
  })
  output$unif <- renderPlot({
    hist(rv$unif, breaks = 30, col = "grey", border = "white",
      main = "500 random draws from a standard uniform distribution")
  })
  output$chisq <- renderPlot({
    hist(rv$chisq, breaks = 30, col = "grey", border = "white",
       main = "500 random draws from a Chi Square distribution with two degree of freedom")
  })

#having the parameters reflect on the text description
  output$senmuO <- renderText({
    as.character(input$senmu)
  })
  output$sensdO <- renderText({
    as.character(input$sensd)
  })
  output$resmuO <- renderText({
    as.character(input$resmu)
  })
  output$ressdO <- renderText({
    as.character(input$ressd)
  })
}

shinyApp(server = server, ui = ui)