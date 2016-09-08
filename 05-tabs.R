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
               selectInput("bStacked","",choices = c("Stacked histogram","Overlapped histogram")),
               selectInput("bDensity","",choices = c("Percentage", "Count"))
        )
      ),
br(),
p("This is the end.")
    ),
    tabPanel(title = "Use the Mixture Model",
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
  
  ####the functions for plotting####
  senmuR <- reactive({log(input$senmu)})
  sensdR <- reactive({log(input$sensd)})
  resmuR <- reactive({log(input$resmu)})
  ressdR <- reactive({log(input$ressd)})
  
  sen_popR <- reactive({rlnorm(input$nn*(1-input$prop_resist),senmuR(),sensdR())})
  res_popR <- reactive({rlnorm(input$nn*input$prop_resist,resmuR(),ressdR())})
  
  
  genData <- reactive({
    #nn <- input$nn
    senmu <- senmuR()
    sensd <- sensdR()
    
    #prop_resist <- input$prop_resist
    resmu <- resmuR()
    ressd <- ressdR()
    
    sen_pop <- sen_popR() #sensitive population
    res_pop <- res_popR() #resistant population
    c(sen_pop,res_pop)
    #total_pop <- c(sen_pop,res_pop)
  })
  
  genData.DF <- reactive({
    cbind(genData(),c(rep(0,length(sen_popR())),rep(1,length(res_popR()))))
  })
  output$ROC <- renderPlot({
    popDF <- genData.DF()
    
    TPR <- sum(res_popR()>=input$cutoff)/length(res_popR())
    FPR <- sum(sen_popR()>=input$cutoff)/length(sen_popR())
    
    true_res <- sum(res_popR()>=input$cutoff)
    fal_res <- sum(sen_popR()>=input$cutoff)
    fal_sen <- sum(res_popR()<input$cutoff)
    true_sen <- sum(sen_popR()<input$cutoff)
    overlay <- paste(true_res," truly resistant, ", fal_res, " falsely resistant \n",fal_sen," falsely sensitive, ",true_sen," truly sensitive")
    
    
    roc(popDF[,2], popDF[,1],  partial.auc.correct=TRUE, partial.auc.focus="sens",ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, show.thres=TRUE, main="Receiver Operating Characteristic (ROC) Curve")
    points((1-FPR),TPR, col="red", pch=19)
    text(.5,.5,overlay, col="red")
  })
  
  output$densityplot <- renderPlot({
    popDF2 <- genData.DF() 
    
    #mxmdl <- normalmixEM(popDF2)
    #plot(mxmdl, which=2)
    popDF2[popDF2[,2]==0,2] <- "Sensitive"
    popDF2[popDF2[,2]==1,2] <- "Resistant"
    
    popDF2 <- as.data.frame(popDF2)
    popDF2[,1] <- as.numeric(as.character(popDF2[,1]))
    names(popDF2) <- c("Half-life (hours)","Sensitivity")
    
    if(input$bDensity == "Percentage"){
      if(input$bStacked=="Stacked histogram"){
        ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour= Sensitivity)) + theme_bw() +
          geom_histogram(aes(y=(..count..)/sum(..count..)), alpha=.8, position="stack",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
          geom_vline(xintercept= input$cutoff, colour="red", size=1) + ylab("Percent") + ggtitle("Stacked Histogram of Simulated Half-Lives")+
          theme(plot.title= element_text(face="bold")) +
          scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
      }
      
      else {
        ggplot(popDF2, aes(x=`Half-life (hours)`)) + theme_bw() +
          geom_histogram(aes(y=(..count..)/sum(..count..), fill=Sensitivity, colour= Sensitivity), alpha=.4, position="identity",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
          geom_vline(xintercept= input$cutoff, colour="red", size=1) + ylab("Percent") + ggtitle("Histogram of Simulated Half-Lives")+
          theme(plot.title= element_text(face="bold")) +
          scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
      }
    }
    
    else{
      if(input$bStacked =="Stacked histogram"){
        ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour= Sensitivity)) + theme_bw() +
          geom_histogram(alpha=.8, position="stack",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
          geom_vline(xintercept= input$cutoff, colour="red", size=1) + ggtitle("Stacked Histogram of Simulated Half-Lives")+
          theme(plot.title= element_text(face="bold")) +
          scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
      }
      
      else {
        ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour= Sensitivity)) + theme_bw() +
          geom_histogram( alpha=.4, position="identity",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
          geom_vline(xintercept= input$cutoff, colour="red", size=1) + ggtitle("Histogram of Simulated Half-Lives")+
          theme(plot.title= element_text(face="bold")) +
          scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
      }
    }
    
    
    
  })
}

shinyApp(server = server, ui = ui)