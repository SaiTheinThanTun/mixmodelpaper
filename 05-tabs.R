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
"then the distribution will look like the following. The proportion of 
resistant population is",
strong(textOutput("prop_resistO",inline=T)),
". The sample size is",
strong(textOutput("nnO",inline=T)),
". The cutoff value is ",
strong(textOutput("cutoffO",inline=T))
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
p("This is the end."),

###test####
p(textOutput("genDataOut")),

############################
###Portions from MixModel###
############################

fileInput(inputId = "file", label = "Select your input file: (simulated_cloneHLdata_SMRUbyyear.csv in this case)"),
plotOutput(outputId = "histoplot2")
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
  output$prop_resistO <- renderText({
    as.character(input$prop_resist)
  })
  output$nnO <- renderText({
    as.character(input$nn)
  })
  output$cutoffO <- renderText({
    as.character(input$cutoff)
  })
  
  ####the functions for plotting####
  # senmuR <- reactive({log(input$senmu)})
  # sensdR <- reactive({log(input$sensd)})
  # resmuR <- reactive({log(input$resmu)})
  # ressdR <- reactive({log(input$ressd)})
  # 
  # sen_popR <- reactive({rlnorm(input$nn*(1-input$prop_resist),senmuR(),sensdR())})
  # res_popR <- reactive({rlnorm(input$nn*input$prop_resist,resmuR(),ressdR())})
  # 
  # 
  # genData <- reactive({
  #   #nn <- input$nn
  #   senmu <- senmuR()
  #   sensd <- sensdR()
  # 
  #   #prop_resist <- input$prop_resist
  #   resmu <- resmuR()
  #   ressd <- ressdR()
  # 
  #   sen_pop <- sen_popR() #sensitive population
  #   res_pop <- res_popR() #resistant population
  #   c(sen_pop,res_pop)
  #   #total_pop <- c(sen_pop,res_pop)
  # })
  # 
  # genData.DF <- reactive({
  #   cbind(genData(),c(rep(0,length(sen_popR())),rep(1,length(res_popR()))))
  # })
  # output$ROC <- renderPlot({
  #   popDF <- genData.DF()
  # 
  #   TPR <- sum(res_popR()>=input$cutoff)/length(res_popR())
  #   FPR <- sum(sen_popR()>=input$cutoff)/length(sen_popR())
  # 
  #   true_res <- sum(res_popR()>=input$cutoff)
  #   fal_res <- sum(sen_popR()>=input$cutoff)
  #   fal_sen <- sum(res_popR()<input$cutoff)
  #   true_sen <- sum(sen_popR()<input$cutoff)
  #   overlay <- paste(true_res," truly resistant, ", fal_res, " falsely resistant \n",fal_sen," falsely sensitive, ",true_sen," truly sensitive")
  # 
  # 
  #   roc(popDF[,2], popDF[,1],  partial.auc.correct=TRUE, partial.auc.focus="sens",ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, show.thres=TRUE, main="Receiver Operating Characteristic (ROC) Curve")
  #   points((1-FPR),TPR, col="red", pch=19)
  #   text(.5,.5,overlay, col="red")
  # })
  # 
  # output$densityplot <- renderPlot({
  #   popDF2 <- genData.DF()
  # 
  #   #mxmdl <- normalmixEM(popDF2)
  #   #plot(mxmdl, which=2)
  #   popDF2[popDF2[,2]==0,2] <- "Sensitive"
  #   popDF2[popDF2[,2]==1,2] <- "Resistant"
  # 
  #   popDF2 <- as.data.frame(popDF2)
  #   popDF2[,1] <- as.numeric(as.character(popDF2[,1]))
  #   names(popDF2) <- c("Half-life (hours)","Sensitivity")
  # 
  #   if(input$bDensity == "Percentage"){
  #     if(input$bStacked=="Stacked histogram"){
  #       ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour= Sensitivity)) + theme_bw() +
  #         geom_histogram(aes(y=(..count..)/sum(..count..)), alpha=.8, position="stack",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
  #         geom_vline(xintercept= input$cutoff, colour="red", size=1) + ylab("Percent") + ggtitle("Stacked Histogram of Simulated Half-Lives")+
  #         theme(plot.title= element_text(face="bold")) +
  #         scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
  #     }
  # 
  #     else {
  #       ggplot(popDF2, aes(x=`Half-life (hours)`)) + theme_bw() +
  #         geom_histogram(aes(y=(..count..)/sum(..count..), fill=Sensitivity, colour= Sensitivity), alpha=.4, position="identity",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
  #         geom_vline(xintercept= input$cutoff, colour="red", size=1) + ylab("Percent") + ggtitle("Histogram of Simulated Half-Lives")+
  #         theme(plot.title= element_text(face="bold")) +
  #         scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
  #     }
  #   }
  # 
  #   else{
  #     if(input$bStacked =="Stacked histogram"){
  #       ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour= Sensitivity)) + theme_bw() +
  #         geom_histogram(alpha=.8, position="stack",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
  #         geom_vline(xintercept= input$cutoff, colour="red", size=1) + ggtitle("Stacked Histogram of Simulated Half-Lives")+
  #         theme(plot.title= element_text(face="bold")) +
  #         scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
  #     }
  # 
  #     else {
  #       ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour= Sensitivity)) + theme_bw() +
  #         geom_histogram( alpha=.4, position="identity",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
  #         geom_vline(xintercept= input$cutoff, colour="red", size=1) + ggtitle("Histogram of Simulated Half-Lives")+
  #         theme(plot.title= element_text(face="bold")) +
  #         scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
  #     }
  #   }
  # 
  # 
  # 
  # })
  
  #######################################################################
  ###For the users data, run the mixture model and draw the histogram####
  #######################################################################
  output$histoplot2 <- renderPlot({
    inFile <- input$file
    mixdat <- read.csv(inFile$datapath)
    
    N<-ncol(mixdat)
    M <- 5
    
    pval<-0.1
    nboot<-100 # number of iterations for bootstrap
    nsim<-1000 # number of iterations for creating probaility of resistance vs HL graphs
    P<-2 # use P or more samples to get geometric means and discard all other samples for permutation analysis
    T<-100 # number of permutations for permutation analysis
    smax=5000
    
    # create output matrices
    output.mu <- matrix(NA,nrow=M,ncol=N)
    output.sigma <- matrix(NA,nrow=M,ncol=N)
    output.lambda <- matrix(NA,nrow=M,ncol=N)
    output.loglik <- matrix(NA,nrow=M,ncol=N)
    output.mu.se <- matrix(NA,nrow=M,ncol=N)
    output.sigma.se <- matrix(NA,nrow=M,ncol=N)
    output.lambda.se <- matrix(NA,nrow=M,ncol=N)
    AIC<-matrix(0,nrow=M,ncol=N)
    AICdelta<-matrix(0,nrow=M,ncol=N)
    
    nb<-na.omit(mixdat[,N])
    
    # fit single component model
    for (i in 1:N){
      # 1 COMPONENT LOG NORMAL
      nmixdat<-na.omit(mixdat[,i])
      lmixdat<- log(nmixdat)
      xll<-fitdistr(lmixdat,"normal")
      output.loglik[1,i]<- xll$loglik
      output.mu[1,i]<-xll$estimate[1]
      output.lambda[1,i]<-1
      output.sigma[1,i]<-xll$estimate[2]
      output.mu.se[1,i]<-xll$sd[1]
      output.sigma.se[1,i]<-xll$sd[2]
      output.lambda.se[1,i]<-0
      AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
      AICdelta[1,i]<-0
    }
    
    # fit multiple component models sequentially
    for (i in 1:N){
      nmixdat<-na.omit(mixdat[,i])
      lmixdat<- log(nmixdat)
      # >=2 COMPONENTS LOG NORMAL
      j<-1
      # stop if j-component model is more parsimonious than (j-1)-compnent model
      while((j<=M-1) && AICdelta[j,i]<=pval){
        j<-j+1
        res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
        resboot <- boot.se(res, B = nboot)
        resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]	
        output.loglik[j,i]<-res$loglik
        AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
        AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
        if(AICdelta[j,i]<=pval){
          output.mu[1:j,i]<-res$mu
          output.sigma[1:j,i]<-res$sigma
          output.lambda[1:j,i]<-res$lambda
          output.mu.se[1:j,i]<-resboot$mu.se
          output.sigma.se[1:j,i]<-resboot$sigma.se
          output.lambda.se[1:j,i]<-resboot$lambda.se		
        }
      }
      
      #Make the table containing the probabilities of each patient belonging to each component distribution. 
      
      
      
      #Calculate the probabilities of an individual patient, typed into the box by the user, belong to each component distribution.
      # Print the results. 
      output$explanation1 <- renderText({"Below are two graphs. The graph on the left
        represents aggregate data from White et al. 2015. The graph on the right is a graph made from your data."})
      
      output$explanation2 <- renderText({"The graph on th left depicts two half life distributions
        with geometric means SOMETHING AND SOMETHING ELSE respectively. 
        The distribution with a geometric mean half life of SOMETHING was intepreted
        as representing patients with parasites sensitive to artemisinin. The distribution
        with a geometric mean half life of SOMETHING ELSE was interpreted as representing
        patients with parasites resistant artemisinin. With this information, you may be able
        to interpret the graph on the rigth, which represents your own data."})
      
      output$explanation3 <- renderText({"Below are some statistics from 
        from the graph representing your data. There is also a list of the probabilities of 
        each patient belonging to each of the component distributions depicted in the graph."})
      
      output$geometric_means_and_proportions <- renderPrint({
        j <- length(na.omit(output.mu))
        cat("The model predicts ", j, " component geometric mean half lives (hours):",
            
            "\n\n")
        
        for (a in 1:j) {
          cat("Distribution",a,"\n",
              
              "Geometric mean = ", exp(output.mu[a]),
              
              "\n", 
              
              "Contribution to composite distribution = ", output.lambda[a],
              
              "\n\n"
              
          )       
        }
      })
      
      
      
      }
    
    for (ds in 1:N){
      nmixdat<-na.omit(mixdat[,ds])
      plam<-na.omit(output.lambda[,ds])
      pmu<-na.omit(output.mu[,ds])
      psig<-na.omit(output.sigma[,ds])
      hist(nmixdat,freq=FALSE,main = paste("Distribution of parasite clearance half lives","\n", "from your data"),xlab = "Clearance half-life (hours)",ylim=c(0,0.6),col="grey",lwd=2,ps=20) #taken out for shiny #,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12)
      x <- seq(0.1, max(nmixdat), length=1000)
      #hx<-plam[1]*dlnorm(x,meanlog=(pmu[1]),sdlog=psig[1])
      hx <- list()
      
      #casting multiple lines for different distributions
      lcolor <- c('blue','red','brown','green','yellow')
      for(k in 1:length(pmu)){
        hx[[k]]<-plam[k]*dlnorm(x,meanlog=(pmu[k]),sdlog=psig[k])
        lines(x,hx[[k]],col=lcolor[k], lwd=5)
      }
    
    
      
      # if(length(plam)>1){
      #   for(k in 2:length(plam)){
      #     hx<-hx+plam[k]*dlnorm(x,meanlog=(pmu[k]),sdlog=psig[k])
      #   }
      # }
      # lines(x,hx,col="red", lwd=5)
    }
    
    means <- reactive(na.omit(output.mu))
    spreads <- reactive(na.omit(output.sigma))
    proportions <- reactive(na.omit(output.lambda))
    nn.mixdat <- reactive(length(nb))#(nrow(mixdat))
    
    #################################
    ####plots from the cutoff app####
    ####the functions for plotting###
    #################################
    
    senmuR <- reactive(means()[1])
    sensdR <- reactive(spreads()[1])    
    resmuR <- reactive(means()[2])
    ressdR <- reactive(spreads()[2])

    ###test####
    output$genDataOut <- renderPrint({as.character(c(nn.mixdat(),proportions(),senmuR(),sensdR()))})
    output$densityLine <- renderPlot({
      
    })

    sen_popR <- reactive({rlnorm(nn.mixdat()*(1-proportions()[2]),senmuR(),sensdR())})
    res_popR <- reactive({rlnorm(nn.mixdat()*proportions()[2],resmuR(),ressdR())})


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
    
    
  })
}

shinyApp(server = server, ui = ui)