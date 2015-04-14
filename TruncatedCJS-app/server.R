# server.R

shinyServer(function(input, output) {

  ## Load package from library
  library(TruncatedCJS)
  logit <- make.link("logit")
  
  get.dataset <- reactive({
  
    inFile <- input$dat
    if (is.null(inFile))
      return(NULL)
    dat <- read.csv(inFile$datapath, header=TRUE)
    return(dat)
  })
  
  output$contents <- renderTable({
    
    dat <- get.dataset()
    if(is.null(dat)){
      return(NULL)
    }

    return(head(dat))
  })
  
  output$k.slider <- renderUI({
    dat <- get.dataset()
    
    if(is.null(dat)){return(NULL)}
    
    sliderInput("k", "Select value", 1, ncol(dat), 2, step = 1)
  })
  
  get.result <- reactive({
    
    dat <- get.dataset()
    T <- ncol(dat)
    if(is.null(dat)){
      return(NULL)
    }
    
    fit <- fitCJS(dat,T,input$k,control=list(trace=1,maxit=10000),hess=TRUE,truth=NULL,debug=FALSE)
    return(fit)
  })
  
  output$survival.fig <- renderPlot({
    
    dat <- get.dataset()
    T <- ncol(dat)
    fit <- get.result()
    
    if(is.null(dat)){
      return(NULL)
    }
    
    plot(NA, NA, ylim=c(0,1), xlim=c(0,T), xlab="Capture Occasion", ylab="Survival Probability", main='Survival Probabilities over Time')
    points(1:(T-1), fit$phi[,1], pch=16)
    for(t in 1:(T-1)){
      lines(rep(t,2), fit$phi[t,2:3])
    }

  })
  
  output$capture.fig <- renderPlot({
    
    dat <- get.dataset()
    T <- ncol(dat)
    fit <- get.result()
    
    if(is.null(dat)){
      return(NULL)
    }
    
    plot(NA, NA, ylim=c(0,1), xlim=c(0,T), xlab="Capture Occasion", ylab="Capture Probability", main='Capture Probabilities over Time')
    
    points(1:(T-1), fit$p[,1], pch=16)
    
    for(t in 1:(T-1)){
      lines(rep(t,2), fit$p[t,2:3])
      
    }
    
  })
  
  output$parameter.estimates <- renderTable({
    
    dat <- get.dataset()
    T <- ncol(dat)
    fit <- get.result()
    
    if(is.null(dat)){
      return(NULL)
    }
    
    result.table <- data.frame(cbind(1:T,rbind(NA,fit$p),rbind(fit$phi,NA)))
    names(result.table) <- c('index','p.hat','2.5% p.hat', '97.5% p.hat', 'phi.hat','2.5% phi.hat', '97.5% phi.hat')
    result.table[,'index'] <- as.integer(result.table[,'index'])
    
    return(result.table)
  }, include.rownames=FALSE)
  
})