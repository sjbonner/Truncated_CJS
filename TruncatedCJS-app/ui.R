# ui.R

shinyUI(fluidPage(
  titlePanel("Truncated CJS Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Upload Data Set:"),
      fileInput("dat", "Select File", accept=c('.csv')),
      helpText("Enter number of occasions to look ahead:"),
      uiOutput('k.slider')
    ),
    
    mainPanel(p("First 6 rows of data Set:"),
              tableOutput('contents'),
              plotOutput('survival.fig'),
              plotOutput('capture.fig'),
              p("Table of Parameter Estimates:"),
              tableOutput('parameter.estimates'))
  )

  
))