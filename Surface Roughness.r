library(shiny)
library(plotly)
library(shinythemes)


# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(theme = shinytheme("flatly"),

  navbarPage("",
   
  tabPanel("Application",
   # Application title
   titlePanel(h1("ISU Surface Roughness Measurement for Metal Castings",align = "center")),

     sidebarPanel(
       fileInput('pointCloud', 'Please upload your point cloud file as a .txt or .csv',
                 accept=c('text/csv', 
                          'text/comma-separated-values,text/plain', 
                          '.csv')),
       helpText(tags$b("Check box if your point cloud file has a header.")),
       checkboxInput('header', 'Header', FALSE),
       radioButtons("sep", "Separator",
                    choices = c(Space = "",
                                Comma = ","
                                ),
                    selected = "")
       #helpText(tags$b("Check box to see point cloud after applying PCA")),
       #checkboxInput("graph","Display Plot",FALSE)
       ),

       fluidRow(
       column(3,
       tableOutput('contents')
       ),
       column(3,
       tableOutput('Plates'))
       
       

     ),
        verbatimTextOutput("points"),
        verbatimTextOutput("roughness"),
   
  plotlyOutput("plott"),

   tags$footer("R Shiny application created by:", tags$b("Hieu Pham"),tags$br() ,
              "Department of Industrial and Manufacturing Systems Engineering",tags$br(),
              "Iowa State University",
              align = "left", style = "
              bottom:0;
              width:100%;
              height:0px;   /* Height of the footer */
              color: black;
              padding: 10px;
              z-index: 1000;")
   
),
tabPanel("Instructions",
         #titlePanel(h3("Recommended setup conditions for optimal surface roughness results",align = "left")),
         titlePanel(h2("Recommended Setup Conditions for Optimal Surface Roughness Results",align = "center")),
         
         
         
         
         p(strong("1)"), "Ensure to perform the calibration test of the laser scanner as per the recommended manufacturer settings"),
         
         p(strong("2)"), "Ensure to perform scans in a controlled environment free from any direct sunlight"),

         p(strong("3)"), "Ensure to perform all scans at an illumination angle that is perpendicular to the casting surface in order to maintain a constant beam spot diameter "),
         
         p(strong("4)"), "Ensure to perform the scan at a standard field of view as per the manufacturer's settings "),
                  
         p(strong("5)"), "On the metal casting surface, sample a planar surface area of 1.5\'' x  1.5\'' free from undercuts, surface bumps and surface curvature "),
                  
         p(strong("6)"), "Apply a 3D scan spray or developer-based spray on the chosen 1.5\'' x 1.5\'' area. These are temporary and water soluble in nature "),
                  
         p(strong("7)"), "Mask the surrounding areas of the chosen area with black masking tape to allow low light transmittance and blackout the incident light outside
                  the area of interest "),
                  
         p(strong("8)"), "It is recommended to achieve 5000 points per 0.75 square inch area "),
         
         img(src='img1.jpg', height=275,width= 275),
         img(src='img2.png', height=275,width= 275),
         img(src='img3.png', height=275,width= 275)
         
         )

)))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  options(shiny.maxRequestSize=30*1024^2) 
  
  
  output$contents <- renderTable({
    data <- input$pointCloud
    if (is.null(data))
      return(NULL)
    data <- read.csv(data$datapath, header=input$header,sep = input$sep)
    #data = data[,c(1,3,2)]
    if(input$header == F){
      colnames(data) = c("x","y","z")
    }
    head(data,10)
  },caption <- "<b> <span style='color:#000000'>First 10 Rows of Point Cloud",
  caption.placement <- getOption("xtable.caption.placement", "top"))
  
  
  output$Plates <-renderTable({
    
    data <- data.frame(matrix(c(599,742,944,1003),ncol = 1))
    data$Plate <- c("A1","A2","A3","A4")
    data <- data[,c(2,1)]
    colnames(data) <- c("Plates","Avg. Sq Value (micro in.)")
    data
  },caption <- "<b> <span style='color:#000000'>ASTM A802 Surface Texture Plates",
  caption.placement <- getOption("xtable.caption.placement", "top"))
  
  #POINTS
  output$points <- renderText({
    data <- input$pointCloud
    if (is.null(data))
      return(NULL)
    data <- read.csv(data$datapath, header=input$header,sep = input$sep)
    data <- data[,c(1,3,2)]
    
    paste("The point cloud has", nrow(data),"points.")
  
  })
  
  
  #PCA PART 
  output$roughness <- renderText({ 
  data <- input$pointCloud
  if (is.null(data))
    return(NULL)
  data <- read.csv(data$datapath, header=input$header,sep = input$sep)
  data <- data[,c(1,3,2)]
  zCords <- round(prcomp(data)$x[,3],10)
  zzz <- quantile(zCords)
  zz <- mean(zCords) - as.numeric(zzz[2]) # ground line i.e. mean - medium which is literally zero + epsilon
  zz <- as.numeric(zz)
  below <- sample(which(zCords <zz),length(which(zCords < zz)))
  above <- sample(which(zCords >zz),length(which(zCords > zz)))
  combine <- matrix(c(zCords[below],zCords[above]),ncol = 1)
  combine <- as.matrix(combine)
  zCords <- combine - zz
  
  #Compute roughness for Ra based on distance from point to plane
  roughRa <- sum(abs(zCords))/nrow(data)
  
  #Compute roughness for Rq based on distance from point to plane
  roughRqALL <- sqrt(mean(zCords^2))
  
  #Convert to microns. This is still PCA values
  RoughRA <- (1/2.54e-5)*roughRa
  RoughRQ <- (1/2.54e-5)*roughRqALL
  
    if(RoughRQ <= exp(3386.3/501.4)){
    paste(" Warning: Measurement may be inaccurate since roughness is lesser than the minimum threshold of 200 micro in.","\n",
          "Note: Sand casting typically produces an average root mean square (RMS) value of 250-900 micro in. ", "\n",
          "The calculated Sq value is",0, "micro in.", "\n",
            "The calculated Sa value is",0, "micro in.", "\n",
            "For reference purpose only, the surface texture is approximately equivalent to or better than ASTM A802 A1-(Level I)")
    }else if(exp(3386.3/501.4)<= RoughRQ && RoughRQ<=exp(3589/501.4)){
    paste(" Warning: Measurement may be inaccurate since roughness is lesser than the minimum threshold of 200 micro in.","\n", "Note: Sand casting typically produces an average root mean square (RMS) value of 250-900 micro in.", "\n",
      "The calculated Sq value is",round(501.4*log(RoughRQ) - 3386.3,0), "micro in. (",round((501.4*log(RoughRQ) - 3386.3)*.025399999,1),"micro m.)", "\n",
          "The calculated Sa value is",round(501.4*log(RoughRA) - 3386.3,0), "micro in. (",round((501.4*log(RoughRA) - 3386.3)*.025399999,1),"micro m.)", "\n",
          "For reference purpose only, the surface texture is approximately equivalent to or better than ASTM A802 A1-(Level I)")
  }else if(exp(3586.3/501.4)<= RoughRQ && RoughRQ <= 3259.8285){
    paste(" The calculated Sq value is",round(501.4*log(RoughRQ) - 3386.3,0), "micro in. (",round((501.4*log(RoughRQ) - 3386.3)*.025399999,1),"micro m.)", "\n",
          "The calculated Sa value is",round(501.4*log(RoughRA) - 3386.3,0), "micro in. (",round((501.4*log(RoughRA) - 3386.3)*.025399999,1),"micro m.)", "\n",
          "For reference purpose only, the surface texture is approximately equivalent to or better than ASTM A802 A1-(Level I)")
  }else if(3259.8285< RoughRQ && RoughRQ <= 4573.182){
      paste(" The calculated Sq value is",round(501.4*log(RoughRQ) - 3386.3,0), "micro in. (",round((501.4*log(RoughRQ) - 3386.3)*.025399999,1),"micro m.)", "\n",
            "The calculated Sa value is",round(501.4*log(RoughRA) - 3386.3,0), "micro in. (",round((501.4*log(RoughRA) - 3386.3)*.025399999,1),"micro m.)", "\n",
            "For reference purpose only, the surface texture is approximately equivalent to or better than ASTM A802 A2-(Level II)")
  }else if(4573.182< RoughRQ && RoughRQ  <= 6169.509){
      paste(" The calculated Sq value is",round(501.4*log(RoughRQ) - 3386.3,0), "micro in. (",round((501.4*log(RoughRQ) - 3386.3)*.025399999,1),"micro m.)", "\n",
            "The calculated Sa value is",round(501.4*log(RoughRA) - 3386.3,0), "micro in. (",round((501.4*log(RoughRA) - 3386.3)*.025399999,1),"micro m.)", "\n",
            "For reference purpose only, the surface texture is approximately equivalent to or better than ASTM A802 A3-(Level III)")
    }else if( RoughRQ > 6169.509){
      paste(" The calculated Sq value is",round(501.4*log(RoughRQ) - 3386.3,0), "micro in. (",round((501.4*log(RoughRQ) - 3386.3)*.025399999,1),"micro m.)", "\n",
            "The calculated Sa value is",round(501.4*log(RoughRA) - 3386.3,0), "micro in. (",round((501.4*log(RoughRA) - 3386.3)*.025399999,1),"micro m.)", "\n",
            "For reference purpose only, the surface texture is approximately equivalent to or better than ASTM A802 A4-(Level IV)")
  }

})
  
  output$plott <- renderPlotly({
      data <- input$pointCloud
      if (is.null(data))
        return(NULL)
      data <- read.csv(data$datapath, header=input$header,sep = input$sep)
      
      dt <- data.frame(prcomp(data)$x)
      axz <- list(
        range <- c(min(dt$PC3)-.65,max(dt$PC3)+.65)
      )
      plot_ly(dt,x = ~PC1,y = ~PC2,z =~PC3,size = ~PC3,marker = 
                list(color = ~PC3, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
        layout(scene = list(xaxis = list(title = 'X-axis (mm)'),
                            yaxis = list(title = 'Y-axis (mm)'),
                            zaxis = list(title = 'Z-axis (mm)'))) %>%
        layout(scene = list(zaxis = axz),
               title = "Point Cloud After Applying PCA",
               annotations = list(text = "Range of z-coordinates \n
                                  after PCA",
                                  showarrow = FALSE,
                                  x =1,
                                  y = .5))
      
    
  })
  
  })

# Run the application 
shinyApp(ui = ui, server = server)

