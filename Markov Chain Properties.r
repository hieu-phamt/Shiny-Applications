#Libraries needed
library(markovchain)
library(shiny)
library(ggplot2)
library(gridExtra)

# Define UI for application that draws a histogram
ui <- shinyServer(
  
  fluidPage(headerPanel("Absorbing Markov Chains"),

#Define user input parameters        
    sidebarPanel(
      
      helpText("Upload a .CSV file with absorbing states in the last columns"),
      helpText("WARNING: No header is needed for your .CSV file"),
      fileInput('file1', 'Choose file to upload',accept = c('text/csv','text/comma-separated-values','text/tab-separated-values','text/plain','.csv','.tsv')),
      textInput("start", "Enter your initial state:",1), 
      withMathJax(),
      helpText('\\(Theorem:\\) Consider an absorbing Markov chain with \\(t\\) transient states. Let \\(F\\) be a \\(t \\times t\\)
               matrix indexed by transient states, where \\(F_{ij}\\) is the expected number of visits to \\(j\\) given that the chain starts in \\(i\\)
               . Then, \\(F=(I-Q)^{-1}\\)'),
      helpText("From the above theorem we can calculate the exact number of steps needed as well as the exact absorption probabilities."),
      helpText("In addition, we know that if we do not have any legitimate absorbing states we can choose one."),
      textInput("absorb", "Enter your absorbing state:"),
      helpText("Exact Expected Number of Steps to Absorbing State"),
      uiOutput("ex1"),
      helpText('Exact Absorption Probability of Transient State \\(i\\) to Absorbing State \\(j\\)'),
      uiOutput("ex2")
      
      ),
    
    mainPanel("Your Transition Matrix", uiOutput("contents"),plotOutput("plot"),plotOutput("hist"),plotOutput("histp")),
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
  )
)
# Define server logic 
server <- shinyServer(function(input, output, session) {
  
# Upload the matrix from CSV file  
  output$contents <- renderTable({
    matrix <- input$file1
    if (is.null(matrix))
      return(NULL)
    matrix <- read.csv(matrix$datapath, header = F)
    matrix <- as.matrix(matrix)
    dim = dim(matrix)[1]
    # Give matrix columns appropriate names ie.e "1","2","3",....   
    states <- c(1:dim)
    states <- as.character(states)
    sub <- states[c(1:dim)]
    colnames(matrix) <- sub
    matrix
  })
  
# Creating out plot  
  output$plot <- renderPlot({
# Upload same matrix as before
      matrix <- input$file1
      if (is.null(matrix))
        return(NULL)
      matrix <- read.csv(matrix$datapath, header = F)
      matrix <- as.matrix(matrix)
      dim = dim(matrix)[1]
# Give matrix columns appropriate names ie.e "1","2","3",....    
      states <- c(1:dim)
      states <- as.character(states)
      sub <- states[c(1:dim)]
      colnames(matrix) <- sub
# This randomly generates our Markov Chain steps   
      chain <- new("markovchain",states = sub,transitionMatrix = matrix)   # markov function
# Initial starting point for Markov Chain     
      start <- as.character(input$start)
# The absorption state for our MC      
      absorb <- as.character(input$absorb)
# Loop for simulate 10 markov chain sequences
      list <- numeric() 
      simlist <- numeric()
# If statement is used if the uploaded CSV matrix has legitimate absorbing states. That is, if the diagonals equals one      
      if(length(absorbingStates(chain) != 0)){
      for(j in 1:2000){
#Pick starting path
        path <- rmarkovchain(ceiling(dim*(dim/2)*15),object = chain,t0 = start)
#keep track of when it happened to reach absorbing
        list[j] <- list(c(start,path))
        simlist[j] <- which(list[[j]] %in% absorbingStates(chain))[1] - 1
      }
# Next bunch of stuff is used ot create our 6 different graphs. I probably should have used a for loop.        
      x1 <- which(list[[1]] %in% absorbingStates(chain))[1] 
      x2 <- which(list[[2]] %in% absorbingStates(chain))[1] 
      x3 <- which(list[[3]] %in% absorbingStates(chain))[1] 
      x4 <- which(list[[4]] %in% absorbingStates(chain))[1] 
      x5 <- which(list[[5]] %in% absorbingStates(chain))[1] 
      x6 <- which(list[[6]] %in% absorbingStates(chain))[1] 
      x <- max(x1,x2,x3,x4,x5,x6) + 2
      y1 <- list[[1]][c(1:x)]
      y2 <- list[[2]][c(1:x)]
      y3 <- list[[3]][c(1:x)]
      y4 <- list[[4]][c(1:x)]
      y5 <- list[[5]][c(1:x)]
      y6 <- list[[6]][c(1:x)]
#Generating output for number of steps for each path
      steplist <- numeric()
      paste1 <- ""
      for(q in 1:6){
        steplist[q] <- which(list[[q]] %in% absorbingStates(chain))[1] - 1
        paste1 <- paste("path",q," =",steplist[q], "steps|| ",paste1 )
      }
#Generating output for confidence interval of steps
      m <- mean(simlist)
      v <- var(simlist)
#interval <- m - 1.96*sqrt(s)/sqrt(6)
      interval <- 1.96*sqrt(v)/sqrt(length(simlist))
      pasteinterval <- paste("95% C.L. for steps from state",start,"to absorption","is", round(m,2),"+-",round(interval,3), "steps")
#Step probability
      listlength <- numeric()
      paste2 <- ""
      problimits <- ""
      listprob <- numeric()
#Generating output for probabilities
      for(q in 1:2000){
        listprob[q] <- list[[q]][which(list[[q]] %in% absorbingStates(chain))[1]]
      }
      for(z in 1:length(unique(listprob))){
        listlength[z] <- length(which(listprob == unique(listprob)[z]))/length(listprob)
        
        if(listlength[z] != 0){
          paste2 = paste("State", absorbingStates(chain)[z]," is", round(listlength[z],3),"||",paste2 )
        }
      }
#Confidence limits for probabilities
      for(y in 1:length(unique(listprob))){
        problimits <- paste("State",absorbingStates(chain)[y]," is", "", round(listlength[y],3),
                           "+-",round(var(listprob)/sqrt(length(listlength)),2)*.1," || ",problimits)
      }
 #Else statement incase the matrix does not have any legitimate absorbing states    
      }else{
      list2 <- numeric()      
        for(k in 1:2000){
#Pick starting path
          path <- rmarkovchain(ceiling(dim*(dim/2)*15),object = chain,t0 = start)
#keep track of when it happened ot reach absorbing
          end_time <- which(path == absorb)[1]
          if(is.na(end_time)){list2[k] = list(c(1:ceiling(dim*(dim/2)*15)))
          simlist[k] <- length(list(list2[k]))}else{
#keep track of when it happened ot reach absorbing
          list2[k] <- list(c(start,path[c(1:end_time)]))
          simlist[k] <- length(path[c(1:end_time)]) 
          }
        }
# Next bunch of stuff is used ot create our 6 different graphs. I probably should have used a for loop.  
        x1 <- length(list2[[1]])
        x2 <- length(list2[[2]])
        x3 <- length(list2[[3]])
        x4 <- length(list2[[4]])
        x5 <- length(list2[[5]])
        x6 <- length(list2[[6]])
        x <- max(x1,x2,x3,x4,x5,x6) + 2 
        p1 <- x - length(list2[[1]]) 
        p2 <- x - length(list2[[2]])
        p3 <- x - length(list2[[3]])
        p4 <- x - length(list2[[4]]) 
        p5 <- x - length(list2[[5]])
        p6 <- x - length(list2[[6]])
        y1 <- append(list2[[1]],rep(list2[[1]][length(list2[[1]])],p1), after = length(list2[[1]]))
        y2 <- append(list2[[2]],rep(list2[[2]][length(list2[[2]])],p2), after = length(list2[[2]]))
        y3 <- append(list2[[3]],rep(list2[[3]][length(list2[[3]])],p3), after = length(list2[[3]])) 
        y4 <- append(list2[[4]],rep(list2[[4]][length(list2[[4]])],p4), after = length(list2[[4]]))
        y5 <- append(list2[[5]],rep(list2[[5]][length(list2[[5]])],p5), after = length(list2[[5]]))
        y6 <- append(list2[[6]],rep(list2[[6]][length(list2[[6]])],p6), after = length(list2[[6]])) 
        steplist <- numeric()
        paste1 <- ""
        for(q in 1:6){
          steplist[q] <- length(list2[[q]]) - 1
          paste1 <- paste("path",q," =",steplist[q], "steps|| ",paste1 )
        }
        if(length(which(is.na(simlist))) > 0 ){
        m <- simlist[!is.na(simlist)]
        m <- mean(simlist)
        v <- simlist[!is.na(simlist)]
        v <- var(simlist)
        interval <- 1.96*sqrt(v)/2000
        pasteinterval <- paste("95% C.L. for steps from state", start,"to state", absorb,"is", round(m,2),"+-",round(interval,3), "steps")
        }else{
        m <- mean(simlist)
        v <- var(simlist)
        interval <- 1.96*sqrt(v)/sqrt(length(simlist))
        pasteinterval <- paste("95% C.L. for steps from state", start,"to state", absorb,"is", round(m,2),"+-",round(interval,3), "steps") 
        }
      }
# When the if statement ends, both the IF and ELSE part meet here.
# To use the GGPLOT2 package we need to make a data frame. This was the hard part.
      df0 <- data.frame(x = 0, y = 0)
      df1 <- data.frame(x = c(1:x),y = y1)
      df2 <- data.frame(x = c(1:x),y = y2)
      df3 <- data.frame(x = c(1:x),y = y3)
      df4 <- data.frame(x = c(1:x),y = y4)
      df5 <- data.frame(x = c(1:x),y = y5)
      df6 <- data.frame(x = c(1:x),y = y6)
# This is our plot in all of its glory.     
      if(start == absorb){
        df_null <- data.frame(x = c(1:x),y = as.numeric(rep(absorb,x)))
        plot <- ggplot() + geom_line(data = df_null, aes(x=x,y=y),group = 1)+
        scale_y_discrete("State",limits = c(1:dim))+
        scale_x_discrete("Step",limits = c(1:x))+labs(title = "Monte Carlo Simulation")
        print(plot)
      }else{
      plot <- ggplot(data = df0, aes(x=x,y=y, group = 1)) +
        geom_line(data = df1, aes(x=x,y=y, group = 1,color = "maroon")) + geom_point(data = df1) +
        geom_line(data = df2,linetype = "longdash",aes(x=x,y=y, group = 1, color ="red"))+geom_point(data = df2) + 
        geom_line(data = df3 ,linetype = "longdash",aes(x=x,y=y, group = 1, color ="blue"))+geom_point(data = df3)+
        geom_line(data = df4,aes(x=x,y=y,linetype = "blank", group = 1, color ="black"))+geom_point(data = df4 )+
        geom_line(data = df5 ,aes(x=x,y=y, group = 1, color ="orange"))+geom_point(data = df5 )+
        geom_line(data = df6,linetype = "longdash",aes(x=x,y=y, group = 1, color ="brown"))+geom_point(data = df6 )+
        scale_x_discrete("Step",limits = c(1:x))+scale_y_discrete("State",limits = c(1:dim))+
        labs(title = "Monte Carlo Simulation of 6 sample paths")+  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
        scale_color_manual(labels = c("path 1","path 2","path 3","path 4","path 5","path 6"),values =c("maroon","red","blue","black","orange","brown"))
      plot = gridExtra::grid.arrange(plot,bottom = paste1)
      plot = gridExtra::grid.arrange(plot,bottom = pasteinterval)
      plot = gridExtra::grid.arrange(plot,bottom = "")
      if(length(absorbingStates(chain) != 0)){
      plot = gridExtra::grid.arrange(plot,bottom = paste("Probability of being absorbed in: ",paste2))
      plot = gridExtra::grid.arrange(plot,bottom = paste("95% C.L. for probability to absorbing :",problimits))}
      print(plot)}
  })
# Creating our histogram #######################################################################################################
  output$hist <-renderPlot({
    # Upload same matrix as before
    # matrix = read.csv("book1.csv", header = F)
    matrix <- input$file1
    if (is.null(matrix))
      return(NULL)
    matrix <- read.csv(matrix$datapath, header = F)
    matrix <- as.matrix(matrix)
    dim <- dim(matrix)[1]
    # Give matrix columns appropriate names ie. "1","2","3",....  
    states <- c(1:dim)
    states <- as.character(states)
    sub <- states[c(1:dim)]
    colnames(matrix) <- sub
    # This randomly generates our Markov Chain steps   
    chain <- new("markovchain",states = sub,transitionMatrix = matrix)   # markov function
    # Initial starting point for Markov Chain     
    start <- as.character(input$start)
    # The absorption state for our MC      
    absorb <- as.character(input$absorb)
    # Loop for simulate 10 markov chain sequences
    list <- numeric() 
    datasteps <- numeric()
    # If statement is used if the uploaded CSV matrix has legitimate absorbing states. That is, if the diagonals equals one      
    if(length(absorbingStates(chain) != 0)){
      for(j in 1:2000){
        #Pick starting path
        path <- rmarkovchain(ceiling(dim*(dim/2)*15),object = chain,t0 = start)
        #keep track of when it happened ot reach absorbing
        list[j] <- list(c(start,path))
        datasteps[j] <-  which(list[[j]] %in% absorbingStates(chain))[1]}
      }else{
      list2 <- numeric()      
      for(k in 1:2000){
        #Pick starting path
        path <- rmarkovchain(ceiling(dim*(dim/2)*15),object = chain,t0 = start)
        #keep track of when it happened ot reach absorbing
        end_time <- which(path == absorb)[1]
        if(is.na(end_time)){list2[k] = ceiling(dim*(dim/2)*15)}else{
        #keep track of when it happened ot reach absorbing
        list2[k] <- list(c(start,path[c(1:end_time)]))
        datasteps[k] <-  length(list2[[k]])
        }
      }
      }
    if(start == absorb){
      datasteps <- as.numeric(rep(0,2000))
      barplot(table(datasteps),main = "Number of Steps to Absorbing State from 2000 Paths",xlab = "steps",col ="green", axis(side=1, at=c(0:5))) 
    }
    barplot(table(datasteps),main = "Number of Steps to Absorbing State from 2000 Paths",xlab = "steps",col ="green")
  })
####Histogram of probabilities###############################################################################################################  
  output$histp <-renderPlot({
    # Upload same matrix as before
    # matrix = read.csv("book1.csv", header = F)
    matrix <- input$file1
    if (is.null(matrix))
      return(NULL)
    matrix <- read.csv(matrix$datapath, header = F)
    matrix <- as.matrix(matrix)
    dim = dim(matrix)[1]
    # Give matrix columns appropriate names ie.e "1","2","3",....   
    states <- c(1:dim)
    states <- as.character(states)
    sub <- states[c(1:dim)]
    colnames(matrix) <- sub
    # This randomly generates our Markov Chain steps   
    chain <- new("markovchain",states = sub,transitionMatrix = matrix)   # markov function
    # Initial starting point for Markov Chain     
    start <- as.character(input$start)
    # The absorption state for our MC      
    absorb <- as.character(input$absorb)
    # Loop for simulate 10 markov chain sequences
    list <- numeric() 
    listprob <- numeric()
    simlist = numeric()
    # If statement is used if the uploaded CSV matrix has legitimate absorbing states. That is, if the diagonals equals one      
    if(length(absorbingStates(chain) != 0)){
      for(j in 1:2000){
        #Pick starting path
        path <- rmarkovchain(dim*25,object = chain,t0 = start)
        #keep track of when it happened ot reach absorbing
        list[j] <- list(c(start,path))
        simlist[j] <- which(list[[j]] %in% absorbingStates(chain))[1] - 1
      }
      #Generating output for probabilities
      for(q in 1:2000){
        listprob[q] <- list[[q]][which(list[[q]] %in% absorbingStates(chain))[1]]
      }
      listprob <- as.numeric(listprob)
      hist(x = listprob,main = "Distribution of Absorbing States from 2000 Paths",xlab = "State", col = "yellow")
       }else{
      return(NULL)
      }
  })
# Creating our first output of expected number of steps to absorbing states
  output$ex1 <- renderTable({
# Upload matrix
    matrix <- input$file1
    if (is.null(matrix))
      return(NULL)
    matrix <- read.csv(matrix$datapath, header = F)
    matrix <- as.matrix(matrix)
    dimm <- dim(matrix)[1]
    # Give matrix columns appropriate names ie.e "1","2","3",....      
    states <- c(1:dimm)
    states <- as.character(states)
    sub <- states[c(1:dimm)]
    colnames(matrix) <- sub
    # This randomly generates our Markov Chain steps   
    chain <- new("markovchain",states = sub,transitionMatrix = matrix)   # markov function
    # Initial starting point for Markov Chain     
# Get length of matrix and the number of legitimate aborbing states
    num <- which(diag(matrix) == 1)
    length <- length(num)
# If we do have legitimate absorbing states we proceed with this if statement
    if(length != 0){
# We define Q like on page 120 in the book
    Q <- matrix[-c(dimm-length+1:length),-c(dimm-length+1:length)]
# We define J like our F matrix on page 125
    J <- solve(diag(dim(Q)[1])-Q)
# We get the row sum of matrix "F". The row sum gives up the expected number of steps. This is the result on page 126
    rowSums(J)
# We are making a table from our row sums
    table <- data.frame(rowSums(J))
# we are giving the column names
    colnames(table) <-c("Steps to an Absorbing State")
# Creates a new column named initial states    
    table$Initial <- c(1:(dimm-length))
    table <- table[c(2,1)]
# Prints table
    table
# ELSE statement if the matrix has no legitimate absorbing states
    }else{start = input$start
#Define starting point and absorbing state
    absorb <- input$absorb
    start <- as.numeric(start)
    absorb <- as.numeric(absorb)
# Define column names    
    states <- c(1:dimm)
    states <- as.character(states)
# Define Q as in the book
    Q <- matrix[-absorb,-absorb]
# Define J like F in the book
    J <- solve(diag(dim(Q)[1])-Q)
# We get the row sum of matrix "F". The row sum gives up the expected number of steps. This is the result on page 126
    rowSums(J)
    sub <- states[1:dim(Q)[1]]
# We are making a table from our row sums    
    table <- data.frame(rowSums(J))
# we are giving the column names
    colnames(table) <- c(paste("Steps to Absorbing State",absorb))
    sub <- states[c(1:dimm)]
    e <- sub[-c(absorb)]
# Creates a new column named initial states  
    table$Initial <- e
    table <- table[c(length(table):1)]
    table}
      })     
# Next output for the absorption probabilities
    output$ex2 <- renderTable({
# Upload matrix
      matrix <- input$file1
      if (is.null(matrix))
        return(NULL)
      matrix <- read.csv(matrix$datapath, header = F)
      matrix <- as.matrix(matrix)
      dimm <- dim(matrix)[1]
      states <- c(1:dimm)
      states <- as.character(states)
      num <- which(diag(matrix)==1)
      length <- length(num)
# If the matrix has legitimate absorbing states we go with this if stastement
      if(length != 0){
# Define Q and R as in the book, J is actually matrix "F" and U is matrix "FR"
        Q <- matrix[-c(dimm-length+1:length),-c(dimm-length+1:length)]
        J <- solve(diag(dim(Q)[1])-Q)
        R <- matrix[-c(dimm-length+1:length),c(dimm-length+1:length)]
        U <- J%*%R
# making a table
        table2 <- data.frame(U)
# this next bunch of stuff is used to create column names.
        list321 <- numeric()
        for(m in 1:length){name = paste("State",num[m]) 
        list321[m] <- c(name)}
        list321 <- list321[c(1:length)]
        names_col <- c(1:length)
        for(n in 1:length(list321)){
        names_col[n]<-list321[n]}
        colnames(table2) <- names_col
        names_row <- states[c(1:(dimm-length))]
# Now we make row names        
        rownames(table2) <- names_row
        table2$Initial <- c(1:(dimm-length))
        table2 <- table2[c(length(table2):1)]
        table2
        }else{
# Make absorption matrix if we choose an absorbing state.
# the explaination is selfexplainatory
        absorb <- input$absorb
        absorb <- as.numeric(absorb)
        Absorbing_State <- matrix(rep(1,dimm-1))
        Absorbing_State <- data.frame(Absorbing_State)
        sub <- states[c(1:dimm)]
        e <- sub[-c(absorb)]
        Absorbing_State$Initial <- e
        Absorbing_State <- Absorbing_State[c(length(Absorbing_State):1)]
        Absorbing_State
        }
  })
})
# Run the application 
shinyApp(ui = ui, server = server)
