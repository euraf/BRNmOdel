library(shiny)
library(igraph)
library(lavaan)
library(purrr)
library(semPlot)
#set.seed(1978)
library(reactlog)

# tell shiny to log all reactivity
reactlog_enable()

#### Global ####


### access functions to circumvent the Error in [: Can't index reactivevalues with `[` error
#' returns value at index
#'
#' @param rv reactiveValues()
#' @param field character vector
#'
#' @return value at index
access_rv = function(rv,field){
  isolate(purrr::pluck(rv, !!!field))
}

set_rv = function(rv,field,value){
  isolate(purrr::pluck(rv, !!!field) <- value )
}


#' Function to run 1 time step of probabilistic transitions
#'
#' @param state logical vector of current state (size = n, number of nodes in the model) OR 2 dimensional array with nodes in rows and repetitions in columns
#' @param transitions list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified, d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
#'
#' @return logical vector of t+1 state (size = number of nodes in the model)
#' @export
#'
#' @examples {transitions<-list(N1=array(c(0.2, 0.4, 0.8, 0.1), dim=c(2,2), dimnames=list(N2=c(0,1),N1=c(0,1))),
#' N2=array(array(c(0.2, 0.4, 0.8, 0.1), dim=c(2,2), dimnames=list(N2=c(0,1),N2=c(0,1)))),
#' N3=array(array(c(0.2, 0.4, 0.8, 0.1), dim=c(2,2), dimnames=list(N2=c(0,1),N3=c(0,1))))
#' )
#' run1step(state=c(0,1,1), transitions=transitions)
#' }
run1step<-function(state, transitions){
  future<-state
  for (ele in names(transitions)){
    transitionarray<-transitions[[ele]]
    contextnodes<-strsplit(names(dimnames(transitionarray))[1], split="_")[[1]]
    codecontext<-paste(state[contextnodes], collapse = "")
    if(! codecontext %in% rownames(transitionarray)) stop(paste("Apparently the transition for", ele, "in state", state[ele], "in context", codecontext, "is not among the provided transition probabilities"))
    proba<-transitionarray[codecontext,as.character(state[ele])]
    #todo: instead of paste/strsplit, use multidimensional arrays, probably more efficient... need to format the truth tables as n-dimensional
    #proba<-transitionarray[matrix(as.character(state[contextnodes]), nrow=1),as.character(state[ele])]
    changes<-rbinom(n=1, size=1, prob=proba)
    if(length(dim(state)>1)){
      future[ele,changes]<- !future[ele,changes] #to do: fix this, it does not work (creates NA)
    } else future[ele][changes]<- !future[ele][changes]
  }
  return(future)
}



#' run simulation with possibility of pesticides (=> the loops are switched between time and reps to all computing mean values to decide pesticied use) and perturbation
#'
#' @param initstate logical vector (named with node names) of initial states
#' @param transitions list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
#' @param timesteps number of timesteps to run
#' @param nbreps number of repetitions to run
#' @param perturbation data.frame with columns (when, what, value) when=timestep at which the perturbation occurs, what=name of the node that is modified by the perturbation and value to impose in all simulations at perturbation time
#' @param pesticides vector (named after the pest nodes) of thresholds ABOVE which to spray pesticide to make the node go back to 0
#'
#' @return an array with dimensions time, nodes, reps
#' @export

runexpe<-function(initstate, transitions, timesteps=52, nbreps=10, perturbation=NULL, pesticides=NULL){
  #keep only non-empty transitions (not done at the beginning because we need to keep the placeholders for static nodes in case we want to make thme dynamic)
  transitions<-transitions[sapply(transitions, function(x) return(names(dimnames(x))[1] !=""))]
  #prepare array to host all results
  allsims<-array(0, dim=c(timesteps, length(initstate), nbreps), dimnames=list(time=NULL,nodes=names(initstate), reps=NULL))
  for(r in 1:nbreps) allsims[1,,r]<-initstate #initialise
  starttime<-2 #initialise start time
  if(!is.null(perturbation)) if(nrow(perturbation)>0) {
    times<-perturbation$when #we extract the perturbation times and order them
    times<-unique(times[order(times)])
    times<-times[times>1 & times<timesteps] #keep only those within the simulations
    for (t in times){ #for each end of interval,
      for(i in starttime:t) { #run until then
        for(r in 1:nbreps) allsims[i,,r]<-run1step(state=allsims[i-1,,r], transitions=transitions) #run  all reps for 1 step
        #now that run1step accepts arrays, do all reps at the same time
        #allsims[i,,]<-run1step(state=allsims[i-1,,], transitions=transitions) #run  all reps for 1 step
        for (p in names(pesticides))  { #apply pesticides
          meanvalue<-mean(allsims[i,p,])
          if(meanvalue>=pesticides[p]) allsims[i,p,]<-0
        }
      }
      #apply perturbation(s)
      todayp<-perturbation[perturbation$when==t,]
      for(r in 1:nbreps) allsims[t,todayp$what,r]<- todayp$value
      starttime<-t+1
    }
  }
  #run the rest
  for(i in starttime:timesteps) {
    for(r in 1:nbreps) allsims[i,,r]<-run1step(state=allsims[i-1,,r], transitions=transitions) #run all reps for 1 step
    #now that run1step accepts arrays, do all reps at the same time
    #allsims[i,,]<-run1step(state=allsims[i-1,,], transitions=transitions) #run  all reps for 1 step
    for (p in names(pesticides))  { #apply pesticides
      meanvalue<-mean(allsims[i,p,])
      if(meanvalue>=pesticides[p]) allsims[i,p,]<-0
    }
  }
  return(allsims)
}
#example use
# transitions<-list(N1=array(c(0.2, 0.4, 0.8, 0.1), dim=c(2,2), dimnames=list(N2=c(0,1),N1=c(0,1))),
#                   N2=array(c(0.2, 0.4, 0.8, 0.1), dim=c(2,2), dimnames=list(N2=c(0,1),list(N2=c(0,1)))),
#                   N3=array(c(0.2, 0.4, 0.8, 0.1), dim=c(2,2), dimnames=list(N2=c(0,1),list(N3=c(0,1))))
# )
#initstate<-c(N1=0,N2=1,N3=1)
#nbreps<-1000
#allsims<-runexpe(initstate=initstate, transitions=transitions, timesteps=52, nbreps=100)
# add perturbation (N2=1 at t=20)
#allsims<-runexpe(initstate=initstate, transitions=transitions, timesteps=52, nbreps=100, perturbation=data.frame(when=c(20, 35), what=c("N1", "N2"), value=c(1, 0)))
#add pesticide against N1 if mean value>0.3
#allsims<-runexpe(initstate=initstate, transitions=transitions, timesteps=52, nbreps=100, pesticides = c(N1=0.3))

plotmeans<-function(allsims,...){
  #take means accross all reps and plot mean values
  means<-apply(allsims,c("time", "nodes"), mean)
  couleurs<-rainbow(n=ncol(means), start=1/6) ; names(couleurs)<-colnames(means)
  plot(0,0, type="n", xlab="time", ylab="mean value", xlim=c(0,nrow(means)), ylim=c(0,1),...)
  for(n in colnames(means)) lines(means[,n], col=couleurs[n],...)
  legend("topright", lty=1, col=couleurs, legend=names(couleurs))
  return(invisible())
}

# allsims<-runexpe(initstate=initstate, transitions=transitions, timesteps=52, nbreps=100,
#                  pesticides = c(N1=0.3),
#                  perturbation=data.frame(when=c(20, 35), what=c("N1", "N2"), value=c(1, 0)))
#
# plotmeans(allsims, main="pesticide and perturbation")

#setwd("~/a_ABSys/DIGITAF/WP2farmers/task2.2_treecropperformance/ColinsPestModel/BooleanPestModel")
library(openxlsx)

#### model definition file read, create template and write

#' Create the dataframes necessary to write the template excel file based on the nodes and edges provided as inputs
#'
#' @param nodes character vector od nodes names was #1 element character vector of node names, separated by commas, semicolons or spaces this comes from the textinput
#' @param listedges data?frame with 2 columns (EffectOn and EffectOf) # was named list (on element per node) of character vectors of causal node names (the names of the list are in the form paste("effectson", i, sep="_") with i the focal node name) this comes from the checkboxgroupinput
#' 
#' @return a list (elements named Nodes, Effects, Transitions, initialisation) of data.frames to be used to write the excel file to be downloaded to be used as template to fill in the transitions of a new model
#' @export
#'
#' @examples
createtemplate<-function(nodes, edges){
  #split the nodes
  Nodes<-nodes
  Effects<-edges
  Transitions<-data.frame()
  for(i in unique(edges$EffectOn)){
    dat<-edges[edges$EffectOn==i,]
    if(nrow(dat)>0){
      EffectOn<-i
      effectnodes<-dat$EffectOf
      #prepare Transitions tab
      EffectOf<-paste(effectnodes, collapse="_")
      context<-as.character(apply(expand.grid(lapply(as.list(effectnodes), function(x) return(c(0,1)))),1, paste, collapse=""))
      toto<-data.frame(EffectOn=EffectOn, EffectOf=EffectOf,context=context, low2high="", high2low="")
      Transitions<-rbind(Transitions, toto)
    }
  }
  initialisation<-data.frame(initAt1=character(0))
  return(list(Nodes=Nodes, Effects=Effects,Transitions=Transitions,initialisation=initialisation))
}

#' Read model definition from excel
#'
#' @param filepath path of the excel file defining a new model
#'
#' @return a named list with elements nodes, edges, transitions and initvalues. nodes=vector of node names, edges = data.frame with 2 columns (EffectOf and EffectOn) with 1 row per edge, transitions=list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified, d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0). initvalues a named logical vector (names = nodes) 
#' @export
#'
#' @examples
readexcel<-function(filepath="DefaultValues.xlsx"){
  nodes<-read.xlsx(filepath, sheet = "Nodes", colNames =FALSE)[,1] #character vector of node names
  names(nodes)<-nodes
  edges<-read.xlsx(filepath, sheet = "Effects") #data.frame with 2 columns: EffectOn	EffectOf (EffectOn is the name of thefocal node and EffectOf is name of the node having an effect on it, (i.e. causing the context of the focal node)
  import<-read.xlsx(filepath, sheet = "Transitions") #data.frame with columns EffectOn	EffectOf	context	low2high	high2low (context is a character value of 010011 codes of context states)
  transitions<-lapply(nodes, function(x){
    toto<-import[import$EffectOn==x,]
    if(nrow(toto)>0){
      dim1<-list(toto$context)
      names(dim1)<-unique(toto$EffectOf)
      dim2<-list(c("0", "1"))
      names(dim2)<-x
      toto<-as.array(as.matrix(toto[,c("low2high", "high2low")]))
      dimnames(toto)<-c(dim1, dim2)
      #order by context in order to find more quickly the values
      toto[order(dim1[[1]], decreasing=FALSE),]
      return(toto)
    } else {
      toto<-array(0, dim=c(0,2))
      dim2<-list(c("0", "1"))
      names(dim2)<-x
      dimnames(toto)<-c(list(NULL), dim2)
      return(toto)
    }
  }) #Transitions is formatted according to the (not very smart) format chosen initially: list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
  initvalues<-rep(FALSE, length(nodes)) ; names(initvalues)<-nodes #named logical vector
  initat1<-read.xlsx(filepath, sheet = "initialisation")[,1] #vector names of nodes initialised at 1
  initvalues[initat1]<-TRUE
  return(list(nodes=nodes, edges=edges, transitions=transitions, initvalues=initvalues))
}
# read in default model
DEFAULTVALUES<-readexcel()

#### UI ####

library(shinythemes)
library(shinyBS)

ui <- fluidPage(
  theme = shinytheme("cosmo"),
  tags$style('#sidebar {border: none; background-color: transparent; padding: 0;}'),
  titlePanel(# app title/description
    "Boolean Regulatory Network mOdel (BRNmOdel) - click on the boxes to expand/collapse them (red boxes are not working properly yet...)"),
  sidebarLayout(
    sidebarPanel(
      id = "sidebar",
      
      bsCollapsePanel(
        "Play with the model",
        uiOutput(outputId="boxinitialisation"),
        numericInput(inputId="nbreps", label="# repetitions", value=100),
        numericInput(inputId="timesteps", label="# time steps", value=52),
        bsCollapsePanel(
          "Outbreaks",
          uiOutput(outputId="boxoutbreaks"),
          style = "success"
        ),
        style = "success"
      ), #end play
      bsCollapsePanel(
        "Define your own model",
        fluidRow(
          fileInput("uploadexcel", "Upload a complete model definition file", accept="xlsx"),
          actionButton("resetdefault", "Reset default values")
        ),
        h2("... OR ..."),
        p("use this interface for step by step model definition:"),
        bsCollapsePanel(
          "Model structure",
          uiOutput(outputId="modeldefinition1"),
          uiOutput(outputId="modeldefinition2"),
          uiOutput(outputId="modeldefinition3"),
          style = "danger"
        ),
        bsCollapsePanel(
          "Model parameterization",
          uiOutput(outputId="boxParameterisation"),
          style = "danger"
        ),
        style = "success"
      ), #end define
      width = 4
    ),
    mainPanel(fluidPage(fluidRow(
      
      plotOutput("plot1"),
      bsCollapsePanel(
        "Model analysis",
        plotOutput("plotGraphmainPanel"),
        p("SEM analysis (only significant effects are included, top row are the values at the previous timestep)"),
        p("To do check the analysis because results are not the same as when I ran the analysis separately (it had a negative effect of diseases_previous on yield)"),
        plotOutput("plotSEM"),
        style = "danger"
      )
    ))
    )
    
  )
)

#### Server ####

server <- function(input, output, session) {
  rvstructure <- reactiveValues(nodes=DEFAULTVALUES$nodes, #vector of node names
                                edges=DEFAULTVALUES$edges) #data.frame with 2 columns: EffectOn	EffectOf (EffectOn is the name of thefocal node and EffectOf is the concatenation (name1_name2_name3 of the nodes causing the contexte of the focal node)
  rvparam<-reactiveValues(transitions=DEFAULTVALUES$transitions) #list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
  rvuserstructure<-reactiveValues(nodes=DEFAULTVALUES$nodes, #vector of node names
                                  edges=DEFAULTVALUES$edges) #data.frame with 2 columns: EffectOn	EffectOf (EffectOn is the name of thefocal node and EffectOf is the concatenation (name1_name2_name3 of the nodes causing the contexte of the focal node)
  rvuserparam<-reactiveValues(transitions=DEFAULTVALUES$transitions) #list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
  # to do: currently, interacting with any control causes the model to re-run, which is annoying because it's time consuming. immediate effect should happen only in the "model settings box. actions on the model definition and parameterization boxes should do nothing until the "update model" buttons are clicked.
  # to try to solve this, separate rvstructure and rvparam from the structure and params comming from the "define your own model box"
  #actually the problem comes from the fact tha we use reactiveValuesToList, which returns all inputs even those we do not want. solution not yet implemented : consstruct first the name of widget we need, and extract them with acces_rv() function 
  
  #### dynamic UI ####
  ####   boxinitialisation ####
  output$boxinitialisation<-renderUI({ #we put it in server side because eventually it is dynamic according to model specification by user
    tagList(
      p("click on the nodes that start at 'high level'"),
      checkboxGroupInput(inputId="initialisation",
                         label="Initialisation",
                         choices = rvstructure$nodes,
                         selected=DEFAULTVALUES$initvalues),
    )
  })
  
  ####   boxoutbreaks ####
  output$boxoutbreaks<-renderUI({
    p("Move the slider to create outbreaks (periods of forced hich value for  given node)")
    lapply(rvstructure$nodes, function(i) {
      return(sliderInput(inputId=paste("outbreakrange", i, sep="_"),
                         label=i,
                         min=-1,
                         max=input$timesteps,
                         step=1,
                         value=c(-1,-1)))
    })
  })
  
  ####   box modeldefinition ####
  output$modeldefinition1<-renderUI({
    ####   boxes effects on ####
    # controls_list <- lapply(rv$transitions, function(i) {
    #   control_type <- "checkboxGroupInput"
    #   focalnode<-names(dimnames(i))[2]
    #   input_id <- paste("effectson", focalnode, sep="_")
    #   #message(paste("creation", input_id))
    #   choices <- names(rv$transitions)
    #   selected<-unlist(strsplit( names(dimnames(i))[1], "_"))
    #   labelinput<-paste("Select the nodes having an effect on", focalnode)
    #   return(checkboxGroupInput(input_id, label = labelinput, choices = choices, selected=selected))
    # })
    tagList(
      textInput(inputId = "nodes", label="Node names", value=paste(rvuserstructure$nodes, collapse=",")),
      actionButton(inputId="buttonvalidatenodes", label="Validate nodes") #because it takes time to type, so it's better to indicate when it is finished
    )
  })
  output$modeldefinition2<-renderUI({ # warning this piece of code is duplicated (in observeEvent(eventExpr=input$nodes, )
    checkboxeseffectof <- lapply(rvuserstructure$nodes, function(i) {
      control_type <- "checkboxGroupInput"
      edges<-rvuserstructure$edges
      effectsof<-edges[edges$EffectOn==i, "EffectOf"] #effectsof is a vector of context nodes
      input_id <- paste("effectson", i, sep="_")
      #message(paste("creation", input_id))
      choices <- setdiff(rvuserstructure$nodes, i)
      selected<-effectsof
      labelinput<-paste("Select the nodes having an effect on", i)
      return(checkboxGroupInput(input_id, label = labelinput, choices = choices, selected=selected))
    })
    tagList(checkboxeseffectof)
  })
  
  output$modeldefinition3<-renderUI({
    tagList(
      plotOutput(outputId = "plotGraphsideMenu"),
      downloadButton(outputId="downloadtemplate", label="Download the excel template"),
      p("save the downloaded file and fill in the transisiotn probabilities in tab 'Transitions', then save the file and upload it with the buton at the top of box 'Define your own model' to load your own model")
    )
  })
  
  
  
  
  ####   boxParameterisation  ####
  output$boxParameterisation<-renderUI({
    controls_list <- lapply(rvuserparam$transitions, function(i) {
      focalnode<-names(dimnames(i))[2]
      effectsof<-names(dimnames(i))[1]
      if(effectsof!=""){
        topinfo<-list(p(strong(paste(focalnode, "depends on", effectsof))))
        controls1node<-list()
        for(context in dimnames(i)[[1]]) {
          input_id1 <- paste("probal2h", focalnode, context, sep="_")
          input_id2 <- paste("probah2l", focalnode, context, sep="_")
          controls1node<-c(controls1node,
                           tagList(fluidRow(column(width=6, #sliderInput(input_id1, label = paste(context, "low to high"),
                                                   #            min=0, max=1, value=i[context,"0"])),
                                                   numericInput(input_id1, label = paste(context, "low to high"),
                                                                value=i[context,"0"], step = 0.005)),
                                            column(width=6, #sliderInput(input_id2, label = paste(context, "high to low"),
                                                   #            min=0, max=1, value=i[context,"1"]))))
                                                   numericInput(input_id2, label = paste(context, "high to low"),
                                                                value=i[context,"1"], step = 0.005))))
                           
          )
        }
        
        return(bsCollapsePanel(title=focalnode,topinfo, controls1node))
      } else {
        topinfo<-list(p(strong(paste(focalnode, "is fixed"))))
        return(bsCollapsePanel(title=focalnode,topinfo))
      }
      
    })
    tagList(
      p("You can further finetune the parameters:"),
      tagList(controls_list),
      actionButton("updatetransitions", label="Update transitions (not coded yet)"),
      downloadButton(outputId="finetunedfile", label="Download fine-tuned file")
      #h2("Once you have completed the transitions, please upload the completed file with the upload button on top of the model definition box"),
      
    )
  })
  
  
  #### server logic ####
  
  #### Perturbations ####
  perturbationData<-reactive({
    perturbation<-data.frame()
    controls <- reactiveValuesToList(input)
    controls<-controls[grepl(x=names(controls), pattern="outbreakrange_", fixed=TRUE)]
    if(length(controls)>0){
      hasbeenset<-sapply(controls, function (x) return(!identical(as.numeric(x), c(-1,-1))))
      controls<-controls[hasbeenset]#we remove those that were not modified
      if(length(controls)>0) for(i in 1:length(controls)){
        what<-gsub(pattern="outbreakrange_", replacement="", names(controls)[i])
        rangetimes<-controls[[i]]
        newperturb<-data.frame(when=max(rangetimes[1],1):max(rangetimes[2],1), what=what, value=1)
        perturbation<-rbind(perturbation, newperturb)
      }
    }
    return(perturbation)
  })#perturbationData() is a data.frame of the perturbations (one line per time-node, with columns when what value)
  
  #### upload excel file and replace model ####
  observeEvent(
    eventExpr=c(input$uploadexcel),
    handlerExpr={
      req(input$uploadexcel)
      newvalues<-readexcel(input$uploadexcel$datapath)
      rvstructure$nodes <- newvalues$nodes #vector of node names
      rvstructure$edges <- newvalues$edges #data.frame with 2 columns: EffectOn	EffectOf (EffectOn is the name of thefocal node and EffectOf is the concatenation (name1_name2_name3 of the nodes causing the contexte of the focal node)
      rvparam$transitions <-newvalues$transitions #list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
      rvuserstructure$nodes <- newvalues$nodes #vector of node names
      rvuserstructure$edges <- newvalues$edges #data.frame with 2 columns: EffectOn	EffectOf (EffectOn is the name of thefocal node and EffectOf is the concatenation (name1_name2_name3 of the nodes causing the contexte of the focal node)
      rvuserparam$transitions <-newvalues$transitions #list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
    },
    ignoreInit =TRUE
  )
  
  #### reset default values ####
  observeEvent(
    eventExpr=c(input$resetdefault),
    handlerExpr={
      rvstructure$nodes <- DEFAULTVALUES$nodes #vector of node names
      rvstructure$edges <- DEFAULTVALUES$edges #data.frame with 2 columns: EffectOn	EffectOf (EffectOn is the name of thefocal node and EffectOf is the concatenation (name1_name2_name3 of the nodes causing the contexte of the focal node)
      rvparam$transitions <-DEFAULTVALUES$transitions #list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
      rvuserstructure$nodes <- DEFAULTVALUES$nodes #vector of node names
      rvuserstructure$edges <- DEFAULTVALUES$edges #data.frame with 2 columns: EffectOn	EffectOf (EffectOn is the name of thefocal node and EffectOf is the concatenation (name1_name2_name3 of the nodes causing the contexte of the focal node)
      rvuserparam$transitions <-DEFAULTVALUES$transitions #list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
      
    },
    ignoreInit =TRUE
  )
  
  #### update the user-defined nodes ####
  observeEvent(eventExpr=input$buttonvalidatenodes,
               handlerExpr={
                 usernodes<-input$nodes
                 usernodes<-gsub(pattern=";", replacement=",", x=usernodes)
                 usernodes<-unname(unlist(strsplit(usernodes, split=",")))
                 if(!identical(usernodes,rvuserstructure$nodes)){
                   rvuserstructure$nodes<-usernodes #this will cause modeldefinition2 UI to update #was <<-
                 }
               },
               ignoreInit = TRUE
  )
  
  #### update the user-defined relationships structure ####
  # icicicic i'm trying to create an observer that reacts to the change of the input$effectson_ checkboxes and only those
  observe({
    edges<-data.frame()
    #print(names(isolate(reactiveValuesToList(input))))
    for(i in rvuserstructure$nodes){
      print(i)
      effects<-access_rv(input, paste("effectson", i, sep="_"))
      print(effects)
      if(!is.null(effects)>0) {
        toto<-data.frame(EffectOn=i, EffectOf=effects)
        edges<-rbind(edges, toto)
      }
    }
    print(edges)
    if(ncol(edges)>0) rvuserstructure$edges<-edges
    #print(effects)
    # controls <- reactiveValuesToList(input)
    # nodes<-gsub(pattern="; ", replacement=",", x=controls[["nodes"]])
    # nodes<-unname(unlist(strsplit(nodes, split=",")))
    # listedges<-controls[grepl(x=names(controls), pattern="effectson_", fixed=TRUE)]
    # edges<-data.frame()
    # if(length(listedges)>0) for(i in names(listedges)){
    #   effects<-listedges[[i]]
    #   if(length(effects)>0) {
    #     toto<-data.frame(EffectOn=gsub(pattern="effectson_", replacement="", x=i), EffectOf=effects)
    #     edges<-rbind(edges, toto)
    #   }
    # } else edges<-data.frame(EffectOn=character(0),EffectOf=character(0) )
    # rvuserstructure$edges<-edges #was <<-
  })
  
  #### create and download template with the user model structure ####
  
  output$downloadtemplate <- downloadHandler(
    filename = function() {
      return("Template_modeldefinition.xlsx")
    },
    content = function(file) {
      nodes<-isolate(rvuserstructure$nodes)
      edges<-isolate(rvuserstructure$edges)
      listdf<-createtemplate(nodes=nodes, 
                             edges=edges)
      write.xlsx(listdf, file = file)
    }
  )
  
  
  
  #### modification of transitions ####
  observeEvent(eventExpr=input$updateTransitions,
               handlerExpr={
                 controls <- reactiveValuesToList(input)
                 controls<-c(controls[grepl(x=names(controls), pattern="^probal2h")],
                             controls[grepl(x=names(controls), pattern="^probah2l")])
                 toto<-strsplit(names(controls), split="_")
                 probatypes<-sapply(toto,"[[", 1)
                 focalnodes<-sapply(toto,"[[", 2)
                 context<-sapply(toto,"[[", 3)
                 #transitions=list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
                 for(n in names(rvparam$Transitions)){
                   subcontext<-context[focalnodes==n]
                   subprobatypes<-context[focalnodes==n]
                   newvalues<-controls[focalnodes==n]
                   rvparam$Transitions[[n]][subcontext, as.character(subprobatypes=="probah2l")]<-newvalues #was <<-
                 }
                 #NB once transitions have been changed to 3dimensionnal arrays, something like this will be easier
                 #rv$Transitions[focalnodes, context, as.character(probatypes=="probah2l")]<-controls #was <<-
               },
               ignoreInit = TRUE
  )
  
  
  
  #### Main plot (and running of model) ####
  # output$plot1 <- renderPlot({
  #   initstate<-logical(length(rvstructure$nodes)); names(initstate)<-rvstructure$nodes
  #   initstate[names(initstate) %in% input$initialisation]<-1 #to do: correct (minor) bug that prevents default init states to be taken account before the user open the model settings box (because then the input$initialisation control does not exist yet)
  #   nbreps<-input$nbreps
  #   timesteps<-input$timesteps
  #   perturbation<-perturbationData()
  #   allsims<-runexpe(initstate=initstate, 
  #                    transitions=rvparam$transitions, 
  #                    timesteps=timesteps, 
  #                    nbreps=nbreps,
  #                    perturbation=perturbation)
  #   plotmeans(allsims)
  #   
  # })
  #do it in an observeEvent so that it is rendered only if inputs from the model settings box are changed, or if the update buttons from the other boxes are clicked
  #to do: fix this, it does not work: it continues redrawing the raph every time anything is changed, but it doesn t take into account the changed probabilities
  observeEvent(handlerExpr=isolate(output$plot1 <- renderPlot({
    initstate<-logical(length(rvstructure$nodes)); names(initstate)<-rvstructure$nodes
    initstate[names(initstate) %in% input$initialisation]<-1 #to do: correct (minor) bug that prevents default init states to be taken account before the user open the model settings box (because then the input$initialisation control does not exist yet)
    nbreps<-input$nbreps
    timesteps<-input$timesteps
    perturbation<-perturbationData()
    allsims<-runexpe(initstate=initstate, 
                     transitions=rvparam$transitions, 
                     timesteps=timesteps, 
                     nbreps=nbreps,
                     perturbation=perturbation)
    plotmeans(allsims)
  })),
  eventExpr={input$nbreps
    input$timesteps
    perturbationData()
    input$updateTransitions
    input$updateedges
  },
  ignoreNULL =FALSE
  )
  
  #### small graph side menu in model definition box ####
  # icicic I don't understand why this does not update when rvuserstructure$nodes is changed !
  output$plotGraphsideMenu<-renderPlot({
    #### plot graph of relationships ####
    #rvuserstructure$nodes #just to update when nodes are modified
    #data<-rvuserstructure$edges
    #does not work let's do it directly from input
    for(i in rvuserstructure$nodes){
      print(i)
      effects<-access_rv(input, paste("effectson", i, sep="_"))
      print(effects)
      if(!is.null(effects)>0) {
        toto<-data.frame(EffectOn=i, EffectOf=effects)
        edges<-rbind(edges, toto)
      }
    }
    print(edges)
    data<-edges
    g <- graph_from_data_frame(data[,c("EffectOf", "EffectOn")], directed=TRUE)
    plot(g, layout=layout.circle, main="Declared effects")
  })
  
  
  #### graphs model analysis ####
  output$plotGraphmainPanel<-renderPlot({
    #### plot graph of relationships ####
    if(TRUE) {
      g <- graph_from_data_frame(rvstructure$edges[,c("EffectOf", "EffectOn")], directed=TRUE)
      plot(g, layout=layout.circle, main="Declared effects")
    } else plot(0,0, type="n", bty="n", xlab="just to test plots in collapsible box", ylab="", xaxt="n", yaxt="n")
    
  })
  
  output$plotSEM<-renderPlot({
    if("Trees" %in% rvstructure$nodes){
      #### run model with trees and without trees (to do: make this generic to run all combinations of fixed nodes states) ####
      initstate<-logical(length(rvstructure$nodes)) ; names(initstate)<-rvstructure$nodes
      initstate["Trees"]<-1
      allsims<-runexpe(initstate=initstate, 
                       transitions=rvparam$transitions, 
                       timesteps=input$timesteps, 
                       nbreps=input$nbreps)
      notrees<-initstate
      notrees["Trees"]<-0
      allsims2<-runexpe(initstate=notrees,
                        transitions=rvparam$transitions,
                        timesteps=52, 
                        nbreps=100#,
                        #                  pesticides = c(N1=0.3),
                        #                  perturbation=data.frame(when=c(20, 35), what=c("N1", "N2"), value=c(1, 0))
      )
      
      #### format simulaed data ####
      #allsimstreenotree<-abind::abind(allsims, allsims2)
      addprevioustimestep<-function(allsims){
        moyennestrees<-as.data.frame(apply(allsims,c(1,2), mean))
        toto<-names(moyennestrees); names(toto)<-toto
        decal<-as.data.frame(lapply(toto, function(x) numeric(nrow(moyennestrees))))
        decal[,toto]<-lapply(moyennestrees[,toto], function(x) {titi<-x ; titi[2:length(x)]<-x[1:(length(x)-1)]; return(titi)})
        names(decal)<-paste(toto,"previous", sep="_")
        moyennestrees$time<-1:nrow(moyennestrees)
        moyennestrees<-cbind(moyennestrees, decal)
        return(moyennestrees)
      }
      
      dataforSEM<-rbind(addprevioustimestep(allsims), addprevioustimestep(allsims2))
      
      #### run the SEM analysis ####
      #install.packages("lavaan", dependencies=T)
      #install.packages("semPlot", dependencies=T)
      
      
      edges<-rvstructure$edges
      model<-paste(sapply(unique(edges$EffectOn), function(x) {
        paste(x, paste(c(paste(x, "previous", sep="_"), 
                         paste(edges[edges$EffectOn==x, "EffectOf"], "previous", sep="_")),
                       collapse="+"), sep="~")
      }), collapse ="\n")
      # model<- '
      # NaturalEnemies~NaturalEnemies_previous+Trees_previous+Weeds_previous+Diseases_previous+Pests_previous
      # Pests~Pests_previous+Trees_previous+NaturalEnemies_previous+Weeds_previous+Diseases_previous
      # Yield~+Yield_previous+Weeds_previous+Diseases_previous+Pests_previous
      # Diseases~Diseases_previous+Pests_previous+Trees_previous+Weeds_previous
      # Weeds~Weeds_previous+Trees_previous+Pests_previous+Yield_previous+Diseases_previous
      # '
      
      fit <- sem(model, data=dataforSEM, estimator = "ML")
      resume<-summary(fit)
      #print(resume)
      # rerun without unsignificant efffects
      testeffets<-resume$pe
      effectstokeep<-testeffets[!is.na(testeffets$pvalue) & testeffets$pvalue<0.05 &testeffets$op=="~",]
      model2<- paste(paste(effectstokeep$lhs, effectstokeep$op, effectstokeep$rhs, sep=""), collapse="\n")
      fit2 <- sem(model2, data=dataforSEM, estimator = "ML")
      resume2<-summary(fit2)
      #print(resume2)
      # same but not taking into account autocorrelation
      #effectstokeepnoauto<-testeffets[!is.na(testeffets$pvalue) & testeffets$pvalue<0.05 &testeffets$op=="~" & paste(testeffets$lhs, "previous", sep="_")!=testeffets$rhs,]
      #model3<- paste(paste(effectstokeepnoauto$lhs, effectstokeepnoauto$op, effectstokeepnoauto$rhs, sep=""), collapse="\n")
      #fit3 <- sem(model3, data=dataforSEM, estimator = "ML")
      #resume<-summary(fit3)
      #plot the significant relationships
      semPaths(fit2, "std", fade = F, layout = "tree", intercepts = F, residuals  = F)
    } else { #Trees in not in the model
      plot(0,0, type="n", bty="n", xlab="This code currently only works if there is a fixed node named Trees", ylab="", xaxt="n", yaxt="n")
      
    }
  })
}

shinyApp(ui = ui, server = server)
