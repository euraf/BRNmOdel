library(shiny)
set.seed(1978)



#### Global ####


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
      future[ele,changes]<- !future[ele,changes]
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

#read the default values from excel
#setwd("~/a_ABSys/DIGITAF/WP2farmers/task2.2_treecropperformance/ColinsPestModel/BooleanPestModel")
library(openxlsx)
nodes<-read.xlsx("DefaultValues.xlsx", sheet = "Nodes", colNames =FALSE)[,1] #character vector of node names
names(nodes)<-nodes
edges<-read.xlsx("DefaultValues.xlsx", sheet = "Effects") #data.frame with 2 columns: EffectOn	EffectOf (EffectOn is the name of thefocal node and EffectOf is the concatenation (name1_name2_name3 of the nodes causing the contexte of the focal node)
import<-read.xlsx("DefaultValues.xlsx", sheet = "Transitions") #data.frame with columns EffectOn	EffectOf	context	low2high	high2low (context is a character value of 010011 codes of context states)
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
initat1<-read.xlsx("DefaultValues.xlsx", sheet = "initialisation")[,1] #vector names of nodes initialised at 1
initvalues[initat1]<-TRUE


#### UI ####

library(shinythemes)
library(shinyBS)

ui <- fluidPage(
  theme = shinytheme("cosmo"),
  tags$style('#sidebar {border: none; background-color: transparent; padding: 0;}'),
  titlePanel(# app title/description
    "Boolean Regulatory Network mOdel (BRNmOdel)"),
  sidebarLayout(
    sidebarPanel(
      id = "sidebar",
      
      bsCollapsePanel(
        "Model settings",
        uiOutput(outputId="boxinitialisation"),
        numericInput(inputId="nbreps", label="# repetitions", value=100),
        numericInput(inputId="timesteps", label="# time steps", value=52),
        bsCollapsePanel(
          "Outbreaks",
          uiOutput(outputId="boxoutbreaks"),
          style = "success"
        ),
        style = "success"
      ),
      bsCollapsePanel(
        "Model definition",
        uiOutput(outputId="modeldefinition"),
        style = "danger"
      ),
      bsCollapsePanel(
        "Parameterisation",
        uiOutput(outputId="boxParameterisation"),
        style = "danger"
      ),
      width = 4
    ),
    mainPanel(plotOutput("plot1"))
    
  )
)

#### Server ####

server <- function(input, output, session) {
  rv <- reactiveValues(nodes=nodes,
                       edges=edges,
                       transitions=transitions) #list (one element per node that can change, size p, named according to node name) of arrays of dimension (F1_F2_F3...=2^d,NX=2) where the dimension name F1_F2_F3... indicates the influencing factors, NX is the name of the node being modified d is the number of influencing factors on the given node, the rows are named by the influencing factors' state and the columns are named 0 (transition from 0 to 1) and 1 (transition from 1 to 0)
  
  #### dynamic UI ####
  ####   boxinitialisation ####
  output$boxinitialisation<-renderUI({ #we put it in server side because eventually it is dynamic according to model specification by user
    tagList(
      checkboxGroupInput(inputId="initialisation",
                         label="Initialisation",
                         choices = rv$nodes,
                         selected=initat1),
    )
  })
  
  ####   boxoutbreaks ####
  output$boxoutbreaks<-renderUI({
    lapply(rv$nodes, function(i) {
      return(sliderInput(inputId=paste("outbreakrange", i, sep="_"),
                         label=i,
                         min=-1,
                         max=input$timesteps,
                         step=1,
                         value=c(-1,-1)))
    })
  })
  
  ####   box modeldefinition ####
  output$modeldefinition<-renderUI({
    ####   boxes effects on ####
    controls_list <- lapply(rv$transitions, function(i) {
      control_type <- "checkboxGroupInput"
      focalnode<-names(dimnames(i))[2]
      input_id <- paste("effectson", focalnode, sep="_")
      #message(paste("creation", input_id))
      choices <- names(rv$transitions)
      selected<-unlist(strsplit( names(dimnames(i))[1], "_"))
      labelinput<-paste("Select the nodes having an effect on", focalnode)
      return(checkboxGroupInput(input_id, label = labelinput, choices = choices, selected=selected))
    })
    tagList(
      textInput(inputId = "nodes", label="Node names", value=paste(names(rv$transitions), collapse=",")),
      tagList(controls_list),
      actionButton(inputId="updatemodel", label="Update model (not yet coded)")
    )
  })
  ####   boxParameterisation  ####
  output$boxParameterisation<-renderUI({
    controls_list <- lapply(rv$transitions, function(i) {
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
      actionButton(inputId="updateTransitions", label="Update model (not yet coded)"),
      tagList(controls_list)
      #put back the button here after debugging
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
                 for(n in names(rv$Transitions)){
                   subcontext<-context[focalnodes==n]
                   subprobatypes<-context[focalnodes==n]
                   newvalues<-controls[focalnodes==n]
                   rv$Transitions[[n]][subcontext, as.character(subprobatypes=="probah2l")]<<-newvalues
                 }
                 #NB once transitions have been changed to 3dimensionnal arrays, something like this will be easier
                 #rv$Transitions[focalnodes, context, as.character(probatypes=="probah2l")]<<-controls
               },
               ignoreInit = TRUE
  )
  
  #### Main plot (and running of model) ####
  output$plot1 <- renderPlot({
    initstate<-logical(length(rv$nodes)); names(initstate)<-rv$nodes
    initstate[names(initstate) %in% input$initialisation]<-1
    nbreps<-input$nbreps
    timesteps<-input$timesteps
    perturbation<-perturbationData()
    allsims<-runexpe(initstate=initstate, 
                     transitions=rv$transitions, 
                     timesteps=timesteps, 
                     nbreps=nbreps,
                     perturbation=perturbation)
    plotmeans(allsims)
    
  })
}

shinyApp(ui = ui, server = server)
