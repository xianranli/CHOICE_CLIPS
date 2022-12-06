library(shiny)
library(Biostrings)
library(shinyjs)
library(brochure)
library(shinyWidgets)
library(DT)
library(data.table)
library(dplyr)
library(seqinr)
library(dendextend)
library(shinyBS)
library(sortable)
library(rjson)
library(shinythemes)


### 12/01/22
########################################################
administrator_path <- '/mnt/bridgecereal/';

database_folder <- paste(administrator_path,"database",'/',sep='');
gff_folder <- paste(administrator_path,"gff",'/',sep='');

script_folder <- paste(administrator_path,"script",'/',sep=''); # http://10.105.85.25/BRIDGEcereal/

User_folder <-paste(administrator_path,"User",'/',sep='');

#All_species<-list.files(database_folder)

source(paste(script_folder,"BRIDGEcereal_Sub.R",sep=''), local = TRUE); 

########################################################


############################################################ Creating a navlink 
nav_links <- tags$ul(


flowLayout(

  tags$li(
    tags$a(href = "/", "Main"), 
  ),

   tags$li(
   tags$a(href = "/Wheat", "Wheat"),
  ),

   tags$li(
   tags$a(href = "/Maize", "Maize"),
  ),

   tags$li(
   tags$a(href = "/Sorghum", "Sorghum"),
  ),

   tags$li(
   tags$a(href = "/Rice", "Rice"),
  ),

   tags$li(
   tags$a(href = "/Barley", "Barley"),
  ),



tags$style(
"li a {font-size: 20px;font-weight: bold;}",
)

)

)


page_0 <- function(){

  page(
    href = "/",


    ui <-  function(request){

      tagList(
        

       # fluidPage(theme = shinytheme("lumen")), ## cerulean, cosmo, cyborg, darkly, flatly, journal, lumen, paper, readable, sandstone, simplex, slate, spacelab, superhero, united, yeti 
        fluidPage(shinythemes::themeSelector()),

        h1("Wellcome! This is BRIDGEcereal app main page",style="text-align:center"),
 #       h2("Subheading?"),
 #       h3("Subheading?"),
   
        nav_links,

# To add ui part for page0
useShinyjs(),
#mainPanel(width = 6, helpText("Add some instructions here?")),
#verbatimTextOutput("info")


sidebarLayout(

sidebarPanel(


#column(8,
#pickerInput(
#  inputId = "Species_Gene", 
#  label = "Lookup your gene model used in BRIDGEcereal:", 
#  choices = c('',list.files(gff_folder)),
#  selected = c(''), ## by default

#  options = list(
#    'actions-box' = TRUE,
#    size = 30,
#    'selected-text-format' = "count > 1"
#  ), 
#  multiple = FALSE,
#)
#),

column(10, textAreaInput("Feedback", "Any question in using BRIDGEcereal to show your gene model ?", "Your gene model and questions... We will figure it out ASAP:", width = "1000px", height ="250px" )),
column(10, actionButton("Submit_Q", label = "Submit Your Questions",class = "btn-warning")),


), # sidebarPanel


mainPanel(

fluidRow(
  
#### progress in middle
tags$head(
    tags$style(
      HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }
           "
      )
    )
  ),
#### progress in middle


#column(10, align="center", imageOutput("Demo_Image")),
#column(12, offset = 0,DT::dataTableOutput("GeneModel1"),style='padding-top:5px; padding-bottom:5px'), 


#column(12, align="center", tags$iframe(width="1200", height="800", src="http://10.105.85.25/BRIDGEcereal/BRIDGEcereal_Instructions.pdf")),

column(6, align="left", tags$iframe(width="550", height="320", src="https://bridgecereal.scinet.usda.gov/BRIDGEcereal_Instructions1.pdf")),
column(6, align="right", tags$iframe(width="550", height="320", src="https://bridgecereal.scinet.usda.gov/BRIDGEcereal_Instructions2.pdf")),
column(6, align="left", tags$iframe(width="550", height="320", src="https://bridgecereal.scinet.usda.gov/BRIDGEcereal_Instructions3.pdf")),
column(6, align="right", tags$iframe(width="550", height="320", src="https://bridgecereal.scinet.usda.gov/BRIDGEcereal_Instructions4.pdf")),


column(12, align="center", tags$iframe(width="1200", height="800", src="https://bridgecereal.scinet.usda.gov/BRIDGEcereal/Demo.mp4", frameborder="0", allowfullscreen=NA)),



#uiOutput("video"),
    


) # fluidRow

) # mainPanel


) # sidebarLayout


      ) # For tagList
    }, # For ui function of page_0
    




# To add server function part for page0

server <- function(input, output, session){


observeEvent(input$Submit_Q, {

########### with progress ...
    progress <- Progress$new(session, min=1, max=5)
    on.exit(progress$close())
    progress$set(message = 'In progress ...',
                 detail = 'This may take a little while...')
    for (i in 1:5) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
####################################################
#output$User_Q <- renderText({ input$Feedback })

num_q <- sample(1:1000000,size=1, replace = FALSE);
# num_q <- as.numeric(read.table(paste(script_folder, "Visit_Times.txt",sep=''),header=F));

writeLines(input$Feedback , paste(script_folder,'User_Q','/',num_q,'_','User_question.txt',sep=''),sep="\n")

shinyjs::disable(id = "Submit_Q")
shinyjs::disable(id = "Feedback")

  })





    } # server function of Page_0


  ) # page for Page_0

} # Page_0 function



############ To combine pages together

 brochureApp(

  page_0(),

#  Species("Wheat","IWGSC",database_folder,gff_folder,script_folder,User_folder), # 'IWGSC' ... defined as default_ref
   Species("Maize","B73",database_folder,gff_folder,script_folder,User_folder)   # 'B73' ... defined as default_ref
#  Species("Sorghum","BTx623",database_folder,gff_folder,script_folder,User_folder), # 'BTx623' ... defined as default_ref
#  Species("Rice","Nipponbare",database_folder,gff_folder,script_folder,User_folder), # 'Nipponbare' ... defined as default_ref
#  Species("Barley","Morex",database_folder,gff_folder,script_folder,User_folder)     # 'Morex' ... defined as default_ref





# To add many other pages

)



