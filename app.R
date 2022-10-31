library(shiny)
library(Biostrings)
library(hutils)
library(shinyjs)
library(brochure)
library(shinyWidgets)
library(DT)
library(data.table)
library(dplyr)
library(msa)
library(seqinr)
library(ape)
library(dendextend)
library(shinyBS)
library(sortable)
library(rjson) 


########################################################
administrator_path <- '/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/';

database_folder <- paste(administrator_path,"database",'/',sep='');
gff_folder <- paste(administrator_path,"gff",'/',sep='');
script_folder <- paste(administrator_path,"script",'/',sep='');
User_folder <-paste(administrator_path,"User",'/',sep='');

#All_species<-list.files(database_folder)

#source(paste(administrator_path,"script/Species.R",sep=''), local = TRUE); 
#source(paste(administrator_path,"script/Page_species_info.R",sep=''), local = TRUE);

########################################################


############################################################ Creating a navlink 
nav_links <- tags$ul(


flowLayout(

  tags$li(
    tags$a(href = "/", "Main"), 
  ),

#   tags$li(
#   tags$a(href = "/Wheat", "Wheat"),
#  ),



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

        h1("Wellcome! This is TestShiny app main page",style="text-align:center"),
 #       h2("Subheading?"),
 #       h3("Subheading?"),
   

        nav_links,

# To add ui part for page0


mainPanel(width = 6, helpText("Add some instructions here?")),

#verbatimTextOutput("info")


      ) # For tagList
    }, # For ui function of page_0
    

# To add server function part for page0

server <- function(input, output, session){



#sometext<-c("References to be added here ???")

#output$info <- renderText({ sometext })


    } # server function of Page_0


  ) # page for Page_0

} # Page_0 function






############ To combine pages together

 brochureApp(

  page_0()

#  Species("Wheat",database_folder,gff_folder,script_folder,User_folder)
 
#  for(sp in All_species){Species(sp,database_folder,gff_folder,script_folder,User_folder)}



# To add many other pages

)



