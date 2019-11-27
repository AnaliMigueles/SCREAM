source("191125_functions_SCREAM.R")

options(shiny.maxRequestSize = 30*1024^2)

### Outside variables that I need in case 
#Zebrafish
background_regen <- read.csv("backgroud.csv") #This DGE table can be found on the supplementary materials of Jiang et. al., 2014
single_fish <- read.csv("single_cell_version2.csv") #This list of markers can be found on the supplementary materials of Lush et. al., 2019
cluster_regen <- read.csv("ClusterA_postcut_manhattan_zebrafish_MASK.csv")


##Planaria
background_irrad <- as.data.frame(read.table("190717_DE_irrad_x_vs_WT.csv", sep = ",", header = TRUE)) #The DGE table can be found on the supplementary materials of SCREAM paper, Migueles-Lozano et.al, 2019
single_plan <- read.csv("180518_CorrectedCategoriesAtlas.csv") #This list of markers table can be found on the supplementary materials of SCREAM paper, Migueles-Lozano et.al, 2019
cluster_irrad <- read.csv(file ="191117_DOWN24HR.csv")


#####################################################################
### Template shynny app
##load library
library(shiny)
library(DT)
library(Rcpp)
library(httpuv)
library(RColorBrewer) 
library("statmod")
library("gplots")
library("markdown")
library("shinythemes")
library("shinycssloaders")
library("shinyjs")
#############

ui <- shinyUI(fluidPage( useShinyjs(),
                         
                         tags$style(HTML("
        input[type=number] {
              -moz-appearance:textfield;
        }
        input[type=number]::{
              -moz-appearance:textfield;
        }
        input[type=number]::-webkit-outer-spin-button,
        input[type=number]::-webkit-inner-spin-button {
              -webkit-appearance: none;
              margin: 0;
        }
    ")),
                         #shinythemes::themeSelector(), 
                         theme = shinytheme("lumen"),
                         h1(id="big-heading","Single Cell RNA-seq Enrichment Analysis Method: SCREAM",  align = "center"),
                         tags$style(HTML('#big-heading{color: darkblue;}', '#big-heading{font-size:60px;}', '#big_heading{font-weight: bold;}')), #"#big-heading{color: black;}", "#big-heading{font-size:60px;}", "#big_heading{font-weight: bold;}",
                         
                         
                         # Sidebar layout with input and output definitions ----
                         sidebarLayout(
                           
                           
                           # Sidebar panel for inputs ----
                           
                           
                           sidebarPanel(
                             
                             p(strong("1. Needed data")), 
                             p("-A single cell markers list"), 
                             p("-A full bulk RNA-seq gene list (background)"),
                             p("-A bulk RNA-seq cluster of interest gene list"),
                             p("To download examples of these tables click in any of these buttons below"),
                             column(3,downloadButton("down1", "Single cell list of marker", style='padding:3px; font-size:70%; width:180px; margin.right:1000px')),
                             column(3,downloadButton("down2", "Full Bulk RNA-seq gene list", style='padding:3px; font-size:70%; width:180px; margin.right:1000px')),
                             column(4,downloadButton("down3", "Bulk RNA-seq cluster", style='padding:3px; font-size:70%; width:180px; margin.right:1000px')),
                             br(),
                             br(),
                             br(),
                            
                              strong("Instruccions", style = "font-si50pt"),
                              p("-In case the single cell marker list have more information such as p value and/or log2FC enrichment columns, name those columns by", strong( "p_val"), "and",strong("avg_logFC")), 
                              p("-In a), b), c) sections, there are 3 options. One is", strong( "User file"),  "option and the remaining two are examples", strong("Example files")),  
                             p("-If you want to use the Example files, just select them in each panel and click GO! "), 
                              p("-After you select all your", strong("3 files, click"), strong("GO!"), "button"), 
                              p("-Every time you want to change the inputs,", strong ("go back to this panel"), ", load the respective files and click 'GO!' button"), 
                             br(),
                             br(),
                             
                            
                             
                             #####
                             fluidRow(
                               
                               column(6,div( selectInput("example1", "a) Upload the single cell markers list you want to use", 
                                                         choices=c("Planaria (Fincher CT et al, 2018)" = "planaria", "Zebrafish lateral line single cell (Lush ME et al., 2019)"="zebrafish","User files"),
                                                         selected=c("zebrafish"))))),
                             
                             # Input: Select a file ----
                             conditionalPanel(
                               condition = "input.example1 == 'User files'",
                               fileInput("file1", "Choose CSV Single Cell Markers File",
                                         multiple = FALSE,
                                         accept = c("text/csv",
                                                    "text/comma-separated-values,text/plain",
                                                    ".csv"))
                             ),
                             
                             conditionalPanel(
                               condition = "input.example1 == 'planaria'",
                               fluidRow(column(6,div( selectInput("tissues", "Tissues", 
                                                                  choices=c("brain","cathepsin+","ciliated neurons", "cluster3 neurons", 
                                                                            "epidermal", "intestine", "muscle", "neural", "neural+brain", "non-ciliated neurons",
                                                                            "parenchymal", "pharynx", "protonephridia", "smedwi-1 high", "smedwi-1+", "sexual"),
                                                                  selected=c("smedwi-1+")))))
                             ),
                             
                             
                             # Horizontal line ----
                             tags$hr(),
                             
                             fluidRow(column(6,div(selectInput("example2", "b) Upload your bulk RNA-seq genes of interests", 
                                                               choices=c("Planaria 24hr post irradiation (Cheng LC, 2018)"="planaria","Lateral line 1hr post neomycin (Jiang, 2014)"="zebrafish","User files"),
                                                               selected=c("zebrafish"))))),
                             
                             # Input: Select a file ----
                             conditionalPanel(
                               condition = "input.example2 == 'User files'", 
                               fileInput("file2", "Choose CSV Genes of interest File",
                                         multiple = FALSE,
                                         accept = c("text/csv",
                                                    "text/comma-separated-values,text/plain",
                                                    ".csv"))
                             ),
                             
                             
                             # Horizontal line ----
                             tags$hr(),
                             
                             fluidRow(column(6,div(selectInput("example3", "c) Upload all the genes in your bulk RNA-seq data (background)", 
                                                               choices=c("Planaria lethal irradiation (Cheng LC, 2018)"="planaria","Lateral line post neomycin (Jiang, 2014)"="zebrafish","User files"),
                                                               selected=c("zebrafish"))))),
                             # Input: Select a file ----
                             conditionalPanel(
                               condition = "input.example3 == 'User files'", 
                               fileInput("file3", "Choose CSV Background gene File",
                                         multiple = FALSE,
                                         accept = c("text/csv",
                                                    "text/comma-separated-values,text/plain",
                                                    ".csv"))
                             ),
                             
                             tags$hr(),
                             p("d) The 'Adjusted p value' is going to be used to filter the significant enriched cell population in 'Download list of genes'"),
                             numericInput(inputId = "pvalue",
                                          label = "Adjusted p value:",
                                          value = 0.05)
                             
                             
                           ),
                           
                           # Main panel for displaying outputs ----
                           mainPanel(
                             tags$head(
                               tags$style(type='text/css', 
                                          ".nav-tabs {font-size: 20px} ")),    
                             tabsetPanel(
                               id = 'dataset', 
                               #tabPanel("Instructions",
                               #        helpText("hola")),
                               
                               tabPanel("2) Enrichment", 
                                        tags$head(
                                          tags$style(HTML("
                    #goButton {
                    display:block;
                    height: 80px;
                    width: 120px;
                    font-size: 10px;
                    font-weight: bold;
                    }"))
                                        ),
                                        ##
                                        br(),
                                        br(),
                                        column(12, align="left",actionButton("prior", "Prior to analysis (Optional)", style = "width:400px")),
                                        hidden(
                                          div(id='text_div2',
                                              htmlOutput("text3"),
                                              column(4,numericInput(inputId = "pval",
                                                                    label = "Adjusted p value:", value=NULL 
                                              )),
                                              
                                              column(4,numericInput(inputId = "posFC",
                                                                    label = "positive log2 FC(enrichment):", value=NULL 
                                              )),
                                              
                                              column(4,numericInput(inputId = "negFC",
                                                                    label = "negative log2 FC(enrichment):", value=NULL 
                                              ))
                                              
                                          )
                                        ),
                                        ##
                                        br(),
                                        br(),
                                        
                                        #column(12, p("Click the", strong("GO!"), "button to start the analysis", style='font-size:200%')),
                                        column(12, align="center", p("Click the", strong ("GO!"), "button to start the analysis", style='font-size:200%')),
                                        
                                        column(12, align="center",actionButton("goButton", "GO!", style='padding:3px; font-size:200%; width:100px')),
                                        #p(strong("Calculate the enrichment (using hypergeometric test)"), style="font-size:200%"),
                                        
                                        #column(5,textOutput("text1")),
                                        #tags$head(tags$style("#text1{color: black;
                                        #               font-size: 25px;
                                        #              font-style: bold;
                                        #             }"
                                        #)
                                        #),
                                        #Action buttons for help (extra explication)        
                                        br(),
                                        br(),
                                        br(),
                                        tags$head(
                                          tags$style(HTML("
      .shiny-output-error-validation {
        color: red;
        font-size: 25px;
        font-style: bold;
      }
    "))
                                        ),
                                        
                                        htmlOutput("error1"),
                                        
                                        br(), 
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        
                                        #column(7,align="left",actionButton("question1", "",icon = icon("question"),style='padding:8px')),
                                        
                                        
                                        div(style = "display:inline-block;vertical-align:top;", 
                                            column(11,textOutput("text1")),
                                            tags$head(tags$style("#text1{color: black;
                                 font-size: 25px;
                                 font-style: bold;
                                 }"
                                            )
                                            ),
                                            
                                            #column(1,actionButton("question1", "",icon = icon("question"),style='padding:8px'))
                                            column(1,uiOutput("ques")),
                                            br(),
                                            br(),
                                            hidden(
                                              div(id='question-answer',
                                                  htmlOutput("answer")
                                              )
                                            )
                                            
                                        ),
                                        br(),
                                        br(),
                                        br(),
                                        
                                        column(width = 12,withSpinner(DT::dataTableOutput('mytable2'), type = 6)), #offset = 0,
                                        
                                        htmlOutput("text4"),
                                        br(),
                                        br(),
                                        uiOutput("down4"),
                                        br(),
                                        br(),
                                        htmlOutput("text5"),
                                        br(),
                                        uiOutput("down5")
                               ),
                               
                               tabPanel("3) Extras",
                                        tags$head(
                                          tags$style(HTML("
                  #up2 {
                    display:block;
                    height: 80px;
                    width: 120px;
                    font-weight: bold;
                    }

                    "))
                                        ),
                                        br(),
                                        br(),
                                        column(12, align="center", p("Click the", strong ("GO!"), "button to Load tables that cointain extra information", style='font-size:200%')),
                                        br(),
                                        br(),
                                        column(12, align="center", actionButton("up2", "GO!", style='padding:4px; font-size:200%')),
                                        br(),
                                        br(),
                                        #Add the dependency 
                                        #h3 ("Number of genes in each single cell cluster"),
                                        
                                        htmlOutput("extra1_text"),   
                                        br(),
                                        br(),
                                        #withSpinner(DT::dataTableOutput("extra1"), type=6),
                                        DT::dataTableOutput("extra1"),
                                        br(),
                                        br(),
                                        #column(12, DT::dataTableOutput("extra1")),
                                        
                                        #downloadButton("downloadData_extra1", "Download extras"),
                                        uiOutput("extra1_down"),
                                        br(),
                                        br(),
                                        #h3 ("Number of genes in each single cell cluster after only staying with genes both present in the single cell and bulk data"),
                                        htmlOutput("extra2_text"),
                                        br(),
                                        br(),
                                        #column(12,DT::dataTableOutput("extra2")),
                                        column(width = 12,withSpinner(DT::dataTableOutput('extra2'), type = 6)),
                                        
                                        #downloadButton("downloadData_extra2", "Download extras"),
                                        uiOutput("extra2_down"),
                                        br(),
                                        br(),
                                        #h3 ("Number of genes in each single cell cluster inside a bulk RNA-seq trend"),
                                        htmlOutput("extra3_text"),
                                        br(),
                                        br(),
                                        column(12,DT::dataTableOutput("extra3")),
                                        #downloadButton("downloadData_extra3", "Download extras"),
                                        uiOutput("extra3_down"),
                                        br(),
                                        br(),
                                        #h3 ("Enrichment Raw p-values"),
                                        htmlOutput("extra4_text"),
                                        br(),
                                        br(),
                                        column(12,DT::dataTableOutput("extra4")),
                                        #downloadButton("downloadData_extra4", "Download extras"),
                                        uiOutput("extra4_down"),
                                        br(),
                                        br(),
                                        #h3 ("Enrichment Adjusted p-values"),
                                        htmlOutput("extra5_text"),
                                        br(),
                                        br(),
                                        column(12,DT::dataTableOutput("extra5")),
                                        #downloadButton("downloadData_extra5", "Download extras"),
                                        uiOutput("extra5_down"),
                                        br(),
                                        br(),
                                        #column(12, p("We can test whether the overlap between two single cell clusters is significant enough to be enriched", style = "font-size:20px")),
                                        htmlOutput("extra6_6_text"),
                                        br(),
                                        br(),
                                        #column(12,h2("Overlap: Being aware of the overlap between single cell clusters can help you to decide whether to merge populations or shrink your makers list by p values or log2FC")),
                                        htmlOutput("heatmap_text"),
                                        br(),
                                        br(),
                                        column(12,plotOutput('overlap_heat')),
                                        br(),
                                        br(),
                                        #h3 ("Overlap: Number of shared genes between single cell clusters (after the filtering)"),
                                        htmlOutput("extra6_6_text2"),
                                        br(),
                                        br(),
                                        column(12,DT::dataTableOutput("extra6_6")),
                                        #downloadButton("downloadData_extra6_6", "Download extras"),
                                        uiOutput("extra6_6_down"),
                                        br(),
                                        br(),
                                        #h3 ("Overlap: Number of shared genes between single cell clusters inside the bulk RNA-seq cluster"),
                                        htmlOutput("extra6_text"),
                                        br(),
                                        br(),
                                        column(12,DT::dataTableOutput("extra6")),
                                        #downloadButton("downloadData_extra6", "Download extras"),
                                        uiOutput("extra6_down"),
                                        br(),
                                        br(),
                                        #h3 ("Overlap: Enrichment raw p-values due to overlap "),
                                        htmlOutput("extra7_text"),
                                        br(),
                                        br(),
                                        column(12,DT::dataTableOutput("extra7")),
                                        #downloadButton("downloadData_extra7", "Download extras"),
                                        uiOutput("extra7_down"),
                                        br(),
                                        br(),
                                        #h3 ("Overlap: Enrichment adjusted p-values due to overlap"),
                                        htmlOutput("extra8_text"),
                                        br(),
                                        br(),
                                        column(12,DT::dataTableOutput("extra8")),
                                        #downloadButton("downloadData_extra8", "Download extras")
                                        uiOutput("extra8_down")
                                        
                               ),
                               
                               tabPanel("Contact",
                                        p(strong("Anali Migueles Lozano")),
                                        p(strong("Alejandro Sánchez Lab")),
                                        p(strong("Stowers Institute for Medical Research")),
                                        p("1000 E 50th St, Kansas City, MO 64110"),
                                        p(strong("Email address:"), "amigueleslozano@stowers.org, anamigueeles@gmail.com")
                               ) 
                             )
                             
                           )#mainPanel brachet
                           
                         )# end sidebarLayout
)#end of fluid page
)#end of shinyUI



server <- function(input, output, session) {
  
  
  single_cell_input <- eventReactive( input$goButton, {
    
    isolate(switch(input$example1,
                   "planaria" = single_plan,
                   "zebrafish" = single_fish,
                   "User files" = read.csv(input$file1$datapath,
                                           sep = ",")))
  } ) #, ignoreNULL = FALSE
  
  
  tissue_input <- eventReactive( input$goButton, {
    
    single<-isolate(single_cell_input())
    single<-single[single$cell_type %in% input$tissues,]
    
  } ) #, ignoreNULL = FALSE
  
  
  clusterx_input <- eventReactive( input$goButton, {
    
    isolate(switch(input$example2,
                   "planaria" = cluster_irrad,
                   "zebrafish" = cluster_regen,
                   "User files" = read.csv(input$file2$datapath,
                                           sep = ",")))
  })
  background_input <- eventReactive( input$goButton, {
    
    isolate(switch(input$example3,
                   "planaria" = background_irrad,
                   "zebrafish" = background_regen,
                   "User files" = read.csv(input$file3$datapath,
                                           sep = ",")))
  })
  
  single_cell_markers<-eventReactive (input$goButton,{
    single_cell<-isolate(single_cell_input())
    a<- input$pval#NULL#pval
    c<- input$negFC #negfc
    d<- input$posFC#2#posfc
    
    #validate(need(is.na(a), "a no es na"))
    #validate(need(is.na(c), "c no es na"))
    #validate(need(is.na(d), "d no es na"))
    
    if (!is.na(a) & !is.na(d) & !is.na(c)) {   ### merge and shrink
      single_cell <- subset(single_cell, (avg_logFC <=c | avg_logFC >= d) & p_val <= a) 
      
    } else if(!is.na(a) & !is.na(d)) {
      single_cell <- subset(single_cell, p_val <=  a & avg_logFC >= d)
      
    } else if(!is.na(a) & !is.na(c)){
      single_cell <- subset(single_cell, p_val <=  a & avg_logFC <=c)
      
    } else if(!is.na(a)){
      single_cell <- subset(single_cell, p_val <=  a )
      
    } else if(!is.na(c) & !is.na(d) ){
      single_cell <- subset(single_cell, avg_logFC <= c | avg_logFC >= d)
      
    } else if(!is.na(c)){
      single_cell <- subset(single_cell, avg_logFC <= c )
      
    } else if(!is.na(d)){
      single_cell <- subset(single_cell, avg_logFC >= d)
    } else{single_cell} #end of oval and fold change
    
  })
  
  ####Create the reactive objects of the constant elements such as background, single cell and cluster. 
  background<-reactive({
    
    input$goButton
    single_cell<-single_cell_markers() #todo
    background<-background_input() #todo
    background<-background[background$ID %in% single_cell$ID, ]
    background
  })
  
  #single_cell_tissue<-tissue_input[tissue_input$ID %in% single_cell$ID, ] 
  observeEvent(input$goButton, {
    #toggle('text_div')
    output$text1 <- renderText({"Calculate the enrichment (using hypergeometric test)"})
    output$text4 <- renderUI({ 
      str1 <- paste("In this summary, it is presented the enrichment of each single cell cluster and which other cluster can be enriched due to the shared genes between.") 
      HTML(paste("<b>",str1,"</b>", sep = '<br/>'))
    }) #close renderUI
    output$text5 <- renderUI({ 
      str1 <- paste("In the gene list it is only going to be considered the single cell clusters that passed the", strong("Adjusted p-value"), "indicated by the user in the left panel.") 
      HTML(paste("<b>",str1,"</b>", sep = '<br/>' ))
    })
    
    output$down4 <- renderUI(downloadButton("downloadData", "Download the summary", style = "width:400px"))
    output$down5 <- renderUI(downloadButton("genes", "Download list of genes"))
    output$ques<-renderUI(actionButton("question1", "",icon = icon("question"),style='padding:8px'))
    #column(1,actionButton("question1", "",icon = icon("question"),style='padding:8px'))
    #column(1,uiOutput("ques"))
    
  })
  
  observeEvent(input$question1, {
    toggle('question-answer')
    output$answer<-renderUI({
      str1<-paste("<font size=3 color=blue> Single cell cluster :</font> In this column are present all single cell clusters inside the single cell marker list")
      str2<-paste("<font size=3 color=blue> Enrichment (adj pvalue):</font> In this column there are all the resulting p values after applying the hypergeometric test to the data sets. Note: If you click <font size=3 color=red> Download list of genes </font>, only the single cell cluster that passed the specified adjusted p value will be considere in the file")
      str3<-paste("<font size=3 color=blue> False positive enrichment:</font> In this column are the single cell clusters which its significant enrichment can be due to markers overlap with other clusters.")
      str4<-paste("<font size=3 color=blue> Note:</font> If cluster A is both in this column and Enrichment (adj pvalue), that suggest no redundancies in the enrichment of this cluster. If cluster A is enriched in Enrichment (adj pvalue) column but is not enriched in this one, that is due to p values correction effect when false postive enrichments are calculated. These are independent from the Enrichment(p value). The significance is determined by the specified adj pvalue. If in this column there are multiple single cell clusters mentioned, it suggests the redundancy of enrichments between clusters") 
      
      HTML(paste("<b>",str1, "<br>",str2, "<br>",str3,"<br>",str4,"</b>"))  
    })
  })
  
  #observeEvent(input$instructions, {
   # toggle('text_div')
    #output$text2 <- renderUI({ 
     # str3 <- paste("-In case the single cell marker list have more information such as p value and/or log2FC enrichment columns, name those by p_val and avg_logFC") 
      #str4 <- paste("-After you select all your 3 files, click 'GO!' button") 
      #str5 <- paste("-Every time you want to change the inputs, go back to this panel, load the respective files and click 'GO!' button") 
      #str6 <- paste("-In a), b), c) sections, there are 3 options. Excluding the User file options, the remaining two are the Example files used in ....") 
      #str7 <- paste("-If you want to use the Example files, just select them and click GO! ") 
      #HTML(paste("<b>",str1,"<br>", str2,"<br>", str3, "<br>",str4, "<br>",str5, "<br>",str6, "<br>", str7, sep = '<br/>'))
      #HTML(paste("-All 3 inputs must have in common 1 column name: <font size=2 color=red> ID </font>","-Single cell marker list MUST have another column named:<font size=2 color=red> cluster </font>",
      #           str3, "<br>", "<b>",str4, str5, str6, str7, "</b>", sep = '<br/>' ))
  #  }) #close renderUI
  #})
  
  observeEvent(input$prior, {
    toggle('text_div2')
    output$text3 <- renderUI({ 
      str1 <- paste("Prior to the analysis, you can extra filter the single cell list of markers based on p value and FC. If not, click GO! to start the analysis") 
      HTML(paste("<b>",str1,"</b>", sep = '<br/>' ))
    }) #close renderUI
  })
  
  
  
  ###  creamos dos single cell (uno antes y otro despues del corte con el bulk-RNA seq)
  single_cell_prior<-reactive({
    input$goButton
    single_cell<-isolate(single_cell_markers()) ##
    if(input$example1=="planaria") {
      tissue<-tissue_input()
      single_cell<-tissue[tissue$ID %in% single_cell$ID,]
    } else {single_cell<-single_cell_markers()}
    
  })
  
  single_cell_post<-reactive({
    input$goButton
    background<-background()
    single_cell<-single_cell_markers() ##
    
    if(input$example1=="planaria") {
      tissue<-tissue_input()
      single_cell<-tissue[tissue$ID %in% single_cell$ID,]
    } else {single_cell<-single_cell_markers()}
    
    single_cell<-single_cell[single_cell$ID %in% background$ID,]
    
  })
  
  clusterx<-reactive({
    input$goButton
    clusterx<-clusterx_input() #todo
    background<-background() #todo
    clusterx<-clusterx[clusterx$ID %in% background$ID, ] #todo (no es tissue dependent) y esto es oslo para reducir el cluster 
    clusterx
  })
  
  
  output$error1<- renderUI({
    validate(need(dim(clusterx())[1]!=0, "No genes in the bulk RNA-seq cluster is present in the single cell marker list with this parameters. Enter another list/parameters"))
  })
  
  #if (nrow(clusterx())==0) {
  ###Aqui poner un mensaje de error si el cluster completo desaparece cuando se compara con el single cell y desaparece. 
  #}
  
  clusterx_single<-reactive({
    input$goButton
    single_cell<-single_cell_post() ##
    clusterx<-clusterx() #todo
    background<-background() #todo
    clusterx_single<-single_cell[single_cell$ID %in% clusterx$ID, ] #todo (no es tissue dependent) y esto es oslo para reducir el cluster 
    
  })
  
  ### New reactive elements
  #background()
  #single_cell()
  #clusterx()
  
  data<-reactive({    
    pvalue<-input$pvalue
    single_cell<-single_cell_post() ## aqui tiene que ser el cambio
    clusterx<- clusterx() #todo
    background<-background() #todo
    clusterx_single<-clusterx_single() #single cell & cluster #aqui se necesita poner el tissue input
    ##########
    
    list_clouds1<-data.frame(table(single_cell$cluster))
    list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
    #sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1]))
    summary_general<-summary_background(single_cell, list_clouds1,output=1)
    
    list_clouds2<-data.frame(table(clusterx_single$cluster))
    list_clouds2<-sort(as.character(list_clouds2$Var1))
    summary_cluster<-summary_background(clusterx_single, list_clouds1,output=1)
    summary_cluster_binary<-summary_background(clusterx_single, list_clouds1,output=2)
    
    enrichment2<-enrichment(background, summary_general, summary_cluster, list_clouds1, clusterx)
    
    ################ 
    #enrichment_logical<-as.matrix((enrichment2))
    enrichment_logical<-enrichment2
    a<-dim(enrichment_logical)[1]
    k<-vector()
    
    for (i in 1:a) {
      k<-enrichment_logical[i,]
      k<-p.adjust(k )
      enrichment_logical[i,]<-k
    }
    
    enrichment_logical_adj<-enrichment_logical
    
    obj<-overlap(binary_m =  summary_cluster_binary, summary_cluster =  summary_cluster, summary_general =  summary_general, background = background, clusterx=clusterx, pvalue = pvalue, enrichment_logical_adj =enrichment_logical_adj, list_clouds = list_clouds1, output=1) #Tabla final
    obj<-obj[,-1]
    colnames(obj)<-c("Single cell cluster", "Enrichment (adj pvalue)", "False Positive enrichment"  )
    obj
    
    
  }) #close of reactive 
  
  list_genes<- reactive({
    input$goButton
    enrichment<-isolate(data())
    pvalue<-isolate(input$pvalue)
    
    single_cell<-single_cell_post() 
    clusterx<-clusterx()
    background<-background()
    
    clusterx_single <- clusterx_single() #single cell & cluster
    
    #list_enrich<-subset(enrichment, as.numeric(enrichment[,3])<=pvalue, select = Population)
    list_enrich<-subset(enrichment, as.numeric(enrichment[,2])<=pvalue, select = c("Single cell cluster"))
    clusterx_single[clusterx_single$cluster %in% list_enrich[,1],]
    
    #if(input$example1 == 'planaria') {
    # #lets merge the new ids con el resto de la tabla
    #trans_table<-trans()
    #data_merge<-merge(clusterx_single, trans(), by="ID")
    #data_merge
    #} else(clusterx_single)
    
  })
  
  
  
  output$mytable3  = DT::renderDataTable({  
    input$goButton
    single_cell<-isolate(single_cell_prior())
    list_clouds1<-data.frame(table(single_cell$cluster))
    list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
    summary_general<-t(summary_background(single_cell, list_clouds1,output=1))
    
  })
  
  output$mytable4  = DT::renderDataTable({   
    input$goButton
    single_cell<-isolate(single_cell_post())
    list_clouds1<-data.frame(table(single_cell$cluster))
    list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
    summary_general<-t(summary_background(single_cell, list_clouds1,output=1))
  })
  
  output$mytable5 = DT::renderDataTable({   
    input$goButton
    pvalue<-isolate(input$pvalue)
    
    single_cell<-isolate(single_cell_post())
    clusterx<-isolate(clusterx())
    background<-isolate(background())
    clusterx_single<-isolate(clusterx_single())  #single cell & cluster
    ###########
    
    list_clouds1<-data.frame(table(single_cell$cluster))
    list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
    summary_general<-summary_background(single_cell, list_clouds1,output=1)
    
    list_clouds2<-data.frame(table(clusterx_single$cluster))
    list_clouds2<-sort(as.character(list_clouds2$Var1))
    summary_cluster<-t(summary_background(clusterx_single, list_clouds1,output=1))
    
    
  })
  #########  
  extra11<-eventReactive(input$up2,{
    #input$up2
    isolate({single_cell<-isolate(single_cell_prior())
    list_clouds1<-data.frame(table(single_cell$cluster))
    list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
    summary_general<-t(summary_background(single_cell, list_clouds1,output=1))
    sum1<-as.data.frame(summary_general)
    sum1$Singlecellcluster<-rownames(sum1)
    sum1<-sum1[c(2,1)]
    colnames(sum1)<-c("Single cell clusters", "Number of genes in each single cell cluster ")
    sum1
    })
  })
  
  extra22<-eventReactive(input$up2,{
    #input$up2
    isolate({single_cell<-isolate(single_cell_post())
    list_clouds1<-data.frame(table(single_cell$cluster))
    list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
    summary_general<-t(summary_background(single_cell, list_clouds1,output=1))
    sum1<-as.data.frame(summary_general)
    sum1$Singlecellcluster<-rownames(sum1)
    sum1<-sum1[c(2,1)]
    colnames(sum1)<-c("Single cell clusters", "Number of genes in each single cell cluster ")
    sum1
    })
  })
  
  extra33<-eventReactive(input$up2,{
    #input$up2
    isolate({single_cell<-isolate(single_cell_post())
    clusterx_single<-isolate(clusterx_single())
    ###########
    list_clouds1<-data.frame(table(single_cell$cluster))
    list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
    summary_cluster<-t(summary_background(clusterx_single, list_clouds1,output=1))
    sum1<-as.data.frame(summary_cluster)
    sum1$Singlecellcluster<-rownames(sum1)
    sum1<-sum1[c(2,1)]
    colnames(sum1)<-c("Single cell clusters", "Number of genes in each single cell cluster ")
    sum1})
  })
  
  extra44<-eventReactive(input$up2,{
    #input$up2
    isolate({ pvalue<-isolate(input$pvalue)
    single_cell<-isolate(single_cell_post())
    clusterx<-isolate(clusterx())
    background<-isolate(background())
    clusterx_single<-isolate(clusterx_single()) #single cell & cluster
    
    ###########
    list_clouds1<-data.frame(table(single_cell$cluster))
    list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
    summary_general<-summary_background(single_cell, list_clouds1,output=1)
    
    list_clouds2<-data.frame(table(clusterx_single$cluster))
    list_clouds2<-sort(as.character(list_clouds2$Var1))
    summary_cluster<-summary_background(clusterx_single, list_clouds1,output=1)
    summary_cluster_binary<-summary_background(clusterx_single, list_clouds1,output=2)
    
    enrichment2<-t(enrichment(background, summary_general, summary_cluster, list_clouds1, clusterx))
    sum1<-as.data.frame(enrichment2)
    sum1$Singlecellcluster<-rownames(sum1)
    sum1<-sum1[c(2,1)]
    colnames(sum1)<-c("Single cell clusters", "Raw p value")
    sum1
    })
  })
  extra55<-eventReactive(input$up2,{
    #input$up2
    isolate({
      pvalue<-isolate(input$pvalue)
      
      single_cell<-isolate(single_cell_post())
      clusterx<-isolate(clusterx())
      background<-isolate(background())
      clusterx_single<-isolate(clusterx_single()) #single cell & cluster
      
      ###########
      list_clouds1<-data.frame(table(single_cell$cluster))
      list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
      summary_general<-summary_background(single_cell, list_clouds1,output=1)
      
      list_clouds2<-data.frame(table(clusterx_single$cluster))
      list_clouds2<-sort(as.character(list_clouds2$Var1))
      summary_cluster<-summary_background(clusterx_single, list_clouds1,output=1)
      summary_cluster_binary<-summary_background(clusterx_single, list_clouds1,output=2)
      
      enrichment2<-enrichment(background, summary_general, summary_cluster, list_clouds1, clusterx)
      
      ################ 
      #enrichment_logical<-as.matrix((enrichment2))
      enrichment_logical<-enrichment2
      a<-dim(enrichment_logical)[1]
      k<-vector()
      
      for (i in 1:a) {
        k<-enrichment_logical[i,]
        k<-p.adjust(k )
        enrichment_logical[i,]<-k
      }
      
      enrichment_logical_adj<-t(enrichment_logical)
      sum1<-as.data.frame(enrichment_logical_adj)
      sum1$Singlecellcluster<-rownames(sum1)
      sum1<-sum1[c(2,1)]
      colnames(sum1)<-c("Single cell clusters", "Adjusted p value")  
      sum1
    })
  })
  
  
  extra66<-eventReactive(input$up2,{
    #input$up2
    isolate({
      pvalue<-isolate(input$pvalue)
      
      single_cell<-isolate(single_cell_post())
      clusterx<-isolate(clusterx())
      background<-isolate(background())
      clusterx_single<-isolate(clusterx_single()) #single cell & cluster
      ###########
      
      list_clouds1<-data.frame(table(single_cell$cluster))
      list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
      summary_general<-summary_background(single_cell, list_clouds1,output=1)
      
      list_clouds2<-data.frame(table(clusterx_single$cluster))
      list_clouds2<-sort(as.character(list_clouds2$Var1))
      summary_cluster<-summary_background(clusterx_single, list_clouds1,output=1)
      summary_cluster_binary<-summary_background(clusterx_single, list_clouds1,output=2)
      
      enrichment2<-enrichment(background, summary_general, summary_cluster, list_clouds1, clusterx)
      
      ################ 
      #enrichment_logical<-as.matrix((enrichment2))
      enrichment_logical<-enrichment2
      a<-dim(enrichment_logical)[1]
      k<-vector()
      
      for (i in 1:a) {
        k<-enrichment_logical[i,]
        k<-p.adjust(k )
        enrichment_logical[i,]<-k
      }
      
      enrichment_logical_adj<-enrichment_logical
      
      overlap(binary_m =  summary_cluster_binary, summary_cluster =  summary_cluster, summary_general =  summary_general, background = background, clusterx=clusterx, pvalue = pvalue, enrichment_logical_adj =enrichment_logical_adj, list_clouds = list_clouds1, output=2)
    })  
  })
  extra6_61<- eventReactive(input$up2,{
    #input$up2
    isolate({
      pvalue<-isolate(input$pvalue)
      
      single_cell<-isolate(single_cell_post())
      clusterx<-isolate(clusterx())
      background<-isolate(background())
      clusterx_single<-isolate(clusterx_single()) #single cell & cluster
      ###########
      
      list_clouds1<-data.frame(table(single_cell$cluster))
      list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
      summary_general<-summary_background(single_cell, list_clouds1,output=1)
      
      list_clouds2<-data.frame(table(clusterx_single$cluster))
      list_clouds2<-sort(as.character(list_clouds2$Var1))
      summary_cluster<-summary_background(clusterx_single, list_clouds1,output=1)
      summary_cluster_binary<-summary_background(clusterx_single, list_clouds1,output=2)
      summary_general_binary<-summary_background(single_cell, list_clouds1,output=2)
      
      enrichment2<-enrichment(background, summary_general, summary_cluster, list_clouds1, clusterx)
      
      ################ 
      #enrichment_logical<-as.matrix((enrichment2))
      enrichment_logical<-enrichment2
      a<-dim(enrichment_logical)[1]
      k<-vector()
      
      for (i in 1:a) {
        k<-enrichment_logical[i,]
        k<-p.adjust(k )
        enrichment_logical[i,]<-k
      }
      
      enrichment_logical_adj<-enrichment_logical
      
      overlap(binary_m =  summary_general_binary, summary_cluster = summary_general, summary_general =  summary_general, background = background, clusterx=background, pvalue = pvalue, enrichment_logical_adj =enrichment_logical_adj, list_clouds = list_clouds1, output=2) #overlap
    })
  })
  
  extra77<-eventReactive(input$up2,{
    #input$up2
    isolate({
      pvalue<-isolate(input$pvalue)
      
      single_cell<-isolate(single_cell_post())
      clusterx<-isolate(clusterx())
      background<-isolate(background())
      clusterx_single<-isolate(clusterx_single()) #single cell & cluster
      ###########
      
      list_clouds1<-data.frame(table(single_cell$cluster))
      list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
      summary_general<-summary_background(single_cell, list_clouds1,output=1)
      
      list_clouds2<-data.frame(table(clusterx_single$cluster))
      list_clouds2<-sort(as.character(list_clouds2$Var1))
      summary_cluster<-summary_background(clusterx_single, list_clouds1,output=1)
      summary_cluster_binary<-summary_background(clusterx_single, list_clouds1,output=2)
      
      enrichment2<-enrichment(background, summary_general, summary_cluster, list_clouds1, clusterx)
      
      ################ 
      #enrichment_logical<-as.matrix((enrichment2))
      enrichment_logical<-enrichment2
      a<-dim(enrichment_logical)[1]
      k<-vector()
      
      for (i in 1:a) {
        k<-enrichment_logical[i,]
        k<-p.adjust(k )
        enrichment_logical[i,]<-k
      }
      
      enrichment_logical_adj<-enrichment_logical
      
      overlap(binary_m =  summary_cluster_binary, summary_cluster =  summary_cluster, summary_general =  summary_general, background = background, clusterx=clusterx, pvalue = pvalue, enrichment_logical_adj =enrichment_logical_adj, list_clouds = list_clouds1, output=3) 
    })
  })
  extra88<-eventReactive(input$up2,{
    #input$up2
    isolate({
      pvalue<-isolate(input$pvalue)
      
      single_cell<-isolate(single_cell_post())
      clusterx<-isolate(clusterx())
      background<-isolate(background())
      clusterx_single<-isolate(clusterx_single()) #single cell & cluster
      ###########
      
      list_clouds1<-data.frame(table(single_cell$cluster))
      list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
      summary_general<-summary_background(single_cell, list_clouds1,output=1)
      
      list_clouds2<-data.frame(table(clusterx_single$cluster))
      list_clouds2<-sort(as.character(list_clouds2$Var1))
      summary_cluster<-summary_background(clusterx_single, list_clouds1,output=1)
      summary_cluster_binary<-summary_background(clusterx_single, list_clouds1,output=2)
      
      enrichment2 <- enrichment(background, summary_general, summary_cluster, list_clouds1, clusterx)
      
      ################ 
      #enrichment_logical<-as.matrix((enrichment2))
      enrichment_logical<-enrichment2
      a<-dim(enrichment_logical)[1]
      k<-vector()
      
      for (i in 1:a) {
        k<-enrichment_logical[i,]
        k<-p.adjust(k )
        enrichment_logical[i,]<-k
      }
      
      enrichment_logical_adj<-enrichment_logical
      
      overlap(binary_m = summary_cluster_binary, summary_cluster =  summary_cluster, summary_general =  summary_general, background = background, clusterx=clusterx, pvalue = pvalue, enrichment_logical_adj =enrichment_logical_adj, list_clouds = list_clouds1, output=4)
    })  
  })
  ############  <font size=3 color=blue> Single cell cluster :</font>
  
  output$extra1  <- DT::renderDataTable ({   
    input$up2
    isolate(extra11())
  }, rownames = FALSE)
  
  output$extra1_text <- renderUI({ 
    input$up2  
    str1 <- eventReactive(input$up2,{paste("<font size=5 >Total number of genes in each single cell cluster (from original data)</font>")}) 
    HTML(paste("<b>",str1(),"</b>"))
  }) #close renderUI
  
  output$extra2_text <- renderUI({ 
    input$up2  
    str1 <- eventReactive(input$up2,{paste("<font size=5 >Number of genes in each single cell cluster after only staying with genes both present in the single cell and bulk data</font>")}) 
    HTML(paste("<b>",str1(),"</b>"))
  }) #close renderUI
  
  output$extra2  <- DT::renderDataTable({   
    input$up2
    isolate(extra22())
  }, rownames = FALSE)
  
  output$extra3 <- DT::renderDataTable({ 
    input$up2
    isolate(extra33())
  }, rownames = FALSE)
  
  output$extra3_text <- renderUI({ 
    input$up2  
    str1 <- eventReactive(input$up2,{paste("<font size=5 >Number of genes in each single cell cluster inside your bulk RNA-seq cluster</font>")}) 
    HTML(paste("<b>",str1(),"</b>"))
  }) #close renderUI
  
  
  output$extra4 <- DT::renderDataTable({ 
    input$up2
    isolate(extra44())
  }, rownames = FALSE)
  
  output$extra4_text <- renderUI({ 
    input$up2  
    str1 <- eventReactive(input$up2,{paste("<font size=5 >Enrichment Raw p-values</font>")}) 
    HTML(paste("<b>",str1(),"</b>"))
  }) #close renderUI
  
  output$extra5 <- DT::renderDataTable({ 
    input$up2
    isolate(extra55())
  }, rownames = FALSE)
  
  output$extra5_text <- renderUI({ 
    input$up2  
    str1 <- eventReactive(input$up2,{paste("<font size=5 > Enrichment Adjusted p-values </font>")}) 
    HTML(paste("<b>",str1(),"</b>"))
  }) #close renderUI
  
  output$extra6 <- DT::renderDataTable({ 
    input$up2
    isolate(extra66())
  }, rownames = FALSE)
  
  output$extra6_text <- renderUI({ 
    input$up2  
    str1 <- eventReactive(input$up2,{paste("<font size=5 > Overlap: Number of shared genes between single cell clusters inside the bulk RNA-seq cluster </font>")}) 
    HTML(paste("<b>",str1(),"</b>"))
  }) #close renderUI
  
  output$extra6_6 <- DT::renderDataTable({ 
    input$up2
    isolate(extra6_61())
  }, rownames = FALSE)
  
  output$extra6_6_text2 <- renderUI({ 
    input$up2  
    str1 <- eventReactive(input$up2,{paste("<font size=5 > Overlap: Number of shared genes between single cell clusters (after the filtering) </font>")}) 
    HTML(paste("<b>",str1(),"</b>"))
  }) #close renderUI
  
  
  output$extra7 <- DT::renderDataTable({ 
    input$up2
    isolate(extra77())
  }, rownames = FALSE)
  
  output$extra7_text <- renderUI({ 
    input$up2  
    str1 <- eventReactive(input$up2,{paste("<font size=5 > Overlap: Enrichment raw p-values due to overlap </font>")}) 
    HTML(paste("<b>",str1(),"</b>"))
  }) #close renderUI
  
  output$extra8 <- DT::renderDataTable({ 
    input$up2
    isolate(extra88())
  }, rownames = FALSE)
  
  output$extra8_text <- renderUI({ 
    input$up2  
    str1 <- eventReactive(input$up2,{paste("<font size=5 > Overlap: Enrichment adjusted p-values due to overlap </font>")}) 
    HTML(paste("<b>",str1(),"</b>"))
  }) #close renderUI
  
  output$heatmap_text <- renderUI({ 
    str1 <- eventReactive(input$up2,{paste("<font size=5 > We can test whether the overlap between two single cell clusters is significant enough to be enriched </font>")}) 
    HTML(paste("<b>",str1(),"</b>"))
  }) #close renderUI
  
  output$overlap_heat <- renderPlot({
    
    map <-eventReactive(input$up2,{
      
      pvalue <-input$pvalue
      single_cell <- single_cell_post() ## aqui tiene que ser el cambio
      clusterx <- clusterx() #todo
      background <- background() #todo
      clusterx_single <- clusterx_single() #single cell & cluster #aqui se necesita poner el tissue input
      
      ##########
      
      list_clouds1<-data.frame(table(single_cell$cluster))
      list_clouds1<-sort(as.character(subset(list_clouds1, Freq!=0, Var1)[,1])) #Listclouds1 is cells_full
      summary_general<-summary_background(single_cell, list_clouds1,output=1)
      
      list_clouds2<-data.frame(table(clusterx_single$cluster))
      list_clouds2<-sort(as.character(list_clouds2$Var1))
      summary_cluster<-summary_background(clusterx_single, list_clouds1,output=1)
      summary_cluster_binary<-summary_background(clusterx_single, list_clouds1,output=2)
      summary_general_binary<-summary_background(single_cell, list_clouds1,output=2)
      
      enrichment2<-enrichment(background, summary_general, summary_cluster, list_clouds1, clusterx)
      
      ################ 
      #enrichment_logical<-as.matrix((enrichment2))
      enrichment_logical<-enrichment2
      a<-dim(enrichment_logical)[1]
      k<-vector()
      
      for (i in 1:a) {
        k<-enrichment_logical[i,]
        k<-p.adjust(k )
        enrichment_logical[i,]<-k
      }
      
      enrichment_logical_adj<-enrichment_logical
      
      #over_heatmap <-overlap(binary_m =  summary_cluster_binary, summary_cluster =  summary_cluster, summary_general =  summary_general, background = background, clusterx=clusterx, pvalue = pvalue, enrichment_logical_adj =enrichment_logical_adj, list_clouds = list_clouds1, output=5)
      
      over_heatmap <-overlap(binary_m =  summary_general_binary, summary_cluster = summary_general, summary_general =  summary_general, background = background, clusterx=background, pvalue = pvalue, enrichment_logical_adj =enrichment_logical_adj, list_clouds = list_clouds1, output=5) #overlap
      
      Colors=colorRampPalette(rev(brewer.pal(n = 7, name ="YlGnBu")))(200)
      Breaks=seq(0,100,0.5)
      heatmap.2(over_heatmap, col=rev(Colors),distfun = function(x) dist(x,method = 'euclidean'), scale = "none", trace = "none", main = "Overlap between single cell clusters", margins = c(10,30), cexCol = 1, cexRow= 1,breaks = Breaks, key.title = "Percentage of overlap")
      
    })
    map()
  })
  
  
  ##Buttons for download
  
  observeEvent (input$up2, {
    output$extra1_down <- renderUI({
      downloadButton("downloadData_extra1", "Download extras")
    })
    
    output$extra2_down <- renderUI({
      downloadButton("downloadData_extra2", "Download extras")
    })
    output$extra3_down <- renderUI({
      downloadButton("downloadData_extra3", "Download extras")
    })
    
    output$extra4_down <- renderUI({
      downloadButton("downloadData_extra4", "Download extras")
    })
    
    output$extra5_down <- renderUI({
      downloadButton("downloadData_extra5", "Download extras")
    })
    
    output$extra6_down <- renderUI({
      downloadButton("downloadData_extra6", "Download extras")
    })
    
    output$extra6_6_down <- renderUI({
      downloadButton("downloadData_extra6_6", "Download extras")
    })
    output$extra7_down <- renderUI({
      downloadButton("downloadData_extra7", "Download extras")
    })
    output$extra8_down <- renderUI({
      downloadButton("downloadData_extra8", "Download extras")
    })
    
  })
  
  
  
  
  
  ############  
  #output$mytable1 = DT::renderDataTable({   
  # input$goButton
  #isolate(background())
  #}) #Closing $mytable1   
  
  
  output$mytable2 = DT::renderDataTable({   
    input$goButton
    isolate(data())
    
  }) #Closing $mytable1   
  
  
  
  ##Download the example tables
  output$down1 <- downloadHandler(
    filename = function() {
      paste("Example_Single_cell_list_of_marker", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(single_fish, file, row.names = FALSE)
    }
  )
  
  output$down2 <- downloadHandler(
    filename = function() {
      paste ("Example_Bulk_RNAseq_genelist", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(background_regen, file, row.names = FALSE)
    }
  )
  
  output$down3 <- downloadHandler(
    filename = function() {
      paste("ExampleBulkRNAseqcluster", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(cluster_regen, file, row.names = FALSE)
    }
  )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("enrichment", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data(), file, row.names = FALSE)
    }
  )     
  
  output$genes <- downloadHandler(
    filename = function() {
      paste("list_of_genes", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(list_genes(), file, row.names = FALSE)
    }
  )
  #####
  output$downloadData_extra1 <- downloadHandler(
    filename = function() {
      paste("extra1", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(extra11(), file, row.names = TRUE)
    }
  )
  
  output$downloadData_extra2 <- downloadHandler(
    filename = function() {
      paste("summary2", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(extra22(), file, row.names = TRUE)
    }
  )
  
  output$downloadData_extra3 <- downloadHandler(
    filename = function() {
      paste("summary3", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(extra33(), file, row.names = FALSE)
    }
  )
  
  output$downloadData_extra4 <- downloadHandler(
    filename = function() {
      paste("summary4", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(extra44(), file, row.names = FALSE)
    }
  )
  
  output$downloadData_extra5 <- downloadHandler(
    filename = function() {
      paste("summary5", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(extra55(), file, row.names = FALSE)
    }
  )
  
  output$downloadData_extra6 <- downloadHandler(
    filename = function() {
      paste("summary6", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(extra66(), file, row.names = TRUE)
    }
  )
  
  output$downloadData_extra6_6 <- downloadHandler(
    filename = function() {
      paste("summary7", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(extra6_61(), file, row.names = TRUE)
    }
  )
  
  
  output$downloadData_extra7 <- downloadHandler(
    filename = function() {
      paste("summary8", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(extra77(), file, row.names = TRUE)
    }
  )
  
  output$downloadData_extra8 <- downloadHandler(
    filename = function() {
      paste("summary9", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(extra88(), file, row.names = TRUE)
    }
  )
  
} #close server
shinyApp(ui, server)
#runApp(shinyApp(ui, server), launch.browser = TRUE)
