## app.R ##
Rlib="/data/boehm/group/shiny_apps/Rlibs3.4.0"
library(shiny,lib.loc=Rlib)
library(shinydashboard,lib.loc=Rlib)
library(rhandsontable,lib.loc=Rlib)
library(DT,lib.loc=Rlib)

#options(shiny.maxRequestSize=5000*1024^2)

ui <- function(request) {dashboardPage(
    dashboardHeader(title = "Dataset selection"),
    ## Sidebar content
    dashboardSidebar(

        textInput(inputId="group", label="Group", value = "", width = NULL, placeholder = NULL),
        selectInput(inputId="projectid", label="ProjectID", choices=c("SELECT A PROJECT","Trancoso512A_B_C")),
        actionButton("submitinput", "Retrieve dataset"),
        textInput(inputId="geneid", label="GeneID", value="",placeholder="TYPE IN GENE ID"),
        actionButton("selectgenes", "Select genes"),
        textOutput("fileDescription"),
        bookmarkButton()
        ),
        
    dashboardBody(
        h2("Single cell RNAseq analysis"),
        uiOutput("resultPanels")
               
    )

 )}


server <- function(input, output, session) {
    
    output$walkThrough<-renderUI(HTML("<ul><li>1.Provide group and project ID information to retrieve a serialized R object containing your dataset. Path will be taken from a database of datasets(3). Click on retrieve dataset. Your data will appear in the InputData tab.</li><li>2.Provide semicolon-separated Gene IDs to plot aggregate expression for.</li><li>3.#This is a currently not implemented## Provide rules for cell assignment to a known class that will be used to facet the plots.Click on Run analysis.#End of not implemented# Your results will appear in the corresponding tabs.</li><li>The order of providing the information matters!</li></ul>"))
    output$FAQ<-renderText("Currently, no uniform gene naming system is prerequisite. You have to provide Gene IDs consistent with the naming used to produce your dataset.\n Merging data from multiple datasets or batch effect removal are currenlty not supported.\n For questions, bug reports or feature requests, contact sikora@ie-freiburg.mpg.de.\n For reporting issues or pull requests on GitHub, go to https://github.com/katsikora/scRNAseq_shiny_app .")
    
    output$fileDescription<-renderText("GeneID: Please provide a semicolon-separated list of Gene IDs you would like to obtain results for.")


################################
    #require(rio,lib.loc=Rlib)
    library(RaceID,lib.loc=Rlib)
    #require(scater,lib.loc=Rlib)
    library(ggplot2,lib.loc=Rlib)
    

    output$sessionInfo <- renderPrint({capture.output(sessionInfo())})

    dbFile<-read.table("/data/manke/group/shiny/sikora/aux_files/scRNAseq.DB.csv",header=TRUE,sep="\t",quote="",as.is=TRUE)


################################
    ##load input data, show head of ndata
    observeEvent(input$submitinput, {
        if((input$group!="")&(input$projectid!="")){
            inGroup<-isolate(input$group)
            inProjectID<-isolate(input$projectid)
  
            datPath<-dbFile$PathToData[(dbFile$Group %in% inGroup) & (dbFile$ProjectID %in% inProjectID)]
            output$debug<-renderPrint({capture.output(print(datPath))})
            }else{output$debug<-renderText("Input file not found.")}
        #sc<-readRDS("/data/manke/group/shiny/sikora/aux_files/sc.minT5000.rds")
        #sc<-readRDS("/data/boehm/sikora/sikora/scRNAseq.lamprey/171201ALL.Lp.Trinity.good.OLD.NEW.cVLR.RaceID.monocle.rBE.CGenes.RPLS.FGenes.cellCycle.workspaceR/sc.minT5000.RDS")
        if(grepl("rds$",datPath,ignore.case=TRUE)){
            sc<-readRDS(datPath)}
        else if (grepl("rdata$",datPath,ignore.case=TRUE)){
           myEnv<-environment()
           sctmp<-load(datPath, envir = myEnv)
           sc <- myEnv[[sctmp]]
            }
        #render the head
        ndata<-as.data.frame(as.matrix(sc@ndata)*5000,stringsAsFactors=FALSE)
        output$datHead<-renderTable({ndata[1:10,1:min(8,ncol(ndata))]},caption="Normalized data",caption.placement = getOption("xtable.caption.placement", "top"),include.rownames=TRUE)
        load("/data/boehm/sikora/trancoso/512/180723ALL.Lp.Trinity.good.OLD.NEW.cVLR_CDA1x2.RaceID3NEW.NOrBE.CGenes.RPLS.FGenes.cellCycle.workspaceR/WGCNA.p6/magClust.RData")
        
        output$configurator<-renderUI({tagList(selectInput(inputId="SHhit",label="CDA/VLR hit:",choices=c("All",unique(as.character(magClust$CDA.VLR.hit)))),
                                               selectInput(inputId="mod",label="Module:",choices=c("All",unique(as.character(magClust$module))))) })
        
        output$magclust<-renderDT({
          data <- magClust  
        if (input$SHhit != "All") {
          data <- data[data$CDA.VLR.hit == input$SHhit,]
        }
        if (input$mod != "All") {
          data <- data[data$module == input$mod,]
        }
          data},options = list(autoWidth = TRUE,scrollX=TRUE), filter = "bottom")

        observeEvent(input$selectgenes, {
            inGenesL<-isolate(input$geneid)
            inGenes<-unlist(strsplit(inGenesL,split=";"))
            output$ingenes<-renderText(inGenes)
            nv<-inGenes[inGenes %in% rownames(ndata)]
            if(length(nv)>0){
                ifelse(length(nv)==1,nt<-nv,nt<-"Selected genes")
            output$tsneAgg<-renderPlot({plotexpmap(sc,nv,n=nt,logsc=as.logical(input$tsnelog))})
            }#fi
            ###produce top correlated genes for aggregated selected gene(s)
            cor.log2<-cor(x=log2(colSums(ndata[rownames(ndata) %in% inGenes,])+0.1),y=t(log2(ndata+0.1))) 
            output$corlog2<-renderPlot({boxplot(t(cor.log2))})
            cor.log2T<-as.data.frame(t(cor.log2),stringsAsFactors=FALSE)
            colnames(cor.log2T)<-"cor"
            cor.log2T$abscor<-abs(cor.log2T$cor)
            cor.log2T<-cor.log2T[order(cor.log2T$abscor,decreasing=TRUE),]
            output$top10cor<-renderTable({head(cor.log2T[,"cor",drop=FALSE],n=10)},caption="Top 10 correlated genes",caption.placement = getOption("xtable.caption.placement", "top"),include.rownames=TRUE)
 


        },ignoreInit=TRUE)#end of observe input$selectgenes

        observeEvent(input$plotpwcor,{
            if((input$pwselX!="")&(input$pwselY!="")){
            inpwselXL<-isolate(input$pwselX)
            inpwselYL<-isolate(input$pwselY)}
            inpwselX<-unlist(strsplit(inpwselXL,split=";"))
            inpwselY<-unlist(strsplit(inpwselYL,split=";"))
            plotdat<-as.data.frame(cbind(colSums(ndata[rownames(ndata) %in% inpwselX,]),colSums(ndata[rownames(ndata) %in% inpwselY,])),stringsAsFactors=FALSE)
            colnames(plotdat)<-c("X","Y")
            output$pwplot<-renderPlot({ggplot(data=plotdat)+geom_point(aes(x=X,y=Y))})

       },ignoreInit=TRUE)#end of observe input$plotpwcor 

    },ignoreInit=TRUE)#end of observe input$submitinput


############################
    output$resultPanels<-renderUI({myTabs<-list(tabPanel(title="WalkThrough",
                                                      fluidPage(
                                                          box(title="Walkthrough",uiOutput("walkThrough")),
                                                          box(title="Miscellaneous information",textOutput("FAQ"))                                                          
                                                               )
                                                          ),
                                                  tabPanel(title="InputData",
                                                      fluidPage(
                                                          fluidRow(
                                                              tableOutput("datHead"),
                                                              tableOutput("countDatHead")
                                                                   ),
                                                          #fluidRow(
                                                              rHandsontableOutput("hot"),
                                                              actionButton(inputId="savetable",label="Run analysis"),
                                                            #      ),
                                                          fluidRow(
                                                              box(renderText("Please complete missing sample information. Group variable will be used for plot faceting.")),
                                                              tableOutput("inTabHead"),
                                                              box(title="Debug",
                                                                  textOutput("debug"),
                                                                  textOutput("ingenes"))
                                                                  ) 
                                                              )
                                                          ),
                                                
                                                tabPanel(title="Annotation.Table",
                                                         fluidRow(
                                                           uiOutput("configurator")
                                                           ),
                                                         #fluidPage(
                                                         DTOutput("magclust"),
                                                         actionButton(inputId="selGenesFromTab",label="Select Gene IDs from table")
                                                         #)
                                                ),
                                                   tabPanel(title="Tsne.Map",
                                                      fluidPage(
                                                          box(plotOutput("tsneAgg"),width=4),
                                                          box(title = "Plot controls",selectInput("tsnelog", "Log scale",choices=c("TRUE","FALSE"),selected="TRUE")),
                                                          box(title="Method Description",renderText("(Log) counts were aggregated over selected genes and the expression levels were colour-coded on the tsne map."))
                                                                )
                                                              ),
                                                   tabPanel(title="Top.Correl.Genes",
                                                      fluidPage(
                                                          box(plotOutput("corlog2"),width=4),
                                                          box(tableOutput("top10cor")),
                                                          box(title="Method Description",renderText("Pearson correlation was calculated between log2-transformed aggregated counts for gene selection and all log2-transformed genes in the ndata slot of the sc object. Top 10 genes are listed."))
                                                               )
                                                          ),
                                                  tabPanel(title="Pairwise.Expression",
                                                      fluidPage(
                                                          box(plotOutput("pwplot")),
                                                          box(title = "Select gene(s) X",textInput(inputId="pwselX", label="Gene symbol(s) for X axis",value="")),
                                                          box(title = "Select gene(s) Y",textInput(inputId="pwselY", label="Gene symbol(s) for Y axis",value="")),
                                                          actionButton(inputId="plotpwcor",label="Plot expression"),
                                                          box(title="Method Description",renderText("Pairwise plot of normalized counts."))
                                                          
                                                               )
                                                          ),
                                                  tabPanel(title="sessionInfo",
                                                      fluidPage(
                                                          verbatimTextOutput("sessionInfo")                                                          
                                                               )
                                                          )

                                                 )

            do.call(tabsetPanel, myTabs)})


}

shinyApp(ui, server,enableBookmarking="url")
