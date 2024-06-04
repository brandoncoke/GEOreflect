library(shiny)
require(openxlsx)
require(ggplot2)
require(DT)
require(plotly)
brandontheme=theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  plot.title = element_text(color="black", size=14,
                            face="plain",hjust = 0.5),
  axis.text = element_text(color="black", size=12,
                           face="plain",hjust = 0.5),
  axis.title = element_text(color="black", size=12, face="plain",
                            hjust = 0.5))
# Define UI for data upload app ----
ui <- fluidPage(
  tags$head(tags$script(src = "message-handler.js")),
  navbarPage("GEOreflect analysis", id="nav",
             tabPanel("DEG list upload",
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("file1", "Choose a .csv, .tsv or .xlsx DEG list",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",
                           ".xlsx",
                           ".txt",
                           ".tsv")),
      
      
      # Horizontal line ----
      tags$hr(),
      selectInput("gene_col", "Gene column",
                  NULL),
      selectInput("logfc_col", "Log fold column",
                  NULL),
      selectInput("pval_col", "p-value column",
                  NULL),
      checkboxInput("unmatched", "Ignore unmatched genes")
      
  
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      DT::dataTableOutput("contents")
      
    )
    
  )
),
tabPanel("Reranking and export",
         fluidRow(
           column(4,
                  sliderInput("minlogfc",
                              "Log fold lower bound",
                              min = -5,  max = 0, value = -1))
           ,
           column(4,
                  sliderInput("maxlogfc",
                              "Log fold upper bound",
                              min = 0,  max = 5,  value = 1))
           ,
           column(4,
                  sliderInput("pvallim",
                              "p-value limit",
                              min = 0,  max = 1,  value = 0.05))
           ),
           column(4,
                  actionButton("georeflect", "Run GEOreflect",
                               class = "btn-primary")
           ),
           column(4,
                  textInput("export_name", "Name for export file",
                            placeholder = "A_comparison")),
         column(4,
                checkboxInput("download", "Export ranked file?\n Check working directory")),
         column(4,
                checkboxInput("tha_plot_button", "Plot ranks- see next tab- it will take some time")),
             # Output: Data file ----
             DT::dataTableOutput("georeflect_frame")),
tabPanel("Plotting ranks",
         sidebarLayout(
           
           # Sidebar panel for inputs ----
           sidebarPanel(
             checkboxInput("extremes", "Highlight extreme differences"),
             tags$hr(),
             sliderInput("minlogfc",
                         "Log fold lower bound",
                         min = -5,  max = 0, value = -1),
             sliderInput("maxlogfc",
                         "Log fold upper bound",
                         min = 0,  max = 5,  value = 1),
             sliderInput("pvallim",
                         "p-value limit",
                         min = 0,  max = 1,  value = 0.05),
             
             actionButton("tha_plot_button", "Replot",
                          class = "btn-primary"),
             checkboxInput("plot_export", "Export plot"),
             sliderInput("export_size",
                         "Plot export size",
                         min = 500,  max = 7500,  value = 2000),
             
             

           ),
           #main plot for rank
           mainPanel(
             
             plotlyOutput("georeflect_plot",
                          width= "700px",
                          height = "700px")
             
           )
         ))
        
)
#This is the end
)



#server code
load("percentile_matrix_p_value_RNAseq.RDS")
get_data= function(input){
  req(input$file1, cancelOutput = T)
  tryCatch(
    {
      #check extension
      extension= "csv"
      if(grepl("[.]xlsx$", input$file1$datapath)){
        extension= "xlsx"
      }
      if(grepl("[.]tsv$", input$file1$datapath)){
        extension= "tsv"
      }
      
      switch(extension,
             csv={
               df= read.csv(input$file1$datapath,
                            header = T,
                            sep = ",")
             },
             tsv={
               df= read.csv(input$file1$datapath,
                            sep= "\t")
             },
             xlsx={
               require(openxlsx)
               df= openxlsx::read.xlsx(input$file1$datapath,sheet = 1,na.strings = T)
             },
             {
               df= read.csv(input$file1$datapath,
                            header = T,
                            sep = ",")
             }
      )
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      stop(safeError(e))
    }
  )
  return(df)
}
#can't reuse the code...
get_cols= function(input){
  req(input$file1, cancelOutput = T)
  tryCatch(
    {
      #check extension
      extension= "csv"
      if(grepl("[.]xlsx$", input$file1$datapath)){
        extension= "xlsx"
      }
      if(grepl("[.]tsv$", input$file1$datapath)){
        extension= "tsv"
      }
      
      switch(extension,
             csv={
               df= read.csv(input$file1$datapath,
                            header = T,
                            sep = ",")
             },
             tsv={
               df= read.csv(input$file1$datapath,
                            sep= "\t")
             },
             xlsx={
               require(openxlsx)
               df= openxlsx::read.xlsx(input$file1$datapath)
             },
             {
               df= read.csv(input$file1$datapath,
                            header = T,
                            sep = ",")
             }
      )
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      stop(safeError(e))
    }
  )
  return(colnames(df))
  
}

min_max_normalisation= function(vector=1:10){ #min max normalisation- 1 = highest, 0= lowest
  (vector-min(vector))/(max(vector)-min(vector))
}
geometric_mean= function(values= c(1,2,1)){
  (prod(values))^(1/length(values))
}
get_platform_percentile_RNA_seq= function(gene_and_pvalue){
  if(as.numeric(gene_and_pvalue[2])== 0 | 
     !(gene_and_pvalue[1] %in% rownames(percentile_matrix_p_value_RNAseq))){
    return(1)
  }else{
    gene= c(t(gene_and_pvalue[1])) #why does this work no clue but gets the string i need- t is transpose so eh#
    which.min(as.numeric(percentile_matrix_p_value_RNAseq[rownames(percentile_matrix_p_value_RNAseq) == gene, ]) <
                as.numeric(gene_and_pvalue[2])) #why do I need to coerce to float- how knows but safer as input can be treated as a string
  }
}
shift_label= function(x= output[1, ], median_shift){
  label= "-"
  if(as.logical(as.numeric(x[4]) > as.numeric(x[6]))){
    label= "↓"
  }
  if(as.logical(as.numeric(x[4]) < as.numeric(x[6]))){
    label= "↑"
  }
  if(abs(as.numeric(x[4]) - as.numeric(x[6])) >
     median_shift){
    label= paste0(label, label)
  }
  

  return(label)
}


GEOreflect_reranking_RNA_seq= function(the_frame, pvalue_indice, gene_indice,
                                       logfc_indice,
                                       minlogfc= -1,
                                       pvallim= 0.05,
                                       maxlogfc=1,
                                       unmatched_bool= T){

  temp= the_frame[, c(gene_indice,
                      pvalue_indice,
                      logfc_indice)]
  
 
  
  
  colnames(temp)= c("genes", "pvalues", "logfc")
  
  temp= temp[temp$logfc < minlogfc | temp$logfc > 
               maxlogfc, ]
  temp= temp[temp$pvalues <= pvallim, ]
  
  
  temp[is.na(temp$pvalues),]= 0
  temp[is.na(temp$logfc),]= 0
  if(unmatched_bool){
    temp= temp[temp$genes %in% rownames(percentile_matrix_p_value_RNAseq), ]
  }
  
  if(!any(temp$genes %in% rownames(percentile_matrix_p_value_RNAseq))){
    return(data.frame(Genes= temp$genes,
                      logFC= temp$logfc,
                      'p-values'= temp$pvalues,
                      pval_rank= rank(temp$pvalues),
                      GEOreflect_rank= NA,
                      Rank_change= NA,
                      Shift= NA,
                      Platform_relative_rank= NA
    ))
  }else{
    temp$pval_normalised= min_max_normalisation(temp$pvalues)
    temp$plat_rank= as.integer(apply(temp, 1,
                          get_platform_percentile_RNA_seq))
    temp$GEOreflect= min_max_normalisation(temp$plat_rank)
    temp$comb_score= apply(temp[,4:5], 1, geometric_mean)
    output= data.frame(Genes= temp$gene,
                       logFC= temp$logfc,
                       'p-values'= temp$pvalues,
                       pval_rank= rank(temp$pvalues),
                       Platform_relative_rank= temp$plat_rank,
                       GEOreflect_rank= rank(temp$comb_score))
    
    
    
    median_shift= as.numeric(
      quantile(abs(output$pval_rank - output$GEOreflect_rank), 
               0.5))
    output= output[order(output$GEOreflect_rank),]
    output$Rank_change= output$GEOreflect_rank - output$pval_rank
    
    output$Shift= apply(output, 1, shift_label, median_shift)
    return(output)
  }
}


server <- function(input, output, session) {
  
  output$contents <- DT::renderDataTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    df= get_data(input)
    
    return(df)
    
  })
  observe({
    colspresent= input$file1
    if (is.null(colspresent))
      colspresent <- ""
    colnames_string <- reactive(get_cols(input))
    colspresent= as.character(colnames_string())
    
    if(any(grepl("gene|hgnc|symbol|ident", as.character(colspresent)))){
      defualt_gene_indice=
        grep("gene|hgnc|symbol|ident", as.character(colspresent), ignore.case = T)[1]
    }else{
      defualt_gene_indice= 1
    }
    if(any(grepl("expr|logfc|fc|fold", as.character(colspresent),
                 ignore.case = T))){
      defualt_logfc= grep("expr|logfc|fc|fold|log_", as.character(colspresent),
                          ignore.case = T)[1]
    }else{
      defualt_logfc= 2
    }
    if(any(grepl("p-val|adj[_|-]p|pval|p[.]val", as.character(colspresent)))){
      defualt_pval=
        grep("p-val|adj[_|-]p|pval|p[.]val", as.character(colspresent), ignore.case = T)[1]
    }else{
      defualt_pval= 3
    }
    
  
    
    updateSelectInput(session,
                      "pval_col", "p-value column",
                      choices= colspresent,
                      selected = as.character(colspresent)[defualt_pval])
    updateSelectInput(session,
                      "logfc_col", "Log fold column",
                      choices= colspresent,
                      selected = as.character(colspresent)[defualt_logfc])
    updateSelectInput(session,
                      "gene_col", "Gene column",
                      choices= colspresent,
                      selected = as.character(colspresent)[defualt_gene_indice])
    
  })
  output$pval_selected <- renderText({
    paste0(input$pval_col)
  })
  observeEvent(input$georeflect, {

    output$georeflect_frame= DT::renderDataTable({
      req(input$file1, cancelOutput = T)
      req(input$gene_col, cancelOutput = T)
      req(input$logfc_col, cancelOutput = T)
      req(input$pval_col, cancelOutput = T)
      df= get_data(input)
      #session$sendCustomMessage(type = 'georeflect_message',
      #                          message = 'Running GEOreflect')
      pvalue_indice= which(colnames(df) == input$pval_col)
      gene_indice= which(colnames(df) == input$gene_col)
      logfc_indice= which(colnames(df) == input$logfc_col)
      df= df[df[, logfc_indice] < input$minlogfc | df[, 2] > 
                                             input$maxlogfc, ]
      df= df[df[, pvalue_indice] <= input$pvallim, ]

      
      
      GEOreflect_output= GEOreflect_reranking_RNA_seq(
        the_frame= df,
        pvalue_indice= pvalue_indice,
        gene_indice= gene_indice,
        logfc_indice= gene_indice,
        unmatched_bool = input$unmatched,
        minlogfc= input$minlogfc,
        pvallim= input$pvallim,
        maxlogfc=input$maxlogfc
      )
      
      
      if (input$download) {
        write.csv(GEOreflect_output,
                  paste0("~/", input$export_name, "_GEOreflect_ranking.csv"),
                  sep="")
      }
      colnames(GEOreflect_output)= gsub(
        "_",
        " ",
        colnames(GEOreflect_output))
      
      
      
      
      return(GEOreflect_output)
      
    })
  })
  
  
  
  
  
  
  observeEvent(input$tha_plot_button, {
    output$georeflect_plot= renderPlotly({
      if (input$tha_plot_button) {
        req(input$file1, cancelOutput = T)
        req(input$gene_col, cancelOutput = T)
        req(input$logfc_col, cancelOutput = T)
        req(input$pval_col, cancelOutput = T)
        df= get_data(input)
        session$sendCustomMessage(type = 'georeflect_message',
                                  message = 'Running GEOreflect and plotting ranks')
        pvalue_indice= which(colnames(df) == input$pval_col)
        gene_indice= which(colnames(df) == input$gene_col)
        logfc_indice= which(colnames(df) == input$logfc_col)
        #df= df[df[,2] < input$minlogfc | df[, 2] > 
        #         input$maxlogfc, ]
        #df= df[df[, 3] <= input$pvallim, ]
        GEOreflect_output= GEOreflect_reranking_RNA_seq(
          the_frame= df,
          pvalue_indice= pvalue_indice,
          gene_indice= gene_indice,
          logfc_indice= gene_indice,
          unmatched_bool = input$unmatched,
          minlogfc= input$minlogfc,
          pvallim= input$pvallim,
          maxlogfc=input$maxlogfc
        )
        
        
        bool_sum= as.integer(nrow(GEOreflect_output) > 15) +
          as.integer(nrow(GEOreflect_output) > 50) +
          as.integer(nrow(GEOreflect_output) > 100) +
          as.integer(nrow(GEOreflect_output) > 250) +
          as.integer(nrow(GEOreflect_output) > 750) +
          as.integer(nrow(GEOreflect_output) > 1000) +
          as.integer(nrow(GEOreflect_output) > 2000) +
          as.integer(nrow(GEOreflect_output) > 5000)
        bool_sum= bool_sum+ 1
        axes_scale=switch(bool_sum,
                          axes_scale= seq(0, nrow(GEOreflect_output), 2),
                          axes_scale= seq(0, nrow(GEOreflect_output), 5),
                          axes_scale= seq(0, nrow(GEOreflect_output), 25),
                          axes_scale= seq(0, nrow(GEOreflect_output), 50),
                          axes_scale= seq(0, nrow(GEOreflect_output), 100),
                          axes_scale= seq(0, nrow(GEOreflect_output), 250),
                          axes_scale= seq(0, nrow(GEOreflect_output), 750),
                          axes_scale= seq(0, nrow(GEOreflect_output), 1000))
        colnames(GEOreflect_output)[colnames(GEOreflect_output) == 
                                      "Genes"]= "Gene"
        
        
        colnames(GEOreflect_output)= gsub(
          "pval_rank",
          "p-value_rank",
          colnames(GEOreflect_output))
        
        colnames(GEOreflect_output)= gsub(
          "_",
          " ",
          colnames(GEOreflect_output))
        #export_frame$opacity= 0.01
        
        
        
        
        tha_plot= ggplot(GEOreflect_output,
                         aes(y=`GEOreflect rank`, x=`p-value rank`, label=Gene)) +
          geom_point(size= 3, 
                     color= "#008080") +
          geom_abline(intercept = 0, slope = 1, color="red",
                      linetype="dashed", linewidth=1.5) +
          labs(x="p-value rank",y="GEOreflect rank",title =
                 "p-value vs GEOreflect rank") +
          scale_y_continuous(expand = c(0, 0), breaks =  axes_scale) +
          scale_x_continuous(expand = c(0, 0), breaks =  axes_scale) +
          brandontheme +
          guides(alpha = "none")
        if(input$plot_export){
          tha_plot
          ggsave(paste0("~/", input$export_name, "_GEOreflect_ranking_plot.tiff")
                 ,width=input$export_size,height= input$export_size,units = "px")
        }
        
        return(ggplotly(tha_plot))
      
      }
      
      
      
      
      
      
      
    })})
  #Ends here   
  }
# Create Shiny app ----
shinyApp(ui, server)