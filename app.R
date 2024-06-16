library(shiny)
require(openxlsx)
require(ggplot2)
require(DT)
require(plotly)
#60MB limit
options(shiny.maxRequestSize=60*1024^2)
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
                                    accept = c(".csv",
                                               ".xlsx",
                                               ".txt",
                                               ".tsv")),
                          checkboxInput("header", "Header present?",
                                        value = T),
                          
                          
                          # Horizontal line ----
                          tags$hr(),
                          selectInput("gene_col", "Gene column- HNGC symbols like USP7",
                                      NULL),
                          selectInput("logfc_col", "Log fold column",
                                      NULL),
                          selectInput("pval_col", "p-value column",
                                      NULL),
                          selectInput("average_exprss_col", "Average expression column",
                                      NULL),
                          checkboxInput("unmatched", "Ignore unmatched genes",
                                        value = T)
                          
                          
                          
                          
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
                            header = input$header,
                            sep = ",")
             },
             tsv={
               df= read.csv(input$file1$datapath,
                            sep= "\t",
                            header= input$header)
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
      if(grepl("[.]tsv$|[.]txt$", input$file1$datapath)){
        extension= "tsv"
      }
      
      switch(extension,
             csv={
               df= read.csv(input$file1$datapath,
                            header = input$header,
                            sep = ",")
             },
             tsv={
               df= read.csv(input$file1$datapath,
                            input$header,
                            sep= "\t")
             },
             xlsx={
               require(openxlsx)
               df= openxlsx::read.xlsx(input$file1$datapath)
             },
             {
               df= read.csv(input$file1$datapath,
                            input$header,
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
  if(as.logical(as.numeric(x[1]) > as.numeric(x[2]))){
    label= "↓"
  }
  if(as.logical(as.numeric(x[1]) < as.numeric(x[2]))){
    label= "↑"
  }
  if(abs(as.numeric(x[1]) - as.numeric(x[2])) >
     median_shift){
    label= paste0(label, label)
  }
  
  
  return(label)
}

GEOreflect_reranking_RNA_seq= function(the_frame,
                                       pvalue_indice,
                                       gene_indice,
                                       logfc_indice, 
                                       average_exprss_indice,
                                       minlogfc= -1,
                                       pvallim= 0.05,
                                       maxlogfc=1,
                                       unmatched_bool= T){
  temp= the_frame[, c(gene_indice,
                      pvalue_indice,
                      logfc_indice,
                      average_exprss_indice)]

  colnames(temp)= c("genes", "pvalues", "logfc", "average_expression")
  temp= temp[temp$genes != "", ]
  if(class(temp$pvalues) != "numeric" |
     class(temp$logfc) != "numeric" |
     class(temp$genes) != "character" |
     sum(is.na(temp$pvalues) == nrow(temp))|
     nrow(temp) == 0){
    return(data.frame(Genes= temp$genes,
                      logFC= temp$logfc,
                      'p-values'= temp$pvalues,
                      pval_rank= rank(temp$pvalues),
                      average_expression_rank= rank(-temp$average_expression),
                      GEOreflect_rank= NA,
                      Rank_change= NA,
                      Shift= NA,
                      Platform_relative_rank= NA
    ))
  }else{
    temp= temp[!is.na(temp$logfc) &
                 !is.na(temp$pvalues), ]
    temp= temp[temp$logfc < minlogfc | temp$logfc > 
                 maxlogfc, ]
    temp= temp[temp$pvalues <= pvallim, ]
    if(unmatched_bool){
      temp= temp[temp$genes %in% rownames(percentile_matrix_p_value_RNAseq), ]
    }
    
    if(!any(temp$genes %in% rownames(percentile_matrix_p_value_RNAseq)) |
       nrow(temp) == 0){
      if(nrow(temp) == 0){
        showModal(modalDialog(
          title = "Filtering step too stringent",
          paste0("Increase the p-value limit or bounds for upper and lower log fold change"),
          easyClose = TRUE,
          footer = NULL
        ))
      }
      
      
      return(data.frame(Genes= temp$genes,
                        logFC= temp$logfc,
                        'p-values'= temp$pvalues,
                        pval_rank= rank(temp$pvalues),
                        average_expression_rank= rank(-temp$average_expression),
                        GEOreflect_rank= NA,
                        Rank_change= NA,
                        Shift= NA,
                        Platform_relative_rank= NA
      ))
    }else{
      temp$pval_normalised= rank(-temp$pvalues)/(length(temp$pvalues))
      temp$express_norm= rank(temp$average_expression)/(length(temp$average_expression))
      
      temp$plat_rank= as.integer(apply(temp, 1,
                                       get_platform_percentile_RNA_seq))
      temp$GEOreflect= ((1000 - temp$plat_rank)+1)/ 1000
      temp$comb_score= apply(temp[,c("GEOreflect", "pval_normalised", 
                                     "express_norm")], 1, geometric_mean)
      output= data.frame(Genes= temp$gene,
                         logFC= temp$logfc,
                         'p-values'= temp$pvalues,
                         pval_rank= rank(temp$pvalues),
                         average_expression_rank= rank(-temp$average_expression),
                         Platform_relative_rank= temp$plat_rank,
                         GEOreflect_rank= rank(-temp$comb_score))
      
      
      
      median_shift= as.numeric(
        quantile(abs(output$pval_rank - output$GEOreflect_rank), 
                 0.75))
      output= output[order(output$GEOreflect_rank),]
      output$Rank_change= output$GEOreflect_rank - output$pval_rank
      
      output$Shift= apply(output[,c('pval_rank', 'GEOreflect_rank')], 1, shift_label, median_shift)
      return(output)
    }
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
    #Symbol defaults
    if(any(grepl("gene|hgnc|symbol|ident", as.character(colspresent)),
           ignore.case = T)
    ){
      defualt_gene_indice=
        grep("gene|hgnc|symbol|ident", as.character(colspresent),
             ignore.case = T)[1]
      #if it actually has the hngc symbols
      if(any(grepl("hgnc|symbol", as.character(colspresent)),
             ignore.case = T)){
        defualt_gene_indice=
          grep("hgnc|symbol", as.character(colspresent), ignore.case = T)[1]
      }
      
    }else{
      defualt_gene_indice= 1
    }
    #logFC defualts
    if(any(grepl("expr|logfc|fc|fold|log_", as.character(colspresent),
                 ignore.case = T))){
      defualt_logfc= grep("expr|logfc|fc|fold|log_", as.character(colspresent),
                          ignore.case = T)[1]
      
      #if it has multiple fcs
      if(any(grepl("log", as.character(colspresent), ignore.case = T)
             &
             grepl("fold|fc", as.character(colspresent), ignore.case = T)
      )
      ){
        bool= grepl("log", as.character(colspresent), ignore.case = T) &
          grepl("fold|fc", as.character(colspresent), ignore.case = T)
        
        defualt_logfc= which(bool)[1]
      }
    }else{
      defualt_logfc= 2
    }
    
    
    defualt_pval= 3
    if(any(grepl("p-val|pval|p[.]val", as.character(colspresent), 
                 ignore.case = T))){
      
      defualt_pval=
        which(grepl("p-val|pval|p[.]val", as.character(colspresent), 
                    ignore.case = T))
      defualt_pval= defualt_pval[1]
      
    }
    #better if using non adjusted p-values
    if(any(grepl("p-val|pval|p[.]val", as.character(colspresent), 
                 ignore.case = T) & 
           !grepl("adj", as.character(colspresent), 
                  ignore.case = T))){
      
      defualt_pval=
        which(grepl("p-val|pval|p[.]val", as.character(colspresent), 
                    ignore.case = T) & 
                !grepl("adj", as.character(colspresent), 
                       ignore.case = T))
      defualt_pval= defualt_pval[1]
      
    }
    default_expression=2
    if(any(grepl("avg_expression|expression", as.character(colspresent), 
                 ignore.case = T))){
      
      default_expression=
        which(grepl("avg[_|.| ]expr|expression|AveExpr|avex", as.character(colspresent), 
                    ignore.case = T))
      default_expression= default_expression[1]
      
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
                      "gene_col", "Gene column- HNGC symbols like USP7",
                      choices= colspresent,
                      selected = as.character(colspresent)[defualt_gene_indice])
    updateSelectInput(session,
                      "average_exprss_col", "Average expression column",
                      choices= colspresent,
                      selected = as.character(colspresent)[default_expression])
    
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
      req(input$average_exprss_col, cancelOutput = T)
      showModal(modalDialog(
        title = "GEOreflect is running. This may take some time",
        paste0("Can take around 3 minutes- needs around 3 GB RAM free."),
        easyClose = TRUE,
        footer = NULL
      ))
      df= get_data(input)
      gene_indice= which(colnames(df) == input$gene_col)
      logfc_indice= which(colnames(df) == input$logfc_col)
      pvalue_indice= which(colnames(df) == input$pval_col)
      average_exprss_indice= which(colnames(df) == input$average_exprss_col)
      GEOreflect_output= GEOreflect_reranking_RNA_seq(
        the_frame= df,
        pvalue_indice= pvalue_indice,
        gene_indice= gene_indice,
        logfc_indice= logfc_indice,
        unmatched_bool = input$unmatched,
        average_exprss_indice= average_exprss_indice,
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
    showModal(modalDialog(
      title = "Plotting data- will take some time",
      paste0("Can take around 3 minutes- needs around 3 GB RAM free."),
      easyClose = TRUE,
      footer = NULL
    ))
    output$georeflect_plot= renderPlotly({
      if (input$tha_plot_button) {
        req(input$file1, cancelOutput = T)
        req(input$gene_col, cancelOutput = T)
        req(input$logfc_col, cancelOutput = T)
        req(input$pval_col, cancelOutput = T)
        df= get_data(input)
        #session$sendCustomMessage(type = 'georeflect_message',
        #                          message = 'Running GEOreflect and plotting ranks')
        pvalue_indice= which(colnames(df) == input$pval_col)
        gene_indice= which(colnames(df) == input$gene_col)
        logfc_indice= which(colnames(df) == input$logfc_col)
        average_exprss_indice= which(colnames(df) == input$average_exprss_col)
        
        #df= df[df[,2] < input$minlogfc | df[, 2] > 
        #         input$maxlogfc, ]
        #df= df[df[, 3] <= input$pvallim, ]
        GEOreflect_output= GEOreflect_reranking_RNA_seq(
          the_frame= df,
          pvalue_indice= pvalue_indice,
          gene_indice= gene_indice,
          logfc_indice= logfc_indice,
          unmatched_bool = input$unmatched,
          average_exprss_indice= average_exprss_indice,
          minlogfc= input$minlogfc,
          pvallim= input$pvallim,
          maxlogfc=input$maxlogfc
        )
        
        
        if(any(as.numeric(GEOreflect_output$GEOreflect_rank))){
          bool_sum= as.integer(nrow(GEOreflect_output) > 15) +
            as.integer(nrow(GEOreflect_output) > 50) +
            as.integer(nrow(GEOreflect_output) > 100) +
            as.integer(nrow(GEOreflect_output) > 250) +
            as.integer(nrow(GEOreflect_output) > 500) +
            as.integer(nrow(GEOreflect_output) > 750) +
            as.integer(nrow(GEOreflect_output) > 1000) +
            as.integer(nrow(GEOreflect_output) > 2500) +
            as.integer(nrow(GEOreflect_output) > 5000)
          bool_sum= bool_sum+ 1
          axes_scale=switch(bool_sum,
                            axes_scale= seq(0, nrow(GEOreflect_output), 2),
                            axes_scale= seq(0, nrow(GEOreflect_output), 5),
                            axes_scale= seq(0, nrow(GEOreflect_output), 25),
                            axes_scale= seq(0, nrow(GEOreflect_output), 50),
                            axes_scale= seq(0, nrow(GEOreflect_output), 50),
                            axes_scale= seq(0, nrow(GEOreflect_output), 100),
                            axes_scale= seq(0, nrow(GEOreflect_output), 100),
                            axes_scale= seq(0, nrow(GEOreflect_output), 150),
                            axes_scale= seq(0, nrow(GEOreflect_output), 200))
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
          
          
          
          
          if(input$extremes){
            GEOreflect_output$opacity= 0.05
            GEOreflect_output$opacity[GEOreflect_output$Shift ==
                                        "↑↑" |
                                        GEOreflect_output$Shift == "↓↓"]= 1
            tha_plot= ggplot(GEOreflect_output,
                             aes(y=`GEOreflect rank`, x=`p-value rank`, label=Gene)) +
              geom_point(size= 3, 
                         color= "#008080", 
                         aes(alpha= opacity)) 
          }else{
            tha_plot= ggplot(GEOreflect_output,
                             aes(y=`GEOreflect rank`, x=`p-value rank`, label=Gene)) +
              geom_point(size= 3, 
                         color= "#008080") 
          }
          
          
          
          tha_plot= tha_plot +
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
        }else{
          showModal(modalDialog(
            title = "Check that your gene/probe column is correct",
            paste0("They should be formatted as HNGC symbols e.g. USP7, ABCA1 and WNT5A"),
            easyClose = TRUE,
            footer = NULL
          ))
          bool_sum= as.integer(nrow(GEOreflect_output) > 15) +
            as.integer(nrow(GEOreflect_output) > 50) +
            as.integer(nrow(GEOreflect_output) > 100) +
            as.integer(nrow(GEOreflect_output) > 250) +
            as.integer(nrow(GEOreflect_output) > 500) +
            as.integer(nrow(GEOreflect_output) > 750) +
            as.integer(nrow(GEOreflect_output) > 1000) +
            as.integer(nrow(GEOreflect_output) > 2500) +
            as.integer(nrow(GEOreflect_output) > 5000)
          bool_sum= bool_sum+ 1
          axes_scale=switch(bool_sum,
                            axes_scale= seq(0, nrow(GEOreflect_output), 2),
                            axes_scale= seq(0, nrow(GEOreflect_output), 5),
                            axes_scale= seq(0, nrow(GEOreflect_output), 25),
                            axes_scale= seq(0, nrow(GEOreflect_output), 50),
                            axes_scale= seq(0, nrow(GEOreflect_output), 50),
                            axes_scale= seq(0, nrow(GEOreflect_output), 100),
                            axes_scale= seq(0, nrow(GEOreflect_output), 100),
                            axes_scale= seq(0, nrow(GEOreflect_output), 150),
                            axes_scale= seq(0, nrow(GEOreflect_output), 200))
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
                           aes(y=`p-value rank`, x=`p-value rank`, label=Gene)) +
            geom_point(size= 3, 
                       color= "#008080")  +
            geom_abline(intercept = 0, slope = 1, color="red",
                        linetype="dashed", linewidth=1.5) +
            labs(x="p-value rank",y="GEOreflect rank",title = 
                   "<a href = 'https://www.nytimes.com/'>The NY TIMES</a>") +
                   #"p-value vs GEOreflect rank") +
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
        
        
        
      }
      
      
      
      
      
      
      
    })})
  #Ends here   
}
# Create Shiny app ----
shinyApp(ui, server)