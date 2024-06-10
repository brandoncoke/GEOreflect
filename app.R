library(shiny)
require(openxlsx)
require(ggplot2)
require(DESeq2)
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
      tags$hr(),
      selectizeInput(
        'CTRL_cols', 'Select control columns', choices = NULL,
        multiple = TRUE),
      selectizeInput(
        'treated_cols', 'Select treated columns',
        choices = NULL, multiple = TRUE),
      checkboxInput("Normalise reads", "Raw reads as input?",
                    value = T),
      
      
      # Horizontal line ----

      selectInput("id_type", "Gene id type- HGNC like CTNNB1, ensembl- ENSG00000168036 or entrez like 1499",
                  c("HGNC",
                    "Entrez id",
                    "Ensembl id")),
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
conversion_table= read.csv("conversion_table.csv")
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

cut_isoforms= function(id= "ENSG00000066827.11"){
  unlist(strsplit(id, "[.]"))[1]
}


GEOreflect_reranking_RNA_seq= function(the_frame,
                                       minlogfc= -1,
                                       pvallim= 0.05,
                                       maxlogfc=1,
                                       unmatched_bool= T){
  temp= the_frame
  colnames(temp)= c("genes", "pvalues", "logfc", "average_expression")
  temp= temp[temp$genes != "", ]
  
  if(class(temp$pvalues) != "numeric" |
     class(temp$logfc) != "numeric" |
     class(temp$genes) != "character" |
     sum(is.na(temp$pvalues) == nrow(temp))){
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
                         Platform_relative_rank= temp$plat_rank,
                         GEOreflect_rank= rank(-temp$comb_score))
      
      median_shift= as.numeric(
        quantile(abs(output$pval_rank - output$GEOreflect_rank), 
                 0.75))
      output= output[order(output$GEOreflect_rank),]
      output$Rank_change= output$GEOreflect_rank - output$pval_rank
      
      output$Shift= apply(output, 1, shift_label, median_shift)
      return(output)
    }
  }
}


server <- function(input, output, session) {
  
  output$contents <- DT::renderDataTable({
    req(input$file1)
    
    df= get_data(input)
    default_control= NULL
    if(any(grepl("CON|CTRL|GFP|scram|WT", as.character(colnames(df)), 
                 ignore.case = T))){
      default_control=
        which(grepl("CON|CTRL|GFP|scram|WT", as.character(colnames(df)), 
                    ignore.case = T))
    }
    default_treated= NULL
    if(any(grepl("treated|sh|KO|KD|CRISPR", as.character(colnames(df)), 
                 ignore.case = T) &
       !grepl("CON|CTRL|GFP|scram|WT", as.character(colnames(df)), 
             ignore.case = T))
       ){
      default_treated=
        which(grepl("treated|sh|KO|KD|CRISPR", as.character(colnames(df)), 
                    ignore.case = T) &
                !grepl("CON|CTRL|GFP|scram|WT", as.character(colnames(df)), 
                       ignore.case = T))
    }
    
    

    updateSelectizeInput(session,
                      "CTRL_cols", 'Select control samples',
                      choices= colnames(df),
                      selected= colnames(df)[default_control])
    updateSelectizeInput(session,
                         "treated_cols", 'Select treated samples',
                         choices= colnames(df),
                         selected= colnames(df)[default_treated])
    colnames(df)[1]= "ID_col"
    
    
    updateSelectInput(session,
                      "id_col", "Gene identifier column",
                      choices= colnames(df))

    ids= df[, 1]
    default_id= "HGNC"
    if(sum(sum(ids %in% conversion_table$entrezgene_id) > sum(ids %in% conversion_table$gene))){
      default_id= "Entrez id"
    }
    if(sum(sum(ids %in% conversion_table$ensembl_gene_id) > sum(ids %in% conversion_table$gene))){
      default_id= "Ensembl id"
    }
    updateSelectInput(session,
                      "id_type", "Gene id type- HGNC like CTNNB1, ensembl- ENSG00000168036 or entrez like 1499",
                      c("HGNC",
                        "Entrez id",
                        "Ensembl id"),
                      selected = default_id)
    
    return(df)
  })
  
  #GEOreflect 
  observeEvent(input$georeflect, {
    req(input$file1, cancelOutput = T)
    req(input$CTRL_cols, cancelOutput = T)
    req(input$treated_cols, cancelOutput = T)
    req(input$id_type, cancelOutput = T)
    output$georeflect_frame= DT::renderDataTable({
    df= get_data(input)

    reads= df[,c(input$CTRL_cols, 
                 input$treated_cols)]
    reads= as.data.frame(reads)
    ids= as.character(df[, 1])
    reads$ids= ids
    if(input$id_type == "Entrez id"){
      colnames(reads)[ncol(reads)]= "entrezgene_id"
      
      reads= merge(conversion_table, reads, all= F)
    }
    if(input$id_type == "Ensembl id"){
      reads$ids= as.character(
        lapply(reads$ids, cut_isoforms))
      colnames(reads)[ncol(reads)]= "ensembl_gene_id"
      reads= merge(conversion_table, reads, all= F)
    }
    if(input$id_type == "HGNC"){
      colnames(reads)[ncol(reads)]= "gene"
    }
    rownames(reads)= paste0(reads$gene,";", 1:nrow(reads))
    reads= reads[,-c(1:3)]
    
    
    gsms <- paste0(
      c(
        rep("0", length(input$CTRL_cols)),
        rep("1", length(input$treated_cols))
      ),
      collapse = ""
    )
    sml <- strsplit(gsms, split="")[[1]]
    
    # filter out excluded samples (marked as "X")
    sel <- which(sml != "X")
    sml <- sml[sel]
    
    # group membership for samples
    gs <- factor(sml)
    groups <- make.names(c("CTRL","Treated"))
    levels(gs) <- groups
    sample_info <- data.frame(Group = gs, row.names = colnames(reads))
    
    # pre-filter low count genes
    # keep genes with at least N counts > 10, where N = size of smallest group
    keep <- rowSums( reads >= 10 ) >= min(table(gs))
    reads <- reads[keep, ]
    
    ds <- DESeqDataSetFromMatrix(countData=reads, colData=sample_info, design= ~Group)
    
    ds <- DESeq(ds, test="Wald", sfType="poscount")
    
    # extract results for top genes table
    r <- results (ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")
    DEG_frame= as.data.frame(r)
    DEG_frame$Gene= rownames(DEG_frame)
    DEG_frame$Gene= as.character(lapply(DEG_frame$Gene, 
                                        function(x){
                                          unlist(strsplit(x, ";"))[1]
                                        }))
    
    DEG_frame= DEG_frame[, c("Gene", "pvalue", "log2FoldChange", "baseMean")]
    output= GEOreflect_reranking_RNA_seq(the_frame= DEG_frame,
                                         minlogfc= input$minlogfc,
                                         pvallim= input$pvallim,
                                         maxlogfc= input$maxlogfc,
                                         unmatched_bool= input$unmatched)

    })
  })
  observeEvent(input$tha_plot_button, {
    showModal(modalDialog(
      title = "Plotting data- will take some time",
      paste0("Can take around 3 minutes- needs around 3 GB RAM free."),
      easyClose = TRUE,
      footer = NULL
    ))
    req(input$file1, cancelOutput = T)
    req(input$id_type, cancelOutput = T)
    req(input$CTRL_cols, cancelOutput = T)
    req(input$treated_cols, cancelOutput = T)
    req(input$id_type, cancelOutput = T)
    output$georeflect_plot= renderPlotly({
      df= get_data(input)
      ids= df[, 1]
      reads= df[,c(input$CTRL_cols, 
                   input$treated_cols)]
      reads= as.data.frame(reads)
      
      reads$ids= ids
      if(input$id_type == "Entrez id"){
        colnames(reads)[ncol(reads)]= "entrezgene_id"
        
        reads= merge(conversion_table, reads, all= F)
      }
      if(input$id_type == "Ensembl id"){
        reads$ids= as.character(
          lapply(reads$ids, cut_isoforms))
        colnames(reads)[ncol(reads)]= "ensembl_gene_id"
        reads= merge(conversion_table, reads, all= F)
      }
      if(input$id_type == "HGNC"){
        colnames(reads)[ncol(reads)]= "gene"
      }
      rownames(reads)= paste0(reads$gene,";", 1:nrow(reads))
      reads= reads[,-c(1:3)]
      
      
      gsms <- paste0(
        c(
          rep("0", length(input$CTRL_cols)),
          rep("1", length(input$treated_cols))
        ),
        collapse = ""
      )
      sml <- strsplit(gsms, split="")[[1]]
      
      # filter out excluded samples (marked as "X")
      sel <- which(sml != "X")
      sml <- sml[sel]
      
      # group membership for samples
      gs <- factor(sml)
      groups <- make.names(c("CTRL","Treated"))
      levels(gs) <- groups
      sample_info <- data.frame(Group = gs, row.names = colnames(reads))
      
      # pre-filter low count genes
      # keep genes with at least N counts > 10, where N = size of smallest group
      keep <- rowSums( reads >= 10 ) >= min(table(gs))
      reads <- reads[keep, ]
      
      ds <- DESeqDataSetFromMatrix(countData=reads, colData=sample_info, design= ~Group)
      
      ds <- DESeq(ds, test="Wald", sfType="poscount")
      
      # extract results for top genes table
      r <- results (ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")
      DEG_frame= as.data.frame(r)
      DEG_frame$Gene= rownames(DEG_frame)
      DEG_frame$Gene= as.character(lapply(DEG_frame$Gene, 
                                          function(x){
                                            unlist(strsplit(x, ";"))[1]
                                          }))
      
      DEG_frame= DEG_frame[, c("Gene", "pvalue", "log2FoldChange", "baseMean")]
      GEOreflect_output= GEOreflect_reranking_RNA_seq(the_frame= DEG_frame,
                                           minlogfc= input$minlogfc,
                                           pvallim= input$pvallim,
                                           maxlogfc= input$maxlogfc,
                                           unmatched_bool= input$unmatched)
      
      
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
      
      
    })
  })
  
  
  
  
  
  
#ends here  
}
# Create Shiny app ----
shinyApp(ui, server)