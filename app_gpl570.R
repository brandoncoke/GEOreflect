library(shiny)
require(openxlsx)
require(ggplot2)
require(DT)
require(plotly)
require(limma)
require(GEOquery)
#60 MB limit
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
                          selectizeInput(
                            'GEO_id', 'Select a GEO id', choices = NULL,
                            multiple = TRUE, options = list(maxItems = 1)
                          ),
                        
                          # Horizontal line ----
                          tags$hr(),
                          selectizeInput(
                            'CTRL_cols', 'Select control columns', choices = NULL,
                            multiple = TRUE),
                          selectizeInput(
                            'treated_cols', 'Select treated columns',
                            choices = NULL, multiple = TRUE),
                          checkboxInput("unmatched", "Ignore unmatched genes",
                                        value= T)
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

load('percentile_matrix.RDS')
get_platform_percentile_GPL570= function(probe_and_pvalue= c("222589_at", 4.43E-09)){
  if(as.numeric(probe_and_pvalue[2])== 0 | 
     !(probe_and_pvalue[1] %in% rownames(percentile_matrix))){
    return(1)
  }else{
    gene= c(t(probe_and_pvalue[1])) #why does this work no clue but gets the string i need- t is transpose so eh#
    which.min(as.numeric(percentile_matrix[rownames(percentile_matrix) == gene, ]) <
                as.numeric(probe_and_pvalue[2])) #why do I need to coerce to float- how knows but safer as input can be treated as a string
  }
}
geometric_mean= function(values= c(1,2,1)){
  (prod(values))^(1/length(values))
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
GEOreflect_reranking_GPL570= function(the_frame,
                                      unmatched_bool= T,
                                      minlogfc= -1,
                                      pvallim= 0.05,
                                      maxlogfc=1,
                                      merge_probe_gene= F){
  temp= the_frame
  colnames(temp)= c("probe", "pvalues", "logfc", "genes", 
                    "average_exprss")
  temp= temp[temp$genes != "", ]
  
  
  if(class(temp$pvalues) != "numeric" |
     class(temp$logfc) != "numeric" |
     class(temp$probe) != "character" |
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
      temp= temp[temp$probe %in% rownames(percentile_matrix), ]
    }
    
    if(!any(temp$probe %in% rownames(percentile_matrix))
       | nrow(temp) == 0){
      
      if(nrow(temp) == 0){
        showModal(modalDialog(
          title = "Filtering step too stringent",
          paste0("Increase the p-value limit or bounds for upper and lower log fold change"),
          easyClose = TRUE,
          footer = NULL
        ))
      }
      
      return(data.frame(Probe= temp$probe,
                        Genes= temp$genes,
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
      temp$aver_exprss_normalised= rank(temp$average_exprss)/(length(temp$average_exprss))
      
      temp$plat_rank= as.integer(apply(temp, 1,
                                       get_platform_percentile_GPL570))
      temp$GEOreflect= ((1000 - temp$plat_rank)+1)/ 1000
      temp$comb_score= apply(temp[,c("pval_normalised",
                                     "aver_exprss_normalised",
                                     "GEOreflect")], 1, geometric_mean)
      
      output= data.frame(Probe= temp$probe,
                         Genes= temp$gene,
                         logFC= temp$logfc,
                         'p-values'= temp$pvalues,
                         pval_rank= rank(temp$pvalues),
                         Platform_relative_rank= temp$plat_rank,
                         Average_expression_rank= temp$average_exprss,
                         GEOreflect_rank= rank(-temp$comb_score))

      
      
      median_shift= as.numeric(
        quantile(abs(output$pval_rank - output$GEOreflect_rank), 
                 0.75))
      output= output[order(output$GEOreflect_rank),]
      output$Rank_change= output$GEOreflect_rank - output$pval_rank
      
      output$Shift= apply(output, 1, shift_label, median_shift)
      if(merge_probe_gene){
        output$Genes= paste0(output$Genes, ": ", output$Probe)
      }
      
      return(output)
    }
  }
  
  
}


get_DEG_table= function(gset, gsms){
  #same procedure as GEO2R
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  sml <- strsplit(gsms, split="")[[1]]
  
  # filter out excluded samples (marked as "X")
  sel <- which(sml != "X")
  sml <- sml[sel]
  gset <- gset[ ,sel]
  
  # log2 transformation
  ex <- exprs(gset)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
  
  # assign samples to groups and set up design matrix
  gs <- factor(sml)
  groups <- make.names(c("1","2"))
  levels(gs) <- groups
  gset$group <- gs
  design <- model.matrix(~group + 0, gset)
  colnames(design) <- levels(gs)
  
  fit <- lmFit(gset, design)  # fit linear model
  
  # set up contrasts of interest and recalculate model coefficients
  cts <- paste(groups[1], groups[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=5000000)
  
  tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC", "Gene.ID",
                            "Gene.symbol","Gene.title"))
  average_exprss= as.data.frame(apply(ex, 1, sum, na.rm= T))
  average_exprss$ID= rownames(average_exprss)
  colnames(average_exprss)[1]= "average_exprss"
  tT= merge(tT, average_exprss)
  
  return(tT)
}
treated_labels = c("overexp", "express", "transgen", "expos", "tg", "induc",
                  "stim", "treated", "transfected", "overexpression",
                  "transformed", "tumor", "tomour", "disease",
                  "infect", "disorder",
                  "knock", "null",
                     "s[hi]rna",
                     "delet",
                     "si[a-zA-z]",
                     "reduc",
                     "kd",
                     "\\-\\/",
                     "\\/\\-",
                     "\\+\\/", "\\/\\+",
                     "cre", "flox",
                     "mut",
                     "defici",
                     "[_| ]ko[_| ]|[_| ]ko$")
control_labels <- c("untreat", "_ns_",
                    "normal",
                    "^wt$|^wt[_| ]|[_| ]wt[_| ]|[_| ]wt$",
                    "gfp",
                    "non[-|.|_| ]smoker",
                    "vehicle",
                    "sensitive",
                    "stable",
                    "ctrl",
                    "non.sense",
                    "nonsense",
                    "baseline",
                    "mock",
                    "_luc_|_luc",
                    "siluc",
                    "wildtype|wild.type",
                    "nontreat",
                    "non.treated",
                    "control",
                    "[ |_]con[ |_]|^con$|^con[_| ]|[_| ]con$", #con sometimes stand in for control- needs to be specific- no futher words before or after
                    "untreated",
                    "no.treat",
                    "undosed",
                    "reference",
                    "standard",
                    "untransfected",
                    "mir.nc",
                    "minc",
                    "_non_|^non_",
                    "scramble",
                    "lucif",
                    "parent",
                    "free.medi",
                    "untransfected",
                    "ntc",
                    "mirNC")
#Platform_meta <- GEOquery::getGEO('GPL570')@header[["series_id"]] 
#otherwise if API changed too much
Platform_meta= read.csv('Platform_meta.csv')[,1]

server <- function(input, output, session) {
  updateSelectizeInput(session, 'GEO_id', 'Select a GEO id',
                       choices = Platform_meta,
                       server = TRUE)
  output$contents <- DT::renderDataTable({
    req(input$GEO_id)
    gset <- getGEO(input$GEO_id, GSEMatrix =T, AnnotGPL=T)
    if(length(gset) > 1 & any(grepl("GPL570", attr(gset, "names")))){
      idx= grep('GPL570', attr(gset, "names")) #bit of a bodge but assumes platform specified is enough to distinguish similar ones
    }else{
      idx= 1}
    gset <- gset[[idx]]
    # make proper column names to match toptable
    fvarLabels(gset) <- make.names(fvarLabels(gset))
    df= data.frame(GEO_id= as.character(gset@phenoData@data$geo_accession),
      Title= as.character(gset@phenoData@data$title),
                   Source= as.character(gset@phenoData@data$source_name_ch1),
                   Characteristics= as.character(
                     gset@phenoData@data$characteristics_ch1))
    the_samples= as.character(apply(df, 1, paste0, collapse= ""))
    
    default_CTRL= NULL
    if(any(grepl(paste0(control_labels, collapse = "|"), the_samples))){
      default_CTRL= 
        df$GEO_id[grep(paste0(control_labels, collapse = "|"), the_samples)]
      
    }
    default_treated= NULL
    if(any(grepl(paste0(treated_labels, collapse = "|"), the_samples))){
      default_treated= 
        df$GEO_id[grep(paste0(treated_labels, collapse = "|"), the_samples)]
      
    }
    
    updateSelectizeInput(session,
                         "CTRL_cols", 'Select control samples',
                         choices= df$GEO_id,
                         selected= default_CTRL)
    updateSelectizeInput(session,
                         "treated_cols", 'Select treated samples',
                         choices= df$GEO_id,
                         selected= default_treated)

    return(df)
    
  })

  observeEvent(input$georeflect, {
    req(input$GEO_id, cancelOutput = T)
    req(input$treated_cols, cancelOutput = T)
    req(input$CTRL_cols, cancelOutput = T)
    showModal(modalDialog(
      title = "GEOreflect is running. This may take some time",
      paste0("Can take around 3 minutes- needs around 3 GB RAM free."),
      easyClose = TRUE,
      footer = NULL
    ))
    output$georeflect_frame= DT::renderDataTable({
      gset <- getGEO(input$GEO_id, GSEMatrix =T, AnnotGPL=T)
      if(length(gset) > 1 & any(grepl("GPL570", attr(gset, "names")))){
        idx= grep('GPL570', attr(gset, "names")) #bit of a bodge but assumes platform specified is enough to distinguish similar ones
      }else{
        idx= 1}
      gset <- gset[[idx]]
      # make proper column names to match toptable
      
      
      fvarLabels(gset) <- make.names(fvarLabels(gset))
      GEO_ids= as.character(gset@phenoData@data$geo_accession)
      gsms= rep("X", length(GEO_ids))
      gsms[which(GEO_ids %in% input$CTRL_cols)]= 0
      gsms[which(GEO_ids %in% input$treated_cols)]= 1
      gsms= paste0(gsms, collapse= "")
      the_frame= get_DEG_table(gset = gset,
                    gsms= gsms)
      the_frame= the_frame[c("ID", "P.Value",
                             "logFC", "Gene.symbol",
                             "average_exprss")]
      output= GEOreflect_reranking_GPL570(the_frame= the_frame,
                                  unmatched_bool= input$unmatched,
                                  minlogfc= input$minlogfc,
                                  pvallim= input$pvallim,
                                  maxlogfc= input$maxlogfc,
                                  merge_probe_gene= F)
      
      return(output)
      
    })
  })
  
  
  
  
  
  
  observeEvent(input$tha_plot_button, {
    output$georeflect_plot= renderPlotly({
      if (input$tha_plot_button) {
        req(input$GEO_id, cancelOutput = T)
        req(input$treated_cols, cancelOutput = T)
        req(input$CTRL_cols, cancelOutput = T)
        
        
        showModal(modalDialog(
          title = "Plotting data- will take some time",
          paste0("Can take around 3 minutes- needs around 3 GB RAM free."),
          easyClose = TRUE,
          footer = NULL
        ))
        gset <- getGEO(input$GEO_id, GSEMatrix =T, AnnotGPL=T)
        if(length(gset) > 1 & any(grepl("GPL570", attr(gset, "names")))){
          idx= grep('GPL570', attr(gset, "names")) #bit of a bodge but assumes platform specified is enough to distinguish similar ones
        }else{
          idx= 1}
        gset <- gset[[idx]]
        # make proper column names to match toptable
        
        
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        GEO_ids= as.character(gset@phenoData@data$geo_accession)
        gsms= rep("X", length(the_samples))
        gsms[which(GEO_ids %in% input$CTRL_cols)]= 0
        gsms[which(GEO_ids %in% input$treated_cols)]= 1
        gsms= paste0(gsms, collapse= "")
        the_frame= get_DEG_table(gset = gset,
                                 gsms= gsms)
        the_frame= the_frame[c("ID", "P.Value",
                               "logFC", "Gene.symbol",
                               "average_exprss")]
        GEOreflect_output= GEOreflect_reranking_GPL570(the_frame= the_frame,
                                            unmatched_bool= input$unmatched,
                                            minlogfc= input$minlogfc,
                                            pvallim= input$pvallim,
                                            maxlogfc= input$maxlogfc,
                                            merge_probe_gene= T)
        
        
        
        
        
        
        
        
        
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
                             aes(y=`p-value rank`, x=`p-value rank`, label=Gene)) +
              geom_point(size= 3, 
                         color= "#008080", 
                         aes(alpha= opacity)) 
          }else{
            tha_plot= ggplot(GEOreflect_output,
                             aes(y=`p-value rank`, x=`p-value rank`, label=Gene)) +
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
        }
        
        
        
        
        
      }
      
      
      
      
      
      
      
    })})
  #Ends here   
}
# Create Shiny app ----
shinyApp(ui, server)