#funcitons to do reranking
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

GEOreflect_reranking_GPL570= function(the_frame, pvalue_indice, gene_indice,
                                      logfc_indice, probe_indice=1,
                                      unmatched_bool= T,
                                      minlogfc= -1,
                                      pvallim= 0.05,
                                      average_exprss_indice,
                                      maxlogfc=1,
                                      merge_probe_gene= F){
  temp= the_frame[, c(probe_indice,
                      pvalue_indice,
                      logfc_indice,
                      gene_indice)]
  
  colnames(temp)= c("probe", "pvalues", "logfc", "genes")
  temp= temp[temp$genes != "", ]
  
  if(class(temp$pvalues) != "numeric" |
     class(temp$logfc) != "numeric" |
     class(temp$probe) != "character" |
     sum(is.na(temp$pvalues) == nrow(temp)) |
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
    
    if(!any(temp$genes %in% rownames(percentile_matrix_p_value_RNAseq))
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
      temp$comb_score= apply(temp[,c("GEOreflect", "pval_normalised", 
                                     "express_norm")], 1, geometric_mean)
      
      output= data.frame(Probe= temp$probe,
                         Genes= temp$gene,
                         logFC= temp$logfc,
                         'p-values'= temp$pvalues,
                         average_expression_rank= rank(-temp$average_expression),
                         pval_rank= rank(temp$pvalues),
                         Platform_relative_rank= temp$plat_rank,
                         GEOreflect_rank= rank(-temp$comb_score))
      
      
      
      median_shift= as.numeric(
        quantile(abs(output$pval_rank - output$GEOreflect_rank), 
                 0.75))
      output= output[order(output$GEOreflect_rank),]
      output$Rank_change= output$GEOreflect_rank - output$pval_rank
      
      output$Shift= apply(output[,c('pval_rank',
                                    'GEOreflect_rank')], 1, shift_label, median_shift)
      if(merge_probe_gene){
        output$Genes= paste0(output$Genes, ": ", output$Probe)
      }
      
      return(output)
    }
  }
  
  
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
      
      output$Shift= apply(output[,c('pval_rank',
                                    'GEOreflect_rank')], 1, shift_label, median_shift)
      return(output)
    }
  }
}

GEOreflect_reranking_GPL570= function(the_frame, pvalue_indice, gene_indice,
                                      logfc_indice, probe_indice=1,
                                      average_exprss_indice,
                                      unmatched_bool= T,
                                      minlogfc= -1,
                                      pvallim= 0.05,
                                      maxlogfc=1){
  temp= the_frame[, c(probe_indice,
                      pvalue_indice,
                      logfc_indice,
                      gene_indice)]
  
  colnames(temp)= c("probe", "pvalues", "logfc", "genes")
  temp= temp[temp$genes != "", ]
  
  if(class(temp$pvalues) != "numeric" |
     class(temp$logfc) != "numeric" |
     class(temp$probe) != "character" |
     sum(is.na(temp$pvalues) == nrow(temp)) |
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
    
    if(!any(temp$genes %in% rownames(percentile_matrix_p_value_RNAseq))
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
      temp$comb_score= apply(temp[,c("GEOreflect", "pval_normalised", 
                                     "express_norm")], 1, geometric_mean)
      
      output= data.frame(Probe= temp$probe,
                         Genes= temp$gene,
                         logFC= temp$logfc,
                         'p-values'= temp$pvalues,
                         average_expression_rank= rank(-temp$average_expression),
                         pval_rank= rank(temp$pvalues),
                         Platform_relative_rank= temp$plat_rank,
                         GEOreflect_rank= rank(-temp$comb_score))
      
      
      
      median_shift= as.numeric(
        quantile(abs(output$pval_rank - output$GEOreflect_rank), 
                 0.75))
      output= output[order(output$GEOreflect_rank),]
      output$Rank_change= output$GEOreflect_rank - output$pval_rank
      
      output$Shift= apply(output[,c('pval_rank',
                                    'GEOreflect_rank')], 1, shift_label, median_shift)
      return(output)
    }
  }
  
  
}






#examples RNA-seq
load("percentile_matrix_p_value_RNAseq.RDS") #For RNAseq data
RNA_seq_DEG_frame= read.csv("example_RNAseq_DEG_frame.csv") #this can be found in Git repo
#get colnames
colspresent= colnames(RNA_seq_DEG_frame);
RNA_seq_reranked= GEOreflect_reranking_RNA_seq(the_frame= RNA_seq_DEG_frame, #DEG frame
                             pvalue_indice= 4, #p-value column indice
                             gene_indice= 7, #Gene column indice
                             logfc_indice= 1, #logFC column indice
                             average_exprss_indice = 2, #average expression column
                             minlogfc= -1, #logfc needs to be lower than this
                             pvallim= 0.05, #p-value needs to be lower than this
                             maxlogfc=1,  #logfc needs to be higher than this
                             unmatched_bool= T) #remove genes not present in percentile frame
RNA_seq_reranked

#GPL570
load('percentile_matrix.RDS') #For GPL570 data
GPL570_DEG_frame= read.csv("example_GPL570_DEG_frame.csv") #this can be found in Git repo
GEOreflect_reranking_GPL570(the_frame= GPL570_DEG_frame, #DEG frame
                             pvalue_indice= 4, #p-value column indice
                             gene_indice= 9, #Gene column indice
                             logfc_indice= 7, #logFC column indice
                            average_exprss_indice= 2,
                            probe_indice = 1,
                             minlogfc= -1, #logfc needs to be lower than this
                             pvallim= 0.05, #p-value needs to be lower than this
                             maxlogfc=1,  #logfc needs to be higher than this
                             unmatched_bool= T) #remove genes not present in percentile frame


#alterntively use this to automate the column indice selection
#logFC defualts
if(any(grepl("expr|logfc|fc|fold|log_", as.character(colspresent),
             ignore.case = T))){
  logfc_column= grep("expr|logfc|fc|fold|log_", as.character(colspresent),
                     ignore.case = T)[1]
  
  #if it has multiple fcs
  if(any(grepl("log", as.character(colspresent), ignore.case = T)
         &
         grepl("fold", as.character(colspresent), ignore.case = T)
  )
  ){
    bool= grepl("log", as.character(colspresent), ignore.case = T) &
      grepl("fold", as.character(colspresent), ignore.case = T)
    
    logfc_column= which(bool)[1]
  }
}else{
  logfc_column= 2
}

#Symbol defaults
if(any(grepl("gene|hgnc|symbol|ident", as.character(colspresent)),
       ignore.case = T)
){
  gene_column=
    grep("gene|hgnc|symbol|ident", as.character(colspresent),
         ignore.case = T)[1]
  #if it actually has the hngc symbols
  if(any(grepl("hgnc|symbol", as.character(colspresent)),
         ignore.case = T)){
    gene_column=
      grep("hgnc|symbol", as.character(colspresent), ignore.case = T)[1]
  }
  
}else{
  gene_column= 1
}

pvalue_column=3
if(any(grepl("p-val|pval|p[.]val", as.character(colspresent), 
             ignore.case = T))){
  
  pvalue_column=
    which(grepl("p-val|pval|p[.]val", as.character(colspresent), 
                ignore.case = T))
  pvalue_column= pvalue_column[1]
  
}
#better if using non adjusted p-values- this switches to it if 
#both p-values and adjusted p-values
if(any(grepl("p-val|pval|p[.]val", as.character(colspresent), 
             ignore.case = T) & 
       !grepl("adj", as.character(colspresent), 
              ignore.case = T))){
  
  pvalue_column=
    which(grepl("p-val|pval|p[.]val", as.character(colspresent), 
                ignore.case = T) & 
            !grepl("adj", as.character(colspresent), 
                   ignore.case = T))
  pvalue_column= pvalue_column[1]
  
}

probe_column= 1
if(any(grepl("ID|_at|probe|affy|hg_u133", as.character(colspresent),
             ignore.case = T))){
  grep("p-val|adj[_|-]p|pval", as.character(colspresent), ignore.case = T)[1]
}else{
  probe_column= 1
}
expression_column=2
if(any(grepl("avg_expression|expression", as.character(colspresent), 
             ignore.case = T))){
  
  expression_column=
    which(grepl("avg[_|.| ]expr|expression|AveExpr|avex", as.character(colspresent), 
                ignore.case = T))
  expression_column= expression_column[1]
  
}


GEOreflect_reranking_RNA_seq(the_frame= RNA_seq_DEG_frame, #DEG frame
                             pvalue_indice= pvalue_column, #p-value column indice
                             gene_indice= gene_column, #Gene column indice
                             logfc_indice= logfc_column, #logFC column indice
                             average_exprss_indice= expression_column,
                             minlogfc= -1, #logfc needs to be lower than this
                             pvallim= 0.05, #p-value needs to be lower than this
                             maxlogfc=1,  #logfc needs to be higher than this
                             unmatched_bool= T) #remove genes not present in percentile frame


