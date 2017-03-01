
#' @importFrom grDevices pdf dev.off

dmDS_plotFit <- function(gene_id, counts, samples, dispersion = numeric(), 
  proportions_full = NULL, proportions_null = NULL, table = NULL, 
  plot_type = "barplot", order = TRUE, 
  plot_full = ifelse(is.null(proportions_full), FALSE, TRUE), 
  plot_null = ifelse(is.null(proportions_null), FALSE, TRUE), 
  plot_main = TRUE, out_dir = NULL){
  
  gene <- gene_id
  counts_gene <- counts[[gene]]
  
  if(nrow(counts_gene) <= 1)
    stop("!Gene has to have at least 2 features! \n")
  
  group <- samples$group
  sample_id <- samples$sample_id
  main <- NULL
  
  if(plot_main){
    
    mean_expression_gene <- mean(colSums(counts_gene), na.rm = TRUE)
    
    main <- paste0(gene, "\n Mean expression = ", round(mean_expression_gene))
    
    if(length(dispersion) > 0){
      
      if(length(dispersion) == 1)
        dispersion_gene <- dispersion
      else
        dispersion_gene <- dispersion[gene]
      
      main <- paste0(main, ", Dispersion = ", round(dispersion_gene, 2))
      
    }
    
    if(!is.null(table)){
      
      table_tmp <- table[table$gene_id == gene, ]
      
      main <- paste0(main, "\n LR = ", round(table_tmp["lr"], 2) , 
        ", P-value = ", sprintf("%.02e", table_tmp["pvalue"]), 
        ", FDR = ", sprintf("%.02e", table_tmp["adj_pvalue"]))    
      
    }
  }
  
  prop_full <- NULL
  prop_null <- NULL
  
  if(plot_full)
    prop_full <- proportions_full[[gene]]
  if(plot_null)
    prop_null <- proportions_null[[gene]]
  
  ggp <- dm_plotProportions(counts = counts_gene, group = group, 
    prop_full = prop_full, prop_null = prop_null, main = main, plot_type = plot_type, 
    order = order)
  
  if(!is.null(out_dir)){
    pdf(paste0(out_dir, "dmfit_", gsub(pattern = "\\.", replacement = "_" , paste0(gene)), ".pdf"), 
      width = 12, height = 7)
    print(ggp)
    dev.off()
  }else{
    return(ggp)  
  }
  
}










