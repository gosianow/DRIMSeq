
#' @importFrom grDevices pdf dev.off

dmSQTL_plotFit <- function(gene_id, snp_id, counts, genotypes, blocks, samples, 
  dispersion = numeric(), fit_full = NULL, fit_null = NULL, table = NULL, 
  plot_type = "boxplot", order = TRUE, 
  plot_full = ifelse(is.null(fit_full), FALSE, TRUE), 
  plot_null = ifelse(is.null(fit_null), FALSE, TRUE), 
  plot_main = TRUE, out_dir = NULL){
  
  gene <- gene_id
  snp <- snp_id
  block <- blocks[[gene]][blocks[[gene]][, "snp_id"] == snp, "block_id"]
  counts_gene <- counts[[gene]]
  
  if(nrow(counts_gene) < 2)
    stop("!Gene has to have at least 2 features! \n")
  
  group <- genotypes[[gene]][block, ]
  
  NAs <- !(is.na(counts_gene[1,]) | is.na(group))
  counts_gene <- counts_gene[, NAs, drop = FALSE]
  group <- factor(group[NAs])
  sample_id <- samples$sample_id[NAs]
  
  main <- NULL
  
  if(plot_main){
    
    mean_expression_gene <- mean(colSums(counts_gene), na.rm = TRUE)
    
    main <- paste0(gene, " : ", snp, " : ", block,
      "\n Mean expression = ", round(mean_expression_gene))
    
    if(length(dispersion) > 0){
      
      if(class(dispersion) == "numeric")
        dispersion_gene <- dispersion
      else
        dispersion_gene <- dispersion[[gene]][block]
      
      main <- paste0(main, ", Dispersion = ", round(dispersion_gene, 2))
      
    }
    
    if(!is.null(table)){
      
      table_tmp <- table[table$gene_id == gene & 
          table$block_id == block & table$snp_id == snp, ]
      
      main <- paste0(main, "\n LR = ", round(table_tmp["lr"], 2) , 
        ", P-value = ", sprintf("%.02e", table_tmp["pvalue"]), 
        ", FDR = ", sprintf("%.02e", table_tmp["adj_pvalue"]))    
      
    }
    
  }
  
  pi_full <- fit_full[[gene]][[block]][, levels(group), drop = FALSE]
  pi_null <- fit_null[[gene]][[block]]
  
  ggp <- dm_plotProportions(counts = counts_gene, group = group, 
    pi_full = pi_full, pi_null = pi_null, main = main, 
    plot_type = plot_type, order = order)
  
  
  if(!is.null(out_dir)){
    pdf(paste0(out_dir, "dmfit_", gsub(pattern = "\\.", replacement = "_" , 
      paste0(gene, "_", snp)), ".pdf"),  width = 12, height = 7)
    print(ggp)
    dev.off()
  }else{
    return(ggp)  
  }
  
  
}










