### Functions that are used in show methods


################################################################################

#' @importFrom utils head tail

show_matrix <- function(object, nhead = 2, ntail = 2){
  # object is a matrix
  
  nr <- nrow(object)
  nc <- ncol(object)
  
  cat(class(object), " with ", nr, ifelse(nr == 1, " row and ", " rows and "), 
    nc, ifelse(nc == 1, " column\n", " columns\n"), sep = "")
  
  if(nr > 0 && nc > 0){
    
    if(is.null(colnames(object))){
      colnames(object) <- paste0("[,", 1:ncol(object), "]")
    }
    if(is.null(rownames(object))){
      rownames(object) <- paste0("[", 1:nrow(object), ",]")
    }
    
    if(nr <= (nhead + ntail)){
      
      out <- object
      
    }else{
      
      out <- do.call(rbind, list(head(object, nhead), matrix(rep.int("...", nc), 
        1, nc, dimnames = list(NULL, colnames(object))), tail(object, ntail)))
      
      nms <- rownames(object)
      if(nhead > 0)
        s1 <- paste0(head(nms, nhead))
      if(ntail > 0)
        s2 <- paste0(tail(nms, ntail))
      
      rownames(out) <- c(s1, "...", s2)
      
    }
    
    if(nc > (nhead + ntail)){
      
      out <- do.call(cbind, list(out[, 1:nhead, drop = FALSE], 
        matrix(rep.int("...", 
          ifelse(nr < (nhead + ntail + 1L), min(nr, nhead + ntail), 
            nhead + ntail + 1L)), ncol = 1, dimnames = list(NULL, "...")), 
        out[, (nc-ntail+1):nc, drop = FALSE]))
      
    }   
    
    ### print adjusted for numeric or character
    if(class(object[1, 1]) == "numeric"){
      print(out, quote = FALSE, right = TRUE, na.print = "NA") 
    }else{
      print(out, quote = TRUE, right = TRUE, na.print = "NA") 
    }
    
  }
  
}


################################################################################

show_numeric <- function(object, nhead = 2, ntail = 2, class = TRUE, print = TRUE){
  
  nl <- length(object)  
  if(class)
    cat(class(object), "of length", length(object), "\n")
  
  if(nl > 0){
    
    if(nl < (nhead + ntail + 1L)) {
      out <- round(object, 2)
    } else {
      dots <- "..."
      if(!is.null(names(object)))
        names(dots) <- "..."
      out <- c(round(head(object, nhead), 2), dots , round(tail(object, ntail), 2))
    }
    if(print)
      print(out, quote = FALSE, right = TRUE)
    else
      return(out)
    
  }else{
    
    if(print)
      print(object)
    else
      return(object)
    
  }
  
}


################################################################################


show_numeric_list <- function(object, nhead = 2){
  
  nl <- length(object)  
  cat(class(object), "of length", nl, "\n")
  
  if(nl > 0){
    np <- min(nl, nhead)
    
    object <- object[1:np]
    
    if(is.null(names(object)))
      print_names <- paste0("[[", 1:np, "]]\n")
    else 
      print_names <- paste0("$", names(object), "\n")
    
    for(i in 1:np){
      
      cat(print_names[i])
      show_numeric(object[[i]])
      cat("\n")
      
    }
    
    if(np < nl){
      
      if(is.null(names(object)))
        cat(paste0("[[...]]\n"))
      else 
        cat(paste0("$...\n"))
      
    }
    
  }else{
    
    print(object)
    
  }
  
}



################################################################################


show_MatrixList_list <- function(object, nhead = 2){
  
  nl <- length(object)  
  cat(class(object), "of length", nl, "\n")
  
  if(nl > 0){
    np <- min(nl, nhead)
    
    object <- object[1:np]
    
    if(is.null(names(object)))
      print_names <- paste0("[[", 1:np, "]]\n")
    else 
      print_names <- paste0("$", names(object), "\n")
    
    for(i in 1:np){
      
      cat(print_names[i])
      print(object[[i]])
      cat("\n")
      
    }
    
    if(np < nl){
      
      if(is.null(names(object)))
        cat(paste0("[[...]]\n"))
      else 
        cat(paste0("$...\n"))
      
    }
    
  }else{
    
    print(object)
    
  }
  
}






