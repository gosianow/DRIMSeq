colorb <- function(n){
  
  clrs <- c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1" ,
    "blueviolet", "firebrick2", "deepskyblue",  "orchid2", "chartreuse3", 
    "gold", "slateblue1", "tomato" , "blue", "magenta", "green3", "yellow", 
    "purple3", "red" ,"darkslategray1", "lightpink1", "lightgreen", "khaki1", 
    "plum3", "salmon")
  
  nc <- length(clrs)
  
  if(n > nc)
    clrs <- rep(clrs, ceiling(n/nc))
  
  clrs[1:n]
  
  # colorRampPalette(clrs)(n)
  
}


# nb <- 24
# barplot(rep(1, nb), col = colorb(nb))
# dev.off()
