
# These codes were borrowed from scRepertoire R package

clonalProportion <- function(clone_data,split = c(10,100,500,1000),cloneCall = "clone_id", 
                             group,group_order,exportTable = FALSE) {

  df <- expression2List(clone_data, split.by=group)
 
  mat <- matrix(0, length(df), length(split), 
                dimnames = list(names(df),paste0('[', c(1, split[-length(split)] + 1), ':', split, ']')))

  df <- lapply(df, '[[', cloneCall)
  df <- lapply(df, na.omit)
  df <- lapply(df, as.data.frame(table))
  
  mat <- mat[lapply(df,nrow) > 0,]
  df <- df[lapply(df,nrow) > 0]
  
  
  for (i in seq_along(df)) {
    df[[i]] <- rev(base::sort(as.numeric(df[[i]][,2])))
    #print(i)
  }

  cut <- c(1, split[-length(split)] + 1)
  for (i in seq_along(split)) {
    mat[,i] <- vapply(df, function (x) sum(na.omit(x[cut[i]:split[i]])), 
                      FUN.VALUE = numeric(1))
  }
  
  if (exportTable == TRUE) {
    return(mat)
  }
  
  mat_melt <- reshape2::melt(mat)
  col <- length(unique(mat_melt$Var2))
  mat_melt$Var1 <- factor(mat_melt$Var1,levels=group_order)
  plot <- ggplot(mat_melt, aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat = "identity", position="fill",color = "gray60", lwd= 0.15) +
    scale_fill_manual(name = "Clonal Indices",values = colorblind_vector(col)) +
    xlab("")+ylab("Cell Fraction (%)") +
    theme_classic()+
    theme(legend.title=element_text(size = 8),
          legend.text=element_text(size = 6),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"),
          axis.title.y = element_text(size = 9),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 8,angle =90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 8))
    
  
  return(plot)
}



expression2List <- function(clone_data, split.by) {
  
  if(inherits(x=clone_data, what ="Seurat")){
    meta <- sc@meta.data
  }
  if(inherits(x=clone_data, what ="data.frame")){
    meta <- clone_data
  }
  
  
  if(is.null(split.by)){
    split.by <- "cluster"
  }
  unique <- str_sort(as.character(unique(meta[,split.by])), numeric = TRUE)
  df <- NULL
  for (i in seq_along(unique)) {
    subset <- subset(meta, meta[,split.by] == unique[i])
    subset <- subset(subset, !is.na(clone_id))
    df[[i]] <- subset
  }
  names(df) <- unique
  return(df)
}


colorblind_vector <- colorRampPalette(c("#d62728","#ffbb78","#98df8a","#2ca02c"))




