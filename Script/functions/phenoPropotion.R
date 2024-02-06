

phenoPropotion <- function(data,x_name,y_name,cols,legend_name){
  
  df_stats <- data %>% 
    dplyr::group_by(get(x_name),get(y_name),.drop = FALSE) %>% 
    dplyr::summarise(n=n())
  df_stats <- as.data.frame(df_stats)
  names(df_stats)[1:2] <- c(x_name,y_name)
  df_stats <- na.omit(df_stats)
  df_stats[,y_name] <-  factor(as.character(df_stats[,y_name]),levels=names(cols))

  p <- ggplot(df_stats, aes(x = get(x_name), y = n, fill = get(y_name))) + 
    geom_bar(position = "fill",stat = "identity",na.rm=T,width=0.8) + #position="stack" gives numbers
    scale_fill_manual(values= cols) +
    theme_classic()+ 
    theme(legend.title=element_text(size = 8),
          legend.text=element_text(size = 6),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"),
          axis.title.y = element_text(size = 9),
          axis.title.x = element_text(size = 8),
          #axis.text.x = element_text(size = 8,angle =30, vjust = 1, hjust=0.75),
          axis.text.x = element_text(size = 7,angle =90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 8))+
    xlab("")+ylab("Cell Fraction (%)")+labs(fill=legend_name)+
    scale_y_continuous(labels = scales::percent)
  
  return(p) 
}

phenoPropotion_flip <- function(data,x_name,y_name,cols,legend_name){
  
  df_stats <- data %>% 
    dplyr::group_by(get(x_name),get(y_name),.drop = FALSE) %>% 
    dplyr::summarise(n=n())
  names(df_stats)[1:2] <- c(x_name,y_name)
  df_stats <- na.omit(df_stats)
  df_stats <- as.data.frame(df_stats)
  df_stats[,y_name] <-  factor(as.character(df_stats[,y_name]),levels=names(cols))
  
  p <- ggplot(df_stats, aes(x = get(x_name), y = n, fill = get(y_name))) + 
    geom_bar(position = "fill",stat = "identity",na.rm=T,width=0.8) + #position="stack" gives numbers
    scale_fill_manual(values= cols) +
    theme_classic()+ 
    theme(legend.position="top",
          legend.title=element_text(size = 8),
          legend.text=element_text(size = 6),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8))+
    xlab("")+ylab("Cell Fraction (%)")+labs(fill=legend_name)+
    scale_y_continuous(labels = scales::percent)+
    coord_flip()
 
   return(p)
}



