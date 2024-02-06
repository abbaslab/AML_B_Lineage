

cloneMutation <- function(clone_data,col_name,colors){
  
  clone_data <- na.omit(clone_data[,c(col_name,"clone_id","mu_freq_seq_r")])
  p <- ggplot(clone_data, aes(x=get(col_name), y=mu_freq_seq_r, fill=get(col_name))) +
    geom_boxplot(width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    ylab("Mutation frequency") +
    guides(fill=F)+labs(fill="")+xlab("")+ylim(0,0.13)+
    theme_classic() +
    theme(legend.position="none",
          axis.title.y = element_text(size = 9),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8,angle =45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 8))+
    scale_fill_manual(values=colors)
  
  return(p)

}


