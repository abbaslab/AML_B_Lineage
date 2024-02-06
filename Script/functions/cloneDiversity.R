
cloneDiversity <- function(clone_data,col_name,colors,legend_name){
  
  clone_data <- na.omit(clone_data[,c(col_name,"clone_id")])
  
  diversity <- alphaDiversity(clone_data, group=col_name, clone="clone_id",uniform=T,
                              min_q=0, max_q=4, step_q=0.1,ci=0.95, nboot=200)
  eval(parse(text=paste0("diversity@diversity$",`col_name`," <- factor(diversity@diversity$",`col_name`,",levels=names(colors))")))
  
  p <- plotDiversityTest(diversity,2, main_title=NULL, legend_title=legend_name,colors=colors)+
    ylab("BCR clonotype diversity\n(1/Simpson's index, mean±SD)")+
    theme(legend.position="none",
          axis.title.y = element_text(size = 9),
          axis.title.x = element_text(size = 9),
          axis.text.x = element_text(size = 8,angle =45, vjust = 1, hjust=1),
          #axis.text.x = element_text(size = 7,angle =90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 8),
          panel.border=element_rect(size=1))
  
  return(p)
}
