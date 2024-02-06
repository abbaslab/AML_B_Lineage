
monocle3_fuction <- function(cds,func,n1=3,n2=5,n3=20){
  
  colour_bk <- c("#f0f0f0",colorRampPalette(c("#006837","#d9ef8b"))(n1),
                          colorRampPalette(c("#d9ef8b","#fee08b"))(n2),
                          colorRampPalette(c("#fee08b","#a50026"))(n3))
  plot_cells(cds, color_cells_by=func, 
             reduction_method="UMAP",show_trajectory_graph=T,
             trajectory_graph_color = "grey60",trajectory_graph_segment_size = 0.75,
             label_roots = F,label_leaves = F,label_cell_groups =F,label_branch_points = F,
             alpha=0.8,cell_size = 0.4,group_label_size = 10)+
    xlab("")+ylab("")+ggtitle("")+del_axis()+
    scale_colour_gradientn(colours = colour_bk)+
    theme(legend.title=element_text(size = 8),
          legend.text=element_text(size = 6),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"))
}



