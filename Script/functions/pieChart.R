pieChart <- function(data,group,title,color,legend_name){
  
  ggpie(data = data, group_key = group, count_type = "full",
        border_size=0.2,border_color = "gray40",
        label_info = "all", label_type = "none",
        label_size = 3, label_pos = "in",label_split = NULL)+
    ggtitle(title)+
    theme(title = element_text(size=8),
          legend.text = element_text(size = 5),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"),
          #legend.position = "none"
          )+
    scale_fill_manual(values=color,name=legend_name)
}
