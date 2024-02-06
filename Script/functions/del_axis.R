

# delete the axises of the ggplot
del_axis <- function(){
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_blank(),
        axis.line.y=element_blank()) #+
    #xlab("")+ylab("")+ggtitle("")
}
  



