
# These codes were borrowed from scRepertoire R package


compareClonotypes <- function(df,cloneCall = "clone_id",samples =c("PT7A","PT7B","PT7C"), 
                              clonotypes = NULL, 
                              numbers = 5, 
                              group = "sample", 
                              graph = "alluvial", 
                              exportTable = FALSE){
  df <- expression2List(df, split.by=group)
 
  Con.df <- NULL
  for (i in seq_along(df)) {
    tbl <- as.data.frame(table(df[[i]][,cloneCall]))
    tbl[,2] <- tbl[,2]/sum(tbl[,2])
    colnames(tbl) <- c("Clonotypes", "Proportion")
    tbl$Sample <- names(df[i])
    Con.df <- rbind.data.frame(Con.df, tbl)
  }
  if (!is.null(samples)) {
    Con.df <- Con.df[Con.df$Sample %in% samples,] }
  if (!is.null(clonotypes)) {
    Con.df <- Con.df[Con.df$Clonotypes %in% clonotypes,] }
  if (!is.null(numbers)) {
    top <- Con.df %>%
      group_by(Con.df[,3]) %>%
      slice_max(n = numbers, order_by = Proportion, with_ties = FALSE)
    Con.df <- Con.df[Con.df$Clonotypes %in% top$Clonotypes,] }
  if (nrow(Con.df) < length(unique(Con.df$Sample))) {
    stop("Reasses the filtering strategies here, there is not 
            enough clonotypes to examine.") }
  if (exportTable == TRUE) { return(Con.df)}
  
  plot <- ggplot(Con.df, aes(x = Sample, fill = Clonotypes, group = Clonotypes,
                             stratum = Clonotypes, alluvium = Clonotypes, 
                             y = Proportion, label = Clonotypes)) +
    theme_classic() +
    theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.text.x=element_text(size=7),
        legend.title=element_text(size = 8),
        legend.text=element_text(size = 5),
        legend.key.height = unit(0.25,"cm"),
        legend.key.width = unit(0.25,"cm"))
  if (graph == "alluvial") {
    plot = plot +  geom_stratum() + geom_flow(stat = "alluvium")
  } else if (graph == "area") {
    plot = plot +
      geom_area(aes(group = Clonotypes), color = "black") }
  return(plot)
}
