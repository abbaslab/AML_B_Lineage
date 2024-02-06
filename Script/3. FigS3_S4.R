
# figS3 and S4, B cells of AML patients by transcriptome features and BCR repertoire & Association between B cells and clinical manifestations

# figS3A
  object_mem <- subset(object_b,subset=celltype %in% c("CD27+ memory B","Atypical memory B"))
  curve <- estimateAbundance(object_mem@meta.data, group="celltype", ci=0.95, nboot=100, clone="clone_id")
  mem_colors <- c(`CD27+ memory B`="#fe9929", `Atypical memory B`="#c49c94")
  plot(curve, colors = mem_colors, legend_title="CellType")
  
  
# figS3B
  cloneDiversity(clone_data=object_mem@meta.data,col_name="celltype",legend_name="CellType",colors=mem_colors)  


# figS3C
  top_clone <- countClones(object_b@meta.data)
  top_clone <- top_clone %>% 
    slice_max(seq_count,n=30,with_ties=T)
  # delete samples with < 40 cells
  sample_use <- names(table(object_b$sample)[table(object_b$sample)>=40])
  object_b_top <- subset(object_b,subset=sample %in% sample_use)
  object_b_top$top_clone <- ""
  object_b_top$top_clone <- object_b_top$clone_id[object_b_top$clone_id %in% top_clone$clone_id]
  DimPlot(object_b_top,label =F,group.by = "top_clone",cols=ColAssign(unique(object_b_top$top_clone)),
          pt.size=0.5,na.value="gray95",order=T)+
    labs(col="Sample")+ggtitle("")+
    theme(legend.title=element_text(size = 8),
          legend.text=element_text(size = 8),
          legend.key.height = unit(0.1,"cm"),
          legend.key.width = unit(0.1,"cm"))+NoAxes()
  
  
# figS3D
  load(file="Data/patient_obs.rda")
  object_allu <- subset(object_b,subset=sample_type == "RelRef")
  object_allu <- subset(object_allu,subset=celltype %in% c("Pro/Pre-B","Immature B"),invert=T)
  object_allu@meta.data <- left_join(object_allu@meta.data,patient_obs[,c("cell_id","cdr3")])
  row.names(object_allu@meta.data) <- object_allu$cell_id
  
  compareClonotypes(df=object_allu@meta.data,samples=c("PT7A","PT7B","PT7C"),cloneCall = "cdr3")
  
  
# figS3E
  family_V <- countGenes(object_b@meta.data[!is.na(object_b$v_call),], 
                         gene="v_call", groups="celltype", mode="family")
  names(family_V)[1] <- "CellType"
  g1 <- ggplot(family_V, aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle("IGHV gene usage") +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous(labels=percent) +
    scale_color_manual(values=celltype.colors) +
    geom_point(aes(color=CellType), size=1, alpha=0.8)
  
  family_D <- countGenes(object_b@meta.data[!is.na(object_b$d_call),], 
                         gene="d_call", groups="celltype", mode="family")
  names(family_D)[1] <- "CellType"
  g2 <- ggplot(family_D[family_D$gene!="",], aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle("IGHD gene usage") +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous(labels=percent) +
    scale_color_manual(values=celltype.colors) +
    geom_point(aes(color=CellType), size=1, alpha=0.8)
  
  family_J <- countGenes(object_b@meta.data[!is.na(object_b$j_call),],
                         gene="j_call", groups="celltype", mode="family")
  names(family_J)[1] <- "CellType"
  g3 <- ggplot(family_J[family_J$gene!="",], aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle("IGHJ gene usage") +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous(labels=percent) +
    scale_color_manual(values=celltype.colors) +
    geom_point(aes(color=CellType), size=1, alpha=0.8)
  g1+g2+g3+plot_layout(nrow=1,guides = "collect")
  
  
  
# figS4A
  plotprop_b <- subset(object_b,is_treatment == "Post-Treat",invert=T)
  # plotprop_b <- subset(plotprop_b,celltype %in% c("IFN-induced memory B","Non-switched memory B"),invert=T)
  cell_percent <- plotprop_b@meta.data
  cell_percent$CellType <- paste0(cell_percent$celltype)
  cell_percent$CellType[cell_percent$celltype %in% c("Pro/Pre-B","Immature B")] <- "Pro/Pre and Immature B"
  cell_percent$CellType[cell_percent$celltype %in% c("Plasma","Plasmablast")] <- "Plasma (blast)"
  
  df <- cell_percent %>% 
    group_by(sample_type,sample,CellType,.drop = FALSE) %>% 
    dplyr::summarise(n=n()) %>% 
    dplyr::mutate(freq = n/sum(n))
  CellType.colors <- c("Pro/Pre and Immature B"="#9467bd",
                       "Naive B"="seagreen3",
                       "CD27+ memory B"="#fe9929",
                       "Atypical memory B"="#c49c94",
                       "Plasma (blast)"="steelblue3")
  df$CellType <- factor(df$CellType,levels=names(CellType.colors))
  
  df_pr <- df %>% filter(CellType == "Pro/Pre and Immature B")
  p1 <- ggplot(df_pr, aes(x = sample_type, y = freq,fill = CellType)) + 
    geom_boxplot(width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
    ylab("Cell Fraction (%)")+ggtitle("Pro/Pre and Immature B")+
    scale_fill_manual(values=CellType.colors)+
    scale_y_continuous(labels = scales::percent)+
    stat_signif(comparisons = list(c("Normal","NewlyDx")),
                textsize=2.5,y_position=c(0.2),map_signif_level=F)+
    theme_classic()+theme(legend.position="none")
    
  df_naive <- df %>% filter(CellType == "Naive B")
  p2 <- ggplot(df_naive, aes(x = sample_type, y = freq,fill = CellType)) + 
    geom_boxplot(width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
    ylab("Cell Fraction (%)")+ggtitle("Naive B")+
    scale_fill_manual(values=CellType.colors)+
    scale_y_continuous(labels = scales::percent)+
    stat_signif(comparisons = list(c("Normal","NewlyDx"),c("NewlyDx","RelRef"),c("Normal","RelRef")),
                textsize=2.5,y_position=c(0.8,0.7,0.9),map_signif_level=F)+
    theme_classic()+theme(legend.position="none")
    
  df_cd27p <- df %>% filter(CellType == "CD27+ memory B")
  p3 <- ggplot(df_cd27p, aes(x = sample_type, y = freq,fill = CellType)) + 
    geom_boxplot(width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
    ylab("Cell Fraction (%)")+ggtitle("CD27+ memory B")+
    scale_fill_manual(values=CellType.colors)+
    scale_y_continuous(labels = scales::percent)+
    stat_signif(comparisons = list(c("Normal","NewlyDx"),c("NewlyDx","RelRef"),c("Normal","RelRef")),
                textsize=2.5,y_position=c(0.6,0.7,0.98),map_signif_level=F)+
    theme_classic()+theme(legend.position="none")
    
  df_atypical <- df %>% filter(CellType == "Atypical memory B")
  p4 <- ggplot(df_atypical, aes(x = sample_type, y = freq,fill = CellType)) + 
    geom_boxplot(width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
    ylab("Cell Fraction (%)")+ggtitle("Atypical memory B")+
    scale_fill_manual(values=CellType.colors)+
    scale_y_continuous(labels = scales::percent)+
    stat_signif(comparisons = list(c("Normal","NewlyDx"),c("NewlyDx","RelRef"),c("Normal","RelRef")),
                textsize=2.5,y_position=c(0.42,0.5,0.6),map_signif_level=F)+
    theme_classic()+theme(legend.position="none")
    
  df_plasma <- df %>% filter(CellType == "Plasma (blast)")
  p5 <- ggplot(df_plasma, aes(x = sample_type, y = freq,fill = CellType)) + 
    geom_boxplot(width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
    ylab("Cell Fraction (%)")+ggtitle("Plasma (blast)")+
    scale_fill_manual(values=CellType.colors)+
    scale_y_continuous(labels = scales::percent)+
    stat_signif(comparisons = list(c("Normal","NewlyDx"),c("NewlyDx","RelRef"),c("Normal","RelRef")),
                textsize=2.5,y_position=c(0.58,0.63,0.69),map_signif_level=F)+
    theme_classic()+theme(legend.position="none")
  p1+p2+p3+p4+p5+plot_layout(nrow=1)
  
  