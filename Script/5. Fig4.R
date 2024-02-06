
# fig4, Pseudotime analysis reveals potential transitional dynamics among B-cell lineage subgroups

#######################
# Data processing
  
  # generated from "4. Fig3_S5.R"
  # naive_mem <- readRDS(file.path("Data","object_naive_mem.rds"))
  
 
  # Monocle3 trajectory analysis
  DEG1 <- FindMarkers(naive_mem,ident.1=c("Atypical MemB AP1_lo","Atypical MemB AP1_hi"),group.by = "celltype_2")   #!!!!!!!!!!!!!!!!!!!
  DEG1$gene <- rownames(DEG1)
  DEG2 <- FindMarkers(naive_mem,ident.1="CD27+ memory B",ident.2="Naive B",group.by = "celltype")
  DEG2$gene <- rownames(DEG2)
  DEGs <- unique(c(DEG1$gene,DEG2$gene))
  
  naive_mem@active.assay <- "RNA"
  cell_metadata <- naive_mem@meta.data
  expression_data <- GetAssayData(naive_mem,slot="counts")
  gene_annotation <- data.frame(gene_short_name=rownames(expression_data),
                                stringsAsFactors = F,row.names = rownames(expression_data))
  naive_mem_cds <- new_cell_data_set(expression_data=expression_data,
                                        cell_metadata = cell_metadata,
                                        gene_metadata = gene_annotation)
  naive_mem_cds <- preprocess_cds(naive_mem_cds, num_dim =50,use_genes = DEGs,method="PCA")
  naive_mem_cds <- reduce_dimension(naive_mem_cds,reduction_method="UMAP",preprocess_method="PCA",
                                       umap.min_dist = 0.45,umap.n_neighbors = 28)
  naive_mem_cds <- cluster_cells(naive_mem_cds)
  naive_mem_cds <- learn_graph(naive_mem_cds,use_partition=T)
  naive_mem_cds <- order_cells(naive_mem_cds)
  
  # saveRDS(naive_mem_cds,file.path("Data","object_naive_mem_cds.rds"))
  
  
  
#######################
# Figures
  
# fig4A
  
  # naive_mem_cds <- readRDS(file.path("Data","object_naive_mem_cds.rds"))
  
  plot_cells(naive_mem_cds, color_cells_by="celltype_2", 
             reduction_method="UMAP",show_trajectory_graph=T,
             trajectory_graph_color = "grey40",trajectory_graph_segment_size = 0.75,
             label_roots = F,label_leaves = F,label_cell_groups =F,label_branch_points = F,
             alpha=1,cell_size = 0.42,group_label_size = 10)+
    xlab("")+ylab("")+ggtitle("")+del_axis()+
    scale_color_manual(values = celltype_2.colors)
          
  plot_cells(naive_mem_cds, color_cells_by="pseudotime", 
             reduction_method="UMAP",show_trajectory_graph=T,
             trajectory_graph_color = "grey40",trajectory_graph_segment_size = 0.75,
             label_roots = F,label_leaves = F,label_cell_groups =F,label_branch_points = F,
             alpha=1,cell_size = 0.45,group_label_size = 10)+
    xlab("")+ylab("")+ggtitle("")+del_axis()

  plot_cells(naive_mem_cds, color_cells_by="sample_type", 
             reduction_method="UMAP",show_trajectory_graph=T,
             trajectory_graph_color = "grey40",trajectory_graph_segment_size = 0.75,
             label_roots = F,label_leaves = F,label_cell_groups =F,label_branch_points = F,
             alpha=1,cell_size = 0.42,group_label_size = 10)+
    xlab("")+ylab("")+ggtitle("")+del_axis()+
    scale_color_manual(values = sample_type.colors)

  plot_cells(naive_mem_cds, color_cells_by="MF_group", 
             reduction_method="UMAP",show_trajectory_graph=T,
             trajectory_graph_color = "grey40",trajectory_graph_segment_size = 0.75,
             label_roots = F,label_leaves = F,label_cell_groups =F,label_branch_points = F,
             alpha=1,cell_size = 0.42,group_label_size = 10)+
    xlab("")+ylab("")+ggtitle("")+del_axis()+
    scale_color_manual(values = MF_group.colors)

  monocle3_fuction(naive_mem_cds,func="AP.1.transcription.factor51")  

  
# fig4B
  
  # naive_mem <- readRDS(file.path("Data","object_naive_mem.rds"))
  mem <- subset(naive_mem,subset=celltype %in% c("CD27+ memory B","Atypical memory B"))
  
  p1 <- ggplot(mem@meta.data,aes(x=celltype,y=CytoTRACE))+
    geom_boxplot(aes(fill=celltype),width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    scale_fill_manual(values=celltype.colors)+
    stat_signif(comparisons = list(c("CD27+ memory B","Atypical memory B")),
                textsize=2.5,y_position=1,map_signif_level=T)+
    theme_classic() +
    theme(legend.position="none",
          axis.title.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(size = 7))
  
  mem$celltype_2 <- factor(mem$celltype_2,
                           levels=c("CD27+ MemB AP1_lo","CD27+ MemB AP1_hi","Atypical MemB AP1_lo","Atypical MemB AP1_hi"))
  p2 <- ggplot(mem@meta.data,aes(x=celltype_2,y=CytoTRACE))+
    geom_boxplot(aes(fill=celltype_2),width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    scale_fill_manual(values=celltype_2.colors)+
    stat_signif(comparisons = list(c("CD27+ MemB AP1_lo","CD27+ MemB AP1_hi"),
                                   c("Atypical MemB AP1_lo","Atypical MemB AP1_hi")),
                textsize=2.5,y_position=c(1,1),map_signif_level=T)+
    theme_classic() +
    theme(legend.position="none",
          axis.title.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(size = 7))
  
  p1+p2
  
  
# fig4C
  
  naive_mem_trans <- naive_mem
  umap_monocle3 <- naive_mem_cds@int_colData$reducedDims$UMAP
  colnames(umap_monocle3) <- c("UMAP_1","UMAP_2")
  naive_mem_trans@reductions$umap@cell.embeddings <- umap_monocle3
  
  ###########
  # SingleR
  cd27mem.query <- subset(naive_mem_trans,subset=celltype == "CD27+ memory B") 
  naive.reference <- subset(object_b,subset=celltype %in% c("Naive B"))   
  #7.802634 mins
  pred <- SingleR(ref=naive.reference@assays$RNA@data,
                  test=cd27mem.query@assays$RNA@data,
                  labels=naive.reference$celltype_2)
  
  pred_singleR <- data.frame(cell_id=row.names(pred),SingleR=pred$labels)
  naive_mem_trans@meta.data <- left_join(naive_mem_trans@meta.data,pred_singleR)
  row.names(naive_mem_trans@meta.data) <- naive_mem_trans$cell_id
  
  DimPlot(naive_mem_trans,group.by="SingleR",
          cols=c("#98df8a","#2ca02c"),na.value="grey95",
          pt.size=0.02)+
    ggtitle("")+
    NoAxes()
  
  p1 <- ggpie(data=naive_mem_trans@meta.data[naive_mem_trans$celltype_2 == "CD27+ MemB AP1_lo",],
              group_key="SingleR",count_type = "full",
              fill_color=c("#98df8a","#2ca02c"),
              border_size=0.2,border_color = "gray40",
              label_info = "all", label_type = "horizon",
              label_size = 2, label_pos = "in",label_split = NULL)+
    ggtitle("CD27+ MemB AP1_lo")+
    theme(title = element_text(size=8),
          legend.text = element_text(size = 5),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"))
  p2 <- ggpie(data=naive_mem_trans@meta.data[naive_mem_trans$celltype_2 == "CD27+ MemB AP1_hi",],
              group_key="SingleR",count_type = "full",
              fill_color=c("#98df8a","#2ca02c"),
              border_size=0.2,border_color = "gray40",
              label_info = "all", label_type = "horizon",
              label_size = 2, label_pos = "in",label_split = NULL)+
    ggtitle("CD27+ MemB AP1_hi")+
    theme(title = element_text(size=8),
          legend.text = element_text(size = 5),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"))
  p1+p2+plot_layout(guides="collect")
 
  
  ###########
  # TransferData
  cd27mem.query <- subset(object_b,subset=celltype == "CD27+ memory B") 
  naive.reference <- subset(object_b,subset=celltype %in% c("Naive B")) 
  
  cd27mem.query <- NormalizeData(cd27mem.query)
  cd27mem.query <- FindVariableFeatures(cd27mem.query)
  cd27mem.query <- ScaleData(cd27mem.query)
  
  naive.reference <- NormalizeData(naive.reference)
  naive.reference <- FindVariableFeatures(naive.reference)
  naive.reference <- ScaleData(naive.reference)
  
  # find anchors
  anchors <- FindTransferAnchors(reference = naive.reference, query = cd27mem.query)
  
  # transfer labels
  predictions <- TransferData(anchorset = anchors, refdata = naive.reference$celltype_2)
  cd27mem.query <- AddMetaData(object = cd27mem.query, metadata = predictions)
  
  
  phenoPropotion_flip(data=cd27mem.query@meta.data,cols=c(`NaiveB AP1_lo`="#98df8a",`NaiveB AP1_hi`="#2ca02c"),
                      x_name="celltype_2",y_name="predicted.id",legend_name="")
  
  
# fig4D
  
  st_cell <- c("NaiveB AP1_lo","NaiveB AP1_hi","CD27+ MemB AP1_lo","CD27+ MemB AP1_hi",
                  "Atypical MemB AP1_lo","Atypical MemB AP1_hi","Plasma","Plasmablast")
  st_b <- subset(object_b,subset=celltype_2 %in% st_cell)
  st_b$celltype_2 <- paste0(st_b$celltype_2) 
  in.dat <- data.frame(Cell_Name=st_b$cell_id, clone.id=st_b$clone_id, patient=st_b$sample,
                       majorCluster=st_b$celltype_2,loc=st_b$is_treatment)
  out <- Startrac.run(in.dat, proj="AML",verbose=F)
  
  
  #####################
  # 各细胞亚型总体值
  Startrac::plot(out,index.type="pairwise.tran",byPatient=T)
  # startrac包这个热图色条不稳定，如果把Atypical分成两类，就只有这两类有颜色，其他的都看不出来
  hp <- out@pIndex.tran[1:length(unique(st_b$celltype_2)),3:ncol(out@pIndex.tran)]
  row.names(hp) <- out@pIndex.tran[1:length(unique(st_b$celltype_2)),"majorCluster"]
  names(hp) <- names(out@pIndex.tran)[3:ncol(out@pIndex.tran)]
  hp <- hp[st_cell,st_cell]
  # 两类Atypical之间数值过高
  hp[hp > 0.03] <- 0.03
  colour_bk <- c("gray70",#colorRampPalette(c("#006837","#d9ef8b"))(40),
                 colorRampPalette(c("#006837","#d9ef8b"))(120),
                 colorRampPalette(c("#d9ef8b","#fee08b"))(30),
                 colorRampPalette(c("#fee08b","#a50026"))(300))
  pheatmap::pheatmap(hp,cluster_rows=F,cluster_cols=F,angle_col=45,fontsize = 10,
           color=colour_bk,na_col="white",
           width=5,height=4)
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  