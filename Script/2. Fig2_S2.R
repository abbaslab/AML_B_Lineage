
# fig2 and S2, B cells of AML patients by transcriptome features and BCR repertoire & Association between B cells and clinical manifestations


# fig2A
  DimPlot(object_b,group.by="celltype",cols=celltype.colors)
  
  p1 <- DimPlot(object_b,cells.highlight = list(IGHM=WhichCells(object_b,expression  = c_call == "IGHM")),
                cols.highlight="#fa9fb5",cols="linen",sizes.highlight=0.05,raster=F)+
    theme(legend.text=element_text(size = 8),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm")) & NoAxes()
  p2 <- DimPlot(object_b,cells.highlight = list(IGHD=WhichCells(object_b,expression  = c_call == "IGHD")),
                cols.highlight="seagreen3",cols="linen",sizes.highlight=0.05,raster=F)+
    theme(legend.text=element_text(size = 8),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm")) & NoAxes()
  p3 <- DimPlot(object_b,cells.highlight = list(IGHG=WhichCells(object_b,expression  = c_call %like% "IGHG")),
                cols.highlight="#2171b5",cols="linen",sizes.highlight=0.05,raster=F) +
    theme(legend.text=element_text(size = 8),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm")) & NoAxes()
  p4 <- DimPlot(object_b,cells.highlight = list(IGHA=WhichCells(object_b,expression  = c_call %like% "IGHA")),
                cols.highlight="#fe9929",cols="linen",sizes.highlight=0.05,raster=F)+
    theme(legend.text=element_text(size = 8),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm")) & NoAxes()
  p1+p2+p3+p4+plot_layout(nrow=1,guides="collect")
    
  
# fig2B
  B_markers <- list(
    `B cell`=c("CD79A","CD79B","MS4A1","CD19","CD37"),
    Plasma=c("JCHAIN","MZB1","SDC1"),
    Proliferation=c("MKI67","TOP2A","TUBB","STMN1","TYMS"),
    Exhausted=c("FCRL5","FCRL3","FCGR2B","CD72"),
    `IFN-indeuced`=c("IFITM1","MX1","IFI44L","IRF7"),
    Naive=c("IL4R","FCER2","IGHD"),
    Immature=c("MME","CD24","CD38"),
    `Pro/Pre`=c("RAG1","RAG2","VPREB1","IGLL1")
  )
  dotplot_b <- subset(object_b,celltype %in% c("CD27+ memory B","Non-switched memory B"),invert=T)
  DotPlot(dotplot_b, features = B_markers,group.by = "celltype",dot.scale=2.5)+
    scale_y_discrete(position="right")+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size=7),
          axis.text.y=element_text(hjust = 1,vjust=0.5,size=7),
          legend.title=element_text(size = 8),
          legend.text=element_text(size = 6),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"),
          panel.border=element_rect(size=0.5),
          strip.text=element_text(size=8),
          #strip.text=element_blank(),
          strip.background=element_blank())+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('chartreuse4',"linen",'chocolate',"chocolate4"))+
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3))


# fig2C
  DimPlot(object_b,label =F,group.by = "MF_group",cols=MF_group.colors,pt.size=0.05,na.value="gray90")+
    labs(col="SHM")+ggtitle("")+NoAxes()
  phenoPropotion_flip(data=object_b@meta.data,x_name="celltype_2",y_name="MF_group",cols=MF_group.colors,legend_name="SHM")


# fig2D
  gdat <- Embeddings(object_b,reduction = 'umap');
  gdat <- as.data.frame(gdat);
  gdat$timepoint <- object_b$sample_type
  ggplot(data = gdat,aes(x = UMAP_1, y = UMAP_2 ))+
    stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)+
    geom_point(color = 'white',size = 0.11)+
    scale_fill_viridis(option="magma")+
    theme_black()+
    facet_wrap(~timepoint,ncol = 1)+NoAxes()
  

# fig2E
  miloplot_b <-subset(object_b,sample_type != "RelRef")
  b_milo <- Milo(as.SingleCellExperiment(miloplot_b))
  b_milo <- buildGraph(b_milo, k = 30, d = 30)
  b_milo <- makeNhoods(b_milo, prop = 0.1, k = 30, d=30, refined = TRUE)
  b_milo <- countCells(b_milo, meta.data = as.data.frame(colData(b_milo)), sample="sample")
  
  #' NewlyDx vs Normal
  b_design <- data.frame(colData(b_milo))[,c("sample", "sample_type")]
  # batch effect is not considered
  b_design <- distinct(b_design)
  rownames(b_design) <- b_design$sample
  b_design$sample_type <- factor(b_design$sample_type,levels=c("Normal","NewlyDx"))
  
  b_milo <- calcNhoodDistance(b_milo, d=30)
  
  da_results <- testNhoods(b_milo, design = ~ sample_type, design.df = b_design)
  head(da_results)
  da_results %>%
    arrange(SpatialFDR) %>%
    head() 

  plot_milo <- annotateNhoods(b_milo, da_results, coldata_col = "celltype")
  head(plot_milo)
  plot_milo$celltype <- ifelse(plot_milo$celltype_fraction < 0.7, "Mixed", plot_milo$celltype)
  plot_milo$celltype <- factor(plot_milo$celltype,levels=c("Mixed",levels(object_b$celltype)))
  
  ggplot(plot_milo,aes(celltype, logFC,color=logFC)) + 
    xlab("") + ylab("Log Fold Change") + 
    geom_jitter(width=0.2,size=0.01)+
    geom_violin(trim=T,colour="black",fill=NA,linetype=1,size=0.3,width=1)+
    #scale_color_gradient2(midpoint=0, low="#1f77b4", mid="gray80",high="#CF3A32")+
    scale_color_gradient2(midpoint=0, low="darkgreen", mid="gray95",high="#A50026")+
    theme_classic()+
    theme(legend.position = "right",
          legend.title=element_text(size = 8),
          legend.text=element_text(size = 6),
          legend.key.height = unit(0.3,"cm"),
          legend.key.width = unit(0.25,"cm"),
          axis.title.y = element_text(size = 9),
          axis.title.x = element_text(size = 9),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 8))+
    coord_flip()
  
  
# fig2F
  phenoPropotion(data=object_b@meta.data,x_name="sample_type",y_name="celltype",cols=celltype.colors,legend_name="cellType")

  
  
# fig2A
  cor.mat <- GetAssayData(object_b,slot = "data")
  cor.mat <- as.matrix(cor.mat)
  cor.mat <- groupMeans(cor.mat,groups = object_b$sample)
  M <-  cor(cor.mat)
  anno_col <- unique(data.frame(object_b$sample,object_b$Cytogenetics,object_b$sample_type))
  row.names(anno_col) <- anno_col$object_b.sample
  anno_col$object_b.sample <- NULL
  names(anno_col) <- c("Cytogenetics","Sample")
  pheatmap(M,clustering_method = "ward.D2",annotation_col=anno_col,treeheight_col = 25,border_color="gray80",
           angle_col="90",fontsize_col=7,annotation_colors=list(Cytogenetics=Cytogenetics.colors,Sample=sample_type.colors))
  
  
# figS2B
  featureplot <- function(gene){
    colour_bk <- c("#f0f0f0",#colorRampPalette(c("#006837","#d9ef8b"))(5),
                   colorRampPalette(c("#d9ef8b","#fee08b"))(5),
                   colorRampPalette(c("#fee08b","#a50026"))(10))
    FeaturePlot(object_b, c(gene),raster=T,order=F,pt.size=2)+
      scale_colour_gradientn(colours = colour_bk)+
      theme_classic()+
      theme(legend.text=element_text(size = 8),
            legend.key.height = unit(0.4,"cm"),
            legend.key.width = unit(0.25,"cm"),
            title=element_text(size = 10))+NoAxes()
  }
  featureplot("IL4R")
  featureplot("TCL1A")
  featureplot("FCER2")
  featureplot("IGHD")
  featureplot("CD27")
  featureplot("CD24")
  
  
# figS2C
  vlnplot_gene <- function(gene){
    VlnPlot(object_b,features=c(gene),group.by="celltype",
            cols=celltype.colors,pt.siz=0)+xlab("")+
      theme_classic() +
      theme(legend.position = "none",
            axis.title.y = element_text(size = 9),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 7,angle=45,vjust=1,hjust=1),
            axis.text.y = element_text(size = 7),
            title=element_text(size = 10))
  }
  vlnplot_gene("FCRL3")
  vlnplot_gene("FCRL5")
  vlnplot_gene("FCGR2B")
  vlnplot_gene("ITGAX")
  

# figS2D
  phenoPropotion_flip(data=object_b@meta.data,x_name="celltype",y_name="c_call",
                      cols=c_call.colors,legend_name="Isotype")
  
# figS2E
  sample_order <- c(paste0("N",1:6),paste0("PT",9:32),sort(str_extract(unique(object_b$sample),pattern="^PT[1-8][A|B|C]")))
  cols_tmp <- ColAssign(sample_order)
  cols_tmp[names(cols_tmp) %in% c("PT7A","PT7B","PT7C")] <- "gray25"  # 改一下颜色！！！！！
  
  DimPlot(object_b,label =F,group.by = "sample",cols=cols_tmp,pt.size=0.05,na.value="gray90")+
    labs(col="Sample")+ggtitle("")+
    theme(legend.title=element_text(size = 8),
          legend.text=element_text(size = 8),
          legend.key.height = unit(0.1,"cm"),
          legend.key.width = unit(0.1,"cm"))+NoAxes()
 
