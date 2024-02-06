
# fig6 and S7, Further differentiation and stimulated status of plasma cells in AML microenvironment

  # 1509 cells
  Plasma <- subset(object_b,subset=celltype_2 %in% c("Plasma","Plasmablast"))

  
  # Cytotrace
  GeneCounts <- as.matrix(Plasma@assays$RNA@counts)
  iOrd <- rowSums(GeneCounts>0)
  GeneCounts <- GeneCounts[iOrd>10,]#only keep genes expressing in more than 10 cell
  Plasma_CytoTRACE <- CytoTRACE(GeneCounts, enableFast = TRUE,subsamplesize = 1000)
  Plasma$CytoTRACE <- Plasma_CytoTRACE$CytoTRACE[colnames(Plasma)]

  
  # Monocle2
  cell_metadata <- Plasma@meta.data
  expression_data <- GetAssayData(Plasma,slot="counts")
  expression_data <- expression_data[rowSums(expression_data>0)>10,]#only keep genes expressing in more than 10 cell
  gene_annotation <- data.frame(gene_short_name=rownames(expression_data),
                                stringsAsFactors = F,row.names = rownames(expression_data))
  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  plasma_cds <- newCellDataSet(cellData=expression_data,phenoData = pd,
                               featureData = fd,expressionFamily=negbinomial.size())
  plasma_cds <- estimateSizeFactors(plasma_cds)
  plasma_cds <- estimateDispersions(plasma_cds)
  plasma_cds <- detectGenes(plasma_cds, min_expr = 0.1)
  
  expressed_genes <- row.names(subset(fData(plasma_cds),num_cells_expressed >= 300))
  # 55.23255 secs
  diff_test_res <- differentialGeneTest(plasma_cds[expressed_genes,],fullModelFormulaStr = "~celltype_2")
  ordering_genes <- unique(row.names(diff_test_res)[order(diff_test_res$qval)][1:1000])
  
  plasma_cds <- setOrderingFilter(plasma_cds, ordering_genes)
  plasma_cds <- reduceDimension(plasma_cds, max_components = 2,reduction_method = 'DDRTree',pseudo_expr=1)
  plasma_cds <- orderCells(plasma_cds,reverse=F)

  
  # Plasma UMAP embeddings 
  monocle_dims <- t(plasma_cds@reducedDimS)
  colnames(monocle_dims) <- c("UMAP_1","UMAP_2")
  Plasma@reductions$umap@cell.embeddings <- monocle_dims

  
  # plasma
  plasma <- subset(Plasma,subset=celltype=="Plasma")
  
  
# fig6A

  plot_cell_trajectory(plasma_cds, color_by = "celltype_2",cell_size=0.3,show_branch_points =F,show_tree =F)+
    scale_color_manual(values=celltype_2.colors[c("Plasmablast","Plasma")])
  
  plot_cell_trajectory(plasma_cds, color_by = "nFeature_RNA",cell_size=0.3,show_branch_points =F,show_tree =F)+
    scale_colour_gradient2(low="seagreen",mid="linen",high="chocolate4",midpoint = 3000)
  df <- pData(plasma_cds)
  df <- df[order(df$Pseudotime),]
  ggplot(df, aes(x=Pseudotime,y=nFeature_RNA))+geom_smooth(fullrange=T,method="loess",span=1,color="steelblue",size=0.5)+
    xlim(0,27)+ theme_classic()
    
  plot_cell_trajectory(plasma_cds, color_by = "CytoTRACE",cell_size=0.1,show_branch_points =F,show_tree =F)+
    scale_colour_gradient2(low="seagreen",mid="linen",high="chocolate4",midpoint = 0.5)
  ggplot(df, aes(x=Pseudotime,y=CytoTRACE))+geom_smooth(fullrange=T,method="loess",span=1,color="steelblue",size=0.5)+
    xlim(0,25)+ theme_classic()
    

# fig6B
  
  df <- plasma@meta.data[plasma$Isotype %in% c("IGHA","IGHG"),]
  df$Isotype <- factor(df$Isotype,levels = c("IGHA","IGHG"))
  ggplot(df,aes(x=Isotype,y=CytoTRACE))+
    geom_boxplot(aes(fill=Isotype),width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    scale_fill_manual(values=Isotype.colors)+
    stat_signif(comparisons = list(c("IGHA","IGHG")),textsize=2.5,y_position=c(0.95),map_signif_level=T)+
    theme_classic()
  

# fig6C
  
  plot_cell_trajectory(plasma_cds, color_by = "c_call",cell_size=0.3,show_branch_points =F,show_tree =F)+
    scale_color_manual(values=c_call.colors)

  df <- pData(plasma_cds)
  # IGHG
  igg_dense <- density(df$Pseudotime[df$Isotype == "IGHG"],adjust=1,na.rm=T)
  #  plot(igg_dense)
  igg_func <- approxfun(igg_dense)
  # points(df$Pseudotime[df$Isotype == "IGHG"],igg_func(df$Pseudotime[df$Isotype == "IGHG"]))
  
  # IGHA
  iga_dense <- density(df$Pseudotime[df$Isotype == "IGHA"],adjust=1,na.rm=T)
  # plot(iga_dense)
  iga_func <- approxfun(iga_dense)
  # points(df$Pseudotime[df$Isotype == "IGHA"],iga_func(df$Pseudotime[df$Isotype == "IGHA"]))
  
  # IGHG vs. IGHA
  igg_iga_ratio <- c()
  for(i in 1:nrow(df)){
    igg_iga_ratio[length(igg_iga_ratio)+1] <- 
      ifelse(df[i,"Isotype"] %in% c("IGHG","IGHA"),
             igg_func(df[i,"Pseudotime"])/iga_func(df[i,"Pseudotime"]),NA)
  }
  
  Plasma@meta.data <- left_join(Plasma@meta.data,
                                data.frame(cell_id=plasma_cds$cell_id,
                                           Pseudotime=plasma_cds$Pseudotime,
                                           IGG_IGA_enrich=igg_iga_ratio))
  rownames(Plasma@meta.data) <- Plasma$cell_id
  FeaturePlot(Plasma,features="IGG_IGA_enrich",pt.size=0.05,order=T)+
    scale_colour_gradient2(low="#006837",mid="#f0f0f0",high="#a50026",midpoint=1,na.value="gray",limits=c(0,2))+
    ggtitle("")
  
  
# fig6D
  
  VlnPlot(plasma,features=c("IGHG1","IGHG2","IGHG3","IGHG4"),
          cols=is_treatment.colors,fill.by="ident",group.by="is_treatment",
          same.y.lims=T,y.max=8,adjust=0.8,stack=T,flip=T)+
    theme(legend.position = "none") + ggtitle("")

  
# fig6E
  
  df <- pData(plasma_cds)
  # AML
  aml_dense <- density(df$Pseudotime[df$sample_labels == "AML"],adjust=1)
  # plot(aml_dense)
  aml_func <- approxfun(aml_dense)
  # points(df$Pseudotime[df$sample_labels == "AML"],aml_func(df$Pseudotime[df$sample_labels == "AML"]))
  
  # Normal
  normal_dense <- density(df$Pseudotime[df$sample_labels == "Normal"],adjust=1)
  # plot(normal_dense)
  normal_func <- approxfun(normal_dense)
  # points(df$Pseudotime[df$sample_labels == "Normal"],normal_func(df$Pseudotime[df$sample_labels == "Normal"]))
  
  # AML vs. Normal
  aml_normal_ratio <- c()
  for(t in df$Pseudotime){
    aml_normal_ratio[length(aml_normal_ratio)+1] <- 
      aml_func(t)/normal_func(t)
  }
  
  Plasma@meta.data <- left_join(Plasma@meta.data,
                                data.frame(cell_id=plasma_cds$cell_id,
                                           Pseudotime=plasma_cds$Pseudotime,
                                           AML_Normal_enrich=aml_normal_ratio))
  rownames(Plasma@meta.data) <- Plasma$cell_id
  FeaturePlot(Plasma,features="AML_Normal_enrich",pt.size=0.1)+
    scale_colour_gradient2(low="#006837",mid="#f0f0f0",high="#a50026",midpoint=1,na.value="#a50026",limits=c(0,2))+
    ggtitle("")
  
  
  df <- pData(plasma_cds)
  df$sample_labels <- factor(df$sample_labels,levels=c("Normal","AML"))
  ggplot(df, aes(x=Pseudotime,y=sample_labels,fill=sample_labels))+
    geom_density_ridges(scale=4) +
    scale_fill_manual(values=sample_labels.colors)+
    #geom_vline(xintercept = c(34.5),linetype=2)+
    theme_ridges()
  

# fig6F
  
  ggplot(plasma@meta.data[plasma$sample_type != "RelRef",],aes(x=sample_type ,y=CytoTRACE))+
    geom_boxplot(aes(fill=sample_type),width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    scale_fill_manual(values=is_treatment.colors)+
    ggtitle("This chohort, Plasma")+
    stat_signif(comparisons = list(c("Normal","NewlyDx")),textsize=2.5,y_position=c(0.95),map_signif_level=T)+
    theme_classic()
  
  
# fig6G
  
  df <- data.frame(t(plasma_cds@reducedDimS))
  df$Pseudotime <- plasma_cds$Pseudotime
  df <- df[order(df$Pseudotime),]
  df$group <- as.numeric(cut2(df$Pseudotime, g=100))
  
  cor.mat <- GetAssayData(Plasma,slot = "data")[,rownames(df)]
  cor.mat <- as.matrix(cor.mat)
  cor.mat <- groupMeans(cor.mat,groups = df$group)
  
  
  # genes variable across pseudotime:
  {
    iOrd <- apply(cor.mat, 1, sd)
    cor.mat <- cor.mat[order(iOrd,decreasing = T),]
    result.x <- cor.mat[1:3000,]
    gene.cor <- data.frame(genes=rownames(result.x),cors=cor(1:100,t(result.x),use="na.or.complete")[1,])
    gene.cor <- gene.cor[order(gene.cor$cors),]
    gene.cor$BCRgenes <- NA
    gene.cor$BCRgenes[grep("^IG",gene.cor$genes)] <- "BCRgenes"
  }
  
  result.x.smth <- result.x[gene.cor$genes,]
  for (i in 1:nrow(result.x)) {
    result.x.smth[i,] <- smth(result.x.smth[i,],window = 0.1,method = "gaussian")
  }
  result.x.smth[,1:5] <- result.x.smth[,6]
  result.x.smth[,96:100] <- result.x.smth[,95]
  bk <- c(seq(-1,-0.1,by=0.02),seq(0,1,by=0.02))
  colour_bk <- c(colorRampPalette(c("steelblue","#d1e5f0"))(38),
                 colorRampPalette(c("#d1e5f0","#f7f7f7"))(10),
                 colorRampPalette(c("#f7f7f7","#fddbc7"))(10),
                 colorRampPalette(c("#fddbc7","darkred"))(39))
  labels_row <- row.names(result.x.smth)
  labels_row[!labels_row %in% c("CD19","CD38","SDC1")] <- ""
  
  pheatmap::pheatmap(result.x.smth, cluster_rows=F, cluster_cols=F,
           scale="row",show_colnames=F,#show_rownames=F,
           breaks = bk,color = colour_bk,angle_col=45,
           clustering_method="ward.D",border_color=NA,
           labels_row=labels_row,fontsize=6,
           annotation_row=gene.cor[,"BCRgenes",drop=F])
 
  
  bp1 <- enrichGO(rownames(result.x.smth)[1:1000],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  p1 <- enrichplot::dotplot(bp1,showCategory=15)
  bp2 <- enrichGO(rownames(result.x.smth)[1001:2000],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  p2 <- enrichplot::dotplot(bp2,showCategory=15)
  bp3 <- enrichGO(rownames(result.x.smth)[2001:3000],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  p3 <- enrichplot::dotplot(bp3,showCategory=15)
  # p1+p2+p3+plot_layout(nrow=1)
  
  write.xlsx(list(BP_top1000=bp1@result[bp1@result$p.adjust < 0.05,],
                  BP_1000_2000=bp2@result[bp2@result$p.adjust < 0.05,],
                  BP_2000_3000=bp3@result[bp3@result$p.adjust < 0.05,]),
             file=file.path("Data/Tables","TableS6_Plasma_top3000BP.xlsx"))
  
  
# fig6H
  
  colour_bk <- c("#f0f0f0",#colorRampPalette(c("#006837","#d9ef8b"))(1),
                 colorRampPalette(c("#d9ef8b","#fee08b"))(5),
                 colorRampPalette(c("#fee08b","#a50026"))(10))
  FeaturePlot(Plasma,features="CD19",pt.size=0.1)+
    scale_colour_gradientn(colours = colour_bk)
  FeaturePlot(Plasma,features="CD38",pt.size=0.1,order=T)+
    scale_colour_gradientn(colours = colour_bk)
  FeaturePlot(Plasma,features="SDC1",pt.size=0.1,order=T)+
    scale_colour_gradientn(colours = colour_bk)
  
  
# fig6I
  
  # SingleR
  plasma.query <- Plasma
  mem.reference <- subset(object_b,subset=celltype %in% c("CD27+ memory B","Atypical memory B"))    
  # 4.182298 mins
  pred <- SingleR(ref=mem.reference@assays$RNA@data,
                  test=plasma.query@assays$RNA@data,
                  labels=mem.reference$celltype)
  pred_singleR <- data.frame(cell_id=row.names(pred),SingleR=pred$labels)
  Plasma@meta.data <- left_join(Plasma@meta.data,pred_singleR)
  row.names(Plasma@meta.data) <- Plasma$cell_id
  
  DimPlot(Plasma,group.by="SingleR",
          cols=c("#c49c94","#fe9929"),na.value="grey95",
          pt.size=0.04)+
    ggtitle("")+
    NoAxes()
  
  
  ggpie(data=Plasma@meta.data,
        group_key="SingleR",count_type = "full",
        fill_color=c("#fe9929","#c49c94"),
        border_size=0.2,border_color = "gray40",
        label_info = "all", label_type = "horizon",
        label_size = 2, label_pos = "in",label_split = NULL)
  
  
# figS7A
  
  ggplot(Plasma@meta.data[!is.na(Plasma$c_call),],aes(x=c_call,y=CytoTRACE))+
    geom_boxplot(aes(fill=c_call),width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    scale_fill_manual(values=c_call.colors)+
    theme_classic()


# figS7C
  
  gseBP <- msigdbr(species = "Homo sapiens",
                   category = "C5", 
                   subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
  gseBP$gs_name <- gsub('GOBP_','',gseBP$gs_name)
  # gseBP$gs_name <- tolower(gseBP$gs_name)
  gseBP$gs_name <- gsub('_',' ',gseBP$gs_name)
  gseBP <- gseBP %>% split(x = .$gene_symbol, f = .$gs_name)
  names(gseBP) <- str_to_title(names(gseBP))
  
  diff <- FindMarkers(plasma,ident.1="NewlyDx",ident.2="Normal",group.by="is_treatment",logfc.threshold =0)
  diff$gene <- row.names(diff)
  diff <- arrange(diff,desc(avg_log2FC))
  geneList <- diff$avg_log2FC
  names(geneList) <- diff$gene
  
  gseaRes <- fgsea(pathways = gseBP, 
                   stats = geneList,
                   minSize=5,
                   maxSize=500)
  Up <- gseaRes[NES > 0 & padj < 0.1][head(order(desc(NES)), n=20), pathway]
  gseaSets_plot <- gseBP[Up]
  plotGseaTable(gseaSets_plot, geneList, gseaRes,gseaParam = 0.5)

  
  write.xlsx(gseaRes[padj < 0.1][order(desc(NES))],file.path("Data/Tables","TableS7_Plasma_NewlyDx_vs_Normal_GSEA.xlsx"))
 
  
  

  
