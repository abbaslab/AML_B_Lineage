
# fig3 and S5, AP-1 activity defines subclusters of Naïve and memory B cells in AML patients


#######################
# Data processing

# 10053 cells
  
  naive_mem <- subset(object_b,subset=celltype %in% c("Naive B","CD27+ memory B","Atypical memory B"))

  
# CytoTRACE
  
  GeneCounts <- as.matrix(naive_mem@assays$RNA@counts)
  iOrd <- rowSums(GeneCounts>0)
  GeneCounts <- GeneCounts[iOrd>10,]#only keep genes expressing in more than 10 cell
  naive_mem_CytoTRACE <- CytoTRACE(GeneCounts, enableFast = TRUE,subsamplesize = 1000)
  naive_mem$CytoTRACE <- naive_mem_CytoTRACE$CytoTRACE[colnames(naive_mem)]

  
# Reclustering of naive_mem
  
  naive_mem <- NormalizeData(naive_mem)  
  # Selection of highly variable genes, 4000 genes, accumulative variance interpretation: 0.737783
  Vars <- apply(naive_mem@assays$RNA@data,1,var)
  Vars <- Vars[order(Vars,decreasing = T)]
  Svars <- sum(Vars)
  df <- data.frame(Ngenes=seq(1,28277,by=50),FracVar=NA)
  for(i in seq_along(df$Ngenes)) {
    df$FracVar[i] <- sum(Vars[1:df$Ngenes[i]])/Svars
  }
  ggplot(df, aes(Ngenes, FracVar)) +   
    geom_point(size = 1, col = "grey")+
    geom_vline(xintercept = 4000)+
    geom_hline(yintercept = df[df$Ngenes==4001,"FracVar"])+theme_classic()
  naive_mem <- FindVariableFeatures(naive_mem,nfeatures=4000)
  # delete TCR,BCR,MT,RP(ribosomal) genes 
  BCR.genes <- c(grep("^IG[KHL]V",rownames(naive_mem),value = T),
                 grep("^IG[KHL]J",rownames(naive_mem),value = T),
                 grep("^IG[KHL]C",rownames(naive_mem),value = T))
  MT.genes <- c(grep("^MT-",rownames(naive_mem),value = T))
  RP.genes <- c(grep("^RP[SL]",rownames(naive_mem),value = T))
  var.genes <- setdiff(VariableFeatures(naive_mem),
                       c(BCR.genes,MT.genes,RP.genes))
  # 3827
  VariableFeatures(object_b) <- var.genes
  naive_mem <- ScaleData(naive_mem,features = var.genes)
  naive_mem <- RunPCA(naive_mem, verbose = FALSE,features = VariableFeatures(object_b),npcs=50)
  # 30
  ElbowPlot(naive_mem,ndims=50)  
  naive_mem <- RunHarmony(naive_mem, dims.use=1:30,group.by.vars="orig.ident")
  naive_mem <- FindNeighbors(naive_mem, dims=1:30,reduction = "harmony",k.param = 30)
  naive_mem <- FindClusters(naive_mem,resolution=0.5,random.seed=46)
  naive_mem <- RunUMAP(naive_mem,reduction = "harmony",dims=1:30,min.dist=0.7,seed.use=1)

    
# functional term (Hallmarks, Hallmarks, KEGG) scores of single cell
  
  # Hallmarks
  Hset <- msigdbr(species = "Homo sapiens",category="H")
  Gset <-  Hset %>% split(x = .$gene_symbol, f = .$gs_name)
  # AP-1
  Gset[[length(Gset)+1]] <- c("FOS","FOSB","FOSL1","FOSL2","JUN","JUNB","JUND")
  names(Gset)[length(Gset)] <- "AP-1 transcription factor"
  # KEGGSets, KEGG pathway gene sets acquired from KEGGREST R package
  load(file.path("Data","KEGGSets.rds"))
  
  # length: 391
  Gset <- c(Gset,KEGGSets)
  
  # 3.037816 mins
  t1 <- Sys.time()
  naive_mem_func <- AddModuleScore(naive_mem,features=Gset,name = names(Gset))
  t2 <- Sys.time()
  t2-t1

  # saveRDS(naive_mem_func,file.path("Data","object_naive_mem.rds"))


    
#######################
# Figures
  
  # naive_mem <- readRDS(file.path("Data","object_naive_mem.rds"))
  
  
# fig3A 
  
  vln_ap1 <- subset(naive_mem,subset=sample_type == "RelRef",invert=T)
  VlnPlot(vln_ap1,features="AP.1.transcription.factor51",
          group.by="sample_type",pt.size=0,cols=sample_type.colors,adjust=1)+
    theme_classic()+ylab("AP-1")+
    stat_signif(comparisons = list(c("Normal","NewlyDx")),
                textsize=2.5,y_position=2.4,map_signif_level=T,tip_length=0)+
    stat_summary(fun.y=mean, geom="point", shape=18,size=1.5, color="red")
 
  colour_bk <- c("#f0f0f0",colorRampPalette(c("#006837","#d9ef8b"))(5),
                 colorRampPalette(c("#d9ef8b","#fee08b"))(5),
                 colorRampPalette(c("#fee08b","#a50026"))(8))
  FeaturePlot(naive_mem, "AP.1.transcription.factor51",order=F,pt.size=0.2)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+
    ggtitle("AP-1")+NoAxes()
  DimPlot(naive_mem,group.by="celltype_2",cols=celltype_2.colors,label=T,pt.size=2,raster=T)+NoAxes()+NoLegend()
  
  gdat <- Embeddings(naive_mem,reduction = 'umap');
  gdat <- as.data.frame(gdat);
  gdat$timepoint <- naive_mem$sample_type
  ggplot(data = gdat,aes(x = UMAP_1, y = UMAP_2 ))+
    stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)+
    geom_point(color = 'white',size = 0.11)+
    scale_fill_viridis(option="magma")+
    theme_black()+
    facet_wrap(~timepoint,ncol = 3)+NoAxes()
  
  
# fig3B
  
  DEGs_naive <- FindMarkers(naive_mem,ident.1="NaiveB AP1_hi",ident.2="NaiveB AP1_lo",group.by="celltype_2",logfc.threshold=0)
  DEGs_memory <- FindMarkers(naive_mem,ident.1="CD27+ MemB AP1_hi",ident.2="CD27+ MemB AP1_lo",group.by="celltype_2",logfc.threshold=0)
  DEGs_aty <- FindMarkers(naive_mem,ident.1="Atypical MemB AP1_hi",ident.2="Atypical MemB AP1_lo",group.by="celltype_2",logfc.threshold=0)
  DEGs_AMLvNormal <- FindMarkers(naive_mem,ident.1="NewlyDx",ident.2="Normal",group.by="sample_type",logfc.threshold=0)
  
  colnames(DEGs_naive) <- paste0(colnames(DEGs_naive),"-","Naive")
  colnames(DEGs_memory) <- paste0(colnames(DEGs_memory),"-","Memory")
  colnames(DEGs_aty) <- paste0(colnames(DEGs_aty),"-","Atypical")
  colnames(DEGs_AMLvNormal) <- paste0(colnames(DEGs_AMLvNormal),"-","AMLvNormal")

  iOrd <- Reduce("intersect",list(rownames(DEGs_naive),rownames(DEGs_memory),rownames(DEGs_aty),rownames(DEGs_AMLvNormal)))
  df <- Reduce("cbind",list(DEGs_naive[iOrd,],DEGs_memory[iOrd,],DEGs_aty[iOrd,],DEGs_AMLvNormal[iOrd,]))
  df$gene <- rownames(df)
  # write.xlsx(df,file.path("Data/Tables","TableS3_Diff_Genes.xlsx"))
  
  #######################
  # scatter plot
  df$label <- rownames(df)
  # down 
  down_label <- c("TXNIP","MS4A1","LTB","CD52","RCSD1","HVCN1","SELL","CD79B")
  # up
  up <- df$label[(df$`avg_log2FC-Naive` > 0.3 & df$`avg_log2FC-Memory` > 0.3)]
  Ap1 <- c("FOS","FOSB","FOSL1","FOSL2","JUN","JUNB","JUND")
  Hset <- msigdbr(species = "Homo sapiens",category="H")
  Gset <-  Hset %>% split(x = .$gene_symbol, f = .$gs_name)
  NFkb <- Gset$HALLMARK_TNFA_SIGNALING_VIA_NFKB
  up_label <- c(intersect(up,Ap1),intersect(up,NFkb))
  df$label[!(df$label %in% c(down_label,up_label))] <- NA
  df$gene <- row.names(df)
  
  ggplot(data=df, aes(x = `avg_log2FC-Naive`, y = `avg_log2FC-Memory`, color=`avg_log2FC-AMLvNormal`,label=label)) +
    geom_point(stroke = 0, alpha = 1,shape = 16,size=1) + 
    theme_classic() +
    geom_text_repel(color="black",size=1.2,segment.size=0.1,point.padding=NA,max.overlaps=30) +
    scale_color_gradient2(midpoint=0, low="steelblue4", mid="#f0f0f0",high="#a50026",space ="Lab") +
    geom_hline(yintercept=0, linetype="dashed", color = "indianred")+
    geom_vline(xintercept=0, linetype="dashed", color = "indianred")+
    ylab("Log2FC of NaiveB AP1_hi vs. NaiveB AP1_lo")+
    xlab("Log2FC of CD27+ MemB AP1_hi vs. CD27+ MemB AP1_lo")+
    labs(color = "Log2FC of NewlyDx vs. Normal")
  
  
  #######################
  # heatmap

  # AML vs Normal 
  DEGs_AMLvNormal_diff <- DEGs_AMLvNormal[abs(DEGs_AMLvNormal$`avg_log2FC-AMLvNormal`) > 0.3,]
  DEGs_AMLvNormal_diff$gene <- rownames(DEGs_AMLvNormal_diff)
  DEGs_AMLvNormal_diff <- DEGs_AMLvNormal_diff %>% arrange(`avg_log2FC-AMLvNormal`)
  DEGs_AMLvNormal_up <- DEGs_AMLvNormal_diff$gene[DEGs_AMLvNormal_diff$`avg_log2FC-AMLvNormal` > 0.3]
  DEGs_AMLvNormal_down <- DEGs_AMLvNormal_diff$gene[DEGs_AMLvNormal_diff$`avg_log2FC-AMLvNormal` < -0.3]
  DEGs_AMLvNormal_hp <- as.matrix(t(data.frame(DEGs_AMLvNormal_diff$`avg_log2FC-AMLvNormal`)))
  row.names(DEGs_NAMLvNormal_hp) <- "Log2FC of NewlyDx vs. Normal"
  names(DEGs_AMLvNormal_hp) <- DEGs_AMLvNormal_diff$gene
  
  colour_bk <- c(#"#f0f0f0",
    colorRampPalette(c("steelblue4","#d9ef8b"))(length(DEGs_NAMLvNormal_down)-1),
    "gray90",
    colorRampPalette(c("#d9ef8b","#a50026"))(length(DEGs_NAMLvNormal_up)))
  
  show_gene <- unique(c(intersect(DEGs_NAMLvNormal_up,Ap1),intersect(DEGs_NAMLvNormal_up,NFkb),
                        intersect(DEGs_NAMLvNormal_down,down_label)))
  show_gene <- setdiff(show_gene,c("BTG2","GADD45B","LITAF","MCL1","PLEK","PPP1R15A",
                                   "RHOB","SGK1","SLC2A3","SQSTM1","ZBTB10","ZFP36"))
  col_anno <-  columnAnnotation(
    show=anno_mark(at= which(names(DEGs_NAMLvNormal_hp) %in% show_gene),
                   labels = names(DEGs_NAMLvNormal_hp)[which(names(DEGs_NAMLvNormal_hp) %in% show_gene)],
                   link_width = unit(1, "mm"),link_height=unit(5, "mm"),labels_gp = gpar(fontsize = 5)))
  Heatmap(DEGs_AMLvNormal_hp,col = colour_bk,border=T,
          show_row_names=F,cluster_columns = F,show_heatmap_legend=F,
          top_annotation =col_anno,
          width = unit(7, "cm"), height = unit(0.8, "cm"))
  # range(DEGs_NAMLvNormal_hp), -1.386976 ~ 1.915153
  
  
  #Atypical MemB AP1_hi vs. Atypical MemB AP1_lo
  DEGs_aty_diff <- DEGs_aty
  DEGs_aty_diff <- DEGs_aty_diff[abs(DEGs_aty_diff$`avg_log2FC-Atypical`) > 0.3,]
  DEGs_aty_diff$gene <- rownames(DEGs_aty_diff)
  DEGs_aty_diff <- DEGs_aty_diff %>% arrange(`avg_log2FC-Atypical`)
  DEGs_aty_up <- DEGs_aty_diff$gene[DEGs_aty_diff$`avg_log2FC-Atypical` > 0.3]
  DEGs_aty_down <- DEGs_aty_diff$gene[DEGs_aty_diff$`avg_log2FC-Atypical` < -0.3]
  DEGs_aty_hp <- as.matrix(t(data.frame(DEGs_aty_diff$`avg_log2FC-Atypical`)))
  row.names(DEGs_aty_hp) <- "Log2FC of Atypical"
  names(DEGs_aty_hp) <- DEGs_aty_diff$gene
  
  colour_bk <- c(#"#f0f0f0",
    colorRampPalette(c("steelblue4","#d9ef8b"))(length(DEGs_aty_down)-1),
    "gray90",
    colorRampPalette(c("#d9ef8b","#a50026"))(length(DEGs_aty_up)))
  
  show_gene <- unique(c(intersect(DEGs_aty_down,down_label),
                        "JUNB","JUN","JUND","FOS","FOSB",
                        "CD69","CD83","DUSP1","DUSP2",
                        "KLF2","KLF4","KLF6",
                        "NFKB1","NFKB2","NFKBIA",
                        "NR4A1","NR4A2","NR4A3",
                        "TNF","TNFSF9"))
  col_anno <-  columnAnnotation(
    show=anno_mark(side="right",at= which(names(DEGs_aty_hp) %in% show_gene),
                   labels = names(DEGs_aty_hp)[which(names(DEGs_aty_hp) %in% show_gene)],
                   link_width = unit(1, "mm"),link_height=unit(5, "mm"),labels_gp = gpar(fontsize = 5)))
  Heatmap(DEGs_aty_hp,col = colour_bk,border=T,
          show_row_names=F,cluster_columns = F,show_heatmap_legend=F,
          bottom_annotation =col_anno,
          width = unit(7, "cm"), height = unit(0.8, "cm"))
  # range(DEGs_aty_hp), -3.908202 ~ 1.731153
  
  
# fig3C
  
  naive_mem_noRR <- subset(naive_mem,subset=sample_type == "RelRef",invert=T)
  ap1 <- naive_mem_noRR@meta.data[,77]
  other <- naive_mem_noRR@meta.data[,c(27:76,78:417)]
  
  cor_ap1_other <- lapply(1:ncol(other),function(i){
    cor_data <- cor.test(ap1,other[,i],method="spearman")
    data.table(description=names(other)[i],cor=cor_data$estimate,p=cor_data$p.value)
  } )
  cor_ap1_other <- list.rbind(cor_ap1_other)
  
  setorder(cor_ap1_other,"cor")
  cor_ap1_other[,rank:=1:nrow(cor_ap1_other)]
  cor_ap1_other1 <- cor_ap1_other
  
  description <- c("HALLMARK_P53_PATHWAY37","HALLMARK_APOPTOSIS7","HALLMARK_HYPOXIA22",
                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB45","HALLMARK_INFLAMMATORY_RESPONSE25","HALLMARK_UNFOLDED_PROTEIN_RESPONSE46",
                   "MAPK.signaling.pathway...Homo.sapiens..human.229")
  term <- c("P53 pathway","Apoptosis","Hypoxia",
            "TNFA signaling via NFKB","Inflammatory Response","Unfolded Protein Response",
            "MAPK signaling pathway")
  term_name <- data.table(description,term)
  

  cor_ap1_other1 <- left_join(cor_ap1_other1,term_name)
  cor_ap1_other1[cor > 0 & p < 0.05,label:="Pos_Cor"] 
  
  ggplot(data=cor_ap1_other1,aes(x=rank,y=cor,color=label,label=term))+
    scale_color_manual(values=c(Pos_Cor=muted("red")))+
    geom_point(alpha = 0.8,shape = 16,size=0.2)+
    theme_classic() +
    geom_label_repel(size=2.5,label.padding=0.08,min.segment.length=0.1)+
    ylab("Correlation coefficients")+
    theme(legend.position = "none",
          axis.title.y = element_text(size = 9),
          axis.title.x = element_text(size = 9),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7))

    
# fig3D
  
  vln_desc <- c("HALLMARK_APOPTOSIS7","HALLMARK_TNFA_SIGNALING_VIA_NFKB45","HALLMARK_INFLAMMATORY_RESPONSE25",
                "HALLMARK_UNFOLDED_PROTEIN_RESPONSE46","MAPK.signaling.pathway...Homo.sapiens..human.229")
  VlnPlot(naive_mem_noRR,features = vln_desc,cols=sample_type.colors,
          group.by="sample_type",fill.by="ident",
          stack = T,flip=T,adjust=1.2)+
    ylab("")+xlab("Activity Scores")+ggtitle("This cohort")+
    theme(title=element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 6))+
    NoLegend()+
    stat_signif(comparisons = list(c("Normal","NewlyDx")),
                textsize=2.5,y_position=0.1,map_signif_level=T,tip_length=0)
  
  
# fig3E
  
  mem <- subset(naive_mem,subset=celltype == "Naive B",invert=T)
  mem$celltype_3 <- paste0(mem$celltype_2)
  mem$celltype_3[mem$celltype_2 %like% "AP1_lo"] <- "MemB AP1_lo"
  mem$celltype_3[mem$celltype_2 %like% "AP1_hi"] <- "MemB AP1_hi"
  mem$celltype_3 <- factor(mem$celltype_3,levels=c("MemB AP1_lo","MemB AP1_hi"))
  celltype_3.colors <- c(`MemB AP1_lo`="gray50", `MemB AP1_hi`=muted("red"))

  cloneDiversity(clone_data=mem@meta.data,col_name="celltype_3",colors=celltype_3.colors,legend_name="CellType")
  
# fig3F
  
  mem_clone <- na.omit(data.table(clone_id=mem$clone_id,AP1=mem$AP.1.transcription.factor51))
  mem_clone$clone_type <- ifelse(duplicated2(mem_clone$clone_id),"Expanded","Unique")
  mem_clone$clone_type <- factor(mem_clone$clone_type,levels=c("Unique","Expanded"))
  ggplot(mem_clone,aes(x=clone_type,y=AP1,fill=clone_type))+
    geom_boxplot(width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    stat_compare_means()+
    scale_fill_manual(values=c("#c49c94","#2ca02c"))+
    theme_classic() +
    ylab("AP-1")+
    theme(legend.position="none",
          axis.title.y = element_text(size = 9),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8,angle=45,hjust=1),
          axis.text.y = element_text(size = 8))

  
  
  
  
# figS5A
  
  DimPlot(naive_mem,cells.highlight = list(PT18=WhichCells(naive_mem,expression = sample=="PT18")),
          cols.highlight=muted("red"),sizes.highlight=3,raster=T)+NoAxes()


# figS5B
  
  phenoPropotion_flip(data=naive_mem@meta.data,x_name="celltype_2",y_name="is_treatment",cols=is_treatment.colors,legend_name="Sample")
  






