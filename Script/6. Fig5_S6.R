
# fig5 and figS6, Atypical memory B cells were specialized for antigen presentation and interacted tightly with AML malignant cells

  naive_mem <- subset(object_b,subset=celltype %in% c("Naive B","CD27+ memory B","IFN-induced memory B","Non-switched memory B","Atypical memory B"))
  
  Hset <- msigdbr(species = "Homo sapiens",category="H")
  Gset <-  Hset %>% split(x = .$gene_symbol, f = .$gs_name)
  load(file.path("Data","KEGGSets.rds"))
  # length: 391
  Gset <- c(Gset,KEGGSets)

  naive_mem <- AddModuleScore(naive_mem,features=Gset,name = names(Gset))

  
# fig5A
  k="B.cell.receptor.signaling.pathway...Homo.sapiens..human.86"
  mem_func <- subset(naive_mem,subset=celltype %in% c("CD27+ memory B","Atypical memory B"))
  VlnPlot(mem_func,features=k,group.by="celltype",pt.size=0,cols=celltype.colors,adjust=1.5)+
    ylab("BCR signaling pathway")+
    theme_classic()+
    theme(legend.position="none",
          title=element_blank(),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(size = 7))+
    stat_summary(fun.y=mean, geom="point", shape=18,size=1.5, color="red")+
    stat_compare_means(size=2,label.x=0.85)

  
# fig5B
  k="Antigen.processing.and.presentation...Homo.sapiens..human.72"
  VlnPlot(naive_mem,features=k,group.by="celltype",pt.size=0,cols=celltype.colors,adjust=1.5)+
    ylab("Antigen processing and presentation")+
    theme_classic()+
    theme(legend.position="none",
          title=element_blank(),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(size = 7))+
    stat_summary(fun.y=mean, geom="point", shape=18,size=1.5, color="red")+
    stat_compare_means(size=2,label.x=1)


# figS6A
  k="Endocytosis...Homo.sapiens..human.140"
  y_intercept <- mean(naive_mem$Endocytosis...Homo.sapiens..human.140[naive_mem$celltype == "Atypical memory B"])
  VlnPlot(naive_mem,features=k,group.by="celltype",pt.size=0,cols=celltype.colors,adjust=1.5)+
    geom_hline(yintercept = y_intercept,linetype="dashed",color="gray50")+
    ylab("Endocytosis")+
    theme_classic()+
    theme(legend.position="none",
          title=element_blank(),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(size = 7))+
    stat_summary(fun.y=mean, geom="point", shape=18,size=1.5, color="red")+
    stat_compare_means(size=2,label.x=1)

  
# figS6A
  VlnPlot(naive_mem,features=c("HLA-A","HLA-B","HLA-C","HLA-F","HLA-E",
                         "HLA-DMA","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5",
                         "B2M","CD74"),
          cols=celltype_2.colors,fill.by="ident",group.by="celltype_2",
          same.y.lims=T,adjust=1,stack=T,flip=T)+
    theme(legend.position = "none") + ggtitle("")+
    theme(title=element_text(size = 8),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 6))
  
# fig5C
  
  # 63184 cells
  object_NewlyDxMyeloid <- readRDS(file.path("Data","object_NewlyDxMyeloid.rds"))
  # 6327 cells
  object_naive_mem <- subset(object_b,subset=sample_type == "NewlyDx" & celltype %in% c("Naive B","CD27+ memory B","Atypical memory B"))
  
  LR <- read.xlsx(file.path("Data/Tables","TableS5_Ligand_receptor pairs.xlsx"))
  
  #1. Ligand or receptor express in at least *% of a B cell subtype:
  {
    cutoff <- 0.1 # percentage of cells expressing genes
    B.genes <- unique(c(LR$Ligand.ApprovedSymbol,LR$Receptor.ApprovedSymbol))
    B.genes <- B.genes[B.genes%in%rownames(object_naive_mem)]
    mat <- GetAssayData(object_naive_mem)[B.genes,]
    mat <- mat>0
    Mat <- groupMeans(mat,groups=object_naive_mem$celltype_2,na.rm = TRUE, sparse = T)
    Mat <- Mat[,c("NaiveB AP1_lo","NaiveB AP1_hi","CD27+ MemB AP1_lo","CD27+ MemB AP1_hi",
                  "Atypical MemB AP1_lo","Atypical MemB AP1_hi")]
    iOrd <- rowMaxs(Mat)
    Mat <- Mat[iOrd>cutoff,]
    iOrd <- apply(Mat,1,which.max)
    # only keep those showing highest expression in Atypical mem
    Mat <- Mat[iOrd %in% c(5,6),]
    B.genes <- data.frame(genes=rownames(Mat),Mat)
    
    #update LR:
    iOrd1 <- LR$Ligand.ApprovedSymbol %in% B.genes$genes
    iOrd2 <- LR$Receptor.ApprovedSymbol %in% B.genes$genes
    LR <- LR[iOrd1 | iOrd2,]
    LR$I.genes1 <- ifelse(LR$Ligand.ApprovedSymbol %in% B.genes$genes,
                          LR$Receptor.ApprovedSymbol,
                          NA) # ligand in B cells
    LR$I.genes2 <- ifelse(LR$Receptor.ApprovedSymbol %in% B.genes$genes,
                          LR$Ligand.ApprovedSymbol,
                          NA) # receptor in B cells
  }
  
  #2. B interacting genes express in at least *% of a myeloid cell subtype:
  {
    #genes left for further test:
    AML.genes <- unique(c(LR$I.genes1,LR$I.genes2))
    AML.genes <- AML.genes[AML.genes%in%rownames(object_NewlyDxMyeloid)]
    mat <- GetAssayData(object_NewlyDxMyeloid)[AML.genes,]
    mat <- mat>0
    Mat <- groupMeans(mat,groups=object_NewlyDxMyeloid$celltype_2 ,na.rm = TRUE, sparse = T)
    Mat <- Mat[,c("Progenitor","Monocytes","cDC","AML")]
    iOrd <- rowMaxs(Mat)
    Mat <- Mat[iOrd>cutoff,]
    iOrd <- apply(Mat,1,which.max)
    Mat <- Mat[iOrd %in% c(4),] 
    AML.genes <- data.frame(genes=rownames(Mat),Mat)
  }
  
  #3. LR pairs with both expressed at *% in >=1 B/Myeloid cell subtype:
  {
    iOrd1 <- (LR$I.genes1 %in% AML.genes$genes)
    iOrd2 <- LR$I.genes2 %in% AML.genes$genes
    LR <- LR[iOrd1 | iOrd2,]
  }
  
  #4. heatmap between receptor and ligand:
  {
    BLigand.LR <- LR[!is.na(LR$I.genes1),]
    BReceptor.LR <- LR[!is.na(LR$I.genes2),]
    bk <- c(seq(-2,-0.1,by=0.02),seq(0,2,by=0.02))
    colour_bk <- c(colorRampPalette(c("#1f77b4","#d1e5f0"))(83),
                   colorRampPalette(c("#d1e5f0","#f7f7f7"))(15),
                   colorRampPalette(c("#f7f7f7","#fddbc7"))(15),
                   colorRampPalette(c("#fddbc7",muted("red")))(84))
    
    #heatmaps showing ligand in B cells and receptor in AML cells
    mat <- GetAssayData(object_naive_mem)[BLigand.LR$Ligand.ApprovedSymbol,]
    mat1 <- groupMeans(mat,groups=object_naive_mem$celltype_2,na.rm = TRUE, sparse = T)
    mat1 <- mat1[,c("NaiveB AP1_lo","NaiveB AP1_hi","CD27+ MemB AP1_lo","CD27+ MemB AP1_hi",
                    "Atypical MemB AP1_lo","Atypical MemB AP1_hi")]
    
    mat <- GetAssayData(object_NewlyDxMyeloid)[BLigand.LR$Receptor.ApprovedSymbol,]
    mat2 <- groupMeans(mat,groups=object_NewlyDxMyeloid$celltype_2,na.rm = TRUE, sparse = T)
    mat2 <- mat2[,c("Progenitor","Monocytes","cDC","AML")]
    
    #heatmaps showing receptor in B cells and ligand in AML cells
    mat <- GetAssayData(object_naive_mem)[BReceptor.LR$Receptor.ApprovedSymbol,]
    mat3 <- groupMeans(mat,groups=object_naive_mem$celltype_2,na.rm = TRUE, sparse = T)
    mat3 <- mat3[,c("NaiveB AP1_lo","NaiveB AP1_hi","CD27+ MemB AP1_lo","CD27+ MemB AP1_hi",
                    "Atypical MemB AP1_lo","Atypical MemB AP1_hi")]
    
    mat <- GetAssayData(object_NewlyDxMyeloid)[BReceptor.LR$Ligand.ApprovedSymbol,]
    mat4 <- groupMeans(mat,groups=object_NewlyDxMyeloid$celltype_2,na.rm = TRUE, sparse = T)
    mat4 <- mat4[,c("Progenitor","Monocytes","cDC","AML")]
    
    del1 <- which(row.names(mat2) %in% c("ITGAV","NCSTN","ITGB2"))
    del2 <- which(row.names(mat4) %in% c("RGMB","SPON2"))
    
    pheatmap::pheatmap(mat1[-del1,], cluster_rows=F, cluster_cols=F,
             scale="row",show_colnames=T,show_rownames=T,
             cellheight=9,cellwidth=14,fontsize=7,border_color = "grey80",
             breaks = bk,color = colour_bk,angle_col=90)
    pheatmap::pheatmap(mat2[-del1,], cluster_rows=F, cluster_cols=F,
             scale="row",
             show_colnames=T,show_rownames=T,
             cellheight=9,cellwidth=14,fontsize=7,border_color = "grey80",
             breaks = bk,color = colour_bk,angle_col=90)
    pheatmap::pheatmap(mat3[-del2,], cluster_rows=F, cluster_cols=F,
             scale="row",show_colnames=T,show_rownames=T,
             cellheight=9,cellwidth=14,fontsize=7,border_color = "grey80",
             breaks = bk,color = colour_bk,angle_col=90)
    pheatmap::pheatmap(mat4[-del2,], cluster_rows=F, cluster_cols=F,
             scale="row",
             show_colnames=T,show_rownames=T,
             cellheight=9,cellwidth=14,fontsize=7,border_color = "grey80",
             breaks = bk,color = colour_bk,angle_col=90)
  }
  
  
# figS6F
  B_ligands <- row.names(mat1[-del1,])
  AML_receptors <- row.names(mat2[-del1,])
  B_receptors <- row.names(mat3[-del2,])
  AML_ligands <- row.names(mat4[-del2,])
  
  object_naive_mem$celltype_2 <- factor(object_naive_mem$celltype_2,
                                        levels=rev(c("NaiveB AP1_lo","NaiveB AP1_hi","CD27+ MemB AP1_lo","CD27+ MemB AP1_hi",
                                                     "Atypical MemB AP1_lo","Atypical MemB AP1_hi")))
  p1 <- DotPlot(object_naive_mem, features = unique(B_ligands),group.by = "celltype_2",dot.scale=4)+
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
    #scale_color_gradientn(values = seq(0,1,0.2),colours = c('chartreuse4',"linen",'chocolate',"chocolate4"))+
    scale_colour_gradient2(low = "chartreuse4", mid = "linen", high = "chocolate", guide = "chocolate4",midpoint = 0)+
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3))
  
  p2 <- DotPlot(object_naive_mem, features = unique(B_receptors),group.by = "celltype_2",dot.scale=4)+
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
    #scale_color_gradientn(values = seq(0,1,0.2),colours = c('chartreuse4',"linen",'chocolate',"chocolate4"))+
    scale_colour_gradient2(low = "chartreuse4", mid = "linen", high = "chocolate", guide = "chocolate4",midpoint = 0)+
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3))
  
  p3 <- DotPlot(object_NewlyDxMyeloid, features = unique(AML_ligands),group.by = "celltype_2",dot.scale=2.5)+
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
    #scale_color_gradientn(values = seq(0,1,0.2),colours = c('chartreuse4',"linen",'chocolate',"chocolate4"))+
    scale_colour_gradient2(low = "chartreuse4", mid = "linen", high = "chocolate", guide = "chocolate4",midpoint = 0)+
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3))
  
  p4 <- DotPlot(object_NewlyDxMyeloid, features = unique(AML_receptors),group.by = "celltype_2",dot.scale=2.5)+
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
    #scale_color_gradientn(values = seq(0,1,0.2),colours = c('chartreuse4',"linen",'chocolate',"chocolate4"))+
    scale_colour_gradient2(low = "chartreuse4", mid = "linen", high = "chocolate", guide = "chocolate4",midpoint = 0)+
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3))
  
  p1+p2+p4+p3+plot_layout(byrow=T,ncol=2,guides="collect",width=c(7,3),height=c(6,4))
  
  
# fig6D
  gseBP <- msigdbr(species = "Homo sapiens",
                   category = "C5", 
                   subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
  gseBP$gs_name <- gsub('GOBP_','',gseBP$gs_name)
  gseBP$gs_name <- gsub('_',' ',gseBP$gs_name)
  gseBP$gs_name <- str_to_title(gseBP$gs_name)
  
  set.seed(123)
  ego <-enricher(row.names(mat2),TERM2GENE=gseBP)
  ego3 <- dplyr::mutate(ego, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  BP_show <- c("Cell Substrate Adhesion","Integrin Mediated Signaling Pathway",
               "Cellular Response To External Stimulus","Regulation Of Bmp Signaling Pathway",
               "Regeneration","Receptor Signaling Pathway Via Stat",
               "Positive Regulation Of Mapk Cascade","Leukocyte Adhesion To Vascular Endothelial Cell")
  ego3 <- ego3[ego3$Description %in% BP_show,] 
  ggplot(ego3,aes(GeneRatio, fct_reorder(Description, GeneRatio))) +
    geom_point(aes(color=p.adjust, size = Count)) +
    geom_segment(aes(xend=0, yend = Description),color="linen",size=0.8) +
    geom_point(aes(color=p.adjust, size = Count)) +
    scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                          trans = "log10",
                          guide=guide_colorbar(reverse=TRUE, order=1))+
    scale_size_continuous(range=c(1,3)) +
    theme_minimal()+
    ggtitle("Biological process")
  
  
# fig5E-G
  
  # lr genes
  B_ligands <- c("ADAM15","ADAM17","ADAM28","APP",
                 "CALM1","CALM1","CALM1","CALM2","CALM2","CALM3","CALM3",
                 "CALR","FLT3LG","GNAI2","GNAI2",
                 "LRPAP1","PSEN1","RGMB","TGFB1","TNFSF13","CD48")
  AML_receptors <- c("ITGA9","NOTCH1","ITGA4","SLC45A3","KCNQ5","PTPRA","SELL",
                     "KCNQ5","SELL","KCNQ5","SELL","SCARF1",
                     "FLT3","F2R","IGF1R","SORL1","NOTCH1","BMPR2","ENG","TNFRSF1A","CD244")
  lr_pairs <- data.table(B_ligands,AML_receptors)
  
  
  # lsc17
  lsc17 <- c("DNMT3B","ZBTB46","NYNRIN","ARHGAP22","LAPTM4B","MMRN1","DPYSL3","KIAA0125","CDK6",
             "CPXM1","SOCS2","SMIM24","EMP1","NGFRAP1","CD34","AKR1C3","GPR56")
  lsc17_coef <- c(0.0874,-0.0347,0.00865,-0.0138,0.00582,0.0258,0.0284,0.0196,
                  -0.0704,-0.0258,0.0271,-0.0226,0.0146,0.0465,0.0338,-0.0402,0.0501)
  names(lsc17_coef) <- lsc17
  
  
  # TCGA
  load(file.path("Data","tcga.rda"))
  expr_tcga <- log2(tcga$tcga_tpm+1)
  expr_tcga <- as.data.frame(t(scale(t(expr_tcga))))
  pheno_tcga <- tcga$tcga_pdata
  pheno_tcga[,Sample:=str_replace_all(Sample,"-",".")]
  
  LR_mean_sample <- c()
  for(s in names(expr_tcga)){
    LR_mean <- c()
    for(lr in 1:nrow(lr_pairs)){
      
      lr_mean <-ifelse(lr_pairs[lr,B_ligands] %in% row.names(expr_tcga) & lr_pairs[lr,AML_receptors] %in% row.names(expr_tcga),
                       mean(expr_tcga[lr_pairs[lr,B_ligands],s],expr_tcga[lr_pairs[lr,AML_receptors],s],na.rm=T),NA)
      LR_mean <- c(LR_mean,lr_mean)
    }
    LR_mean_sample  <- c(LR_mean_sample,mean(LR_mean,na.rm=T))
  }
  
  lsc_tcga <- apply(expr_tcga[intersect(row.names(expr_tcga),lsc17),],2,function(x) x %*% lsc17_coef[intersect(row.names(expr_tcga),lsc17)])
  
  signatures_tcga <- data.table(Sample=names(expr_tcga),
                                LR_mean=LR_mean_sample,
                                LSC17=lsc_tcga)
  signatures_tcga[,LR_mean:=scale(LR_mean)]
  pheno_tcga <- left_join(pheno_tcga,signatures_tcga)
  
  # HR (+1SD) = 1.3 [1.081-1.563]     Cox P = 0.00535
  # L: 57, M: 61, H: 55
  summary(coxph(Surv(OS_days,OS) ~ LR_mean, data = pheno_tcga))
  smoothCoxph(pheno_tcga$OS_days,pheno_tcga$OS,pheno_tcga$LR_mean,xlab="Expr. ligand-receptor (SD)")

  
  pheno_tcga$Group <- "M"
  pheno_tcga$Group[pheno_tcga$LR_mean > 0.5] <- "H"
  pheno_tcga$Group[pheno_tcga$LR_mean <= -0.5] <- "L"
  pheno_tcga$Group <- factor(pheno_tcga$Group,levels=c("L","M","H"))
  fit=survfit(Surv(OS_days,OS) ~ Group, data = pheno_tcga)
  p2 <- ggsurvplot(fit, data = pheno_tcga, risk.table = T, conf.int = F,pval=T,
                   pval.size=2.5,size=0.3,censor.size = 3,
                   surv.scale="percent", xscale=365,break.x.by=2*365,xlim=c(0,365*8),
                   palette=c("#2ca02c","gray","#d62728"),legend=c(0.7,0.8),xlab="Survival Time (years)")
  p2 <- p2$plot+
    theme(title=element_text(size=9),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7))
  
  # H vs L, p = 0.013
  pheno_tcga2 <- pheno_tcga[pheno_tcga$Group != "M",]
  surv_pvalue(surv_fit(Surv(OS_days,OS) ~ Group, data = pheno_tcga2))
  
  
  p3 <- ggplot(pheno_tcga2,aes(x=Group,y=LSC17,fill=Group))+
    geom_boxplot(color='gray40',outlier.shape=NA,size=0.1,width=0.6)+
    geom_point(position = position_jitterdodge(),color='gray50',size=0.8,shape=21,stroke=0.2)+
    stat_compare_means(comparisons = list(c("L","H")),method="wilcox.test",size=2.5)+
    scale_fill_manual(values=c("#2ca02c","#d62728"))+
    theme_classic() +
    ylab("LSC17 Signature score")+
    theme(legend.position="none",
          axis.title.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 8))
  p2+p3+plot_layout(widths = c(1.5,1))
  
  
  
  
  
  
  
 