

# Packages (Versions)
  
  # 4.1.1
  library(Seurat)
  # 0.1.0
  library(harmony)
  # 1.4.0
  library(miloR)
  # 1.10.0
  library(SingleR)
  
  # 1.14.2
  library(data.table)
  # 0.4.6.2
  library(rlist)
  # 7.5.1
  library(msigdbr)
  # 1.0.10
  library(dplyr)
  # 1.3.2
  library(tidyverse)
  # 4.4.4
  library(clusterProfiler)
  # 1.22.0
  library(fgsea)
  # 0.4.9
  library(survminer)
  # 1.2.0
  library(scales)
  # 0.7.8
  library(questionr)
  # 4.7-0
  library(Hmisc)
  # 4.2.5
  library(openxlsx)
  # 1.44.0
  library(phenoTest)
  
  # 3.4.0
  library(ggplot2)
  # 0.6.3
  library(ggsignif)
  # 0.5.2
  library(ggpmisc)
  # 4.2.4
  library(ggthemes)
  # 0.4.0
  library(ggpubr)
  # 0.9.1
  library(ggrepel)
  # 0.5.3
  library(ggridges)
  # 0.2.5
  library(ggpie)
  # 0.6.2
  library(viridis)
  # 1.0.12
  library(pheatmap)
  # 2.12.0
  library(ComplexHeatmap)
  # 1.6.1
  library(simplifyEnrichment)
  # 1.1.1
  library(patchwork)
  # 1.1
  library(smoother)
  # 0.12.3
  library(ggalluvial)
  
  # 1.2.0
  library(alakazam)
  # 1.1.1
  library(shazam)
  # 0.1.0
  library(Startrac)
  # 1.7.0
  library(scRepertoire)
  
  # 0.3.3
  library(CytoTRACE)
  # 2.24.1
  library(monocle)
  # 1.2.9
  library(monocle3)
  # 1.0.3
  library(SCENT)



# Functions
  
  fucs <- list.files("Script/functions",full.names=T,pattern=".R$")
  # 注：有的程序如果读不进来 用notepad改成utf-8编码
  lapply(fucs,function(x) {print(x);source(x,encoding="utf-8");return("Yes")})
  


# Data
  
  #############################################
  # scBCR data, patient_obs: 22317 cells
  load(file="Data/patient_obs.rda")
  
  
  
  #############################################
  # scRNA object combined with scBCR data, obejct_b: 12423 cells
  
  # object_b_count
  load(file.path("Data","object_b_count.rda"))
  # object_b_metadata
  load(file.path("Data","object_b_metadata.rda"))
  
  
  object_b <- CreateSeuratObject(counts=object_b_count,meta.data=object_b_metadata)
  object_b <- NormalizeData(object_b)  
  
  # Find highly variable genes, Selection of highly variable genes, 4000 genes, accumulative variance interpretation: 0.7473641
  Vars <- apply(object_b@assays$RNA@data,1,var)
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
  object_b <- FindVariableFeatures(object_b,nfeatures=4000)
  # delete TCR,BCR,MT,RP(ribosomal) genes 
  TCR.genes <- grep("^TR[ABGD][VJ]",rownames(object_b),value = T)
  BCR.genes <- c(grep("^IG[KHL]V",rownames(object_b),value = T),
                 grep("^IG[KHL]J",rownames(object_b),value = T),
                 grep("^IG[KHL]C",rownames(object_b),value = T))
  MT.genes <- c(grep("^MT-",rownames(object_b),value = T))
  RP.genes <- c(grep("^RP[SL]",rownames(object_b),value = T))
  
  var.genes <- setdiff(VariableFeatures(object_b),
                       c(TCR.genes,BCR.genes,MT.genes,RP.genes))
  # 3769
  VariableFeatures(object_b) <- var.genes
  
  object_b <- ScaleData(object_b,features = var.genes)
  object_b <- RunPCA(object_b, verbose = FALSE,features = VariableFeatures(object_b),npcs=50)
  ElbowPlot(object_b,ndims=50)  
  object_b$orig.ident <- as.character(object_b$orig.ident)
  object_b <- RunHarmony(object_b, dims.use=1:30,group.by.vars="orig.ident")
  object_b <- FindNeighbors(object_b, dims=1:30,reduction = "harmony",k.param = 30)
  object_b <- FindClusters(object_b,resolution=1.5,random.seed=1)
  object_b <- RunUMAP(object_b,reduction = "harmony",dims=1:30,min.dist=0.5,seed.use=1)
  
  # To maintain consistency with the UMAP plot in the article, it is recommended to replace the UMAP embedding.
  umap_embeddings <- as.matrix(object_b_metadata[,c("UMAP_1","UMAP_2")])
  colnames(umap_embeddings) <- c("umap_1","umap_2")
  object_b@reductions$umap@cell.embeddings <- umap_embeddings
  
  
  
  #############################################
  # colors
  sample_labels.colors <- c(Normal="#1f77b4",AML="#aec7e8")
  sample_type.colors <- c(Normal="#2ca02c",NewlyDx="#1f77b4",RelRef="#CF3A32")
  is_treatment.colors <- c(Normal="#1f77b4",NewlyDx="#aec7e8",`Prior-Treat`="#98df8a",`Post-Treat`="#2ca02c")
  celltype.colors <- c(`Pro/Pre-B`="#9467bd",`Immature B`="#c5b0d5",`Naive B`="seagreen3",
                       `CD27+ memory B`="#fe9929",`IFN-induced memory B`="#d62728",`Non-switched memory B`="#EC6564",
                       `Atypical memory B`="#c49c94",
                       Plasma="#aec7e8",Plasmablast="#1f77b4")
  celltype_2.colors <- c(
    `Pro/Pre-B`="#9467bd",`Immature B`="#c5b0d5",`NaiveB AP1_lo`="#98df8a",`NaiveB AP1_hi`="#2ca02c",
    `CD27+ MemB AP1_lo`="#ffbb78",`CD27+ MemB AP1_hi`="#ff7f0e",`IFN-induced MemB`="#d62728",`Non-switched MemB`="#EC6564",
    `Atypical MemB AP1_lo`="#c49c94",`Atypical MemB AP1_hi`="#8c564b",
    Plasma="#aec7e8",`Plasmablast`="#1f77b4")
  c_call.colors <- c(
    IGHM="#fa9fb5",IGHD="seagreen3",
    IGHG1="#c6dbef",IGHG2="#6baed6",IGHG3="#2171b5",IGHG4="steelblue4",
    IGHA1="#fec44f",IGHA2="#fe9929")
  Isotype.colors <- c(IGHM="#fa9fb5",IGHD="seagreen3",IGHG="#2171b5",IGHA="#fe9929")
  MF_group.colors <- c(Germline="#74add1",Low="#ffffbf",Medium="#fdae61",High="#d73027")
  Cytogenetics.colors <- c(diploid="linen",del5="#6BC463",del7="#FC8115",double="#CF3A32",other="#1F77B4")
  Best_response.colors <- c(CR="#2ca02c",PR="#1f77b4",SD="#ff7f0e",NR="#d62728")
  
  
  
  
  
  
  
  




