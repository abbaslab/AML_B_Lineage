
# fig1 and S1 - Single-cell BCR repertoire in the bone marrow of AML patients

  # 22317 cells
  patient_obs  
  # 20560 cells
  patient_obs <- patient_obs[!patient_obs$sample  %like% "PT[1-8][B|C]",]

  
# fig1B
  sample_order <- c(paste0("N",1:6),paste0("PT",9:32),paste0("PT",1:8,"A"))
  clonalProportion(clone_data=patient_obs,split = c(3,50,100,200),
                 group="sample",group_order=sample_order)

  
# fig1C
  cloneDiversity(clone_data=patient_obs,col_name="sample_type",legend_name="Sample",colors=sample_type.colors)
  

# fig1D
  cloneMutation(clone_data=patient_obs,col_name="sample_type",colors=sample_type.colors)
  phenoPropotion(data=patient_obs,x_name="sample_type",y_name="MF_group",cols=MF_group.colors,legend_name="SHM")
  

# fig1E
  count_clone <- data.table(countClones(patient_obs,group=c("sample","sample_type","Isotype"))) 
  isotype_clone <- data.table(table(count_clone$sample,count_clone$Isotype))
  names(isotype_clone) <- c("sample","Isotype","num")
  isotype_clone <- left_join(isotype_clone,data.table(table(count_clone$sample)),by=c("sample"="V1"))
  isotype_clone[,percent:=num/N]
  isotype_clone <- left_join(isotype_clone,unique(patient_obs[,c("sample","sample_type")]))
  isotype_clone <- isotype_clone[Isotype %like% "IGH[A|G]"]
  isotype_clone[,Isotype:=factor(Isotype,levels=c("IGHG","IGHA"))]
  ggplot(isotype_clone,aes(x=Isotype,y=percent,fill=factor(sample_type,levels=names(sample_type.colors))))+
    geom_boxplot(color='gray40',outlier.shape=NA,size=0.1)+
    geom_point(position = position_jitterdodge(),color='black',size=0.9,shape=21)+
    labs(fill="Sample")+ylab("Fraction of B cells (%)")+
    scale_fill_manual(values=sample_type.colors)+
    scale_y_continuous(labels = scales::percent,limits=c(0,0.8))+
    theme_classic() +
    theme(legend.title=element_text(size = 8),
          legend.text=element_text(size = 6),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"),
          axis.title.y = element_text(size = 9),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8))
  # 0.03619 
  wilcox.test(isotype_clone[Isotype=="IGHG" & sample_type=="RelRef",percent],
              isotype_clone[Isotype=="IGHG" & sample_type %in% c("NewlyDx","Normal"),percent])
  # 0.003923
  wilcox.test(isotype_clone[Isotype=="IGHA" & sample_type=="RelRef",percent],
              isotype_clone[Isotype=="IGHA" & sample_type %in% c("NewlyDx","Normal"),percent])
  
  phenoPropotion(data=patient_obs,x_name="sample_type",y_name="Isotype",cols=Isotype.colors,legend_name="Isotype")
  

  
# figS1B
  gene <- countGenes(patient_obs, gene="v_call", groups="sample_type", mode="gene")
  ighv <- gene %>%
    mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name")))
  ggplot(ighv, aes(x=gene, y=seq_freq)) +
    theme_bw() +
    theme(legend.title=element_text(size = 8),
          legend.text=element_text(size = 6),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 7,angle =45, vjust = 1, hjust=1),
          # axis.text.x = element_text(size = 7,angle =90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 6),
          panel.border=element_rect(size=1))+
    ylab("Percent of repertoire") +xlab("") +labs(color="Sample")+
    scale_y_continuous(labels=percent) +
    geom_point(aes(color=sample_type), size=1,alpha=1)+
    scale_color_manual(values= sample_type.colors)
 
  
# figS1C
  family <- countGenes(patient_obs, gene="v_call", groups=c("sample_type", "Isotype"),
                       clone="clone_id", mode="family") %>% filter(!is.na(Isotype))
  ggplot(family, aes(x=gene, y=clone_freq)) +
    theme_bw() +
    #ggtitle("Clonal Usage") +
    theme(legend.title=element_text(size = 8),
          legend.text=element_text(size = 6),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          #axis.text.x = element_text(size = 7,angle =30, vjust = 1, hjust=0.75),
          axis.text.x = element_text(size = 7,angle =90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 6),
          panel.border=element_rect(size=1))+
    ylab("Percent of repertoire") + xlab("") +
    scale_y_continuous(labels=percent) +
    scale_color_manual(values= Isotype.colors)+
    geom_point(aes(color=Isotype), size=1.5, alpha=1) +
    facet_grid(. ~ sample_type)
  
  
# figS1D
  cloneDiversity(clone_data=patient_obs,col_name="Isotype",legend_name="Isotype",colors=Isotype.colors)

  
# figS1E
  phenoPropotion(data=patient_obs,x_name="c_call",y_name="MF_group",cols=MF_group.colors,legend_name="SHM")
  

# figS1F
  cloneMutation(clone_data=patient_obs,col_name="c_call",colors=c_call.colors)
  

