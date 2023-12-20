## AATD scRNAseq
library(Seurat)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(dplyr)
library(presto)
library(msigdbr)
library(fgsea)
library(dorothea)
library(pheatmap)
library(tidyr)
library(viper)

AATD <- readRDS(file = "AATD.Rds")

##Figure 1
#1A UMAP labeled by donor disease state
Idents(AATD) <-("annot2")
levels(AATD) <- c('Control', 'COPD', 'AATD')
DimPlot(AATD, reduction = "umap", cols=c("darkslategray4","cornflowerblue","brown2")) +ggtitle("Donor Disease State")  +theme(plot.title = element_text(hjust = 0.5))

#1B UMAP labeled by cell identity
Idents(AATD) <-"annot3"
levels(AATD)<-c('AT1','AT2','RASC','Secretory','Basal','Ciliated','PNEC','AlvMacrophage','Monocyte','Mast','T&NK','B','Endothelial','Lymphatic','SmoothMuscle','FibroblastsMyofibroblasts','NA')
DimPlot(subset(AATD,idents =c('AT1','AT2','RASC','Secretory','Basal','Ciliated','PNEC','AlvMacrophage','Monocyte','Mast','T&NK','B','Endothelial','Lymphatic','SmoothMuscle','FibroblastsMyofibroblasts')))+
  ggtitle("Cell Identity")  +theme(plot.title = element_text(hjust = 0.5))

#1C Dot plot of pulmonary epithelial marker gene expression
epithelial_markers <- c("AGER","CAV1",
                        "SFTPC","SFTPA1","SFTPA2",
                        "MGP","SOX4","SCGB3A2", #RASC
                        "SCGB1A1","SCGB3A1","LCN2","TSPAN8", #secretory
                        "KRT5","KRT17","S100A2", #basal
                        "FOXJ1","TMEM190","C1orf194", #ciliated
                        "CHGA","CPE","GRP") #pnec
Idents(AATD)<-"annot3"
levels(AATD)<-rev(c('AT1','AT2','RASC','Secretory','Basal','Ciliated','PNEC','AlvMacrophage','Myeloid','Mast','T&NK','B','Endothelial','Lymphatic','SmoothMuscle','FibroblastsMyofibroblasts','NA'))
DotPlot(AATD,features=epithelial_markers)+ RotatedAxis()
DotPlot(subset(AATD, idents = c('AT1','AT2','RASC','Secretory','Basal','Ciliated','PNEC')),
        features=epithelial_markers) +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  scale_colour_gradient(low = "#e6eeff", high = "#002b80") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

#1D Dot plot of immune cell marker gene expression among immune cell clusters
immune_markers <- c("FABP4","SERPING1","APOC1","CD52", #aegeter alv mac markers
                    "LGMN","MARCKS", #aegeter interstitial macrophages
                    "FCN1","S100A8", #VB Neutrophil
                    #"CLEC10A","PKIB","FCER1A","CD1C","CD1E","CCL17", #VB Dendritic cell - removed due to no expression
                    "CPA3","TPSAB1","TPSB2","MS4A2","HPGDS", #VB Mast cells
                    "CD3E","CD3D","TRAC","CD7","CCL5","TRBC1", #T and NK cells
                    "CD79A","JCHAIN","MZB1","IGKC") #B cells 
Idents(AATD)<-"annot3"
levels(AATD)<-rev(c('AT1','AT2','RASC','Secretory','Basal','Ciliated','PNEC','AlvMacrophage','Neutrophil','Mast','T&NK','B','Endothelial','Lymphatic','SmoothMuscle','FibroblastsMyofibroblasts','NA'))
DotPlot(subset(AATD, idents = c('AlvMacrophage','Neutrophil','Mast','T&NK','B')),
        features=immune_markers) +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  scale_colour_gradient(low = "#e6eeff", high = "#002b80") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

#1E SERPINA1 feature plot
FeaturePlot(AATD,features="SERPINA1")+scale_colour_gradient(low = "#e6eeff", high = "#002b80")

#1F SERPINA1 expression in AT2s and AMs
AT2_alvmacs <- subset(AATD,idents=c("AT2","AlvMacrophage"))
Idents(AT2_alvmacs) <- "annot2"
AT2_alvmacs$annot2 <- factor (x=AT2_alvmacs$annot2, levels=c('Control','COPD','AATD'))
AT2_alvmacs$annot3 <- factor (x=AT2_alvmacs$annot3, levels=c('AT2','AlvMacrophage'))
VlnPlot(AT2_alvmacs,
        features = "SERPINA1",
        assay = "RNA",
        group.by="annot3",
        split.by="annot2",
        pt.size=0.001,
        cols=c("darkslategray4","cornflowerblue","brown2")) +geom_boxplot(width=0.9, color="black",alpha=0) +theme(axis.title.x = element_blank())

##Figure 2
#2A, 2E FGSEA
Idents(AATD) <- "annot3"
AT2_alvmacs <- subset(AATD,idents=c("AT2","AlvMacrophage"))
AT2<- subset(AT2_alvmacs, idents="AT2")
AM <- subset(AT2_alvmacs, idents="AlvMacrophage")
Idents(AT2) <- "annot2"
Idents(AM) <- "annot2"

sc<-AT2  #run for AT2, AM
list.genes <- wilcoxauc(sc, 'annot2')
head(list.genes)
dplyr::count(list.genes, group)
m_df<- msigdbr(species = "Homo sapiens", category = "H") #hallmark gene sets
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

AATD.genes<- list.genes %>%
  dplyr::filter(group == "AATD") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

ranks<- deframe(AATD.genes)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

ggplot(fgseaResTidy %>% filter(padj < 0.008) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


#2B, 2E Violin plot: Hallmark NFkB in AT2s, AM
Hallmark_NFkB <- list(c("ABCA1","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40","BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCND1","CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6","CXCR7","CYR61","DDX58","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1","EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM","GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2","IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1","IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF","LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5","NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1","PANX1","PDE4B","PDLIM5","PER1","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR","PLEK","PLK2","PMEPA1","PNRC1","PPAP2B","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1","REL","RELA","RELB","RHOB","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2","TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10","TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36"))
Idents(AATD) <- "annot3"
AT2_alvmac <- AddModuleScore(AT2_alvmac,Hallmark_NFkB, name = "Hallmark_NFkB")
VlnPlot(subset(AT2_alvmacs, idents="AT2"), 
        features="Hallmark_NFkB1",
        assay = "RNA",
        split.by = "annot2",
        pt.size=0.001,
        cols=c("darkslategray4","cornflowerblue","brown2")) +
  ggtitle("TNFa Signaling via NFkB (Hallmark)")+
  geom_boxplot(width=0.9, color="black",alpha=0) +
  theme(axis.title.x = element_blank())
VlnPlot(subset(AT2_alvmac, idents="AlvMacrophage"),
        features="Hallmark_NFkB1",
        assay = "RNA",
        split.by = "annot2",
        pt.size=0.001,
        cols=c("darkslategray4","cornflowerblue","brown2")) +
  ggtitle("TNFa Signaling via NFkB (Hallmark)")+
  geom_boxplot(width=0.9, color="black",alpha=0) +
  theme(axis.title.x = element_blank())

#2C, 2G Violin plot: Hallmark UPR in AT2s, AM
Hallmark_UPR <- list(c("ALDH18A1","ARFGAP1","ASNS","ATF3","ATF4","ATF6","ATP6V0D1","BAG3","BANF1","CALR","CCL2","CEBPB","CEBPG","CHAC1","CKS1B","CNOT2","CNOT4","CNOT6","CXXC1","DCP1A","DCP2","DCTN1","DDIT4","DDX10","DKC1","DNAJA4","DNAJB9","DNAJC3","EDC4","EDEM1","EEF2","EIF2AK3","EIF2S1","EIF4A1","EIF4A2","EIF4A3","EIF4E","EIF4EBP1","EIF4G1","ERN1","ERO1A","EXOC2","EXOSC1","EXOSC10","EXOSC2","EXOSC4","EXOSC5","EXOSC9","FKBP14","FUS","GEMIN4","GOSR2","H2AX","HERPUD1","HSP90B1","HSPA5","HSPA9","HYOU1","IARS1","IFIT1","IGFBP1","IMP3","KDELR3","KHSRP","KIF5B","LSM1","LSM4","MTHFD2","MTREX","NABP1","NFYA","NFYB","NHP2","NOLC1","NOP14","NOP56","NPM1","PAIP1","PARN","PDIA5","PDIA6","POP4","PREB","PSAT1","RPS14","RRP9","SDAD1","SEC11A","SEC31A","SERP1","SHC1","SLC1A4","SLC30A5","SLC7A5","SPCS1","SPCS3","SRPRA","SRPRB","SSR1","STC2","TARS1","TATDN2","TSPYL2","TTC37","TUBB2A","VEGFA","WFS1","WIPI1","XBP1","XPOT","YIF1A","YWHAZ","ZBTB17"))
Idents(AATD) <- "annot3"
AT2_alvmacs <- AddModuleScore(AT2_alvmacs,Hallmark_UPR, name = "Hallmark_UPR")
AT2_alvmacs$annot2 <- factor (x=AT2_alvmac$annot2, levels=c('Control','COPD','AATD'))
FeaturePlot(AT2_alvmac, features = "Hallmark_UPR1")
VlnPlot(subset(AT2_alvmacs, idents="AT2"),
        features="Hallmark_UPR1",
        assay = "RNA",
        split.by = "annot2",
        pt.size=0.001,
        cols=c("darkslategray4","cornflowerblue","brown2")) +
  ggtitle("Unfolded Protein Response (Hallmark)")+
  geom_boxplot(width=0.9, color="black",alpha=0) +
  theme(axis.title.x = element_blank())
VlnPlot(subset(AT2_alvmacs, idents="AlvMacrophage"),
        features="Hallmark_UPR1",
        assay = "RNA",
        split.by = "annot2",
        pt.size=0.001,
        cols=c("darkslategray4","cornflowerblue","brown2")) +
  ggtitle("Unfolded Protein Response (Hallmark)")+
  geom_boxplot(width=0.9, color="black",alpha=0) +
  theme(axis.title.x = element_blank())

#2D,2H Regulon analysis
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

##obtain the regulons based on interactions with confidence level A 
regulon <- dorothea_regulon_human %>%
  dplyr::filter(tf %in% c("AR","ARNTL","ATF1","ATF2","ATF4","ATF6","CDX2","CEBPA","CEBPB","CREB1","CTCF","E2F1","E2F2","E2F3","E2F4","EGR1","ELK1","EPAS1","ERG","ESR1","ESR2","ETS1","ETS2","ETV4","FLI1","FOS","FOSL1","FOSL2","FOXA1","FOXL2","FOXM1","FOXO1","FOXO3","FOXO4","GATA1","GATA2","GATA3","GLI2","HIF1A","HNF1A","HNF4A","IRF1","JUN","JUND","KLF4","KMT2A","LEF1","MITF","MYB","MYC","MYCN","NFATC2","NFIC","NFKB1","NR2F2","NR3C1","NR5A1","PAX6","PAX8","PGR","POU2F1","PPARA","PPARG","RARA","REL","RELA","RFX5","RUNX1","RXRA","SMAD3","SMAD4","SOX10","SOX2","SOX9","SP1","SP3","SPI1","SREBF1","SREBF2","SRF","STAT1","STAT3","STAT5A","STAT5B","STAT6","TAL1","TCF7L2","TFAP2A","TFAP2C","TP53","TWIST1","USF1","USF2","VDR","WT1","YY1"))
AT2 <- ScaleData(AT2)
Idents(AT2) <- "annot3"

AT2 <- run_viper(AT2, regulon,
                 options = list(method = "scale", minsize = 1, 
                                eset.filter = FALSE, cores = 1, 
                                verbose = FALSE))
DefaultAssay(object = AT2) <- "dorothea"
AT2 <- ScaleData(AT2)
viper_scores_df <- GetAssayData(AT2, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame() %>%
  t()

CellsClusters <- data.frame(cell = names(Idents(AT2)),
                            cell_type = as.character(AT2$annot2),
                            stringsAsFactors = FALSE,
                            check.names=FALSE) %>% dplyr::mutate(., cell = sapply(.$cell, function(x) sub("-", ".", x)) %>% unname())

viper_scores_clusters <- viper_scores_df  %>% 
  data.frame(check.names=FALSE) %>% 
  rownames_to_column("cell") %>% 
  gather(tf, activity, -cell) %>% 
  inner_join(CellsClusters)

summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(200, var) %>%
  distinct(tf)


summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

#reorder by disease state
summarized_viper_scores_df <- summarized_viper_scores_df[c("Control","COPD","AATD"), ]

#reorder by difference (AATD-Control)
summarized_viper_scores_df <- summarized_viper_scores_df[,c("RELA","NFKB1","REL","FOXL2","RARA","CEBPB","JUN","SMAD3","CEBPA","FOXO3","FOXO4","CREB1","ETS2","HIF1A","FLI1","ESR1","EPAS1","ETS1","LEF1","E2F3","PPARA","EGR1","HNF4A","FOS","GATA3","FOXM1","POU2F1","ATF6","ELK1","MITF","FOSL2","FOXO1","CDX2","ETV4","PGR","MYB","KLF4","PAX6","ATF4","AR","FOXA1","FOSL1","NR5A1","ATF2","KMT2A","MYC","IRF1","RXRA","RFX5","NR3C1")]

palette_length = 100
my_color = colorRampPalette(c("#0517BA", "white","#BA1505"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "AT2 Regulon Expression", angle_col = 0,
                       treeheight_col = 0,  border_color = NA,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE) 


##repeated above for alveolar macrophages
AM <- ScaleData(AM)
Idents(AM) <- "annot3"

