#This script reproduces the single cell RNAseq analysis of KF and KFCA mice described in Menjivar et al., 2023
#Data was processed in line with the Seurat workflow:
#Website: https://satijalab.org/seurat/index.html
#
#Reference: Stuart et al., Comprehensive Integration of Single-Cell Data. 
#Cell, 2019;177(7):1888-1902.e21. PMID: 31178118 PMCID: PMC6687398. 
#
#Figures 2G-J, 2S2B-E, 2S3C, 2S3F-H, 3B-F, and 3I are represented in this script. 
#Preprocessing of data for Figure 1C&D and 1S1A-B are described in https://github.com/PascaDiMagliano-Lab/MultimodalMappingPDA-scRNASeq. 
#Preprocessing of data for Figure 1F-G and 1S1C-E are described in https://github.com/PascaDiMagliano-Lab/-Inactivation-of-WNT-signaling-sensitizes-pancreatic-cancer-to-immunotherapy-in-mice.
#
#---------------------------------------------------------------#
#R version 4.1.1 (2021-08-10) -- "Kick Things"
#Seurat Version 4.0.4
#SeuratObject Version 4.0.2
#clusterProfiler Version 4.0.5

#Load Required Packages:
library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)
library(tibble)
library(RColorBrewer)
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db) #Mouse data set
library(DOSE)
library(ReactomePA)

#### Load Required Functions: -------------------------------------------------------------------------------------------------- ####
cleanDE <- function(DE_name){
  DE_name <- subset(DE_name, p_val_adj <= 0.05)
  DE_name$gene <- row.names(DE_name)
  
  rpl <- grep("Rpl", row.names(DE_name), value = T)
  rps <- grep("Rps", row.names(DE_name), value = T)
  mt <- grep("mt-", row.names(DE_name), value = T)
  hbb <- grep("Hbb-", row.names(DE_name), value = T)
  hba <- grep("Hba-", row.names(DE_name), value = T)
  
  x <- DE_name[ !(DE_name$gene %in% rpl), ]
  x <- x[ !(x$gene %in% rps), ]
  x <- x[ !(x$gene %in% hbb), ]
  x <- x[ !(x$gene %in% hba), ]
  DE_name <- x[ !(x$gene %in% mt), ]
  return(DE_name)
}
#### Preprocessing: ------------------------------------------------------------------------------------------------------------ ####

#Load in Datasets:
KF_data <- Read10X("~raw_feature_bc_matrix")
KF_Arg1KO_data <- Read10X("~raw_feature_bc_matrix")

#Create Seurat Objects:
KF <- CreateSeuratObject(KF_data, min.cells = 3, min.features = 100)
KF_Arg1KO <- CreateSeuratObject(KF_Arg1KO_data, min.cells = 3, min.features = 100)

#Add Desired Metadata:
KF[["Group"]] <- "KF"
KF_Arg1KO[["Group"]] <- "KFCA"

#Merge Seurat Objects:
KF_KFArg1KO_merged <- merge(x = KF, y = KF_Arg1KO, add.cell.ids = (c("c", "f")))

#Normalize Data:
KF_KFArg1KO_merged  <- NormalizeData(object = KF_KFArg1KO_merged, normalization.method = "LogNormalize", scale.factor = 10000)

#Apply Unbiased QC Cutoffs:
KF_KFArg1KO_merged [["percent.mt"]] <- PercentageFeatureSet(object = KF_KFArg1KO_merged, pattern = "^mt-")
KF_KFArg1KO_merged  <- subset(x = KF_KFArg1KO_merged, subset = nCount_RNA > 1000 & nCount_RNA < 60000 & percent.mt < 15)

#Find Variable Genes:
KF_KFArg1KO_merged  <- FindVariableFeatures(object = KF_KFArg1KO_merged, selection.method = "vst", nfeatures = 2000)

#Scale Data:
KF_KFArg1KO_merged  <- ScaleData(object = KF_KFArg1KO_merged , vars.to.regress = c("nCount_RNA"), features = rownames(KF_KFArg1KO_merged))

#Run PCA and Determine Dimensions for 90% Variance:
KF_KFArg1KO_merged  <- RunPCA(KF_KFArg1KO_merged, verbose = FALSE)
st_dev <- KF_KFArg1KO_merged@reductions$pca@stdev
var <- st_dev^2
sum(var[1:28])/ sum(var) 

#Find Neighbors and Cluster Cells:
KF_KFArg1KO_merged<- FindNeighbors(object = KF_KFArg1KO_merged, dims = 1:28) 
KF_KFArg1KO_merged<- FindClusters(object = KF_KFArg1KO_merged, resolution = 1.7)

#Run UMAP:  
KF_KFArg1KO_merged<- RunUMAP(object = KF_KFArg1KO_merged, dims = 1:28) 
DimPlot(object = KF_KFArg1KO_merged, reduction = "umap", label = T, pt.size = 1)

#Identify Red Blood Cell Clusters:
FeaturePlot(iKras, features = c("Hbb-bt", "Ptprc", "Krt19", "Try4", "Col1a2"), cols = c("gainsboro", "firebrick1"))
DotPlot(iKras, features = rev(c("Hbb-bt", "Ptprc", "Krt19", "Try4", "Col1a2")), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

#Label red blood cell clusters:
KF_KFArg1KO_merged <- RenameIdents(KF_KFArg1KO_merged,
                                   "0" = "RBC", 
                                   "1" = "RBC", 
                                   "2" = "Keep", 
                                   "3" = "Keep", 
                                   "4" = "Keep", 
                                   "5" = "Keep", 
                                   "6" = "Keep",  
                                   "7" = "Keep", 
                                   "8" = "Keep", 
                                   "9" = "Keep", 
                                   "10" = "Keep", 
                                   "11" = "Keep", 
                                   "12" = "Keep", 
                                   "13" = "Keep", 
                                   "14" = "Keep", 
                                   "15" = "Keep", 
                                   "16" = "Keep", 
                                   "17" = "Keep", 
                                   "18" = "RBC", 
                                   "19" = "Keep", 
                                   "20" = "Keep", 
                                   "21" = "Keep", 
                                   "22" = "Keep", 
                                   "23" = "Keep", 
                                   "24" = "Keep", 
                                   "25" = "Keep", 
                                   "26" = "Keep", 
                                   "27" = "Keep", 
                                   "28" = "Keep", 
                                   "29" = "Keep")

#Remove red blood cells:
KF_KFArg1KO_merged <- subset(KF_KFArg1KO_merged, idents = "Keep")

#### Creating Annotated Global Object: ----------------------------------------------------------------------------------------- ####

#Find Variable Genes:
KF_KFArg1KO_merged <- FindVariableFeatures(object = KF_KFArg1KO_merged, selection.method = "vst", nfeatures = 2000)

#Scale Data:
KF_KFArg1KO_merged <- ScaleData(object = KF_KFArg1KO_merged, vars.to.regress = c("nCount_RNA"), features = rownames(KF_KFArg1KO_merged))

#Run PCA and Determine Dimensions for 90% Variance:
KF_KFArg1KO_merged <- RunPCA(KF_KFArg1KO_merged, verbose = FALSE)
st_dev <- KF_KFArg1KO_merged@reductions$pca@stdev
var <- st_dev^2
sum(var[1:27])/ sum(var)

#Find Neighbors and Cluster Cells:
KF_KFArg1KO_merged <- FindNeighbors(object = KF_KFArg1KO_merged, dims = 1:27)
KF_KFArg1KO_merged <- FindClusters(object = KF_KFArg1KO_merged, resolution = 1.7)

#Run UMAP:
KF_KFArg1KO_merged <- RunUMAP(object = KF_KFArg1KO_merged, dims = 1:27)
DimPlot(object = KF_KFArg1KO_merged, reduction = "umap", label = T, pt.size = 1)

#Identify Clusters Based On Marker Expression:
DotPlot(KF_KFArg1KO_merged, features = c("Krt19", "Krt18","Try4", "Col1a2", "Acta2", "Clec3b", "Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3", "Ccr7", "Cd8a", "Nkg7","Il2ra","Il1rl1","Arg1","Arg2", "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                                         "Itgae","Clec9a","Batf3", "Ly6c2", "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Ccna2", "Pecam1", "Cdh5", "Mki67", "Hbb-bt"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

#Annotate Clusters:
KF_KFArg1KO_merged <- RenameIdents(KF_KFArg1KO_merged,"0" = "Fibroblast", 
                                   "1" = "Fibroblast", 
                                   "2" = "Fibroblast", 
                                   "3" = "Fibroblast",
                                   "4" = "Fibroblast",
                                   "5" = "Cd4", 
                                   "6" = "Cd8", 
                                   "7" = "Treg", 
                                   "8" = "Th2", 
                                   "9" = "Fibroblast",
                                   "10" = "Acinar", 
                                   "11" = "Macrophage", 
                                   "12" = "Fibroblast",
                                   "13" = "Acinar", 
                                   "14" = "NKT",
                                   "15" = "Fibroblast",
                                   "16" = "Mesothelial",
                                   "17" = "B cell", 
                                   "18" = "Endothelial 1",
                                   "19" = "Macrophage", 
                                   "20" = "Plasma cell",
                                   "21" = "DC",
                                   "22" = "Macrophage", 
                                   "23" = "Proliferating",
                                   "24" = "gdT", 
                                   "25" = "MDSC",
                                   "26" = "NK",
                                   "27" = "Mast cells",
                                   "28" = "Endothelial 2",
                                   "29" = "Myeloid",
                                   "30" = "Epithelial")
KF_KFArg1KO_merged[["manual_clusters"]] <- KF_KFArg1KO_merged@active.ident

KF_KFArg1KO_merged <- RenameIdents(KF_KFArg1KO_merged,"Fibroblast" = "Fibroblasts",
                                   "Cd4" = "CD4+ T Cells",
                                   "Cd8" = "CD8+ T Cells",
                                   "Treg" = "Tregs",
                                   "Th2" = "CD4+ T Cells",
                                   "gdT" = "gd T Cells",
                                   "Acinar" = "Acinar Cells",
                                   "Macrophage" = "Macrophages",
                                   "NKT" = "NK T Cells",
                                   "Mesothelial" = "Mesothelial Cells",
                                   "B cell" = "B Cells",
                                   "Endothelial 1" = "Endothelial Cells",
                                   "Plasma cell" = "Plasma Cells",
                                   "DC" = "Dendritic Cells",
                                   "Proliferating" = "Proliferating Cells",
                                   "MDSC" = "MDSCs",
                                   "NK" = "NK Cells",
                                   "Mast cells" = "Mast Cells",
                                   "Endothelial 2" = "Endothelial Cells",
                                   "Myeloid" = "Dendritic Cells",
                                   "Epithelial" = "Epithelial Cells")
new_order <- c("Epithelial Cells", "Acinar Cells", "Mesothelial Cells", "Fibroblasts", "Endothelial Cells", "Macrophages",  "MDSCs", 
               "Dendritic Cells", "CD4+ T Cells","CD8+ T Cells", "Tregs", "gd T Cells", "NK T Cells", "NK Cells", "B Cells", 
               "Plasma Cells", "Mast Cells", "Proliferating Cells")
KF_KFArg1KO_merged@active.ident <- factor(KF_KFArg1KO_merged@active.ident, levels = new_order)
KF_KFArg1KO_merged[["collapsed_clusters"]] <- KF_KFArg1KO_merged@active.ident

#Save Seurat Object:
save(KF_KFArg1KO_merged, file = 'KF_KFArg1KO_merged.RData')

# Identify Switches - Run these to switch the identity of your object.
Idents(object = KF_KFArg1KO_merged) <- "seurat_clusters"
Idents(object = KF_KFArg1KO_merged) <- "Group"
Idents(object = KF_KFArg1KO_merged) <- "manual_clusters"
Idents(object = KF_KFArg1KO_merged) <- "collapsed_clusters"

#### Annotating Myeloid Cells: ------------------------------------------------------------------------------------------------- ####

#Subset Myeloid Cells:
Myeloid <- subset(x = KF_KFArg1KO_merged, idents = c( "Macrophages", "Dendritic Cells", "Proliferating Cells"))

#Find Variable Genes:
Myeloid <- FindVariableFeatures(object = Myeloid, selection.method = "vst", nfeatures = 2000)

#Scale Data:
Myeloid <- ScaleData(object =Myeloid, vars.to.regress = c("nCount_RNA"), features = rownames(Myeloid))

#Run PCA and Determine Dimensions for 90% Variance:
Myeloid <- RunPCA(Myeloid, verbose = FALSE)
st_dev <- Myeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:35])/ sum(var) 

#Find Neighbors and Cluster Cells:
Myeloid <- FindNeighbors(object = Myeloid, dims = 1:35)
Myeloid <- FindClusters(object = Myeloid, resolution = 1.2)

#Run UMAP:
Myeloid <- RunUMAP(object = Myeloid, dims = 1:35) 
DimPlot(object = Myeloid, reduction = "umap", label = TRUE, pt.size = 2)

#Identify Clusters Based On Marker Expression:
DotPlot(Myeloid, features = c("Il1rl1","Arg1","Arg2", "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                              "Itgae","Clec9a","Batf3", "Ly6c2", "Mki67", "Col1a1", "Clec3b", "Acta2", "Ptprc"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

#Annotate Clusters:
Idents(Myeloid) <- "seurat_clusters"
Myeloid <- RenameIdents(Myeloid,
                        "0" = "Macrophage", 
                        "1" = "Macrophage", 
                        "2" = "Macrophage", 
                        "3" = "DC", 
                        "4" = "Macrophage", 
                        "5" = "Macrophage", 
                        "6" = "Fibroblast", 
                        "7" = "Macrophage", 
                        "8" = "Cycling", 
                        "9" = "Cycling DC", 
                        "10" = "DC", 
                        "11" = "DC")

#Subset Only Macrophages:
Myeloid <- subset(x = Myeloid, idents = "Macrophage")

#Find Variable Genes:
Myeloid <- FindVariableFeatures(object = Myeloid, selection.method = "vst", nfeatures = 2000)

#Scale Data:
Myeloid <- ScaleData(object =Myeloid, vars.to.regress = c("nCount_RNA"), features = rownames(Myeloid))

#Run PCA and Determine Dimensions for 90% Variance:
Myeloid <- RunPCA(Myeloid, verbose = FALSE)
st_dev <- Myeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:38])/ sum(var) 

#Find Neighbors and Cluster Cells:
Myeloid <- FindNeighbors(object = Myeloid, dims = 1:38) 
Myeloid <- FindClusters(object = Myeloid, resolution = 1.2) 

#Run UMAP:
Myeloid <- RunUMAP(object = Myeloid, dims = 1:38)
DimPlot(object = Myeloid, reduction = "umap", label = TRUE, pt.size = 2)

#Annotate Clusters:
Idents(Myeloid) <- "seurat_clusters"
Myeloid <- RenameIdents(Myeloid,
                        "0" = "Macrophage 1",
                        "5" = "Macrophage 1",
                        "1" = "Macrophage 2",
                        "2" = "Macrophage 3", 
                        "4" = "Macrophage 4",
                        "3" = "Macrophage 5")
Myeloid[["manual_clusters"]] <- Myeloid@active.ident

#Save Seurat Object:
save(Myeloid, file = 'Myeloid.RData')

#### Annotating T Cells: ------------------------------------------------------------------------------------------------------- ####

#Subset NK and T Cells:
NKTcells <- subset(KF_KFArg1KO_merged, idents = c("NK","NKT","Cd8","Cd4","Treg","Th2","gdT" ))

#Find Variable Genes:
NKTcells <- FindVariableFeatures(object = NKTcells, selection.method = "vst", nfeatures = 2000)

#Scale Data:
NKTcells <- ScaleData(object =NKTcells, vars.to.regress = c("nCount_RNA"), features = rownames(NKTcells))

#Run PCA and Determine Dimensions for 90% Variance:
NKTcells <- RunPCA(NKTcells, verbose = FALSE)
st_dev <- NKTcells@reductions$pca@stdev
var <- st_dev^2
sum(var[1:37])/ sum(var) 

#Find Neighbors and Cluster Cells:
NKTcells <- FindNeighbors(object = NKTcells, dims = 1:37)
NKTcells <- FindClusters(object = NKTcells, resolution = 1.2) 

#Run UMAP:
NKTcells <- RunUMAP(object = NKTcells, dims = 1:37) 
DimPlot(object = NKTcells, reduction = "umap", label = T, pt.size = 2)

#Annotate Clusters:
Idents(NKTcells) <- "seurat_clusters"
NKTcells <- RenameIdents(NKTcells,
                         "0" = "Exhausted CD8+ T Cell", 
                         "1" = "Treg", 
                         "2" = "CD4 T Helper", 
                         "3" = "gd T Cell", 
                         "4" = "CD4 T Helper", 
                         "5" = "Treg", 
                         "6" = "Naive CD8+ T Cell", 
                         "7" = "Myeloid Cont", 
                         "8" = "Not NKT cells", 
                         "9" = "Cytotoxic CD8+ T Cell", 
                         "10" = "NKT Cell", 
                         "11" = "Not NKT cells",
                         "12" = "CD4 T Helper",
                         "13" = "NK Cell" )
new_order <- c("Naive CD8+ T Cell","Cytotoxic CD8+ T Cell", "Exhausted CD8+ T Cell", "CD4 T Helper","Treg", "gd T Cell","NKT Cell","NK Cell","Not NKT cells","Myeloid Cont")
NKTcells@active.ident <- factor(NKTcells@active.ident, levels = new_order)
NKTcells[["manual_clusters"]] <- NKTcells@active.ident

#Remove Contaminating Cells:
NKTcellslabeled <- subset(NKTcells, idents = c("Not NKT cells","Myeloid Cont"), invert = TRUE)

#Save Seurat Object:
save(NKTcellslabeled,file="NKTcellslabeled.RData")

#### Figures: ------------------------------------------------------------------------------------------------------------------ ####
#### Figure 2G, Global UMAP ####
DimPlot(object = KF_KFArg1KO_merged, reduction = "umap", label = F, pt.size = 1, cols = c("Fibroblasts" = "#a6cee3",
                                                                                          "CD4+ T Cells" = "#cab2d6",
                                                                                          "CD8+ T Cells"= "#e31a1c",
                                                                                          "Tregs" = "#1b7837",
                                                                                          "Acinar Cells" = "#de77ae",
                                                                                          "Macrophages" = "#7fbc41",
                                                                                          "NK T Cells" = "#ff7f00",
                                                                                          "Mesothelial Cells" = "#80cdc1",
                                                                                          "B Cells" = "#b15928",
                                                                                          "Endothelial Cells" = "#192CB5",
                                                                                          "Plasma Cells" = "#bf812d",
                                                                                          "Dendritic Cells" = "#1f78b4",
                                                                                          "Proliferating Cells" = "#c51b7d",
                                                                                          "gd T Cells" = "#ffd700",
                                                                                          "MDSCs"= "#800026",
                                                                                          "NK Cells" = "#434343",
                                                                                          "Mast Cells" = "#fb9a99",
                                                                                          "Epithelial Cells"= "#6a3d9a"))

#### Figure 2H, Macrophage DE Heatmap ####
mac_markers <- FindAllMarkers(Myeloid, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
mac_markers <- cleanDE(mac_markers)
top10 <- mac_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
mac_colors <- c("Macrophage 1" = "#FE361E", "Macrophage 2" = "#243372","Macrophage 3" = "#83938D", "Macrophage 4" = "#28AFEE", "Macrophage 5" = "#fb9600")
DoHeatmap(Myeloid, features = top10$gene, group.colors = mac_colors)

#### Figure 2I, Macrophage Arg1 and Arg2 Dotplot ####
DotPlot(Myeloid, features = c("Arg1","Arg2"), cols = "RdBu", dot.scale = 8, split.by = "Group") + RotatedAxis() 

#### Figure 2J and Figure 2 S3G, Macrophage Pathway Analysis via Reactome ####

#Create DE List of Macrophages in KF vs KFCA:
mac_de <- FindMarkers(Myeloid, ident.1 = "KF", ident.2 = "KF_Arg1KO", group.by = "Group", test.use = "MAST") 
mac_de <- cleanDE(mac_de)

#Isolate Genes Enriched in Each Group:
up_in_KF <- subset(mac_de, avg_log2FC > 0)
up_in_KO <- subset(mac_de, avg_log2FC < 0)

#Make Gene List for KF:
gene_list_KF <- up_in_KF$gene_names
gene_list_KF <- replace(gene_list_KF, gene_list_KF=="AC149090.1", "Pisd")
gene_convert_KF <- bitr(gene_list_KF, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Mm.eg.db", drop = F)

#Make Gene List for KO:
gene_list_KO <- up_in_KO$gene_names
gene_convert_KO <- bitr(gene_list_KO, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Mm.eg.db", drop = F)

#Run Gene Lists Through Reactome and Visualize (Figure 2J):
Reactome_KO <- enrichPathway(gene=gene_convert_KO$ENTREZID, pvalueCutoff = .05, organism = 'mouse', readable = T)
dotplot(Reactome_KO, showCategory=5)

#Run Gene Lists Through Reactome and Visualize (Figure 2 S3G):
Reactome_KF <- enrichPathway(gene=gene_convert_KF$ENTREZID, pvalueCutoff = .05, organism = 'mouse', readable = T)
dotplot(Reactome_KF, showCategory=5)

#### Figure 2 S2B, Macrophage Global UMAP Split by Condition ####
DimPlot(object = KF_KFArg1KO_merged, reduction = "umap", label = F, pt.size = 1, cols = c("Fibroblasts" = "#a6cee3",
                                                                                          "CD4+ T Cells" = "#cab2d6",
                                                                                          "CD8+ T Cells"= "#e31a1c",
                                                                                          "Tregs" = "#1b7837",
                                                                                          "Acinar Cells" = "#de77ae",
                                                                                          "Macrophages" = "#7fbc41",
                                                                                          "NK T Cells" = "#ff7f00",
                                                                                          "Mesothelial Cells" = "#80cdc1",
                                                                                          "B Cells" = "#b15928",
                                                                                          "Endothelial Cells" = "#192CB5",
                                                                                          "Plasma Cells" = "#bf812d",
                                                                                          "Dendritic Cells" = "#1f78b4",
                                                                                          "Proliferating Cells" = "#c51b7d",
                                                                                          "gd T Cells" = "#ffd700",
                                                                                          "MDSCs"= "#800026",
                                                                                          "NK Cells" = "#434343",
                                                                                          "Mast Cells" = "#fb9a99",
                                                                                          "Epithelial Cells"= "#6a3d9a"), split.by = "Group")

#### Figure 2 S2C, Global Dotplot ####
DotPlot(KF_KFArg1KO_merged, features = c("Krt19", "Cdh1","Try4", "Amy2a2","Msln", "Col1a2", "Pdpn","Pecam1", "Cdh5", "Ptprc", "Cd68","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1","Itgae","Clec9a","Batf3", "Cd3e", "Cd4", "Cd8a", "Foxp3", "Il2ra", "Trdc", "Nkg7","Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Il1rl1","Ccna2", "Mki67"), cols = "RdBu", dot.scale = 8) + RotatedAxis()

#### Figure 2 S2D&E, Population Numbers ####
Idents(KF_KFArg1KO_merged) <- "Group"
KF <- subset(KF_KFArg1KO_merged, idents = "KF")
KF_ARG1KO <- subset(KF_KFArg1KO_merged, idents = "KF_Arg1KO")
Idents(KF) <- "collapsed_clusters"
Idents(KF_ARG1KO) <- "collapsed_clusters"
table(KF@active.idents)
table(KF_ARG1KO@active.idents)

#### Figure 2 S3C, Split Macrophage UMAP ####
DimPlot(object = Myeloid, reduction = "umap", label = FALSE, pt.size = 2, split.by= "Group")
#### Figure 2 S3F, Split Macrophage Apoe Dotplot ####
DotPlot(Myeloid, features = c("Apoe"), cols = "RdBu", dot.scale = 8, split.by = "Group") + RotatedAxis()
#### Figure 2 S3H, Antigen Presentation DE Dotplot ####
DotPlot(Myeloid, features = c("H2-K1","Psmb8","H2-Q7", "Ubb", "H2-T23", "Psmb9", "H2-T22", "Uba52", "H2-Q6", "Mrc1", "Rbx1", "Psmd8", "Psmb5", "Znrf2"), cols = "RdBu", dot.scale = 8, split.by = "Group") + RotatedAxis()
#### Figure 3B, Split Lymphocyte UMAP ####
DimPlot(object = NKTcellslabeled, reduction = "umap", label = F, pt.size = 2, cols = c(  "Naive CD8+ T Cell"= "#f4de13",
                                                                                         "Cytotoxic CD8+ T Cell"= "#f21914",
                                                                                         "Exhausted CD8+ T Cell"= "#1c2bf4",
                                                                                         "CD4 T Helper"= "#b3d1de",
                                                                                         "Treg"= "#cebbea",
                                                                                         "gd T Cell"= "#ff99cc",
                                                                                         "NKT Cell"="#F2BE8F",
                                                                                         "NK Cell"="#434343"), split.by = "Group")
#### Figure 3C, Lymphocyte Marker Dimplot ####
DotPlot(NKTcells, features = c('Cd3e','Cd4','Foxp3','Cxcr3','Tbx21','Ifng','Il17a','Ccr6','Rorc','Eomes','Pdcd1','Ctla4','Gzmk','Cd8a','Ccr7','Sell','Nkg7','Gzmb','Trdc'), cols = "PuBuGn", dot.scale = 6) + RotatedAxis()

#### Figure 3D&E, Population Numbers ####
Idents(object = NKTcellslabeled) <- 'Collapsed_labels'
table (Idents(NKTcellslabeled))
prop.table(table(Idents(NKTcellslabeled), NKTcellslabeled$Group), margin = 2)
table(Idents(NKTcellslabeled), NKTcellslabeled$Group)
CellNumberTable.collapsed <- table(Idents(NKTcellslabeled), NKTcellslabeled$Group)
write.csv(CellNumberTable.collapsed, file = 'CellNumberTable.collapsed.csv')

#### Figure 3F&I, Violin Plots ####

#Subset CD8 T Cells:
CD8 <- subset(NKTcells, idents = c("Naive CD8+ T Cell","Cytotoxic CD8+ T Cell", "Exhausted CD8+ T Cell"))

#Figure 3F:
VlnPlot(CD8, features = c("Gzmb"),group.by = "Group")
VlnPlot(CD8, features = c("Prf1"),group.by = "Group")
VlnPlot(CD8, features = c("Ifng"),group.by = "Group")

#Figure 3I:
VlnPlot(CD8, features = c("Ctla4"),group.by = "Group")
VlnPlot(CD8, features = c("Furin"),group.by = "Group")
VlnPlot(CD8, features = c("Lag3"),group.by = "Group")
VlnPlot(CD8, features = c("Pdcd1"),group.by = "Group")