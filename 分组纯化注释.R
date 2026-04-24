library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggplot2)
setwd("D:/R/R/PROJECT/IMNM/MergeMC")


#1.定义屏蔽marker
muscle_genes <- c(
  # 原有核心肌肉基因
  "MYH1", "MYH2", "MYH3", "MYH7",
  "ACTA1", "ACTN2",
  "TTN", "NEB",
  "TNNT1", "TNNT2", "TNNT3",
  "TNNI1", "TNNI2",
  "TPM1", "TPM2", "TPM3",
  "MYBPC1", "MYBPC2",
  "MYL1", "MYL2", "MYL3", "MYLPF",
  "CKM", "DES", "DMD", "XIRP2", "CMYA5", "MYOT",
  
  # 新增：基于近期分析发现的高表达污染标记
  "TNNC1", "TNNC2",       # 肌钙蛋白C亚型
  "TCAP",                 # 肌节相关蛋白
  "CSRP3",                # 肌肉 LIM 蛋白
  "TRDN",                 # 三联体蛋白 (Triadin)
  "MYL2",                 # 肌球蛋白轻链（已包含，再次确认）
  "CAPN3",                # 钙蛋白酶3（常见肌肉污染）
  "MYOG",                 # 肌生成素（成肌标记）
  "ACTC1",                 # 心肌/骨骼肌肌动蛋白
  "CA3",
  "CHAMP",
  "GGTA1P"
)
rbc_genes <- c("HBB", "HBA1", "HBA2", "HBM", "HBQ1", "HBZ") # 游离血红蛋白污染
ig_genes <- grep("^IG[HKL][VDJC]", rownames(srat), value = TRUE)
mt_genes <- grep("^MT-", rownames(srat), value = TRUE)
rrna_genes <- c("MT-RNR1", "MT-RNR2", "RNA18SN5", "RNA28SN5")
neural_genes <- c(
  "NOVA1","ROBO2","LRRTM4",
  "BDNF","DLGAP2","ELFN1"
)
y_chromosome_genes <- c(
  "DDX3Y", "PRKY", "USP9Y", "UTY", "TTTY14", "TTTY15", "NLGN4Y"
)
x_genes <- c("XIST", "TSIX", "FTX")

low_quality_genes <- c(
  "MALAT1", "NEAT1", "XIST", "HSPA1A", "HSPA1B", "FOS", "JUN", "JUNB",
  "MTRNR2L1", "MTRNR2L8", "MTRNR2L12"
)
fibroblast_genes <- c(
  "CTGF", "COL1A1", "COL1A2", "COL3A1", "COL6A2", "LUM", "APOD", "MFAP5",
  "PCOLCE2", "PRG4", "ADAMTS5", "DCN", "FAP", "PDGFRA"
)
smc_pericyte_genes <- c(
  "MYOCD", "CARMN", "RGS6", "TBX2", "GJC1", "EGFLAM", "CDH6", "RIMS1", "NFASC"
)
endothelial_genes <- c(
  "PTPRB", "MECOM", "TEK", "BTNL9", "ANO2", "NOTCH4", "ANXA3", "DIPK2B",
  "CNTNAP3B", "CYYR1", "VWF", "FLT1", "KDR"
)
satellite_genes <- c(
  "PAX7", "CALCR", "CDH4", "GNA14", "MEGF10", "MUSK"
)
adipocyte_genes <- c(
  "CIDEA", "CIDEC", "PLIN1", "PLIN5", "ADIPOQ", "GYG2", "GLYAT", "KLB", "DGAT2"
)

#合并
exclude_genes <- unique(c(
  muscle_genes,
  rbc_genes,
  ig_genes,
  mt_genes,
  rrna_genes,
  y_chromosome_genes,
  neural_genes
  #x_genes
  #fibroblast_genes,
  #smc_pericyte_genes,
  #endothelial_genes,
  #satellite_genes,
  #adipocyte_genes,
  #low_quality_genes
))
#二、分组处理
srat<-readRDS("srat_myeloid_subset.rds")#change rds name to switch dataset
srat <- FindVariableFeatures(srat, nfeatures = 3000)
features_use <- setdiff(VariableFeatures(srat), exclude_genes)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(srat)
library(scDblFinder)

sce <- scDblFinder(
  sce,
  samples = sce$sample_id  # ❗非常关键
)
table(sce$scDblFinder.class)
sce <- sce[, sce$scDblFinder.class == "singlet"]
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat)
srat <- ScaleData(srat)
srat <- RunPCA(srat)
ElbowPlot(srat)
# 1. 提取 PCA
pca_embed <- Embeddings(srat, "pca")

# 2. metadata
meta <- srat@meta.data

# 3. 直接跑 harmony（核心函数，不走Seurat封装）
harmony_embed <- harmony::HarmonyMatrix(
  data_mat = pca_embed,
  meta_data = meta,
  vars_use = "sample_id",
  do_pca = FALSE
)

# 4. 放回 Seurat
srat[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embed,
  key = "harmony_"
)
srat <- RunUMAP(srat, reduction = "harmony", dims = 1:30,
                umap.method = "umap-learn", metric = "correlation")
srat <- FindNeighbors(srat, reduction = "harmony", dims = 1:30)
srat <- FindClusters(srat, resolution = 0.5)
# 查看可用的降维结果
Reductions(srat)

# 可视化聚类结果
DimPlot(srat,reduction = "umap",label = T,group.by = "seurat_clusters")+
  ggtitle("Merged Myeloid cells group")
# 3. 过滤feature
Idents(srat) <- "seurat_clusters"
all_genes_clean <- setdiff(rownames(srat), exclude_genes)
srat_markers <- FindAllMarkers(srat, 
                              features = all_genes_clean, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)
top10_srat <- srat_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
print(n=200,top10_srat[ c("cluster", "gene")])
DoHeatmap(srat,features = top10_srat$gene)+
  ggtitle("macrophages subgroup clusters top10 markers")
DimPlot(srat, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
FeaturePlot(srat,reduction = "umap",features = c("CD14","CD68","CSF1R","FCGR3A","ITGAM","TYROBP"),ncol = 3)
VlnPlot(srat,group.by = "seurat_clusters",features = c("CD14","CD68","CSF1R","FCGR3A","ITGAM","TYROBP"),ncol = 3)
# 计算每个簇的平均表达（仅限髓系相关基因）
myeloid_genes <- c("CD14", "C1QC", "C1QA", "C1QB", "CSF1R", "FCER1G", "CD68", "MS4A7", "LYVE1", "FOLR2")
avg_exp <- AverageExpression(srat, features = myeloid_genes, group.by = "seurat_clusters")$RNA
t(avg_exp)[,"12"]  # 查看簇12的平均表达
#查看可疑簇的top50
cluster12_markers <- FindMarkers(srat, ident.1 = 12, only.pos = TRUE, logfc.threshold = 0.1)
# 查看基因名（前50个）
head(rownames(cluster12_markers), 50)
# #由于ibm_ctrl_mac的cluster2的Y chromosome genes高表达，屏蔽后重新聚类一次
# srat <- NormalizeData(srat) %>%
#   FindVariableFeatures(nfeatures = 3000) %>%
#   ScaleData(features = features_use) %>%
#   RunPCA(features = features_use) %>%
#   RunUMAP(dims = 1:20) %>%
#   FindNeighbors(dims = 1:20) %>%
#   FindClusters(resolution = 0.4)
# DimPlot(srat,reduction = "umap",label = T,repel = T)
# features_use <- setdiff(VariableFeatures(srat), exclude_genes)
# Idents(srat) <- "seurat_clusters"
# srat_markers <- FindAllMarkers(srat, 
#                               only.pos = TRUE, 
#                               min.pct = 0.25, 
#                               logfc.threshold = 0.25)
# top10_srat <- srat_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
# print(n=200,top10_srat[ c("cluster", "gene")])
# new.cluster.ids <- c("NAMs-like mφ","TLR4+ mφ","pollution","pollution",
#                      "pollution","M2-like","pollution","Migrating mφ",
#                      "Resident mφ","LYVE1+ mφ","Plasma cell","FAPs",
#                      "cDC1","mDC","Proliferation mφ","pDC")
# 
# names(new.cluster.ids) <- levels(srat)
# srat <- RenameIdents(srat, new.cluster.ids)
# DimPlot(srat,reduction = "umap",label = T,repel = T)
# saveRDS(srat,"ibm_ctrl_mac.rds")
# 

#二、剔除污染
# #1.剔除污染簇
# #imnm
# #mac_clean <- subset(srat, idents = c(3, 4), invert = TRUE)
# #ibm#
# srat@active.ident<-srat$seurat_clusters
# bad_clusters <- c("0","2","3","4","6","10","11")
# clean<- subset(srat, idents = bad_clusters, invert = TRUE)
#asys
# bad_clusters <- c(1, 4, 7)   # 肌核/代谢污染、T细胞、B细胞
# mac_clean <- subset(srat, idents = bad_clusters, invert = TRUE)
#dm 
# bad_clusters <- c(1, 3, 6, 7)
# mac_clean <- subset(srat, idents = bad_clusters, invert = TRUE)
# cluster_anno <- c(
#   "0" = "Resident/Neuro-associated Mφ (LYVE1+)",
#   "2" = "Resident/AP Mφ (F13A1+ HLA+)",
#   "4" = "ISG-high (pDC-like)",
#   "5" = "SLC11A1+ Inflammatory Myeloid",
#   "8" = "pDC/Early DC"
# )
# names(cluster_anno) <- levels(mac_clean)
# mac_clean <- RenameIdents(mac_clean, cluster_anno)
# mac_clean$cell_type <- Idents(mac_clean)
# DimPlot(mac_clean,reduction = "umap")+
#   ggtitle("ASyS CTRL macrophages group")
srat@active.ident<-srat$seurat_clusters
bad_clusters <- c("1","11","16")
mac_clean <- subset(srat, idents = bad_clusters, invert = TRUE)
#2.重新定义 variable features
mac_clean<- FindVariableFeatures(mac_clean, nfeatures = 3000)
#重新降维、整合、聚类
mac_clean <- NormalizeData(mac_clean)
mac_clean <- FindVariableFeatures(mac_clean)
mac_clean <- ScaleData(mac_clean, features = VariableFeatures(mac_clean))
mac_clean <- RunPCA(mac_clean, npcs = 30)
ElbowPlot(mac_clean, ndims =30)
# 1. 提取 PCA
pca_embed <- Embeddings(mac_clean, "pca")

# 2. metadata
meta <- mac_clean@meta.data

# 3. 直接跑 harmony（核心函数，不走Seurat封装）
harmony_embed <- harmony::HarmonyMatrix(
  data_mat = pca_embed,
  meta_data = meta,
  vars_use = "sample_id",
  do_pca = FALSE
)

# 4. 放回 Seurat
mac_clean[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embed,
  key = "harmony_"
)
# mac_clean <- RunUMAP(mac_clean, reduction = "harmony", dims = 1:20,
#                 umap.method = "umap-learn", metric = "correlation")
# 使用 features 参数（不能同时指定 reduction 和 dims）
mac_clean <- RunUMAP(mac_clean,
                     reduction = "harmony",
                     dims=1:20,
                     umap.method = "umap-learn",
                     metric = "correlation")
mac_clean <- FindNeighbors(mac_clean, reduction = "harmony", dims = 1:20)
mac_clean <- FindClusters(mac_clean, resolution = 0.5)
DimPlot(mac_clean,reduction = "umap",group.by = "seurat_clusters",
        label = T)+
  ggtitle("Merged myeloids clean1")


features_use <- setdiff(
  VariableFeatures(mac_clean),
  exclude_genes   # muscle / neuronal / IG / MT
)
# #4.Neg_Score精确过滤
# # 假设您已定义了以下标记基因列表（请替换为您实际使用的列表）
# neuron_markers <- c("NOVA1", "ROBO2", "LRRTM4", "DSCAML1", "FGF13")
# macro_markers  <- c("CSF1R", "C1QA", "CD14", "F13A1", "LYVE1")
# muscle_markers <- c("DMD", "TTN", "MYH7", "TNNC2", "ACTC1")
# t_markers      <- c("CD3D", "BCL11B", "LCK", "SKAP1", "THEMIS")
# plasma_markers <- c("IGHG1", "IGKC", "JCHAIN", "IGLC2", "IGHGP")
# 
# # 过滤出在对象中实际存在的基因
# neuron_markers <- intersect(neuron_markers, rownames(clean))
# macro_markers  <- intersect(macro_markers,  rownames(clean))
# muscle_markers <- intersect(muscle_markers, rownames(clean))
# t_markers      <- intersect(t_markers,      rownames(clean))
# plasma_markers <- intersect(plasma_markers, rownames(clean))
# 
# # 检查过滤后各列表长度
# sapply(list(Neuron = neuron_markers, Macro = macro_markers, 
#             Muscle = muscle_markers, Tcell = t_markers, Plasma = plasma_markers), length)
# #安全计算
# # 仅当列表非空时才计算
# if (length(neuron_markers) > 0) {
#   clean <- AddModuleScore(clean, list(neuron_markers), name = "Neuron")
# }
# if (length(macro_markers) > 0) {
#   clean <- AddModuleScore(clean, list(macro_markers), name = "Macro")
# }
# if (length(muscle_markers) > 0) {
#   clean <- AddModuleScore(clean, list(muscle_markers), name = "Muscle")
# }
# if (length(t_markers) > 0) {
#   clean <- AddModuleScore(clean, list(t_markers), name = "Tcell")
# }
# if (length(plasma_markers) > 0) {
#   clean <- AddModuleScore(clean, list(plasma_markers), name = "Plasma")
# }
# # 获取 CSF1R 表达值作为强制保留依据
# csf1r_expr <- FetchData(clean, vars = "CSF1R")[,1]
# #设定阈值
# # 查看各评分的分布
# VlnPlot(clean, features = c("Neuron1", "Muscle1", "Tcell1", "Plasma1"), 
#         group.by = "seurat_clusters", pt.size = 0, ncol = 2)
# 
# # 统计各评分的分位数
# quantile(clean$Neuron1, probs = c(0.5, 0.75, 0.9, 0.95))
# quantile(clean$Muscle1, probs = c(0.5, 0.75, 0.9, 0.95,0.99))
# quantile(clean$Tcell1,  probs = c(0.5, 0.75, 0.9, 0.95))
# quantile(clean$Plasma1, probs = c(0.5, 0.75, 0.9, 0.95))
# #联合过滤
# # 注意：AddModuleScore 默认将评分列命名为 "Neuron1"、"Macro1" 等
# clean_2 <- subset(
#   clean,
#   subset = Neuron1 < 0.105 &
#     Muscle1 < 1.607 &
#     Tcell1  < 0.3 &
#     Plasma1 < 0.105
# )
# # # 重新降维聚类（使用屏蔽肌肉基因后的纯净特征）
# clean <- NormalizeData(clean) %>%
#   FindVariableFeatures(nfeatures = 3000) %>%
#   ScaleData(features = features_use) %>%
#   RunPCA(features = features_use) %>%
#   RunUMAP(dims = 1:20) %>%
#   FindNeighbors(dims = 1:20) %>%
#   FindClusters(resolution = 0.8)
# # mac_clean <- RunPCA(mac_clean)
# # mac_clean <- RunUMAP(mac_clean, dims = 1:30)
# # mac_clean <- FindNeighbors(mac_clean, dims = 1:30)
# # mac_clean <- FindClusters(mac_clean, resolution = 0.5)
DimPlot(mac_clean,reduction = "umap",label = T,repel = T)
Idents(mac_clean) <- "seurat_clusters"
mac_clean_markers <- FindAllMarkers(mac_clean, 
                                    features = features_use,
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
top10_clean <- mac_clean_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
print(n=200,top10_clean[ c("cluster", "gene")])
# 典型单核/巨噬细胞核心基因集（确认身份）
core_macrophage_genes <- c("CD14", "CD68", "CSF1R", "FCER1G", "C1QA", "C1QB", "C1QC",
                           "ITGAM", "ITGAX", "MRC1", "CD163", "FOLR2", "LYVE1")

# 炎症（M1-like）基因集
m1_genes <- c("IL1B", "IL6", "TNF", "CXCL10", "CXCL9", "IDO1", "NOS2", "FCGR1A", "CD80", "CD86")

# 抗炎/修复（M2-like）基因集
m2_genes <- c("ARG1", "CD163", "MRC1", "IL10", "TGFB1", "CD209", "PPARG", "CLEC10A", "VSIG4", "MS4A7")

# 溶酶体/吞噬功能基因集
phagocytosis_genes <- c("CD36", "SCARB2", "MSR1", "MERTK", "TREM2", "LAMP1", "LAMP2", "CTSK", "CTSD")
# 添加核心巨噬细胞模块评分
mac_clean <- AddModuleScore(mac_clean,
                                     features = list(core_macrophage_genes),
                                     name = "Core_Macro_Score")

# 添加 M1 评分
mac_clean <- AddModuleScore(mac_clean,
                                     features = list(m1_genes),
                                     name = "M1_Score")

# 添加 M2 评分
mac_clean <- AddModuleScore(mac_clean,
                                     features = list(m2_genes),
                                     name = "M2_Score")

# 添加吞噬评分
mac_clean <- AddModuleScore(mac_clean,
                                     features = list(phagocytosis_genes),
                                     name = "Phagocytosis_Score")
library(ggplot2)
library(patchwork)

# 小提琴图（按 cell_type）
p1 <- VlnPlot(mac_clean, features = "Core_Macro_Score1", group.by = "seurat_clusters") + 
  ggtitle("Core Macrophage Score")
p2 <- VlnPlot(mac_clean, features = "M1_Score1", group.by = "seurat_clusters") + 
  ggtitle("M1 Score")
p3 <- VlnPlot(mac_clean, features = "M2_Score1", group.by = "seurat_clusters") + 
  ggtitle("M2 Score")
p4 <- VlnPlot(mac_clean, features = "Phagocytosis_Score1", group.by = "seurat_clusters") + 
  ggtitle("Phagocytosis Score")

p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
DoHeatmap(mac_clean,features = top10_clean$gene)+
  ggtitle("Myeloids subgroup clusters top10 markers")
# #查看可疑簇top50
cluster12_markers <- FindMarkers(mac_clean, ident.1 = 12, only.pos = TRUE, logfc.threshold = 0.1)
head(rownames(cluster12_markers), 50)
# 注释
mac_clean$cell_type<-"Uknown"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(0)]<-"Pollution"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(1)]<-"resident"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(2)]<-"resident"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(3)]<-"Mono"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(4)]<-"Mono"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(5)]<-"Mono"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(6)]<-"Activated"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(7)]<-"Activated"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(8)]<-"Polltuion"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(9)]<-"DC"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(10)]<-"Activated"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(11)]<-"resident"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(12)]<-"DC"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(13)]<-"Pollution"
mac_clean$cell_type[mac_clean$seurat_clusters%in%c(14)]<-"Pollution"
DimPlot(mac_clean,reduction = "umap",group.by = "cell_type",label = T,
        repel = T)+ggtitle("Merged Myeloid cells groups")
saveRDS(mac_clean,"Purification_myeloids.rds")
# 或者使用 plyr::revalue 或 dplyr::recode

# # 查看同时表达巨噬细胞标志和 T 细胞标志的细胞
# FeatureScatter(mac_clean, feature1 = "CSF1R", feature2 = "LCK", 
#                    cells = WhichCells(mac_clean, idents = "2"))
#cluster2混杂了巨噬细胞、T细胞、以及doublet，下面进行簇内亚聚类
#3.剔除污染簇
# 剔除明确污染和非目标簇
# bad_clusters <- c(0, 1, 3, 5, 9, 10, 11, 12)   # 17 肥大细胞若需剔除可加入
# mac<- subset(clean, idents = bad_clusters, invert = TRUE)
# 


#第二轮剔除
# 剔除明确污染簇
# 1. 剔除明确污染簇和DC
bad_clusters <- c("0", "8", "9", "12", "13", "14")
clean2 <- subset(mac_clean, idents = bad_clusters, invert = TRUE)
#2.重新定义 variable features
clean2<- FindVariableFeatures(clean2, nfeatures = 3000)
#重新降维、整合、聚类
clean2 <- NormalizeData(clean2)
clean2 <- FindVariableFeatures(clean2)
clean2 <- ScaleData(clean2, features = VariableFeatures(clean2))
clean2 <- RunPCA(clean2, npcs = 30)
ElbowPlot(clean2)
# 1. 提取 PCA
pca_embed <- Embeddings(clean2, "pca")

# 2. metadata
meta <- clean2@meta.data

# 3. 直接跑 harmony（核心函数，不走Seurat封装）
harmony_embed <- harmony::HarmonyMatrix(
  data_mat = pca_embed,
  meta_data = meta,
  vars_use = "sample_id",
  do_pca = FALSE
)

# 4. 放回 Seurat
clean2[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embed,
  key = "harmony_"
)
clean2 <- RunUMAP(clean2, reduction = "harmony", dims = 1:15,
                     umap.method = "umap-learn", metric = "correlation")
clean2 <- FindNeighbors(clean2, reduction = "harmony", dims = 1:15)
clean2 <- FindClusters(clean2, resolution = 0.2)
DimPlot(clean2,reduction = "umap",group.by = "seurat_clusters",
        label = T)+
  ggtitle("Merged Macrophages group")


features_use <- setdiff(
  VariableFeatures(clean2),
  exclude_genes   # muscle / neuronal / IG / MT
)
Idents(clean2) <- "seurat_clusters"
clean2_markers <- FindAllMarkers(clean2, 
                                 features = features_use,
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.25)
top10 <- clean2_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
print(n=200,top10[ c("cluster", "gene")])
pdf("macrophage_heatmap.pdf", width = 12, height = 10)
DoHeatmap(clean2,features = top10$gene)+
  ggtitle("macrophages subgroup clusters top10 markers")
dev.off()
#第三轮剔除
# 1. 剔除明确污染簇和DC
bad_clusters <- c("9")
mac_final <- subset(clean2, idents = bad_clusters, invert = TRUE)
#2.重新定义 variable features
mac_final<- FindVariableFeatures(mac_final, nfeatures = 3000)
#重新降维、整合、聚类
mac_final <- NormalizeData(mac_final)
mac_final <- FindVariableFeatures(mac_final)
mac_final <- ScaleData(mac_final, features = VariableFeatures(mac_final))
mac_final <- RunPCA(mac_final, npcs = 30)
ElbowPlot(mac_final)
# 1. 提取 PCA
pca_embed <- Embeddings(mac_final, "pca")

# 2. metadata
meta <- mac_final@meta.data

# 3. 直接跑 harmony（核心函数，不走Seurat封装）
harmony_embed <- harmony::HarmonyMatrix(
  data_mat = pca_embed,
  meta_data = meta,
  vars_use = "sample_id",
  do_pca = FALSE
)

# 4. 放回 Seurat
mac_final[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embed,
  key = "harmony_"
)
mac_final <- RunUMAP(mac_final, reduction = "harmony", dims = 1:20,
                     umap.method = "umap-learn", metric = "correlation")
mac_final <- FindNeighbors(mac_final, reduction = "harmony", dims = 1:20)
mac_final <- FindClusters(mac_final, resolution = 0.6)
DimPlot(mac_final,reduction = "umap",group.by = "seurat_clusters",
        label = T)+
  ggtitle("Merged Macrophages group")


features_use <- setdiff(
  VariableFeatures(mac_final),
  exclude_genes   # muscle / neuronal / IG / MT
)
Idents(mac_final) <- "seurat_clusters"
mac_final_markers <- FindAllMarkers(mac_final, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.25)
top10 <- mac_final_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
print(n=200,top10[ c("cluster", "gene")])
DoHeatmap(mac_final,features = top10$gene)+
  ggtitle("macrophages subgroup clusters top10 markers")

 saveRDS(mac_final, file = "Purification_macrophages.rds")
cluster0_markers <- FindMarkers(mac_final, ident.1 = 0, only.pos = TRUE, logfc.threshold = 0.1)
head(rownames(cluster0_markers), 50)
# #2.重新定义 variable features
# mac_final<- FindVariableFeatures(mac_final, nfeatures = 3000)
# 
# features_use <- setdiff(
#   VariableFeatures(mac_final),
#   exclude_genes   # muscle / neuronal / IG / MT
# )
# #3.重新降维去除污染簇影响
# mac_final <- ScaleData(mac_final, features = features_use)
# mac_final <- RunPCA(mac_final, features = features_use)
# mac_final <- RunUMAP(mac_clean, dims = 1:20)
# mac_clean <- FindNeighbors(mac_clean, dims = 1:20)
# mac_clean <- FindClusters(mac_clean, resolution = 0.4)
# DimPlot(mac_clean,reduction = "umap",label = T,repel = T)+
#   ggtitle("IBM pure myeloids")
# Idents(mac_clean) <- "seurat_clusters"
# mac_clean_markers <- FindAllMarkers(mac_clean, 
#                                 only.pos = TRUE, 
#                                 min.pct = 0.25, 
#                                 logfc.threshold = 0.25)
# top10_mac_clean <- mac_clean_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
# print(n=200,top10_mac_clean[ c("cluster", "gene")])
# # #查看可疑簇top50
# cluster2_markers <- FindMarkers(mac_clean, ident.1 = 2, only.pos = TRUE, logfc.threshold = 0.1)
# print(head(cluster2_markers, 50))
# # 2. 定义簇注释向量（按剔除后保留簇的顺序）
# # 仅剔除簇6
# mac_clean <- subset(mac_clean, idents = 6, invert = TRUE)

# # 更新注释向量（包含簇2）
# cluster_anno <- c(
#   "2" = "ISG-high (IFN-I)",
#   "6" = "pDC/Early DC"
# )
# 
# names(cluster_anno) <- levels(mac_clean)
# mac_clean <- RenameIdents(mac_clean, cluster_anno)
# mac_clean$cell_type <- Idents(mac_clean)
# 
# # 3. 确保顺序一致
# DimPlot(mac_clean,reduction = "umap",pt.size=0.5,label = T,repel = T)+
#   ggtitle("DM Ctrl macrophages groups")
# DimPlot(mac_clean, reduction = "umap", group.by = "cell_type", 
#         split.by = "group", ncol = 2,pt.size=0.5) +
#   ggtitle("DM vs CTRL: Macrophage Subsets")
# # 设置顺序，使 CTRL 在前
# # 查看当前 antibody 列的水平
# table(mac_clean$antibody)
# 
# # 将 "Others" 替换为 "seronegative"
# mac_clean$antibody[mac_clean$antibody == "Others"] <- "seronegative"
# 
# # 验证修改结果
# table(mac_clean$antibody, mac_clean$group)
# mac_clean$antibody <- factor(mac_clean$antibody, levels = c("CTRL", "HMGCR", "SRP"))
# DimPlot(mac_clean, reduction = "umap", group.by = "cell_type", 
#         split.by = "antibody", ncol = 3, label = FALSE,pt.size = 0.5) +
#   ggtitle("Macrophage Subsets by Antibody Serotype")
#剔除5, 9, 14, 15
bad_clusters <- c( 5,9, 14, 15)
mac_final<-subset(mac_final,idents = bad_clusters, invert = TRUE)
mac_final$cell_type<-"Uknown"
mac_final$cell_type[mac_final$seurat_clusters%in%c(0)]<-"Mature Resident Mφ"
mac_final$cell_type[mac_final$seurat_clusters%in%c(1)]<-"FOLR2+ Resident Mφ"
mac_final$cell_type[mac_final$seurat_clusters%in%c(2)]<-"Intermediate Mφ"
mac_final$cell_type[mac_final$seurat_clusters%in%c(3)]<-"M1-like"
mac_final$cell_type[mac_final$seurat_clusters%in%c(4)]<-"PPARγ drived angiogenesis Mφ"
mac_final$cell_type[mac_final$seurat_clusters%in%c(5)]<-"cDC2"
mac_final$cell_type[mac_final$seurat_clusters%in%c(6)]<-"Patrolling Monocyte"
mac_final$cell_type[mac_final$seurat_clusters%in%c(7)]<-"Classical Monocyte"
mac_final$cell_type[mac_final$seurat_clusters%in%c(8)]<-"Mix Mφ"
mac_final$cell_type[mac_final$seurat_clusters%in%c(10)]<-"LYVE1+ Resident Mφ"
mac_final$cell_type[mac_final$seurat_clusters%in%c(11)]<-"LAMs/TREM2-like Mφ"
mac_final$cell_type[mac_final$seurat_clusters%in%c(12)]<-"MHC-II high Monocyte"
mac_final$cell_type[mac_final$seurat_clusters%in%c(13)]<-"Epigenetically Primed Monocyte"
mac_final$cell_type[mac_final$seurat_clusters%in%c(14)]<-"Antiviral Mφ"
mac_final$cell_type[mac_final$seurat_clusters%in%c(16)]<-"Early-activated Mφ"
#再剔除5
bad_clusters <- c(5)
mac_final<-subset(mac_final,idents = bad_clusters, invert = TRUE)
DimPlot(mac_final,reduction = "umap",group.by = "cell_type",label = T,
        repel = T)+ggtitle("Merged Monocyte-Macrophages cells groups")
#开始组间比较
cell_counts<-table(mac_final$cell_type, mac_final$group)
df_long <- as.data.frame(cell_counts)
colnames(df_long) <- c("cell_type", "group", "n")
# 查看转换结果
head(df_long)
# 列百分比（每组内细胞类型的比例）
df_col_pct <- df_long %>%
  group_by(group) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

# 行百分比（每种细胞类型在不同组中的分布比例）
df_row_pct <- df_long %>%
  group_by(cell_type) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()
#可视化
library(ggplot2)
#纵向百分比
ggplot(df_col_pct, aes(x = group, y = percent, fill = cell_type)) +
  geom_col(position = "stack") +
  labs(y = "Percentage (%)", title = "Cell type composition by group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")
#横向百分比
ggplot(df_row_pct, aes(x = cell_type, y = percent, fill = group)) +
  geom_col(position = "fill") +   # position = "fill" 使每个柱子的高度都为1，反映比例
  labs(y = "Proportion", x = "Cell type", 
       title = "Distribution of each cell type across groups") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +                   # 横向放置，避免标签重叠
  theme_minimal() +
  theme(legend.position = "bottom")

# 4. 保存最终对象
#saveRDS(mac_clean, file = "IMNM_CTRL_Macrophages_Cleaned.rds")
#saveRDS(mac_clean, file = "IBM_CTRL_Macrophages_Cleaned.rds")
#saveRDS(mac_clean, file = "ASyS_CTRL_Macrophages_Cleaned.rds")
# saveRDS(mac_clean, file = "DM_CTRL_Macrophages_Cleaned.rds")
saveRDS(mac_final, "Purification_macrophages.rds")
# 5.再次降维整合聚类，用于跨样本比较
#重新计算 PCA (使用你的屏蔽基因列表)
# 确保 features_use 不包含任何已剔除或污染基因
#round6
bad_clusters<-c(0,9,12,14)
mac_final<-subset(mac_final,idents = bad_clusters, invert = TRUE)
#循环
mac_final <- RunPCA(mac_final, features = features_use)
ElbowPlot(mac_final)
library(harmony)
# 1. 提取 PCA
pca_embed <- Embeddings(mac_final, "pca")

# 2. metadata
meta <- mac_final@meta.data

# 3. 直接跑 harmony（核心函数，不走Seurat封装）
harmony_embed <- harmony::HarmonyMatrix(
  data_mat = pca_embed,
  meta_data = meta,
  vars_use = "sample_id",
  do_pca = FALSE
)

# 4. 放回 Seurat
mac_final[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embed,
  key = "harmony_"
)
mac_final <- RunUMAP(mac_final, reduction = "harmony", dims = 1:15,
                     umap.method = "umap-learn", metric = "correlation")
mac_final <- FindNeighbors(mac_final, reduction = "harmony", dims = 1:15)
mac_final <- FindClusters(mac_final, resolution = 0.5)
DimPlot(mac_final,reduction = "umap",group.by = "seurat_clusters",
        label = T)+
  ggtitle("Merged Macrophages group")
Idents(mac_final) <- "seurat_clusters"
mac_final_markers <- FindAllMarkers(mac_final, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.25)
top10 <- mac_final_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
print(n=200,top10[ c("cluster", "gene")])
DoHeatmap(mac_final,features = top10$gene)+
  ggtitle("macrophages subgroup clusters top10 markers")

# round7 跳出循环
# 仅剔除 Cluster 11
mac_ready <- subset(mac_final, idents = "11", invert = TRUE)
DimPlot(mac_ready, reduction = "umap", group.by = "seurat_clusters",label = TRUE)+
  ggtitle("Macrophages ready for mrVI analysis")
# 直接将 mac_ready 保存为 RDS，进入下一步分析！
saveRDS(mac_ready, "macrophages_ready_for_miloR.rds")
