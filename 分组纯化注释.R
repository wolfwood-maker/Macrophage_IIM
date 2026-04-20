library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggplot2)
setwd("D:/R/R/PROJECT/IMNM/MergeMC")


#1.定义屏蔽marker
muscle_genes <- c(
  "MYH1","MYH2","MYH3","MYH7",
  "ACTA1","ACTN2",
  "TTN","NEB",
  "TNNT1","TNNT2","TNNT3",
  "TNNI1","TNNI2",
  "TPM1","TPM2","TPM3",
  "MYBPC1","MYBPC2",
  "MYL1","MYL2","MYL3","MYLPF",
  "CKM","DES","DMD","XIRP2","CMYA5","MYOT"
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
  rrna_genes
  # y_chromosome_genes,
  # x_genes
  #neural_genes,
  #fibroblast_genes,
  #smc_pericyte_genes,
  #endothelial_genes,
  #satellite_genes,
  #adipocyte_genes,
  #low_quality_genes
))
#二、分组处理
srat<-readRDS("dm_ctrl_mac.rds")#change rds name to switch dataset
DimPlot(srat,reduction = "umap",label = T,repel = T)
srat <- FindVariableFeatures(srat, nfeatures = 3000)

# 3. 过滤feature
features_use <- setdiff(VariableFeatures(srat), exclude_genes)
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

#查看可疑簇的top50
cluster0_markers <- FindMarkers(srat, ident.1 = 0, only.pos = TRUE, logfc.threshold = 0.1)
print(head(cluster0_markers, 50))
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
bad_clusters <- c(1, 3, 6, 7)
mac_clean <- subset(srat, idents = bad_clusters, invert = TRUE)
cluster_anno <- c(
  "0" = "Resident/Neuro-associated Mφ (LYVE1+)",
  "2" = "Resident/AP Mφ (F13A1+ HLA+)",
  "4" = "ISG-high (pDC-like)",
  "5" = "SLC11A1+ Inflammatory Myeloid",
  "8" = "pDC/Early DC"
)


names(cluster_anno) <- levels(mac_clean)
mac_clean <- RenameIdents(mac_clean, cluster_anno)
mac_clean$cell_type <- Idents(mac_clean)
DimPlot(mac_clean,reduction = "umap")+
  ggtitle("ASyS CTRL macrophages group")
#2.重新定义 variable features
mac_clean<- FindVariableFeatures(mac_clean, nfeatures = 3000)

features_use <- setdiff(
  VariableFeatures(mac_clean),
  exclude_genes   # muscle / neuronal / IG / MT
)
#3.重新降维去除污染簇影响
mac_clean <- ScaleData(mac_clean, features = features_use)
mac_clean <- RunPCA(mac_clean, features = features_use)
mac_clean <- RunUMAP(mac_clean, dims = 1:30)
mac_clean <- FindNeighbors(mac_clean, dims = 1:30)
mac_clean <- FindClusters(mac_clean, resolution = 0.4)


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
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
top10_clean <- mac_clean_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
print(n=200,top10_clean[ c("cluster", "gene")])
# #查看可疑簇top50
cluster5_markers <- FindMarkers(mac_clean, ident.1 = 5, only.pos = TRUE, logfc.threshold = 0.1)
print(head(cluster5_markers, 50))
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
# 1. 剔除明确污染簇（8）
bad_clusters <- c(0,1,3, 4,5)
mac_final <- subset(mac_clean, idents = bad_clusters, invert = TRUE)
# saveRDS(mac_clean, file = "ASyS_CTRL_Macrophages_Cleaned.rds")
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
# mac_final <- RunUMAP(mac_final, dims = 1:20)
# mac_final <- FindNeighbors(mac_final, dims = 1:20)
# mac_final <- FindClusters(mac_final, resolution = 0.4)
# DimPlot(mac_final,reduction = "umap",label = T,repel = T)+
#   ggtitle("IBM pure myeloids")
# Idents(mac_final) <- "seurat_clusters"
# mac_final_markers <- FindAllMarkers(mac_final, 
#                                 only.pos = TRUE, 
#                                 min.pct = 0.25, 
#                                 logfc.threshold = 0.25)
# top10_mac_final <- mac_final_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
# print(n=200,top10_mac_final[ c("cluster", "gene")])
# # #查看可疑簇top50
# cluster2_markers <- FindMarkers(mac_final, ident.1 = 2, only.pos = TRUE, logfc.threshold = 0.1)
# print(head(cluster2_markers, 50))
# # 2. 定义簇注释向量（按剔除后保留簇的顺序）
# # 仅剔除簇6
# mac_final <- subset(mac_final, idents = 6, invert = TRUE)

# 更新注释向量（包含簇2）
cluster_anno <- c(
  "2" = "ISG-high (IFN-I)",
  "6" = "pDC/Early DC"
)

names(cluster_anno) <- levels(mac_final)
mac_final <- RenameIdents(mac_final, cluster_anno)
mac_final$cell_type <- Idents(mac_final)

# 3. 确保顺序一致
DimPlot(mac_final,reduction = "umap",pt.size=0.5,label = T,repel = T)+
  ggtitle("DM Ctrl macrophages groups")
DimPlot(mac_final, reduction = "umap", group.by = "cell_type", 
        split.by = "group", ncol = 2,pt.size=0.5) +
  ggtitle("DM vs CTRL: Macrophage Subsets")
# 设置顺序，使 CTRL 在前
# 查看当前 antibody 列的水平
table(mac_final$antibody)

# 将 "Others" 替换为 "seronegative"
mac_final$antibody[mac_final$antibody == "Others"] <- "seronegative"

# 验证修改结果
table(mac_final$antibody, mac_final$group)
mac_final$antibody <- factor(mac_final$antibody, levels = c("CTRL", "HMGCR", "SRP"))
DimPlot(mac_final, reduction = "umap", group.by = "cell_type", 
        split.by = "antibody", ncol = 3, label = FALSE,pt.size = 0.5) +
  ggtitle("Macrophage Subsets by Antibody Serotype")
# 4. 保存最终对象
#saveRDS(mac_final, file = "IMNM_CTRL_Macrophages_Cleaned.rds")
#saveRDS(mac_final, file = "IBM_CTRL_Macrophages_Cleaned.rds")
#saveRDS(mac_final, file = "ASyS_CTRL_Macrophages_Cleaned.rds")
saveRDS(mac_final, file = "DM_CTRL_Macrophages_Cleaned.rds")
