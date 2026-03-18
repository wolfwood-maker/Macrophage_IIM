library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(harmony)
library(tidyr)
setwd("D:/R/R/PROJECT/IMNM/filterMC")
load("C:/Users/13067/Downloads/srat.RData")
#1. Basic QC
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
srat <- subset(
  srat,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    nCount_RNA > 500 &
    percent.mt < 5
)
#2. Normalization
srat <- NormalizeData(srat)

srat <- FindVariableFeatures(
  srat,
  selection.method = "vst",
  nfeatures = 2000
)

srat <- ScaleData(srat)

srat <- RunPCA(srat, npcs = 30)
saveRDS(srat,"Pipeline_srat.rds")
#3. Doublet delete
#关键修正：Seurat v5 必须执行 JoinLayers 才能兼容老插件
#srat <- readRDS("Pipeline_srat.rds")
# 假设细胞名是用 "_" 分隔的，第一部分是样本名
srat$orig.ident <- sapply(strsplit(colnames(srat), "_"), `[`, 1)

# 检查一下现在每个样本有多少细胞
table(srat$orig.ident)
# 确保层已合并
srat <- JoinLayers(srat)

# 获取所有独立的样本 ID（你刚才 table 出来的那些：ASY2, IBM1, IBM2...）
sample_ids <- unique(srat$orig.ident)

# 1. 预先初始化结果容器（在循环外）
srat$doublet_finder <- "Singlet"
sample_ids <- unique(srat$orig.ident)

for (sid in sample_ids) {
  cat("\n-----------------------------------\n")
  cat("Processing Sample:", sid, "\n")
  
  # 2. 【核心修复】不使用 subset，而是重新创建一个“干净”的 Seurat 对象
  # 这样可以彻底切断 v5 对象复杂的元数据关联
  tmp_counts <- GetAssayData(srat, layer = "counts")[, srat$orig.ident == sid]
  srat_grp <- CreateSeuratObject(counts = tmp_counts)
  
  # 3. 标准预处理（DoubletFinder 必须在局部空间运行）
  srat_grp <- NormalizeData(srat_grp, verbose = FALSE)
  srat_grp <- FindVariableFeatures(srat_grp, verbose = FALSE)
  srat_grp <- ScaleData(srat_grp, verbose = FALSE)
  srat_grp <- RunPCA(srat_grp, verbose = FALSE)
  
  # 4. 参数扫描 (pK)
  sweep.res <- paramSweep(srat_grp, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # 【手动提取 pK，防止数据框报错】
  pk_metric <- as.numeric(as.character(bcmvn$BCmetric))
  pk_list <- as.numeric(as.character(bcmvn$pK))
  pK_opt <- pk_list[which.max(pk_metric)]
  cat("Optimal pK found:", pK_opt, "\n")
  
  # 5. 计算预期双细胞数
  n_cells <- ncol(srat_grp)
  rate <- (n_cells / 1000) * 0.008
  nExp <- round(n_cells * rate)
  
  # 6. 运行 DoubletFinder
  srat_grp <- doubletFinder(srat_grp, PCs = 1:20, pN = 0.25, pK = pK_opt, nExp = nExp, sct = FALSE)
  
  # 7. 提取结果并写回原对象
  df_col <- grep("DF.classifications", colnames(srat_grp@meta.data), value = TRUE)[1]
  # 转化为纯字符向量，确保没有属性污染
  res_vector <- as.character(srat_grp@meta.data[[df_col]])
  names(res_vector) <- colnames(srat_grp)
  
  # 使用 match 安全赋值
  srat$doublet_finder[match(names(res_vector), colnames(srat))] <- res_vector
  
  cat("Sample", sid, "done. Found", sum(res_vector == "Doublet"), "doublets.\n")
}

# 最后检查
table(srat$doublet_finder)
#去除doublet
srat_clean<-subset(srat, subset = doublet_finder == "Singlet")
saveRDS(srat_clean, "FAPclean.rds")
srat<-readRDS("D:/R/R/PROJECT/IMNM/filterMC/FAPclean.rds")

# #4.再次降维聚类
# srat <- RunUMAP(srat, dims = 1:20)
# 
# srat <- FindNeighbors(srat, dims = 1:20)
# 
# srat <- FindClusters(srat, resolution = 0.5)
# 
# DimPlot(srat, reduction = "umap",label = TRUE)
# 
# #5.提取CD45+免疫细胞
# FeaturePlot(srat, features = "PTPRC",reduction = "umap")
# immune <- subset(srat, subset = PTPRC > 1)
# saveRDS(immune, "immune.rds")
# ###6.去除muscle contamination
# #定义markers
# ##muscle.genes <- c(
# # "MYH1","MYH2","MYH7",
# #"ACTA1","TTN","CKM","TNNT3"
# #)
# #计算muscle score
# ##immune <- AddModuleScore(
#   #immune,
#   #features = list(muscle.genes),
#   #name = "muscle_score"
# #)
# #FeaturePlot(immune, features = "muscle_score1",reduction = "umap")
# 
# #过滤高muscle score细胞
# #immune <- subset(
#  # immune,
#   #subset = muscle_score1 < 0.3
# #)
# #saveRDS(immune, "Immune_cells.rds")
# #8.harmony聚类
# immune<-readRDS("immune.rds")
# 
# srat <- RunPCA(srat, npcs = 30)
# 
# immune <- RunHarmony(
#   immune,
#   group.by.vars = "Sample"
# )
# #harmony聚类
# ElbowPlot(immune,reduction = "harmony")
# immune <- RunUMAP(
#   immune,
#   reduction = "harmony",
#   dims = 1:15
# )
# 
# immune <- FindNeighbors(
#   immune,
#   reduction = "harmony",
#   dims = 1:15
# )
# 
# immune <- FindClusters(
#   immune,
#   resolution = 0.5
# )
# 
# DimPlot(immune,reduction = "umap",label = T)
# saveRDS(immune,"immune_harmony.rds")
# immune<-readRDS("immune_harmony.rds")
# #9.提取巨噬细胞
# exclude_genes<- c("CD3D","CD19","MS4A1","IGHA1")
# FeaturePlot(immune,reduction = "umap",features = exclude_genes,ncol = 2)
# T_genes<-c("CD3D","SKAP1","THEMIS","ITK")
# FeaturePlot(immune,reduction = "umap",features = T_genes,ncol = 2)
# 
# endothelial_markers <- c("PECAM1", "VWF", "CLDN5")
# FeaturePlot(immune, reduction="umap",features = endothelial_markers, order = TRUE)
# 
# APC_genes <- c(
#   "HLA-DRA","HLA-DRB1","CD74","CIITA"
# )
# FeaturePlot(immune,reduction = "umap",features = APC_genes,ncol = 2)
# FeaturePlot(immune,reduction = "umap",features = c("ITGAM","ITGAX"),pt.size = 1.0,ncol = 2)
# FeaturePlot(immune,reduction = "umap",features = c("MKI67"),pt.size = 1.0,ncol = 1)
# #得到单核-巨噬细胞subgroup
# #分离MKI67+ cells
# Idents(immune) <- "seurat_clusters"
# mki67_cells<- subset(immune,
#                      idents=c(9))
# mki67_cells <- NormalizeData(mki67_cells)
# 
# mki67_cells <- FindVariableFeatures(mki67_cells)
# 
# mki67_cells <- ScaleData(mki67_cells)
# 
# mki67_cells <- RunPCA(mki67_cells)
# 
# ElbowPlot(mki67_cells,reduction = "pca")
# 
# mki67_cells <- RunUMAP(mki67_cells, dims = 1:10)
# 
# mki67_cells <- RunTSNE(mki67_cells, dims = 1:10)
# 
# mki67_cells <- FindNeighbors(mki67_cells, dims = 1:10)
# 
# mki67_cells <- FindClusters(mki67_cells, resolution = 0.5)
# 
# DimPlot(mki67_cells, label = TRUE, reduction="umap",pt.size = 1.0)
# 
# FeaturePlot(mki67_cells,reduction = "umap",features = T_genes)
# Idents(mki67_cells)<- "seurat_clusters"
# mki_mac <- subset(mki67_cells,idents=c(0))
# macs_mki<- colnames(mki_mac)
# #
# Idents(immune) <- "seurat_clusters"
# mon_mac <- subset(
#   immune,
#   idents = c(1,2,3,4,13)
# )
# cells_other<-colnames(mon_mac)
# mac_all<- c(macs_mki, cells_other)
# mon_mac<-subset(immune, cells= mac_all)
# saveRDS(mon_mac, "Macrophages_filter.rds")

#8.直接harmony
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(harmony)
library(tidyr)
setwd("D:/R/R/PROJECT/IMNM/filterMC")
load("C:/Users/13067/Downloads/srat.RData")
#1. Basic QC
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
srat <- subset(
  srat,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    nCount_RNA > 500 &
    percent.mt < 5
)
#2. Normalization
srat <- NormalizeData(srat)

srat <- FindVariableFeatures(
  srat,
  selection.method = "vst",
  nfeatures = 2000
)

srat <- ScaleData(srat)

srat <- RunPCA(srat, npcs = 30)
saveRDS(srat,"Pipeline_srat.rds")
#3. Doublet delete
#关键修正：Seurat v5 必须执行 JoinLayers 才能兼容老插件
#srat <- readRDS("Pipeline_srat.rds")
# 假设细胞名是用 "_" 分隔的，第一部分是样本名
srat$orig.ident <- sapply(strsplit(colnames(srat), "_"), `[`, 1)

# 检查一下现在每个样本有多少细胞
table(srat$orig.ident)
# 确保层已合并
srat <- JoinLayers(srat)

# 获取所有独立的样本 ID（你刚才 table 出来的那些：ASY2, IBM1, IBM2...）
sample_ids <- unique(srat$orig.ident)

# 1. 预先初始化结果容器（在循环外）
srat$doublet_finder <- "Singlet"
sample_ids <- unique(srat$orig.ident)

for (sid in sample_ids) {
  cat("\n-----------------------------------\n")
  cat("Processing Sample:", sid, "\n")
  
  # 2. 【核心修复】不使用 subset，而是重新创建一个“干净”的 Seurat 对象
  # 这样可以彻底切断 v5 对象复杂的元数据关联
  tmp_counts <- GetAssayData(srat, layer = "counts")[, srat$orig.ident == sid]
  srat_grp <- CreateSeuratObject(counts = tmp_counts)
  
  # 3. 标准预处理（DoubletFinder 必须在局部空间运行）
  srat_grp <- NormalizeData(srat_grp, verbose = FALSE)
  srat_grp <- FindVariableFeatures(srat_grp, verbose = FALSE)
  srat_grp <- ScaleData(srat_grp, verbose = FALSE)
  srat_grp <- RunPCA(srat_grp, verbose = FALSE)
  
  # 4. 参数扫描 (pK)
  sweep.res <- paramSweep(srat_grp, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # 【手动提取 pK，防止数据框报错】
  pk_metric <- as.numeric(as.character(bcmvn$BCmetric))
  pk_list <- as.numeric(as.character(bcmvn$pK))
  pK_opt <- pk_list[which.max(pk_metric)]
  cat("Optimal pK found:", pK_opt, "\n")
  
  # 5. 计算预期双细胞数
  n_cells <- ncol(srat_grp)
  rate <- (n_cells / 1000) * 0.008
  nExp <- round(n_cells * rate)
  
  # 6. 运行 DoubletFinder
  srat_grp <- doubletFinder(srat_grp, PCs = 1:20, pN = 0.25, pK = pK_opt, nExp = nExp, sct = FALSE)
  
  # 7. 提取结果并写回原对象
  df_col <- grep("DF.classifications", colnames(srat_grp@meta.data), value = TRUE)[1]
  # 转化为纯字符向量，确保没有属性污染
  res_vector <- as.character(srat_grp@meta.data[[df_col]])
  names(res_vector) <- colnames(srat_grp)
  
  # 使用 match 安全赋值
  srat$doublet_finder[match(names(res_vector), colnames(srat))] <- res_vector
  
  cat("Sample", sid, "done. Found", sum(res_vector == "Doublet"), "doublets.\n")
}

# 最后检查
table(srat$doublet_finder)
#去除doublet
srat_clean<-subset(srat, subset = doublet_finder == "Singlet")
saveRDS(srat_clean, "FAPclean.rds")
srat<-readRDS("D:/R/R/PROJECT/IMNM/filterMC/FAPclean.rds")
#4.harmony
srat <- RunPCA(srat, features = VariableFeatures(srat), npcs = 50)
ElbowPlot(srat, ndims = 50)  # 观察拐点位置
dims_use=1:20
srat <- RunHarmony(srat, group.by.vars = "orig.ident",
                   dims.use = dims_use, plot_convergence = FALSE)
#5.降维聚类可视化
# UMAP（基于 harmony）
srat <- RunUMAP(srat, reduction = "harmony", dims = dims_use, n.neighbors = 30, min.dist = 0.3)

# 构建 KNN 图并聚类
srat <- FindNeighbors(srat, reduction = "harmony", dims = dims_use)
srat <- FindClusters(srat, resolution = 0.5)  # 可调整分辨率
#可视化
DimPlot(srat, reduction = "umap", 
        group.by = "orig.ident", label = FALSE) + ggtitle("After Harmony")
DimPlot(srat, reduction = "umap", 
        group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
#初步注释
srat.markers<-FindAllMarkers(srat,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
srat.markers%>%group_by(cluster) %>% top_n(n=2,wt=avg_log2FC)
top10<-srat.markers%>%group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
top10_table <- top10 %>%
  filter(cluster %in% 0:15) %>%
  group_by(cluster) %>%
  mutate(rank=row_number()) %>%
  pivot_wider(names_from = cluster,
              values_from = gene)
write.csv(top10_table, "cluster_top10_table.csv", row.names=FALSE)
DoHeatmap(srat,features = top10$gene)+NoLegend()
srat$cell_type<-"Uknown"
srat$cell_type[srat$seurat_clusters%in%c(0)]<-"Type II muscle"
srat$cell_type[srat$seurat_clusters%in%c(1)]<-"Type I muscle"
srat$cell_type[srat$seurat_clusters%in%c(2)]<-"Hybrid muscle"
srat$cell_type[srat$seurat_clusters%in%c(3)]<-"Muscle progenitor"
srat$cell_type[srat$seurat_clusters%in%c(4)]<-"Activated satellite"
srat$cell_type[srat$seurat_clusters%in%c(5)]<-"FAP"
srat$cell_type[srat$seurat_clusters%in%c(6)]<-"M2-like"
srat$cell_type[srat$seurat_clusters%in%c(7)]<-"T cells"
srat$cell_type[srat$seurat_clusters%in%c(8)]<-"Satellite cells"
srat$cell_type[srat$seurat_clusters%in%c(9)]<-"Endothelial"
srat$cell_type[srat$seurat_clusters%in%c(10)]<-"Pericyte/smooth muscle"
srat$cell_type[srat$seurat_clusters%in%c(11)]<-"Cycling cells"
srat$cell_type[srat$seurat_clusters%in%c(12)]<-"M1-like"
srat$cell_type[srat$seurat_clusters%in%c(13)]<-"Unclear"
srat$cell_type[srat$seurat_clusters%in%c(14)]<-"Fibroblast"
srat$cell_type[srat$seurat_clusters%in%c(15)]<-"Adipocyte"
DimPlot(srat, reduction = "umap", group.by = "cell_type",label = T, repel = T)
#6.筛选myeloid

myeloid<-subset(srat, idents=c(6,11,12))

 #positive scale
myeloid.core <- c(
  "C1QA","C1QB","C1QC",
  "TYROBP","LST1",
  "LYZ","AIF1"
)
mac.markers <- c(
  "CD68","MS4A7",
  "SPP1","GPNMB","APOE",
  "MARCO","MRC1",
  "FCGR3A","ITGAM"
)
dc.markers <- c(
  "ITGAX",        # CD11c
  "CD1C","CLEC10A","FCER1A",   # cDC2
  "XCR1","CLEC9A","BATF3",     # cDC1
  "LAMP3","CCR7","MARCH1"      # migratory/mature DC
)
 #negative scale
t.markers <- c("CD3D","CD3E","TRAC","CD2","BCL11B")
b.markers <- c(
  "MS4A1","CD79A","CD79B",
  "IGHG1","IGHA1","IGHM","JCHAIN"
)
nk.markers <- c("NKG7","GNLY","KLRD1","PRF1")
muscle.markers <- c(
  "TTN","MYH1","MYH2","ACTA1",
  "CKM","DES","NEB"
)
fibro.markers <- c("FAP","COL1A1","COL3A1","DCN","PDGFRA")
endo.markers <- c("PECAM1","VWF","KDR","CLDN5")
saveRDS(myeloid,"myeloid.rds")
#7.降维聚类流程
# 1. 确保当前活跃 ident 为聚类（可选，但方便检查）
Idents(myeloid) <- "seurat_clusters"

# 2. 重新归一化（如果尚未进行，或想基于子集重新计算）
myeloid <- NormalizeData(myeloid)

# 3. 重新识别高变基因（子集的高变基因可能与总数据集不同）
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000)

# 4. 数据缩放（回归技术变量可选）
myeloid <- ScaleData(myeloid, vars.to.regress = c("percent.mt", "nCount_RNA"))  # 根据您的 QC 变量调整

# 5. PCA 降维
myeloid <- RunPCA(myeloid, features = VariableFeatures(myeloid), npcs = 50)

# 6. 确定使用的 PCs（可运行 ElbowPlot 辅助选择）
ElbowPlot(myeloid, ndims = 50)
# 假设选择前 20 个 PCs
dims_use <- 1:20

# 7. UMAP 降维
myeloid <- RunUMAP(myeloid, dims = dims_use)

# 8. 聚类（构建 KNN 图并寻找聚类）
myeloid <- FindNeighbors(myeloid, dims = dims_use)
myeloid <- FindClusters(myeloid, resolution = 0.5)  # 分辨率可根据需要调整

# 9. 可视化检查
DimPlot(myeloid, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Myeloid subclusters")

# 10. 保存新对象
saveRDS(myeloid, file = "myeloid_processed.rds")

#8.myeloid annotation
mye.markers<-FindAllMarkers(myeloid,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
mye.markers%>%group_by(cluster) %>% top_n(n=2,wt=avg_log2FC)
top10<-mye.markers%>%group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
top10_table <- top10 %>%
  filter(cluster %in% 0:12) %>%
  group_by(cluster) %>%
  mutate(rank=row_number()) %>%
  pivot_wider(names_from = cluster,
              values_from = gene)
write.csv(top10_table, "cluster_top10_table.csv", row.names=FALSE)
DoHeatmap(myeloid,features = top10$gene)+NoLegend()
 #去除污染
#first round

myeloid<-subset(myeloid,idents=c(0,1,3,4,5,6,11))
#second round
myeloid_clean <- subset(myeloid, idents = c(4,7), invert = TRUE)
tcell_cells <- WhichCells(myeloid_clean, expression = TRAC > 0 | CD3D > 0 | CD2 > 0)
myeloid_clean <- subset(myeloid_clean, cells = tcell_cells, invert = TRUE)
myeloid<-myeloid_clean
myeloid$cell_type<-"Uknown"
myeloid$cell_type[myeloid$seurat_clusters%in%c(0)]<-"LYVE1+ Resident Macrophages"
myeloid$cell_type[myeloid$seurat_clusters%in%c(1)]<-"GPNMB+ Remodeling Macs"
myeloid$cell_type[myeloid$seurat_clusters%in%c(2)]<-"Monocyte-derived Macrophages"
myeloid$cell_type[myeloid$seurat_clusters%in%c(3)]<-"Stroma-associated Macrophages"
myeloid$cell_type[myeloid$seurat_clusters%in%c(4)]<-"THBS1+ Stress-responsive Macs"
myeloid$cell_type[myeloid$seurat_clusters%in%c(5)]<-"Lipid-associated Macrophages (LAMs)"
myeloid$cell_type[myeloid$seurat_clusters%in%c(6)]<-"Cycling Myeloid cells"
myeloid$cell_type[myeloid$seurat_clusters%in%c(7)]<-"cDC1"
myeloid$cell_type[myeloid$seurat_clusters%in%c(8)]<-"cDC2"
myeloid$cell_type[myeloid$seurat_clusters%in%c(9)]<-"Mast cells"
DimPlot(myeloid, reduction = "umap",group.by = "cell_type", 
        pt.size=1.0,label = T, repel = T)
saveRDS(myeloid,"Myeloid_annotation.rds")
# #4.再次降维聚类
# srat <- RunUMAP(srat, dims = 1:20)
# 
# srat <- FindNeighbors(srat, dims = 1:20)
# 
# srat <- FindClusters(srat, resolution = 0.5)
# 
# DimPlot(srat, reduction = "umap",label = TRUE)

# #5.提取CD45+免疫细胞
# FeaturePlot(srat, features = "PTPRC",reduction = "umap")
# immune <- subset(srat, subset = PTPRC > 1)
# saveRDS(immune, "immune.rds")
# ###6.去除muscle contamination
# #定义markers
# ##muscle.genes <- c(
# # "MYH1","MYH2","MYH7",
# #"ACTA1","TTN","CKM","TNNT3"
# #)
# #计算muscle score
# ##immune <- AddModuleScore(
#   #immune,
#   #features = list(muscle.genes),
#   #name = "muscle_score"
# #)
# #FeaturePlot(immune, features = "muscle_score1",reduction = "umap")
# 
# #过滤高muscle score细胞
# #immune <- subset(
#  # immune,
#   #subset = muscle_score1 < 0.3
# #)
# #saveRDS(immune, "Immune_cells.rds")
# #8.harmony聚类
# immune<-readRDS("immune.rds")
# 
# srat <- RunPCA(srat, npcs = 30)
# 
# immune <- RunHarmony(
#   immune,
#   group.by.vars = "Sample"
# )
# #harmony聚类
# ElbowPlot(immune,reduction = "harmony")
# immune <- RunUMAP(
#   immune,
#   reduction = "harmony",
#   dims = 1:15
# )
# 
# immune <- FindNeighbors(
#   immune,
#   reduction = "harmony",
#   dims = 1:15
# )
# 
# immune <- FindClusters(
#   immune,
#   resolution = 0.5
# )
# 
# DimPlot(immune,reduction = "umap",label = T)
# saveRDS(immune,"immune_harmony.rds")
# immune<-readRDS("immune_harmony.rds")
# #9.提取巨噬细胞
# exclude_genes<- c("CD3D","CD19","MS4A1","IGHA1")
# FeaturePlot(immune,reduction = "umap",features = exclude_genes,ncol = 2)
# T_genes<-c("CD3D","SKAP1","THEMIS","ITK")
# FeaturePlot(immune,reduction = "umap",features = T_genes,ncol = 2)
# 
# endothelial_markers <- c("PECAM1", "VWF", "CLDN5")
# FeaturePlot(immune, reduction="umap",features = endothelial_markers, order = TRUE)
# 
# APC_genes <- c(
#   "HLA-DRA","HLA-DRB1","CD74","CIITA"
# )
# FeaturePlot(immune,reduction = "umap",features = APC_genes,ncol = 2)
# FeaturePlot(immune,reduction = "umap",features = c("ITGAM","ITGAX"),pt.size = 1.0,ncol = 2)
# FeaturePlot(immune,reduction = "umap",features = c("MKI67"),pt.size = 1.0,ncol = 1)
# #得到单核-巨噬细胞subgroup
# #分离MKI67+ cells
# Idents(immune) <- "seurat_clusters"
# mki67_cells<- subset(immune,
#                      idents=c(9))
# mki67_cells <- NormalizeData(mki67_cells)
# 
# mki67_cells <- FindVariableFeatures(mki67_cells)
# 
# mki67_cells <- ScaleData(mki67_cells)
# 
# mki67_cells <- RunPCA(mki67_cells)
# 
# ElbowPlot(mki67_cells,reduction = "pca")
# 
# mki67_cells <- RunUMAP(mki67_cells, dims = 1:10)
# 
# mki67_cells <- RunTSNE(mki67_cells, dims = 1:10)
# 
# mki67_cells <- FindNeighbors(mki67_cells, dims = 1:10)
# 
# mki67_cells <- FindClusters(mki67_cells, resolution = 0.5)
# 
# DimPlot(mki67_cells, label = TRUE, reduction="umap",pt.size = 1.0)
# 
# FeaturePlot(mki67_cells,reduction = "umap",features = T_genes)
# Idents(mki67_cells)<- "seurat_clusters"
# mki_mac <- subset(mki67_cells,idents=c(0))
# macs_mki<- colnames(mki_mac)
# #
# Idents(immune) <- "seurat_clusters"
# mon_mac <- subset(
#   immune,
#   idents = c(1,2,3,4,13)
# )
# cells_other<-colnames(mon_mac)
# mac_all<- c(macs_mki, cells_other)
# mon_mac<-subset(immune, cells= mac_all)
# saveRDS(mon_mac, "Macrophages_filter.rds")
# #10.巨噬细胞亚组分析聚类
# mac<-readRDS("D:/R/R/PROJECT/IMNM/filterMC/Macrophages_filter.rds")
# mac <- NormalizeData(mac)
# 
# mac <- FindVariableFeatures(mac)
# 
# mac <- ScaleData(mac)
# 
# mac <- RunPCA(mac)
# 
# ElbowPlot(mac,reduction = "pca")
# 
# mac <- RunUMAP(mac, dims = 1:17)
# 
# mac <- RunTSNE(mac, dims = 1:17)
# 
# mac <- FindNeighbors(mac, dims = 1:17)
# 
# mac <- FindClusters(mac, resolution = 0.5)
# 
# DimPlot(mac, label = TRUE, reduction="umap",pt.size = 1.0)
# DimPlot(mac, label = TRUE, reduction = "tsne", pt.size = 1.0)
# #去掉污染细胞
# Idents(mac)<- "seurat_clusters"
# mac<-subset(mac,idents=c(0,1,2,3,4,6,9))
# saveRDS(mac,"Macrophages_filter.rds")
# mac<-readRDS("Macrophages_filter.rds")
# mac <- NormalizeData(mac)
# 
# mac <- FindVariableFeatures(mac)
# 
# mac <- ScaleData(mac)
# 
# mac <- RunPCA(mac)
# 
# ElbowPlot(mac,reduction = "pca")
# 
# mac <- RunUMAP(mac, dims = 1:17)
# 
# mac <- RunTSNE(mac, dims = 1:17)
# 
# mac <- FindNeighbors(mac, dims = 1:17)
# 
# mac <- FindClusters(mac, resolution = 0.5)
# 
# DimPlot(mac, label = TRUE, reduction="umap",pt.size = 1.0)
# #11. Cells annotation
# mac.markers<-FindAllMarkers(mac,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
# mac.markers%>%group_by(cluster) %>% top_n(n=2,wt=avg_log2FC)
# top10<-mac.markers%>%group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
# top10_table <- top10 %>%
#   filter(cluster %in% 0:7) %>%
#   group_by(cluster) %>%
#   mutate(rank=row_number()) %>%
#   pivot_wider(names_from = cluster,
#   values_from = gene)
# write.csv(top10_table, "cluster_top10_table.csv", row.names=FALSE)
# DoHeatmap(mac,features = top10$gene)+NoLegend()
#病种比较
DimPlot(myeloid, 
        reduction = "umap", 
        group.by = "cell_type", 
        split.by = "Group", 
        label = F, 
        repel = TRUE,
        pt.size = 0.6,
        ncol = 2)+
  ggtitle("IIM groups")
#血清型比较
table(myeloid$Group, myeloid$Antibody, useNA = "ifany")
DimPlot(myeloid, 
        reduction = "umap", 
        group.by = "cell_type", 
        split.by = "Antibody", 
        label = F, 
        repel = TRUE,
        ncol = 3)+
  ggtitle("serotype")
#IMNM内部比较
imnm_cells <- subset(myeloid, subset = Group == "IMNM")
DimPlot(imnm_cells, 
        repel = T,
        label = F,
        reduction = "umap", 
        group.by = "cell_type", 
        split.by = "Antibody",
        ncol = 2) +
  ggtitle("IMNM serotype")
#IBM内部比较
ibm_cells <- subset(myeloid, subset = Group == "IBM")
DimPlot(ibm_cells, 
        repel = T,
        label = F,
        reduction = "umap", 
        group.by = "cell_type", 
        split.by = "Antibody",
        ncol = 2) +
  ggtitle("IBM serotype")
# #dendritic cells
# dendritic.genes<-c(
#   "ZBTB46",
#   "CLEC9A",
#   "CD1C",
#   "FCER1A",
#   "CLEC10A",
#   "ITGAX"
# )
# FeaturePlot(mac,reduction = "umap",pt.size = 1.0,
#             features = dendritic.genes, ncol=3)
# #resident_macrophages
# FeaturePlot(mac, reduction = "umap", pt.size = 1.0, 
#             features = c("C1QA","C1QB","C1QC",
#                          "FOLR2","LYVE1","CD163"),
#             ncol = 3)
# #monocyte_derived_macrophages
# FeaturePlot(mac, reduction = "umap", pt.size = 1.0, 
#             features = c("S100A8","S100A9","VCAN",
#                          "LST1","FCN1","CCR2"),
#             ncol = 3)
# #Inflammatory macrophage
# FeaturePlot(mac,features = c("IL1B","NLRP3","CXCL8",
#                              "CCL3","CCL4","PTGS2"),
#             reduction = "umap", pt.size = 1.0,ncol = 3)
# #Interferon-activated macrophage
# FeaturePlot(mac,features = c("MX1","OAS1","OAS2",
#                              "STAT1","CXCL10","GBP1"),
#             reduction = "umap", pt.size = 1.0,ncol = 3)
# #Antigen-presenting macrophage                             
# mhc_2_genes<-c("HLA-DRA","HLA-DRB1","HLA-DPA1",
#                "HLA-DPB1","CD74","CIITA")
# FeaturePlot(mac,features = mhc_2_genes,
#             reduction = "umap", pt.size = 1.0,ncol = 3)
# #比较病种间差异
# DotPlot(mac, 
#         features = mhc_2_genes, 
#         group.by = "Group") + 
#   RotatedAxis() + # 旋转X轴标签，避免重叠
#   scale_color_gradient2(low = "blue", mid = "white", high = "red") + # 自定义颜色
#   labs(x = "Gene", y = "Disease Group", title = "MHC-II gene expression across IIM subtypes")
# #比较血清型间MHC-II
# DotPlot(mac, 
#         features = mhc_2_genes, 
#         group.by = "Antibody") + 
#   RotatedAxis() + # 旋转X轴标签，避免重叠
#   scale_color_gradient2(low = "blue", mid = "white", high = "red") + # 自定义颜色
#   labs(x = "Gene", y = "Antibody Group", title = "MHC-II gene expression across Antibody subtypes")
# # 查看各组的细胞数
# table(mac$Antibody)
# #Phagocytic / debris-clearing macrophage
# FeaturePlot(mac,features = c("GPNMB","SPP1","TREM2",
#                              "APOE","CTSB","LGALS3"),
#             reduction = "umap", pt.size = 1.0,ncol = 3)
# #Pro-repair macrophage
# FeaturePlot(mac,features = c("MRC1","CD163","IGF1",
#                              "TGFB1","SPP1"),
#             reduction = "umap", pt.size = 1.0,ncol = 3)
# #注释annotation
# mac$cell_type<-"Uknown"
# mac$cell_type[mac$seurat_clusters%in%c(0)]<-"LYVE1+ Res Mac"
# mac$cell_type[mac$seurat_clusters%in%c(1)]<-"SPP1+ Mac"
# mac$cell_type[mac$seurat_clusters%in%c(2)]<-"Monocyte-derived Mac"
# mac$cell_type[mac$seurat_clusters%in%c(3)]<-"cDC2"
# mac$cell_type[mac$seurat_clusters%in%c(4)]<-"LAMs"
# mac$cell_type[mac$seurat_clusters%in%c(5)]<-"Cycling myeloid cells"
# mac$cell_type[mac$seurat_clusters%in%c(6)]<-"cDC1"
# mac$cell_type[mac$seurat_clusters%in%c(7)]<-"Interstital Cells"
# DimPlot(mac, reduction = "umap", group.by = "cell_type",label = T, repel = T)

#12.Go analysis
# 按 cluster 拆分基因列表
gene_list<-split(mye.markers$gene,mye.markers$cluster)
# 将每个列表中的 SYMBOL 转换为 ENTREZ ID
gene_entrez <- lapply(gene_list, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
})
# 注意：compareCluster 要求输入是命名列表（cluster 名为元素名）
comp_ego <- compareCluster(geneCluster = gene_entrez,
                           fun = enrichGO,
                           OrgDb = org.Hs.eg.db,
                           ont = "BP",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           readable = TRUE)

# 查看结果
head(comp_ego)
go_simplified <- simplify(comp_ego, 
                          cutoff = 0.7, 
                          by = "p.adjust", 
                          select_fun = min)
# 保存结果（可选）
saveRDS(comp_ego, "D:/R/R/PROJECT/IMNM/filterMC/mac_GO_BP.rds")
saveRDS(go_simplified,"D:/R/R/PROJECT/IMNM/filterMC/Sim_GO_BP.rds")
comp_ego<- readRDS("D:/R/R/PROJECT/IMNM/filterMC/mac_GO_BP.rds")
go_simplified<- readRDS("D:/R/R/PROJECT/IMNM/filterMC/Sim_GO_BP.rds")
# 必须先计算相似度
go_results_sim <- pairwise_termsim(comp_ego)
go_results_reduced <- simplify(go_results_sim, 
                               cutoff = 0.5, 
                               by = "p.adjust", 
                               select_fun = min)
go_df <- as.data.frame(go_results_reduced)
go_top5 <- go_df %>%
  group_by(Cluster) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  ungroup()

# 检查一下现在剩下多少个，如果剩下 20-30 个，图就会非常漂亮
nrow(go_results_reduced@compareClusterResult)
#可视化 GO 结果
#dotplot
p<-dotplot(go_simplified, showCategory = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
library(stringr)
p + scale_y_discrete(labels = function(x) str_wrap(x, width = 200))
#简化到top5
dotplot(go_results_reduced, showCategory = NULL) + 
  aes(y = factor(Description, levels = rev(unique(go_top5$Description)))) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) + # 自动折行
  theme(axis.text.y = element_text(size = 8))

# cnetplot
# 只选前 30 个最显著的画网络图，避免线条过密
# 1. 确保计算了相似度（如果之前跑过可以跳过）
gene_list <- mac.markers$avg_log2FC
names(gene_list) <- rownames(mac.markers)

cnetplot(go_results_sim, 
         showCategory = 5, 
         foldChange = gene_list) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)



#heatmap
# 提取结果表格
result_df <- comp_ego@compareClusterResult

# 选择需要的列：Cluster, Description, p.adjust
df_sub <- result_df[, c("Cluster", "Description", "p.adjust")]

# 计算 -log10(p.adjust)，便于热图着色
df_sub$logP <- -log10(df_sub$p.adjust)

# 转换为宽格式：每一行是一个通路（Description），每一列是一个 Cluster
library(tidyr)
wide_df <- df_sub %>%
  dplyr::select(Cluster, Description, logP) %>%
  pivot_wider(names_from = Cluster, values_from = logP, values_fill = 0)

# 将通路名作为行名，删除通路名列
mat <- as.matrix(wide_df[, -1])
rownames(mat) <- wide_df$Description

# 检查矩阵维度
dim(mat)
# 提取每个 cluster 的 top 10 通路（按 p.adjust 升序）
top10_terms <- result_df %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  pull(Description) %>%
  unique()   # 去重

# 用这些通路筛选矩阵
mat_filt <- mat[rownames(mat) %in% top10_terms, ]

# 查看筛选后的维度
dim(mat_filt)
saveRDS(mat_filt, "mat_filt.rds")
#pheatmap
library(pheatmap)
mat_filt<- readRDS("mat_filt.rds")
pheatmap(mat_filt,
         color = colorRampPalette(c("grey", "yellow", "firebrick3"))(50),
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 7,
         fontsize_col = 10,
         main = "GO enrichment (top 10 per cluster)",
         filename = "GO_heatmap.pdf",
         width = 9,
         height = 10)


