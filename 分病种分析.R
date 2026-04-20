library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggplot2)
setwd("D:/R/R/PROJECT/IMNM/MergeMC")
#一、瘦身
# 1. 读取下载好的对象
srat<-readRDS("D:/snRNA-seq data to share/dm_final.rds")
# 2. 检查当前内存占用 (你会发现它很大)
format(object.size(srat), units = "Gb")

# 3. 核心步骤：删掉 scale.data 层
# 这不会影响你已经算好的 UMAP 和 Harmony 结果
# 这个函数会保留必要的 data 和 counts，但切掉 scale.data
srat <- DietSeurat(srat, 
                   counts = TRUE, 
                   data = TRUE, 
                   scale.data = FALSE, # 删掉这个大的
                   dimreducs = c("pca", "harmony", "umap"))

# 运行完后检查内存，你会发现它真的瘦了
gc()
format(object.size(srat), units = "Gb")

#二、检查效果
DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("DM+CTRL group")
# 看看不同 sample_id 是不是均匀分布在各个群里
DimPlot(srat, reduction = "umap", group.by = "sample_id")

#三、寻找高变基因
# 基于整合后的 harmony 空间构建近邻图
srat <- FindNeighbors(srat, reduction = "harmony", dims = 1:30)

# 寻找聚类
srat <- FindClusters(srat, resolution = 0.5)
# 确保 Idents 设置为聚类标签
Idents(srat) <- "seurat_clusters"

# 运行全群 Marker 寻找
# parallel 包可以加速，但本地建议先不加，稳字当头
all_markers <- FindAllMarkers(srat, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

# 保存结果，防止 RStudio 意外崩溃
#write.csv(all_markers, "D:/snRNA-seq data to share/all_clusters_markers.csv")
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
print(n=210,top10_markers[c("cluster", "gene")])

#注释
new.cluster.ids <- c("RBC", "SMC/Neural", "Slow Muscle","FAPs/Neural",
                     "Endothelial", "Regenerative Muscle", "Resident mφ","noise",
                     "Satellite cells", "Pericyte/SMC", "ISG-high","noise",
                     "T cells", "Lymphatic Endothelial", "noise","Endothelial/Mesothelial",
                     "MT pollution", "SMC", "Neural-assosiated","Adipocyte",
                     "Muscle")

names(new.cluster.ids) <- levels(srat)
srat <- RenameIdents(srat, new.cluster.ids)
saveRDS(srat,"E:/IIM/dm_ctrl.rds")
# 画一张带名字的 UMAP 图
DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE)+
  ggtitle("DM+CTRL Group (87581 cells)")

#四、巨噬细胞亚组分析
#1.分离亚组
mac<-subset(srat, idents=c("Resident mφ"))
#2.定义屏蔽marker
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
rbc_genes <- c(
  "HBB","HBA1","HBA2","HBQ1",
  "ALAS2","AHSP"
)
ig_genes <- grep("^IG[HKL][VDJC]", rownames(mac), value = TRUE)
mt_genes <- grep("^MT-", rownames(mac), value = TRUE)
rrna_genes <- c("MT-RNR1","MT-RNR2")
neural_genes <- c(
  "NOVA1","ROBO2","LRRTM4",
  "BDNF","DLGAP2","ELFN1"
)
y_chromosome_genes <- c(
  "DDX3Y", "PRKY", "USP9Y", "UTY", "TTTY14", "TTTY15", "NLGN4Y"
)
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
  neural_genes,
  fibroblast_genes,
  smc_pericyte_genes,
  endothelial_genes,
  satellite_genes,
  adipocyte_genes,
  low_quality_genes
))
# 3. variable genes
mac <- FindVariableFeatures(mac, nfeatures = 3000)

# 4. 过滤feature
features_use <- setdiff(VariableFeatures(mac), exclude_genes)

# 5. 降维
mac <- ScaleData(mac, features = features_use)
mac <- RunPCA(mac, features = features_use)
ElbowPlot(mac)
# 5. 聚类
mac <- RunUMAP(mac, dims = 1:10)
mac <- FindNeighbors(mac, dims = 1:10)
mac <- FindClusters(mac, resolution = 0.3)
# 6. 可视化检查
DimPlot(mac, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("Macrophages subclusters")
# 7.高变基因
Idents(mac) <- "seurat_clusters"
all_genes_clean <- setdiff(rownames(mac), exclude_genes)
mac_markers <- FindAllMarkers(mac, 
                              features = all_genes_clean, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)
top10_mac <- mac_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
print(n=200,top10_mac[ c("cluster", "gene")])
DoHeatmap(mac,features = top10_mac$gene)+ggtitle("macrophages subgroup clusters top10 markers")
saveRDS(integrated,"integrated_for_annotated.rds")
# 8.初步注释
mac.cluster.ids <- c("Neural-assosiated","pollution","Resident-like mφ","inflammatory mφ",
                     "IFN mφ","Remodeling mφ","pollution","pollution",
                     "DCs")

names(mac.cluster.ids) <- levels(mac)
mac <- RenameIdents(mac, mac.cluster.ids)
saveRDS(mac,"dm_ctrl_mac.rds")
# 画一张带名字的 UMAP 图
DimPlot(mac, reduction = "umap", label = TRUE, repel = TRUE)+
  ggtitle("DM+CTRL macrophages (3223 cells)")
#轻度purification
# 示例 marker
obj<-readRDS("asys_ctrl_mac.rds")
FeaturePlot(obj, features = c("PTPRC","CD68","C1QA","LYZ","COL1A1","VWF"))
obj <- subset(obj, subset = 
                !COL1A1 & !VWF & !ACTA1
)
saveRDS(obj,"asys_ctrl_mac_p.rds")
