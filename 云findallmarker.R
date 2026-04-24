# ======================================================
#  单细胞数据分析完整流程（内存充裕版，90GB RAM）
#  数据来源：final_analysis_ready.rds (8.8GB)
#  作者：AI Assistant
#  日期：2026-04-23
# ======================================================

# ------------------- 0. 环境设置 ----------------------
# 清空环境，避免残留
rm(list = ls())
gc()

# 设置临时文件目录（可选，若磁盘空间紧张可放到 /tmp）
# 默认 R 临时目录在 /tmp，一般有足够空间，但若担心可改到有大空间的地方
# tempdir <- "/root/Rtmp"
# dir.create(tempdir, showWarnings = FALSE)
# Sys.setenv(TMPDIR = tempdir)

# 设置内存上限（90GB 内存，给 R 分配 80GB 上限，留下给系统）
options(future.globals.maxSize = 80 * 1024^3)

# 加载必要的包（如果未安装，请先 install.packages）
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)

# 如果需要 Harmony，加载之
library(harmony)

# ------------------- 1. 加载数据 ----------------------
# 修改为你的实际文件路径
rds_path <- "/root/final_analysis_ready.rds"   # 请按实际修改

cat("正在加载 RDS 文件，文件大小：", file.info(rds_path)$size / 1e9, "GB\n")
srat <- readRDS(rds_path)
cat("加载完成。Seurat 对象包含", ncol(srat), "个细胞，", nrow(srat), "个基因\n")
cat("当前内存使用量：", format(object.size(srat), units = "GB"), "\n")

# 查看当前有哪些 reductions（如果有现成的 PCA/UMAP 可以保留）
print(names(srat@reductions))
#UMAP图
# 如果已经有 umap 降维结果
if ("umap" %in% names(srat@reductions)) {
  p_raw <- DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE) + 
    ggtitle("Original UMAP (before removal)")
  ggsave("umap_original.png", plot = p_raw, width = 10, height = 8, dpi = 300)
  print(p_raw)  # 在服务器上显示（如果使用 X11 或 RStudio）
} else {
  cat("未找到现成的 UMAP，需要先运行 RunUMAP...\n")
  # 如果只有 PCA，先基于 PCA 跑 UMAP
  if (!"pca" %in% names(srat@reductions)) {
    srat <- RunPCA(srat, npcs = 50, verbose = FALSE)
  }
  srat <- RunUMAP(srat, reduction = "pca", dims = 1:30, verbose = FALSE)
  p_raw <- DimPlot(srat, reduction = "umap", label = TRUE) + 
    ggtitle("Original UMAP (computed now)")
  ggsave("umap_original.png", plot = p_raw, width = 10, height = 8, dpi = 300)
}
# ------------------- 2. 定义屏蔽基因列表 ----------------------
# 肌肉基因（你的列表）
muscle_genes <- c(
  "MYH1","MYH2","MYH3","MYH7","ACTA1","ACTN2","TTN","NEB",
  "TNNT1","TNNT2","TNNT3","TNNI1","TNNI2","TPM1","TPM2","TPM3",
  "MYBPC1","MYBPC2","MYL1","MYL2","MYL3","MYLPF","CKM","DES",
  "DMD","XIRP2","CMYA5","MYOT"
)

# 红细胞污染基因
rbc_genes <- c("HBB", "HBA1", "HBA2", "HBM", "HBQ1", "HBZ")

# 线粒体基因 (MT- 开头)
mt_genes <- grep("^MT-", rownames(srat), value = TRUE)

# rRNA 基因
rrna_genes <- c("MT-RNR1", "MT-RNR2", "RNA18SN5", "RNA28SN5")

# 神经相关基因（若组织不是脑，可能是污染）
neural_genes <- c("NOVA1","ROBO2","LRRTM4","BDNF","DLGAP2","ELFN1")

# Y 染色体基因（若样本均为女性则可屏蔽）
y_chromosome_genes <- c("DDX3Y", "PRKY", "USP9Y", "UTY", "TTTY14", "TTTY15", "NLGN4Y")

# IG 基因（免疫球蛋白污染）
ig_genes <- grep("^IG[KHL]", rownames(srat), value = TRUE)

# 合并所有要屏蔽的基因（在差异表达中排除）
exclude_genes <- unique(c(
  muscle_genes, rbc_genes, mt_genes, rrna_genes,
  neural_genes, y_chromosome_genes, ig_genes
))
exclude_genes <- intersect(exclude_genes, rownames(srat))
cat("屏蔽基因数量：", length(exclude_genes), "\n")

# 所有可测试的基因（排除屏蔽基因）
test_genes <- setdiff(rownames(srat), exclude_genes)
cat("参与差异分析的基因数量：", length(test_genes), "\n")

# ------------------- 3. 第一次 FindAllMarkers（识别污染簇） ----------------------
# 注意：这一步可能会比较耗时（取决于细胞数和簇数），但内存安全
# 确保 layers 合并（Seurat v5 需要）
srat <- JoinLayers(srat)

# 如果还没有标准化，需要先标准化（假设原始 rds 可能已经标准化过，但为了保险）
if (!"data" %in% Layers(srat, assay = "RNA")) {
  srat <- NormalizeData(srat)
}

cat("开始 FindAllMarkers，预计需要一段时间...\n")
markers_all <- FindAllMarkers(
  srat,
  features = test_genes,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  assay = "RNA"
)

# 提取每个簇的 top10 标记基因
top10 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

# 保存结果到 CSV 文件，方便查看
write.csv(top10, "top10_markers_all_clusters.csv", row.names = FALSE)
cat("Top10 标记基因已保存至 top10_markers_all_clusters.csv\n")

# 打印前几条，方便快速浏览
print(head(top10, 20))
DoHeatmap(srat,features = top10$gene)+
  ggtitle("macrophages subgroup clusters top10 markers")
# ------------------- 4. 手动识别污染簇并剔除 ----------------------
# ！！！重要！！！
# 请打开 top10_markers_all_clusters.csv，根据每个簇的前10个标记基因，判断哪些簇是污染。
# 常见的污染信号：
#   - 红细胞污染：HBB, HBA1 等高表达
#   - IG 基因高：可能是 B 细胞/浆细胞污染
#   - 线粒体基因比例异常高（但 snRNA 中线粒体比例通常较低）
#   - 神经基因高（如果研究不是脑组织）
#   - Y 染色体基因高（如果样本都是女性）
#   - 肌肉基因高（如果研究不是肌肉组织）

# 根据你的判断，设置要删除的簇编号。例如：
# bad_clusters <- c("2", "5", "8")   # 注意：簇编号是字符型或数字型，取决于数据
# 下面是一个示例，你需要手动修改！！！
bad_clusters <- c("3", "4", "8", "9", "11", "12", "14")   # 请根据实际结果修改

# 剔除污染簇
srat_clean <- subset(srat, subset = !(seurat_clusters %in% bad_clusters))
cat("剔除污染簇后剩余细胞数：", ncol(srat_clean), "\n")

# 可选：释放原始对象内存（如果你不再需要）
rm(srat); gc()

# ------------------- 5. 对纯净子集进行标准化、降维、聚类 ----------------------
# 5.1 标准化和寻找高变基因
srat_clean <- NormalizeData(srat_clean)
srat_clean <- FindVariableFeatures(srat_clean, nfeatures = 3000)

# 5.2 缩放数据 (ScaleData)
# 这一步可能会消耗较多内存，但 90GB 应该足够
cat("正在进行 ScaleData...\n")
srat_clean <- ScaleData(srat_clean)

# 5.3 PCA
srat_clean <- RunPCA(srat_clean, npcs = 50, verbose = FALSE)

# 5.4 确定是否需要进行 Harmony 批次校正
# 假设你的 metadata 中有一列叫做 "batch" 或 "sample"，代表不同样本/批次。
# 如果没有，则跳过 Harmony。请根据你的实际情况修改列名。
batch_column <- "batch"   # 改成你 metadata 中的批次列名，如果没有则保持 NULL

if (!is.null(batch_column) && batch_column %in% colnames(srat_clean[[]])) {
  cat("检测到批次列 ", batch_column, "，运行 Harmony 进行整合...\n")
  srat_clean <- RunHarmony(srat_clean, group.by.vars = batch_column)
  reduction_use <- "harmony"
  dims_use <- 1:30
} else {
  cat("未检测到批次列，跳过 Harmony，直接使用 PCA 降维\n")
  reduction_use <- "pca"
  dims_use <- 1:30
}

# 5.5 聚类和 UMAP
srat_clean <- FindNeighbors(srat_clean, reduction = reduction_use, dims = dims_use)
srat_clean <- FindClusters(srat_clean, resolution = 0.5)
srat_clean <- RunUMAP(srat_clean, reduction = reduction_use, dims = dims_use)

# 5.6 可视化
p1 <- DimPlot(srat_clean, label = TRUE) + ggtitle("Clean Subset UMAP")
p2 <- DimPlot(srat_clean, group.by = batch_column) + ggtitle("Batch Distribution")
p_combined <- p1 + p2
print(p_combined)

# 保存图片
ggsave("umap_clean_subset.png", plot = p_combined, width = 12, height = 6)

# ------------------- 6. 最终注释 ----------------------
# 再次运行 FindAllMarkers（排除屏蔽基因）以获得最终聚类群的标记基因
cat("正在执行最终差异分析以注释细胞类型...\n")
markers_final <- FindAllMarkers(
  srat_clean,
  features = setdiff(rownames(srat_clean), exclude_genes),
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top10_final <- markers_final %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(top10_final, "top10_final_clusters.csv", row.names = FALSE)
cat("最终聚类 top10 标记已保存至 top10_final_clusters.csv\n")

# 基于 top10_final 表格中的标记基因，手动确定每个簇的细胞类型。
# 例如：
new.cluster.ids <- c("Myonuclei", "Fibroblasts", "Endothelial", "Immune", "Other")
# 注意：向量的长度必须等于簇的数量，按照 cluster 的数值顺序排列。
# 你可以先查看有哪些簇：
cat("当前簇编号：", levels(srat_clean), "\n")

# 请根据实际填充下面的向量
# new.cluster.ids <- c("Cluster0_type", "Cluster1_type", ...)

# 重命名
if (length(new.cluster.ids) == length(levels(srat_clean))) {
  names(new.cluster.ids) <- levels(srat_clean)
  srat_clean <- RenameIdents(srat_clean, new.cluster.ids)
  cat("细胞类型重命名完成\n")
} else {
  cat("警告：new.cluster.ids 长度与簇数不一致，请手动修正后再运行 RenameIdents\n")
}

# 保存最终对象
saveRDS(srat_clean, "srat_clean_final.rds")
cat("最终 Seurat 对象已保存至 srat_clean_final.rds\n")

# 导出 UMAP 坐标和 metadata 为 CSV，便于外部查看
umap_coords <- Embeddings(srat_clean, "umap")
write.csv(umap_coords, "umap_coordinates.csv")
write.csv(srat_clean[[]], "metadata_clean.csv")

# ------------------- 7. 清理临时文件（可选） ----------------------
# 如果磁盘空间紧张，可以删除一些中间文件
unlink("top10_markers_all_clusters.csv")
# unlink("srat_clean_final.rds")   # 请勿删除最终结果

cat("分析完成！\n")