#1.准备数据集
setwd("D:/R/R/PROJECT/IMNM/MergeMC")
mac_clean_filtered<- readRDS("Purification_myeloids.rds")

bad_clusters <- c("0", "11", "12")
mac_clean_filtered<-subset(mac_clean_filtered,idents = bad_clusters, invert = TRUE)
# 加载所需包
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(Seurat)

# 下载 KEGG 基因集
# c2 <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG")
# term2gene <- c2[, c("gs_name", "gene_symbol")]
# colnames(term2gene) <- c("term", "gene")
# 方法一：加载全部 C2 后过滤（推荐，最保险）
c2_all <- msigdbr(species = "Homo sapiens", category = "C2")
kegg <- c2_all %>% filter(gs_subcat == "CP:KEGG_MEDICUS")

# 查看有多少 KEGG 通路
length(unique(kegg$gs_name))

# 构建 term2gene
term2gene <- kegg[, c("gs_name", "gene_symbol")]
colnames(term2gene) <- c("term", "gene")
#2.GESA分析
# 获取所有簇编号（确保不包含已剔除的簇）
all_clusters <- as.character(unique(Idents(mac_clean_filtered)))
gsea_results_c2 <- list()

for (cl in all_clusters) {
  cat("Processing cluster", cl, "\n")
  
  # 计算该簇 vs 其他所有细胞的差异（全基因）
  markers <- FindMarkers(mac_clean_filtered,
                         ident.1 = cl,
                         logfc.threshold = 0,
                         min.pct = 0.05,
                         assay = "RNA")
  
  # 构建排序列表
  gene_list <- markers$avg_log2FC
  names(gene_list) <- rownames(markers)
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!is.na(gene_list)]
  
  # GSEA
  gsea_res <- GSEA(gene_list,
                   TERM2GENE = term2gene,
                   pvalueCutoff = 0.05,
                   verbose = FALSE)
  
  gsea_results_c2[[cl]] <- gsea_res
}
#3.可视化
# 提取各簇显著通路（按 NES 绝对值前5）
library(ggplot2)
library(patchwork)

# 创建空列表存储非空结果
gsea_list <- list()

for (cl in names(gsea_results_c2)) {
  # 检查结果是否为 NULL 或结果数据框为空
  if (!is.null(gsea_results_c2[[cl]])) {
    res <- gsea_results_c2[[cl]]@result
    if (nrow(res) > 0) {
      res <- res[order(abs(res$NES), decreasing = TRUE), ]
      res$Cluster <- cl
      gsea_list[[cl]] <- head(res, 5)
    } else {
      cat("Cluster", cl, "has no significant terms.\n")
    }
  } else {
    cat("Cluster", cl, "result is NULL.\n")
  }
}

# 合并所有非空结果
if (length(gsea_list) > 0) {
  gsea_df <- do.call(rbind, gsea_list)
  print(head(gsea_df))
} else {
  cat("No significant GSEA results for any cluster.\n")
}

# 绘制 NES 热图（每个簇的前5通路）
library(reshape2)
library(pheatmap)

nes_mat <- dcast(gsea_df, ID ~ Cluster, value.var = "NES", fill = 0)
rownames(nes_mat) <- nes_mat$ID
nes_mat <- as.matrix(nes_mat[, -1])

pheatmap(nes_mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top KEGG Pathways (|NES|) per cluster")
if (length(gsea_list) > 0) {
  library(ggplot2)
  ggplot(gsea_df, aes(x = Cluster, y = ID)) +
    geom_point(aes(size = -log10(p.adjust), color = NES)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    theme_bw() +
    labs(title = "KEGG Pathway Enrichment (Top 5 |NES| per cluster)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}