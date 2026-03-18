setwd("D:/R/R/PROJECT/IMNM/Milo")
#1.数据整理
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(scran)
library(qs)
library(BiocParallel)
library(stringr)
myeloid<-readRDS("D:/R/R/PROJECT/IMNM/filterMC/Myeloid_annotation.rds")
# 将 Group 重命名为 group
myeloid$group <- myeloid$Group
# 将 Antibody 重命名为 antibody
myeloid$antibody <- myeloid$Antibody
# 将Sample 重命名为sample
myeloid$sample <- myeloid$Sample
# 可选：删除旧列以避免混淆（不是必须）
myeloid$Group <- NULL
myeloid$Antibody <- NULL
myeloid$Sample <- NULL
# 检查新列名
colnames(myeloid[[]])
# 检查数据
table(myeloid$group, useNA = "ifany")
table(myeloid$antibody, useNA = "ifany")
table(myeloid$sample)  # 确保每个样本有足够细胞



#4.运行疾病组间比较
# 1. 提取 IMNM 和 IBM 的细胞
cells_imnm_ibm <- WhichCells(myeloid, expression = group %in% c("IMNM", "IBM"))
sub_imnm_ibm <- subset(myeloid, cells = cells_imnm_ibm)

# 2. 创建新的 Milo 对象（必须重新构建邻域）
sce_sub <- as.SingleCellExperiment(sub_imnm_ibm)
scRNA <- Milo(sce_sub)
reducedDim(scRNA, "PCA") <- reducedDim(sce_sub, "PCA")
umap_coords <- Embeddings(sub_imnm_ibm, "umap")[colnames(scRNA), ]
reducedDim(scRNA, "UMAP") <- umap_coords

scRNA <- buildGraph(scRNA, k = 30, d = 20, reduced.dim = "PCA")
scRNA <- makeNhoods(scRNA, prop = 0.1, k = 30, d = 20, refined = TRUE)
plotNhoodSizeHist(scRNA)
# 3. 添加样本信息并计数
colData(scRNA)$sample <- sub_imnm_ibm$sample
scRNA <- countCells(scRNA, meta.data = data.frame(colData(scRNA)), samples = "sample")

# 4. 创建设计数据框（样本级别）
design_sub <- data.frame(colData(scRNA))[, c("sample", "group")]
design_sub <- distinct(design_sub)
rownames(design_sub) <- design_sub$sample
design_sub$group <- factor(design_sub$group, levels = c("IBM", "IMNM"))

# 5. 运行 testNhoods
da_result <- testNhoods(scRNA,
                      design = ~ group,
                      design.df = design_sub,
                      fdr.weighting = "graph-overlap")

# 6. 查看结果列名（应包含单列 logFC）
colnames(da_result)
da_result %>%
  arrange(SpatialFDR) %>%
  head() 
ggplot(da_result, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_result, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

# 7. 可视化
scRNA <- buildNhoodGraph(scRNA)

plotNhoodGraphDA(scRNA, da_result, 
                 alpha = 0.1, 
                 color_by = "logFC", 
                 reduced.dim = "UMAP") +
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
    values = rescale(c(min(da_result$logFC), 0, max(da_result$logFC))),
    name = "logFC (IMNM vs IBM)"
  ) +
  labs(title = "IMNM vs IBM")
#8.数据注释
da_result<- annotateNhoods(scRNA,
                             da_result,
                             coldata_col = "cell_type")
head(da_result)

ggplot(da_result, aes(cell_type_fraction))+geom_histogram(bins = 50)

str(da_result)
table(da_result$cell_type)
range(da_result$SpatialFDR)
plotDAbeeswarm(da_result, group.by = "cell_type",alpha = 0.9) # alpha默认0.1
#9.find DA markers
## Add log normalized count to Milo object
scRNA <- logNormCounts(scRNA)
# alpha默认0.1,所以SpatialFDR < 0.9是修改过的
da_result$NhoodGroup <- as.numeric(da_result$SpatialFDR < 0.9 & da_result$logFC < 0)
da_nhood_markers <- findNhoodGroupMarkers(scRNA, da_result, 
                                          subset.row = rownames(scRNA), 
                                          aggregate.samples = TRUE, # 建议设置为TRUE
                                          sample_col = "orig.ident")
head(da_nhood_markers)

#10. Group the Nhood
## Run buildNhoodGraph to store nhood adjacency matrix
scRNA <- buildNhoodGraph(scRNA)
## Find groups
da_result <- groupNhoods(scRNA, da_result, 
                          max.lfc.delta = 2,da.fdr = 0.9)#da.fdr默认0，1
head(da_result)

plotNhoodGroups(scRNA, da_result, layout="UMAP") 
plotDAbeeswarm(da_result, "NhoodGroup",alpha = 0.1) #alpha默认是0.1

#11. Find HVGs
## Exclude zero counts genes
keep.rows <- rowSums(logcounts(scRNA)) != 0
scRNA <- scRNA[keep.rows, ]

## Find HVGs
dec <- modelGeneVar(scRNA)
hvgs <- getTopHVGs(dec, n=2000)
head(hvgs)
# [1] "F13A1"  "SPP1"   "PDE4D"  "RNASE1" "CD36"   "MAMDC2"

# 找到每个邻里群的one-vs-all差异基因表达
nhood_markers <- findNhoodGroupMarkers(scRNA, da_results, subset.row = hvgs, 
                                       aggregate.samples = TRUE, 
                                       sample_col = "orig.ident")

head(nhood_markers)[1:4,1:4]
# GeneID     logFC_1 adj.P.Val_1      logFC_2
# 1    A2M -0.16734295  0.56307711  0.275780760
# 2  AAGAB  0.01489963  0.93378744 -0.031319948
# 3  ABCA1 -0.45589406  0.09844721  0.506785395
# 4 ABCA10  0.04148363  0.08259512 -0.007768381

# 找到group1的差异基因
gr1_markers <- nhood_markers[c("logFC_1", "adj.P.Val_1")] 
colnames(gr1_markers) <- c("logFC", "adj.P.Val")
head(gr1_markers[order(gr1_markers$adj.P.Val), ])
#      logFC    adj.P.Val
# 695  0.33708668 7.359700e-07
# 1903 0.64765859 9.437441e-06
# 1739 0.23144226 2.951303e-04
# 1914 0.17622927 3.979074e-04
# 730  0.20825211 1.886975e-03
# 194  0.07300323 1.941041e-03


# 找到group2的差异基因
gr2_markers <- nhood_markers[c("logFC_2", "adj.P.Val_2")] 
colnames(gr2_markers) <- c("logFC", "adj.P.Val")
head(gr2_markers[order(gr2_markers$adj.P.Val), ])

# 如果明确了对第二组细胞感兴趣
nhood_markers <- findNhoodGroupMarkers(scRNA, da_results, subset.row = hvgs, 
                                       aggregate.samples = TRUE, 
                                       sample_col = "orig.ident",
                                       subset.groups = c("2")
)

head(nhood_markers)
#        logFC_2  adj.P.Val_2 GeneID
# TPRG1  1.4356223 2.594446e-14  TPRG1
# FNIP2  0.9308552 1.519445e-11  FNIP2
# GPNMB  1.1801120 2.123099e-11  GPNMB
# SGK1   0.9691252 3.824231e-10   SGK1
# MYO1E  1.1277269 2.812210e-09  MYO1E
# PLA2G7 0.3732118 4.372560e-08 PLA2G7

# 如果明确了相比较的邻域
nhood_markers <- findNhoodGroupMarkers(scRNA, da_results, subset.row = hvgs,
                                       subset.nhoods = da_results$NhoodGroup %in% c('1','2'),
                                       aggregate.samples = TRUE, sample_col = "orig.ident")

head(nhood_markers)
#   GeneID     logFC_1 adj.P.Val_1     logFC_2 adj.P.Val_2
# 1    A2M -0.25470801   0.4606486  0.25470801   0.4606486
# 2  AAGAB -0.01037831   0.9540445  0.01037831   0.9540445
# 3  ABCA1 -0.50051349   0.1681579  0.50051349   0.1681579
# 4 ABCA10  0.01652818   0.7669192 -0.01652818   0.7669192
# 5  ABCA3 -0.01405609   0.6581034  0.01405609   0.6581034
# 6  ABCA6 -0.02016020   0.9820740  0.02016020   0.9820740

#12. HVGs plot
ggplot(nhood_markers, aes(logFC_1,-log10(adj.P.Val_1 ))) + 
  geom_point(alpha=0.5, size=0.5) +
  geom_hline(yintercept = 2)

# 筛选标准可以自行修改
markers <- rownames(nhood_markers)[nhood_markers$adj.P.Val_2 < 0.01 & nhood_markers$logFC_2 > 0]

plotNhoodExpressionGroups(scRNA, da_results, features=markers,
                          subset.nhoods = da_results$NhoodGroup %in% c('11','2'), 
                          scale=TRUE,
                          grid.space = "fixed")
#13.比较同一邻域的不同疾病组细胞
dge_9 <- testDiffExp(scRNA, da_results, design = ~ group, 
                     meta.data = data.frame(colData(scRNA)),
                     subset.row = rownames(scRNA),
                     subset.nhoods=da_results$NhoodGroup=="4")
dge_9

#14.GO analysis
# 假设您已通过 testDiffExp 或 findNhoodGroupMarkers 获得了每个组的差异基因列表
# 例如，nhood_markers 包含 logFC_1, adj.P.Val_1 等列

library(clusterProfiler)
library(org.Hs.eg.db)

# 以 group1 为例，筛选显著上调基因（logFC > 0.5, adj.P.Val < 0.05）
sig_genes_group1 <- nhood_markers %>%
  filter(logFC_1 > 0.5 & adj.P.Val_1 < 0.05) %>%
  pull(GeneID)
sig_genes_group2 <- nhood_markers %>%
  filter(logFC_2 > 0.5 & adj.P.Val_2 < 0.05) %>%
  pull(GeneID)
gene_list <- list(Group1 = sig_genes_group1, Group2 = sig_genes_group2)
comp_ego <- compareCluster(geneCluster = gene_list,
                           fun = enrichGO,
                           OrgDb = org.Hs.eg.db,
                           keyType = "SYMBOL",   # 指定输入为 SYMBOL
                           ont = "BP",            # 生物学过程，可选 "CC", "MF" 或 "ALL"
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           readable = TRUE)       # 将结果中的 ENTREZ ID 转回 SYMBOL
dotplot(comp_ego, showCategory = 10) +          # 每个组显示前10条通路
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("GO enrichment comparison: Group1 vs Group2")
#15.尝试注释细胞
# 预先定义的巨噬细胞亚型标记
mac_markers <- list(
  Resident = c("LYVE1", "MAMDC2", "RNASE1"),
  SPP1_Fibrotic = c("SPP1", "GPNMB", "TREM2"),
  Inflammatory = c("FCN1", "VCAN", "IL1B"),
  Proliferating = c("MKI67", "TOP2A", "STMN1")
)

# 对每个组，计算这些标记的富集分数（例如使用 AddModuleScore）
for (group_name in names(mac_markers)) {
  features <- list(mac_markers[[group_name]])
  myeloid <- AddModuleScore(myeloid, features = features, name = group_name)
}
# 将得分映射到 Milo 对象（假设细胞名一致）
common_cells <- intersect(colnames(myeloid), colnames(scRNA))
for (group_name in names(mac_markers)) {
  colData(scRNA)[[paste0(group_name, "_score")]] <- 
    myeloid[[paste0(group_name, "1")]][common_cells, ]
}
# 绘制小提琴图或箱线图，比较这些分数在各组中的分布
VlnPlot(scRNA, features = c("Resident1", "SPP1_Fibrotic1", "Inflammatory1"), 
        group.by = "nhood_group", pt.size = 0) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white")
