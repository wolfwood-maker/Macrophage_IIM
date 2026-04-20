library(monocle3)
library(Seurat)
library(SeuratWrappers) # 推荐使用这个包进行格式转换
library(magrittr)

# 1. 转换格式
cds <- as.cell_data_set(mac_final)

# 2. 获取 UMAP 坐标和分群信息（Monocle3 必须建立在已有的 UMAP 上）
cds <- estimate_size_factors(cds)
# 1. 运行 Monocle 3 特有的聚类和分区
# 这一步会计算 partition
cds <- cluster_cells(cds, resolution=1e-3) 
# 第一步：学习图结构（这一步会计算轨迹线）
cds <- learn_graph(cds, use_partition = TRUE)

# 第二步：可视化（此时就不会报错了）
plot_cells(cds, 
           color_cells_by = "seurat_clusters", # 看看你的 Seurat 分群在轨迹上的位置
           label_groups_by_cluster = FALSE,
           label_leaves = TRUE, 
           label_branch_points = TRUE)

# 第三步：指定起点（假设 Cluster 5 是起点）
# 运行后点击 Cluster 5 所在的根部节点，然后关闭窗口
cds <- order_cells(cds)
#4. 结果检查
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_cell_groups = FALSE, 
           label_leaves = FALSE, 
           label_branch_points = FALSE)
