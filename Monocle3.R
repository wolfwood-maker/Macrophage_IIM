library(monocle3)
library(Seurat)
library(SeuratWrappers) # 推荐使用这个包进行格式转换
library(magrittr)
mac_final<-readRDS("Purification_macrophages.rds")
bad_clusters <- c( 9, 14, 15)
mac_final<-subset(mac_final,idents = bad_clusters, invert = TRUE)
#对肌肉污染的cluster 1/6 使用decontX assay
#mac_final<-subset(mac_final,idents = setdiff(levels(mac_final), c("10", "12", "13", "14", "15")))#踢掉污染
# 1. 转换格式
cds <- as.cell_data_set(mac_final)

# 2. 获取 UMAP 坐标和分群信息（Monocle3 必须建立在已有的 UMAP 上）
cds <- estimate_size_factors(cds)
# 1. 运行 Monocle 3 特有的聚类和分区
# 这一步会计算 partition
cds <- cluster_cells(cds, resolution=1e-3) 
# 第一步：学习图结构（这一步会计算轨迹线）
# 核心参数：增加 ncenter 增加图密度，将 close_loop 设为 FALSE 允许长距离路径
# 重新计算轨迹，禁用分区以强制跨区域连接
 cds <- learn_graph(cds, use_partition = FALSE, 
                    learn_graph_control = list(ncenter = 700))
# 调整 learn_graph 时的参数，增加分组的独立性
#cds <- learn_graph(cds, close_loop = FALSE, learn_graph_control = list(ncenter = 500))
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
           label_branch_points = FALSE)+
  ggtitle("Monocyte-Macropahges trajectory")
