#!/usr/bin/env Rscript

# ==============================================================================
# 脚本说明：Seurat V5 + Harmony 内存优化版（带分步检查点）
# 适用环境：Linode 48GB 或类似内存受限服务器
# ==============================================================================

library(Seurat)
library(harmony)
library(dplyr)

# --- 核心：检查点辅助函数 ---
# 如果文件存在则读取；如果不存在则执行代码块并保存结果
load_or_run <- function(file_path, run_code_block) {
  if (file.exists(file_path)) {
    cat(paste0("\n[CHECKPOINT] 发现备份文件：", file_path, "，直接载入...\n"))
    return(readRDS(file_path))
  } else {
    cat(paste0("\n[COMPUTING] 未发现备份：", file_path, "，开始计算阶段任务...\n"))
    result <- run_code_block
    cat(paste0("[SAVING] 阶段任务完成，保存检查点至：", file_path, "\n"))
    saveRDS(result, file = file_path, compress = FALSE) # 不压缩以换取保存速度
    return(result)
  }
}

cat("=== Start Optimized Harmony pipeline with Checkpoints ===\n")

# ===========================
# 阶段 1：数据加载、拆分与预处理
# ===========================
# 目标：生成处理好的 final_sample_list
final_sample_list <- load_or_run("cp1_preprocessed_list.rds", {
  files <- c("IMNM" = "imnm_ctrl_mac.rds", 
             "IBM"  = "ibm_ctrl_mac.rds", 
             "ASyS" = "asys_ctrl_mac.rds", 
             "DM"   = "dm_ctrl_mac.rds")
  
  tmp_list <- list()
  for (nm in names(files)) {
    cat("  Reading", nm, "...\n")
    obj <- readRDS(files[nm])
    obj$disease <- nm
    obj <- DietSeurat(obj, counts = TRUE, data = TRUE, scale.data = FALSE)
    
    # 拆分并存入 list
    samples <- SplitObject(obj, split.by = "sample_id")
    tmp_list <- c(tmp_list, samples)
    rm(obj, samples); gc()
  }
  
  cat("  Normalizing and Finding Variable Features per sample...\n")
  min_cells <- 20
  tmp_list <- lapply(tmp_list, function(x) {
    if (ncol(x) < min_cells) return(NULL)
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, nfeatures = 2000, verbose = FALSE)
    return(x)
  })
  Filter(Negate(is.null), tmp_list)
})

# ===========================
# 阶段 2：高变基因选择与最终合并
# ===========================
# 目标：生成 merge 后的对象（解决你目前卡住 4 小时的步骤）
merged <- load_or_run("cp2_merged_obj.rds", {
  cat("Selecting integration features...\n")
  select_features <- SelectIntegrationFeatures(object.list = final_sample_list, nfeatures = 2000)
  
  cat("Final merging (This may take a while)...\n")
  # 记录高变基因到第一个对象中，方便合并后继承
  merged_obj <- merge(x = final_sample_list[[1]], y = final_sample_list[-1])
  VariableFeatures(merged_obj) <- select_features
  
  rm(final_sample_list); gc()
  merged_obj
})

# ===========================
# 阶段 3：ScaleData 与 PCA
# ===========================
merged <- load_or_run("cp3_pca_obj.rds", {
  cat("Scaling and Running PCA...\n")
  # 关键优化：仅 Scale 高变基因
  merged <- ScaleData(merged, features = VariableFeatures(merged), verbose = FALSE)
  merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
  gc()
  merged
})

# ===========================
# 阶段 4：Harmony 整合
# ===========================
merged <- load_or_run("cp4_harmony_obj.rds", {
  cat("Running Harmony...\n")
  merged <- RunHarmony(merged, group.by.vars = "sample_id", reduction = "pca", assay.use = "RNA")
  gc()
  merged
})

# ===========================
# 阶段 5：UMAP 与 聚类（最终产物）
# ===========================
cat("Final UMAP and Clustering...\n")
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.4, verbose = FALSE)

cat("Saving final integrated object...\n")
saveRDS(merged, file = "Harmony_integrated_macrophage_final.rds", compress = FALSE)

cat("=== ALL DONE ===\n")