

lapply(c("TNK/group_FL/slinshot/"), FUN = dir.create)
options(future.globals.maxSize = 100 * 1024^3)

CD8_T<-subset(scRNA_harmony, subset = (celltype_treg == "CD8_T"))

######
#BiocManager::install("kstreet13/slingshot")
library(slingshot)
library(Seurat)
library(devtools)
library(Seurat)
library(cowplot)
library(ggplot2)
library(Matrix)
library(dplyr)
#BiocManager::install("tradeSeq")
library(tradeSeq)
library(RColorBrewer)

#github：https://github.com/stellar/slingshot
#官网使用手册：https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html

#slingshot轨迹分析需要准备两个文件
#1、具有降维信息的cell Matrix
#2、celltype_CXCL13 cluster label

#============================================================================
#                            step1:准备分析文件run slingshot
#============================================================================
#CD8_T <- readRDS("D:/KS项目/公众号文章/monocle3拟时分析/CD8_T.rds")

#我们大多使用Seurat单细胞对象，这里准备文件的后见一般有两种方式

###method1，最简单的办法，直接将Seurat对象转为SingleCellExperiment，因为slingshot的input支持：matrix, SingleCellExperiment, SlingshotDataSet, and PseudotimeOrdering.
sce <- as.SingleCellExperiment(CD8_T, assay = "RNA")#至于assay选啥按照自己实际数据，可能你的数据是integrated or SCT
#slingshot可以选择轨迹起点和终点，这里我们设置起点，起点的选择一般依据数据实际意义，生物学背景或者可以使用
#我们之前讲过的CytoTRACE去辅助推断轨迹root，不过还是要结合具体的生物学背景，不可盲目
table(CD8_T$celltype_CXCL13)
sce_slingshot1 <- slingshot(sce, #这里是SingleCellExperiment单细胞对象
                            reducedDim = 'UMAP', #降维方式
                            clusterLabels = sce$celltype_CXCL13, #celltype
                            start.clus = 'CD8_T_negative', #轨迹起点,也可以不定义
                            approx_points = 150)


###method2，一步一步提取
#提取表达矩阵并构建SingleCellExperiment单细胞对象
expdata <- as.matrix(GetAssayData(object = CD8_T, slot = "data"))
sce <- SingleCellExperiment(assays = List(counts = expdata))
#添加降维信息到SingleCellExperiment对象
reduction_sce <- Embeddings(object = CD8_T, reduction = "umap")
reducedDims(sce) <- SimpleList(UMAP = reduction_sce)
colData(sce)$celltype <- CD8_T@meta.data$celltype
colData(sce)$orig.ident <- CD8_T@meta.data$orig.ident #这里的orig.ident在我这个数据里面相当于分组

#run slingshot
sce_slingshot2 <- slingshot(sce, 
                            clusterLabels = 'celltype_CXCL13', 
                            reducedDim = 'UMAP', 
                            start.clus = 'CD8_T_negative', 
                            approx_points = 150)


# SlingshotDataSet(sce_slingshot2) 查看轨迹信息

#可以看出，我们这个数据推断出3条轨迹
# class: SlingshotDataSet 
# 
# Samples Dimensions
# 6025          2
# 
# lineages: 3 
# Lineage1: CD8_T_negative  CD8_T_positive  NK  T_CD4conv_negative  T_CD4conv_positive  
# Lineage2: CD8_T_negative  CD8_T_positive  NK  PMN(3)  
# Lineage3: CD8_T_negative  CD8_T_positive  PMN(6)  PMN(7)  
# 
# curves: 3 
# Curve1: Length: 18.39	Samples: 3918.34
# Curve2: Length: 19.222	Samples: 3543.35
# Curve3: Length: 20.293	Samples: 2525.62


#============================================================================
#                            step2:可视化轨迹
#============================================================================
#定义颜色
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

#----------------------------------------------------------------------------

#设置颜色并可视化。这里分组是celltype_CXCL13。在cell type上可视化轨迹
library(scales)
library(RColorBrewer)
pdf("TNK/group_FL/slinshot/CD8_T_CXCL13_umap图_拟时序.pdf", width =5, height =5)
cell_colors <- cell_pal(sce_slingshot2$celltype_CXCL13, brewer_pal("qual", "Set2"))

plot(reducedDims(sce_slingshot2)$UMAP, col = cell_colors, pch=16, asp = 1, cex = 0.4)
lines(SlingshotDataSet(sce_slingshot2), lwd=2, col='black')

#计算细胞celltype_CXCL13中心坐标，用于标记
celltype_CXCL13_label <- CD8_T@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(celltype_CXCL13 = CD8_T@meta.data$celltype_CXCL13) %>%
  group_by(celltype_CXCL13) %>%
  summarise(UMAP1 = median(UMAP_1),
            UMAP2 = median(UMAP_2))

for (i in 1:8) {
  text(celltype_CXCL13_label$celltype_CXCL13[i], x=celltype_CXCL13_label$UMAP1[i]-1, y=celltype_CXCL13_label$UMAP2[i])
}
dev.off()


#还可以将分组设置为样本
group_colors <- cell_pal(sce_slingshot2$group_FL, brewer_pal("qual", "Set2"))
plot(reducedDims(sce_slingshot2)$UMAP, col = group_colors, pch=16, asp = 1, cex = 0.4)
lines(SlingshotDataSet(sce_slingshot2), lwd=2, col='black')
#添加legend即可
legend("top",
       pch=15,#形状
       legend=unique(sce$group_FL),#文字
       bty="n",
       col =unique(group_colors),#颜色
       border="black",
       horiz=T) 


#----------------------------------------------------------------------------

##从前面的图上可以看出，我们的数据有三条轨迹,分别可视化三条轨迹拟时顺序
pdf("TNK/group_FL/slinshot/CD8_T_CXCL13_umap图_拟时序_热图2.pdf", width =5, height =5)
plot(reducedDims(sce_slingshot2)$UMAP, asp = 1, pch = 16,cex = 0.4, xlab = "UMAP_1", ylab = "UMAP_2",
     col = hcl.colors(100, alpha = 0.5)[cut(sce_slingshot2$slingPseudotime_1, breaks = 100)])
lines(SlingshotDataSet(sce_slingshot2), linInd = 1, lwd = 2, col = 'black')
dev.off()

plot(reducedDims(sce_slingshot2)$UMAP, asp = 1, pch = 16, xlab = "UMAP_1", ylab = "UMAP_2",
     col = hcl.colors(100, alpha = 0.5)[cut(sce_slingshot2$slingPseudotime_2, breaks = 100)])
lines(SlingshotDataSet(sce_slingshot2), linInd = 2, lwd = 2, col = 'black')

plot(reducedDims(sce_slingshot2)$UMAP, asp = 1, pch = 16, xlab = "UMAP_1", ylab = "UMAP_2",
     col = hcl.colors(100, alpha = 0.5)[cut(sce_slingshot2$slingPseudotime_3, breaks = 100)])
lines(SlingshotDataSet(sce_slingshot2), linInd = 3, lwd = 2, col = 'black')

#----------------------------------------------------------------------------

##拟时密度图，类似于monocle2中的图形，绘制细胞或者样本分组在拟时轨迹的分布
#这里我们以第1条轨迹为示例
table(CD8_T$celltype_CXCL13)
SlingshotDataSet(sce_slingshot2)#这里查看有哪些轨迹，轨迹中包含哪些细胞
density1 <- list(a_1 = density(slingPseudotime(sce_slingshot2)[colData(sce_slingshot2)$celltype_CXCL13 == "CD8_T_negative", 1], na.rm = T), #这里的1表示取第一条轨迹
                 a_2 = density(slingPseudotime(sce_slingshot2)[colData(sce_slingshot2)$celltype_CXCL13 == "CD8_T_positive", 1], na.rm = T))

#作图范围
xlim <- range(c(density1$a_1$x, density1$a_2$x))
ylim <- range(c(density1$a_1$y, density1$a_2$y))

par(mar=c(6, 6, 6, 6), xpd = TRUE)
plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")#基本作图
#添加密度图需要修改
for(i in 1:2) {
  
  polygon(c(min(density1[[i]]$x),density1[[i]]$x, max(density1[[i]]$x)), c(0, density1[[i]]$y, 0),
          col = alpha(brewer.pal(2, "Set1")[i], alpha = 0.5))#有几个样本颜色就填几
  
}

legend("right",
       inset= -0.2,#legend位于画框右侧正中
       pch=15,#形状
       legend=c("CD8_T_negative","CD8_T_positive"),#文字
       bty="n",
       col = alpha(brewer.pal(2, "Set1"), alpha = 0.5),#修改下颜色
       border="black",
       title = "celltype_CXCL13",
       cex = 0.8) 


###也可以看样本的密度图
density1_s <- list(a_1 = density(slingPseudotime(sce_slingshot2)[colData(sce_slingshot2)$group_FL == "nTLS", 1], na.rm = T), #这里的1表示取第一条轨迹
                   a_2 = density(slingPseudotime(sce_slingshot2)[colData(sce_slingshot2)$group_FL == "FL", 1], na.rm = T))

#作图范围
xlim <- range(c(density1_s$a_1$x, density1_s$a_2$x))
ylim <- range(c(density1_s$a_1$y, density1_s$a_2$y))

par(mar=c(6, 6, 6, 6), xpd = TRUE)
plot(xlim, ylim, col = "white", xlab = "", ylab = "")#基本作图
#添加密度图
polygon(c(min(density1_s[[1]]$x),density1_s[[1]]$x, max(density1_s[[1]]$x)), c(0, density1_s[[1]]$y, 0),
        col = alpha(brewer.pal(7, "Set1")[6], alpha = 0.5))#有几个样本颜色就填几
polygon(c(min(density1_s[[2]]$x),density1_s[[2]]$x, max(density1_s[[2]]$x)), c(0, density1_s[[2]]$y, 0),
        col = alpha(brewer.pal(7, "Set1")[7], alpha = 0.5))#有几个样本颜色就填几

legend("right",
       inset= -0.3,#legend位于画框右侧正中
       pch=15,#形状
       legend=c("nTLS","FL"),#文字
       bty="n",
       col = brewer.pal(8, "Set1")[7:8],#颜色
       border="black",
       title = "group",
       cex = 0.8) 
#做一下检验
ks.test(slingPseudotime(sce_slingshot2)[colData(sce_slingshot2)$group_FL == "nTLS", 1], na.rm = T,
        slingPseudotime(sce_slingshot2)[colData(sce_slingshot2)$group_FL == "FL", 1], na.rm = T)

# Asymptotic two-sample Kolmogorov-Smirnov test
# 
# data:  slingPseudotime(sce_slingshot2)[colData(sce_slingshot2)$group_FL == "nTLS", 1] and slingPseudotime(sce_slingshot2)[colData(sce_slingshot2)$group_FL == "FL", 1]
# D = 0.05001, p-value = 0.009179
# alternative hypothesis: two-sided

####其他轨迹的可视化参照上面的代码


#----------------------------------------------------------------------------
#当然了，slingshot结果可以与最先的单细胞seurat对象联动
#这里我们将三个轨迹的pseudotime添加到单细胞对象中（Metadata）
pseudotime = slingPseudotime(sce_slingshot2)%>% as.data.frame() 
Lineages = colnames(pseudotime)

for(i in 1:2){
  pseudotime_sub <- pseudotime[,i]
  CD8_T <- AddMetaData(object = CD8_T,
                       metadata = pseudotime_sub,
                       col.name = Lineages[i])
}

#这样我们就可以使用ggplot做密度图了，可调整性也就变宽了
#定义一个函数，作图更加方便
pseudotime_density <- function(seurat_obj,#已经添加了轨迹pseudotime
                               Lineage,#轨迹
                               cluster_label,#分群的列名，例如我注释的celltype_CXCL13列名是celltype_CXCL13，也可以是sample或者分组
                               color){
  df <- data.frame(seurat_obj[[cluster_label]], seurat_obj[[Lineage]]) 
  colnames(df) <- c("celltype_CXCL13", "Lineage")
  na.omit(df) -> df
  p <- ggplot(df, aes(x=Lineage, fill=celltype_CXCL13)) +
    geom_density(alpha=0.5) + 
    theme_bw()+
    scale_fill_manual(values = color)
  print(p)
  
}


library(dittoSeq)
p1 <- pseudotime_density(CD8_T, Lineage="Lineage1", cluster_label = "celltype_CXCL13", color = dittoColors())+
  theme(axis.title.x  = element_blank())
p2 <- pseudotime_density(CD8_T, Lineage="Lineage1", cluster_label = "group_FL", color = c("#DC050C","#8DD3C7"))

#拼图
pdf("TNK/group_FL/slinshot/CD8_T_celltype_CXCL13_group_FL_umap图_拟时序_火山图.pdf", width =5, height =5)
library(patchwork)
p1+p2+plot_layout(ncol = 1,heights=c(1,1),guides = 'collect')
dev.off()


#----------------------------------------------------------------------------
#基因在不同轨迹随拟时表达变化
slingsce<-SlingshotDataSet(sce_slingshot2)
pseudotimeED <- slingPseudotime(slingsce, na = FALSE)
cellWeightsED <- slingCurveWeights(slingsce)

#这里分析我纳入的是虽有的基因
counts<-sce_slingshot2@assays@data@listData$counts
sce_slinghot <- fitGAM(counts = counts, pseudotime = pseudotimeED, cellWeights = cellWeightsED, nknots = 5, verbose = T)

#这里分析我纳入部分基因
library(Matrix)
# 提取感兴趣的基因
genes_plot <- c("CXCL13", "ITGAE", "IL7R", "KLF2")
# 确保 counts 是一个稀疏矩阵
counts <- sce_slingshot2@assays@data@listData$counts
# 检查 counts 的行名是否包含在 genes_plot 中
genes_available <- rownames(counts) %in% genes_plot
# 提取包含这些基因的行
subset_counts <- counts[genes_available, ]
# 确保新的矩阵仍然是稀疏矩阵
subset_counts <- as(subset_counts, "dgCMatrix")
# 运行 fitGAM
sce_slinghot <- fitGAM(counts = subset_counts, pseudotime = pseudotimeED, cellWeights = cellWeightsED, nknots = 5, verbose = TRUE)
####
mean(rowData(sce_slinghot)$tradeSeq$converged)
rowData(sce_slinghot)$assocRes <- associationTest(sce_slinghot, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sce_slinghot)$assocRes
#作图
gene_dynamic <- list()
genes_plot <- c("CXCL13","ITGAE","IL7R","KLF2")
for(i in 1:length(genes_plot)){
  p = plotSmoothers(sce_slinghot, assays(sce_slinghot)$counts,
                    gene =genes_plot[i], alpha = 0.6, border = T, lwd = 2)+
    ggtitle(genes_plot[i])
  gene_dynamic[[i]] <- p
}
pdf("TNK/group_FL/slinshot/CD8_T_celltype_4个基因_拟时序.pdf", width =7, height =5)
Seurat::CombinePlots(gene_dynamic, ncol = 2)
dev.off()



library(ggplot2)
library(SingleCellExperiment)
library(tradeSeq)

# 设定基因列表
genes_plot <- c("CXCL13", "ITGAE", "IL7R", "KLF2")

# 创建一个空的 data.frame 来存储所有基因的数据
combined_data <- data.frame()

# 提取并整理每个基因的数据
for (gene in genes_plot) {
  # 提取平滑曲线数据
  p <- plotSmoothers(sce_slinghot, assays(sce_slinghot)$counts, gene = gene)
  
  # 提取数据点
  smooth_points <- as.data.frame(ggplot_build(p)$data[[2]])  # 只提取平滑曲线的数据
  
  # 添加基因信息
  smooth_points$gene <- gene
  
  # 合并数据
  combined_data <- rbind(combined_data, smooth_points)
}

# 修改列名以匹配 ggplot2 期望的输入
colnames(combined_data)[colnames(combined_data) == "y"] <- "expression"
colnames(combined_data)[colnames(combined_data) == "x"] <- "pseudotime"

# 自定义颜色
colors <- c("CXCL13" = "#E41A1C", "IL7R" = "#377EB8", "ITGAE" = "#4DAF4A", "KLF2" = "#984EA3")

# 作图
p <- ggplot(combined_data, aes(x = pseudotime, y = expression, color = gene)) +
  geom_line(size = 1.5) + # 绘制平滑线
  scale_color_manual(values = colors) + # 使用自定义颜色
  labs(title = "Gene Expression Over Pseudotime", x = "Pseudotime", y = "Expression (log(CP10k + 1))") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, size = 1) # 添加边框
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # 去除网格底

# 显示图形
print(p)

# 保存图形
ggsave("TNK/group_FL/slinshot/CD8_T_celltype_4_genes_pseudotime.pdf", plot = p, width = 7, height = 5)

sce_slinghot
assays(sce_slinghot)$counts





















#============================================================================
#                            step3:拟时基因热图分析
#============================================================================
#前面我们依据推断了轨迹，并进行了一些简单的可视化，那么后续其实最主要的就是关注
#随着拟时轨迹的基因变化。我们可以借助tradeSeq包进行分析。

#这里我们示例分析一下轨迹1,将轨迹1的拟时添加到Seurat对象
sce_slingshot2$sling_pseudotime = sce_slingshot2[[paste0("slingPseudotime_1")]]
CD8_T$sling_pseudotime = sce_slingshot2$sling_pseudotime

#仅分析在轨迹拟时中的细胞，去除NA
sce_slingshot2_l1 = sce_slingshot2[,!is.na(sce_slingshot2$sling_pseudotime)]
seur = CD8_T[,!is.na(CD8_T$sling_pseudotime)]

#这里分析比较慢，会花费数小时时间。
sce_slingshot2_l1 <- fitGAM(counts(sce_slingshot2_l1), 
                            cellWeights = rep(1, ncol(sce_slingshot2_l1)), 
                            pseudotime = sce_slingshot2_l1$sling_pseudotime)

ATres <- associationTest(sce_slingshot2_l1)
association_test_tab = as_tibble(cbind(gene = rownames(ATres), ATres))




#基因随拟时的分析完成后，我们需要利用Seurat对象和上面的分析结果
#构建拟时热图matrix，一步一步做很复杂，这里我们直接提供一个函数，用来获取matrix
#至于分散的步骤可自行查看函数

slingshot_for_plotMatrix <- function(seurat_obj,#添加了拟时并去除NA的seurat对象
                                     n_bins,#拟时您需要分割的区间，将相似的拟时区间合并，这类似于我们monocle3中的方式
                                     min_exp)#需要过滤的基因（这里是表达scale后的数值,选择>0的数值)
{
  
  #提取metadata，并按照拟时先后排序
  seurat_meta = seurat_obj@meta.data
  seurat_meta = as_tibble(cbind(cell.id = as.character(rownames(seurat_meta)), seurat_meta))
  seurat_meta = seurat_meta[order(seurat_meta$sling_pseudotime),]
  
  #排序好的cell id
  pl_cells = as.character(seurat_meta$cell.id)
  
  #提取表达矩阵,并将cell id的exp排序与前面排序好的cell id一致
  exp = seurat_obj@assays$RNA@data
  exp = exp[,colnames(exp) %in% pl_cells]
  expr_mat = exp[,order(match(colnames(exp), pl_cells))]
  
  expr_mat = as.matrix(expr_mat[rownames(expr_mat) %in% association_test_tab$gene,])
  
  clust_expr_mat = matrix(nrow = nrow(expr_mat), 
                          ncol = n_bins, dimnames = list(rownames(expr_mat), 1:n_bins))
  
  max_pseudotime = max(seurat_meta$sling_pseudotime)
  pseudotime_bin_size = max_pseudotime/n_bins
  
  pseudotime_cluster_stat = NULL
  seurat_obj$pseudotime_bin = NA_integer_
  
  
  for (i in 1 : n_bins){
    
    bin_cells = seurat_meta$cell.id[(seurat_meta$sling_pseudotime > (i-1)*pseudotime_bin_size & 
                                       seurat_meta$sling_pseudotime <= i*pseudotime_bin_size)]
    
    
    seurat_obj$pseudotime_bin[colnames(seurat_obj) %in% bin_cells] = i
    
    #计算基因平均表达量
    if (length(bin_cells)>10){
      m2 = expr_mat[,colnames(expr_mat) %in% bin_cells]
      clust_expr_mat[,i] = apply(m2, 1, mean, na.rm = TRUE)
    }
    
  }
  
  #数据缩放一下，为了更好的展现热图，并删除低表达基因
  mm1 = clust_expr_mat - apply(clust_expr_mat, 1, mean, na.rm = TRUE)
  mm2 = mm1[apply(abs(mm1),1, max, na.rm = TRUE)>min_exp,]
  
  return(mm2)
  
}


#一键获取拟时热图矩阵
mm = slingshot_for_plotMatrix(seurat_obj = seur, n_bins = 20, min_exp = 0.2)

#pheatmap作图

max_range = max(range(is.finite(mm)))
lim = c(-max_range, max_range)

library(pheatmap)
heatmap1 = pheatmap(mm, show_rownames=F, cluster_rows = TRUE,
                    cluster_cols = FALSE, show_colnames = FALSE, 
                    clustering_distance_rows = "euclidean",
                    clustering_method = "ward.D2",
                    treeheight_row = 50,
                    cutree_rows = 5, #这里聚类module数自行选择
                    color = colorRampPalette(c("blue", "white", "red"))(250),
                    breaks = seq(lim[1]/4, lim[2]/4, length.out = 251),
                    border_color = NA)

#提取上面的热图聚类信息，提取module基因，并重新作图
tree_module = cutree(heatmap1$tree_row, k=5)
tree_module = tibble(gene = names(tree_module), module = as.character(tree_module))
tree_module = tree_module[heatmap1$tree_row$order,]
#这里得到的tree_module后续可进行富集分析

#后续热图的修饰和GO富集添加等参考monocle2





