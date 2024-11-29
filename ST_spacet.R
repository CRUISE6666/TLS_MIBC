





load("slide1_normal_flow.Rdata")
#######################SpaCET###################去卷积######
lapply(c("SpaCET_TNKmajor"), FUN = dir.create)
options(future.globals.maxSize = 100 * 1024^3)
# install.packages("devtools")
#devtools::install_github("data2intelligence/SpaCET",force = T)
library(SpaCET)
library(ggplot2)
library(Seurat)
#分到足够细的亚群，然后进行mia分析 
# 设置文件路径
file_path <- "C:/Users/linjiaxing/Desktop/TLSandB/BLCA_TLS/BLCA_scRNA/all_subtype_harmony_res2_markers_de.Rdata"
# 加载RData文件
load(file_path)
table(scRNA_harmony@meta.data$patient_id)
#scRNA<-subset(x=scRNA_harmony, subset = (patient_id == "slide1"))
scRNA<-scRNA_harmony
stRNA <- slide1
#构建SpaCET对象
counts <- as.matrix(stRNA@assays$Spatial@counts)
counts[1:6,1:5]
#提取坐标信息
# img <- P1T@images$image_P11_T@coordinates[,4:5]
# img <- img/27.8933
spotCoordinates <- as.data.frame(stRNA@images$T4463206@coordinates)
spotCoordinates <- spotCoordinates[,2:3]
#spotCoordinates <- spotCoordinates/27.8933
colnames(spotCoordinates) <- c("X","Y")
spotCoordinates$X <- -spotCoordinates$X
#spotCoordinates$Y <- -spotCoordinates$Y
#spotCoordinates$Y <- 600-spotCoordinates$Y
spotCoordinates[1:5,]

SpaCET_obj <- create.SpaCET.object(
  counts=counts,
  spotCoordinates=spotCoordinates,
  imagePath=NA,
  platform = "P1"
)
str(SpaCET_obj)
# calculate the QC metrics
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj)

# plot the QC metrics
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "QualityControl", 
  spatialFeatures=c("UMI","Gene")
)


#提取单细胞表达矩阵
sc_counts <- as.matrix(scRNA@assays$RNA@counts)
sc_counts[1:6,1:5]

table(scRNA$celltype)
table(scRNA$celltype_major)
table(scRNA$celltype_subtype)
#提取细胞信息

sc_annotation <- scRNA@meta.data
sc_annotation <- data.frame(cellID=rownames(sc_annotation),bio_celltype=sc_annotation[,"celltype_subtype"])
rownames(sc_annotation) <- sc_annotation[,1]
####这个所有细胞都有
sc_lineageTree <- list(TNK=c("CD4_exhausted", "CD4_Memory","CD4_Naive", "CD4_proliferation","CD4_Treg","CD8_cytotoxic", 
                             "CD8_exhausted","CD8_NKlike", "CD8_proliferation","NK_CD16n", "NK_CD16p"),
                       B_cell=c("memory_B", "Naïve_B"),
                       Plasma=c("GC_B", "Pan_plasma","SPC"),
                       Macrophage=c("Macro_C1QC", "Macro_MERTK","Macro_MKI67", "Macro_SPP1"),
                       Monocyte=c("Mono_CD14", "Mono_CD14CD16","Mono_CD16"),
                       Dendritic=c("CDC1_CLEC9A", "CDC2_CD1C"),
                       Mast=c("Mast_KIT"),
                       Neutrophil=c("Neutrophil"),
                       Epithelial=c("Epithelial"),
                       Endothelial=c("ACKR1_EC", "PLVAP_EC","neg_EC"),
                       Fibroblast=c("iCAF", "mCAF")
)
SpaCET_obj <- SpaCET.deconvolution.matched.scRNAseq(
  SpaCET_obj, 
  sc_counts=sc_counts, 
  sc_annotation=sc_annotation, 
  sc_lineageTree=sc_lineageTree, 
  coreNo=1
)

####可视化所有亚群
# show the spatial distribution of all cell types.
p1<-SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures="All", 
  sameScaleForFraction = TRUE,
  pointSize = 0.1, 
  nrow=6
)
p1
ggsave("SpaCET_TNKmajor/SpaCET_卷积结果_分细胞出图.pdf", p1,width = 12, height = 10)


#######映射回去#####
data<-SpaCET_obj@results$deconvolution$propMat
slide1[["SpaCET_TNKmajor"]] <- CreateAssayObject(data)
DefaultAssay(slide1) <- "SpaCET_TNKmajor"
celltypes <- rownames(slide1)
celltypes

######
saveRDS(SpaCET_obj, "SpaCET_TNKmajor/SpaCET_TNKmajor_slide1.rds")
save(slide1,file="slide1_normal_flow.Rdata")


########对每个基因画图##############
dir.create("SpaCET_TNKmajor/celltype", showWarnings = FALSE)
# 定义要绘制的基因列表
genes <- celltypes
# 循环绘制多个基因的图并保存
for (gene in genes) {
  # 使用 SpatialFeaturePlot 绘制单个基因的图
  P <- SpatialFeaturePlot(slide1, features = gene)
  img <- P$data[, 1:3]
  plot <- ggplot(img, aes(x = imagecol, y = 600 - imagerow, color = .data[[gene]])) +
    geom_point(size = 1.5) +
    scale_color_gradientn(colors = c("#040103", "#0D0B1D","#49006a", "#ae017e", "#f768a1", "#fcc5c0", "#fff7f3")) +
    theme_minimal() +
    theme(panel.grid = element_blank(), aspect.ratio = 1) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
    theme(panel.background = element_rect(fill = 'black', color = 'black')) +
    labs(x = "", y = "", title = "") +
    theme(aspect.ratio = 1) +
    theme(legend.background = element_rect(fill = 'white', color = 'white'),
          legend.direction = "horizontal",
          legend.position = "top")
  
  # 保存图像到文件夹中
  ggsave(paste0("SpaCET_TNKmajor/celltype/", gene, "_expr.pdf"), plot, height = 6, width = 6)
}

######针对TLS分组画小提琴图和箱型图
library(ggplot2)
library(ggpubr)
###########提取打分以及分组
dataMCP<-slide1[["SpaCET_TNKmajor"]]###修改提取的打分
data<-as.data.frame(dataMCP@counts)
data<-t(data)
data<-as.data.frame(data)
colnames(data) <- gsub("-", "_", colnames(data))
data$id<-row.names(data)
dataarea<-as.data.frame(slide1@meta.data)####可以只提取想要的分组，也可以都提取
dataarea$id<-row.names(dataarea)
scoredata<-inner_join(data,dataarea,by="id")
colnames(scoredata) <- gsub("-", "_", colnames(scoredata))


lapply(c("SpaCET_TNKmajor/TLS_vlin","SpaCET_TNKmajor/TLS_box"), FUN = dir.create)
options(future.globals.maxSize = 100 * 1024^3)
####
colnames(scoredata)
table(scoredata$area)
#########################
#cell_types <- c("T_cells", "Cytotoxic_lymphocytes", "B_lineage", "NK_cells",
#               "Monocytic_lineage", "Myeloid_dendritic_cells", "Neutrophils", 
#               "Endothelial_cells", "Fibroblasts")###修改想要画图的细胞
cell_types<-colnames(data[,-ncol(data)])###全部细胞,减去id名
my_comparisons <- list(c("noTLS_area", "TLS_area"))
for (cell_type in cell_types) {
  # 为每种细胞类型生成箱线图
  p_box <- ggboxplot(scoredata, x = "area", y = cell_type, color = "area") + ##修改分组名称
    stat_compare_means(comparisons = my_comparisons, method = "t.test",
                       label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("***", "**", "*", "ns")))
  
  ggsave(paste0("SpaCET_TNKmajor/TLS_box/", cell_type, "_box_合起来打分.pdf"), plot = p_box, width = 4, height = 4)
  
  # 为每种细胞类型生成小提琴图
  p_violin <- ggplot(scoredata, aes(x = area, y = scoredata[[cell_type]], fill = area)) + ###修改分组名称
    geom_violin() +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                       label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("***", "**", "*", "ns"))) +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.2) +
    scale_fill_manual(values = c(noTLS_area = "#66a61e", TLS_area = "#a6761d")) +  ###有分组，修改对应颜色
    theme(panel.grid.major = element_line(colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.minor = element_blank()) +
    theme_classic() +
    labs(y = cell_type)  # 添加y轴标签
  
  ggsave(paste0("SpaCET_TNKmajor/TLS_vlin/", cell_type, "_vlin_合起来打分.pdf"), plot = p_violin, width = 5, height = 4)
}




######################取占比最高的点作为点代表整个点
dataMCP<-slide1[["SpaCET_TNKmajor"]]###修改提取的打分
data<-as.data.frame(dataMCP@counts)
test<-data[1:6,1:6]
data1<-data[-c(1:8,12:13),]
test<-data1[1:6,1:6]
test

# 我们将使用which.max函数来确定每一列中最大值的位置
# 找到每一列最大值的行名，即含量最高的细胞类型
max_cell_types <- apply(data1, 2, function(column) rownames(data1)[which.max(column)])
# 创建一个新的数据框来存储结果
result <- data.frame(Cell_Type = max_cell_types)
# 打印结果数据框
print(result)
colnames(result)<-"SpaCET_TNKmajor_celltype"
table(result$Cell_Type)

###合并信息到正常的seurat对象
slide1 <- AddMetaData(slide1, metadata = result)
save(slide1,file="slide1_normal_flow.Rdata")
######################
# 设定Seurat对象的身份为"area"并绘制空间图
table(slide1$SpaCET_TNKmajor_celltype)
Idents(slide1) <- "SpaCET_TNKmajor_celltype"
plot <- SpatialDimPlot(slide1, label = TRUE, label.size = 3, crop = TRUE, pt.size.factor = 4.25, stroke = 0.0)
plot
# 保存图像
ggsave("SpaCET_TNKmajor/SpaCET_TNKmajor_celltype_map.pdf", plot, width = 7, height = 5)

####为限制颜色画图####
p1 <- SpatialDimPlot(slide1, group.by = "SpaCET_TNKmajor_celltype", label.size = 3, pt.size.factor = 4.4) +
  plot_layout(guides = "collect") + plot_layout() +
  ggsci::scale_fill_igv()
ggsave("SpaCET_TNKmajor/SpaCET_TNKmajor_celltype_map_bayes.pdf", p1, width = 7, height = 5)


#####使用ggplot绘制散点图定制颜色####
library(ggsci)
library(ggplot2)
# 创建一个示例数据集
df <- data.frame(value = 1:9, group = as.factor(1:9))
# 生成图表
p <- ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_igv()
# 显示图表
print(p)
# 尝试提取颜色值
igv_colors <- ggsci::pal_igv()(10)
print(igv_colors)
##
p1 <- SpatialDimPlot(slide1, group.by = "SpaCET_TNKmajor_celltype", label.size = 3, pt.size.factor = 4.4,stroke = 0) +
  plot_layout(guides = "collect") + plot_layout() +
  scale_fill_manual(values=c("B-cell" = "#5050FFFF",
                             "Dendritic" = "#CE3D32FF",
                             "Endothelial" = "#749B58FF",
                             "Epithelial" = "#F0E685FF",
                             "Fibroblast" = "#466983FF",
                             "Macrophage" = "#BA6338FF",
                             "Monocyte" = "#5DB1DDFF",
                             "Plasma" = "#802268FF",
                             "TNK" = "#6BD76BFF"
  ))
p1
ggsave("SpaCET_TNKmajor/SpaCET_TNKmajor_celltype_map_bayes.pdf", p1, width = 7, height = 5)


data<-slide1@assays$SpaCET_TNKmajor@counts@Dimnames
write.csv(data,"SpaCET_TNKmajor/data.csv")
#########################################################画图###############################
library(ggplot2)
library(ggpubr)
###########提取打分以及分组
DefaultAssay(slide1) <- "SpaCET_TNKmajor"
dataSpaCET_TNKmajor<-slide1[["SpaCET_TNKmajor"]]###修改提取的打分
data<-as.data.frame(dataSpaCET_TNKmajor@counts)
data<-t(data)
data<-as.data.frame(data)
colnames(data) <- gsub("-", "_", colnames(data))
data$id<-row.names(data)
dataCottrazm_TLS<-as.data.frame(slide1@meta.data)####可以只提取想要的分组，也可以都提取
dataCottrazm_TLS$id<-row.names(dataCottrazm_TLS)
scoredata<-inner_join(data,dataCottrazm_TLS,by="id")
colnames(scoredata) <- gsub("-", "_", colnames(scoredata))

######针对Cottrazm_TLS分组画小提琴图和箱型图
table(slide1$Cottrazm_TLS)
lapply(c("SpaCET_TNKmajor/Cottrazm_TLS_vlin","SpaCET_TNKmajor/Cottrazm_TLS_box"), FUN = dir.create)
options(future.globals.maxSize = 100 * 1024^3)
####
colnames(scoredata)
table(scoredata$Cottrazm_TLS)
#########################
#cell_types <- c("T_cells", "Cytotoxic_lymphocytes", "B_lineage", "NK_cells",
#               "Monocytic_lineage", "Myeloid_dendritic_cells", "Neutrophils", 
#               "Endothelial_cells", "Fibroblasts")###修改想要画图的细胞
cell_types<-colnames(data[,-ncol(data)])###全部细胞,减去id名
my_comparisons <- list(c("Bdy","Mal"),c("Bdy","nMal"),c("Bdy","TLS"),c("Mal","nMal"),c("Mal","TLS"),c("nMal","TLS"))
for (cell_type in cell_types) {
  # 为每种细胞类型生成箱线图
  p_box <- ggboxplot(scoredata, x = "Cottrazm_TLS", y = cell_type, color = "Cottrazm_TLS") + ##修改分组名称
    stat_compare_means(comparisons = my_comparisons, method = "t.test",
                       label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("***", "**", "*", "ns")))
  
  ggsave(paste0("SpaCET_TNKmajor/Cottrazm_TLS_box/", cell_type, "_box_合起来打分.pdf"), plot = p_box, width = 4, height = 4)
  
  # 为每种细胞类型生成小提琴图
  p_violin <- ggplot(scoredata, aes(x = Cottrazm_TLS, y = scoredata[[cell_type]], fill = Cottrazm_TLS)) + ###修改分组名称
    geom_violin() +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                       label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("***", "**", "*", "ns"))) +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.2) +
    scale_fill_manual(values = c(Bdy = "#66c2a5", Mal = "#fc8d62", nMal = "#8da0cb", TLS = "#e78ac3")) +  ###有分组，修改对应颜色
    theme(panel.grid.major = element_line(colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.minor = element_blank()) +
    theme_classic() +
    labs(y = cell_type)  # 添加y轴标签
  
  ggsave(paste0("SpaCET_TNKmajor/Cottrazm_TLS_vlin/", cell_type, "_vlin_合起来打分.pdf"), plot = p_violin, width = 5, height = 4)
}







