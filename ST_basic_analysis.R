

lapply(c("Normol_flow"), FUN = dir.create)
options(future.globals.maxSize = 100 * 1024^3)

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)

data_dir <-"D:\\TLSandB\\lianchuan_ST\\outs"
list.files(data_dir)
file_name <-"filtered_feature_bc_matrix.h5"
#读入数据并创建Seurat对象；
slide1<- Load10X_Spatial(data.dir = data_dir, filename = file_name,slice ="slide1")
slide1@project.name <-"slide1"
Idents(slide1) <-"slide1"
slide1$orig.ident <-"slide1"

# 使用小提琴图可视化统计指标
VlnPlot(slide1, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
ggsave('Normol_flow/ncount_nfeature.pdf', width = 5, height = 5)

#找不能用的spot####
# 计算线粒体基因含量
slide1<- PercentageFeatureSet(slide1, "^mt-", col.name = "percent_mito")
# 计算血红蛋白基因含量
slide1<- PercentageFeatureSet(slide1, "^Hb.*-", col.name = "percent_hb")
#table(slide1$percent_hb)
#colnames(slide1)[which(!is.na(slide1$percent_hb))]
#colnames(slide1)[which(is.na(slide1$percent_hb))]
#slide1<- subset(slide1, cells = colnames(slide1)[which(!is.na(slide1$percent_hb))])

#数据预处理
#原始HE图
library(ggplot2)
library(ggspatial)
SpatialFeaturePlot(slide1, features = "nCount_Spatial", alpha = 1, pt.size.factor = 0) +
  theme(legend.position = "right")
ggsave('Normol_flow/原始HE图.pdf', width = 5, height = 5)

p1 <- VlnPlot(slide1, features ="nCount_Spatial",
              pt.size = 0,cols ="tomato") + NoLegend()
p2 <- SpatialFeaturePlot(slide1, features ="nCount_Spatial",alpha = c(0.1, 1), pt.size.factor = 4) +
  theme(legend.position ="right")
#p1+p2
p1 | p2
ggsave('Normol_flow/ncount_spatial.pdf', width = 5, height = 5)

p3 <- VlnPlot(slide1, features ="nFeature_Spatial",
              pt.size = 0,cols ="tomato") + NoLegend()
p4 <- FeatureScatter(slide1, feature1 ="nCount_Spatial",
                     feature2 ="nFeature_Spatial")+ NoLegend()
p3 | p4
ggsave('Normol_flow/nFeature_pca.pdf', width = 5, height = 5)

#用SCTransform()对数据进行标准化, 同时检测高变基因, 输出结果储存在 SCT assay中；
#耗时约：2min
slide1<- SCTransform(slide1, assay = "Spatial", verbose = FALSE)
#提取表达量变变化最高的10个基因;
top10 <- head(VariableFeatures(slide1),10)
top10
#[1] "IGKC"   "IGHG1"  "IGHG3"  "HMGCS2" "IGLC1"  "IGHA1"  "IGHM"   "MYH11"  "DES"    "KRT19" 

#类似常规转录组，我们也可以展示高变基因；
p5 <- VariableFeaturePlot(slide1,cols = c("gray60", "red"))+NoLegend()
p5
p6 <- LabelPoints(plot = p5,points = top10, repel = TRUE,xnudge=0,ynudge=0)
p6
ggsave('Normol_flow/高变基因前10.pdf', width = 5, height = 5)

#展示特定基因的表达水平#####
#SpatialFeaturePlot(slide1, features = c("DES", "WFDC2","MYH11","IGHG1","PLA2G2A"))
SpatialFeaturePlot(slide1, features = top10,ncol = 5,pt.size.factor = 2,alpha = c(0.1, 1.5))+geom_point(stroke = 0)
ggsave('Normol_flow/高变基因前10Spatial图.pdf', width = 15, height = 7)
#也可以对“点”进行个性化设置：点的大小和透明度设置；
#p7 <- SpatialFeaturePlot(slide1, features ="JCHAIN", pt.size.factor = 4,stroke = 0.0)
#p8 <- SpatialFeaturePlot(slide1, features ="COL1", pt.size.factor = 4,stroke = 0.0)
#p8 <- SpatialFeaturePlot(slide1, features ="COL1", pt.size.factor = 4,alpha = c(0.1, 1.5))#尝试将其设置为alphac（0.1，1），以降低具有较低表达式的点的透明度
#p7 + p8
#SpatialFeaturePlot(slide1, features ="IGHG1",pt.size.factor = 4, alpha = c(0.1, 1.5))
#目的基因表达情况精修####

###########################已经画好了,开始画图
P<-SpatialFeaturePlot(slide1, features ="PECAM1")
img <- P$data[,1:3]
#图2 整个大背景黑色
ggplot(img, aes(x=imagecol, y=600-imagerow,color = PECAM1)) +
  geom_point(size = 1.5)+
  scale_color_gradient2(low="black",high="white",mid="red",midpoint = median(img$PECAM1)
  )+
  theme_minimal() + #去掉背景
  theme(panel.grid=element_blank(),aspect.ratio = 1)+ #去掉网格线
  theme(axis.text.x=element_blank(), #去掉y轴坐标
        axis.text.y=element_blank()) +#去掉x轴坐标
  theme(plot.background = element_rect(fill = "black" ))+ 
  #theme(panel.background = element_rect(fill = 'black', color = 'black'))
  theme(legend.background = element_rect(fill = 'white', color = 'white'),
        legend.direction = "horizontal",#图例方向
        legend.position = "top"#图例位置
  )
#图3 得去掉x轴y轴坐标

#colnames(img)[3] <- "gene"
ggplot(img, aes(x=imagecol, y=600-imagerow,color = PECAM1)) +
  geom_point(size = 1.5)+
  scale_color_gradient2(low="#07091E",
                        high="white",
                        mid="red",
                        midpoint = median(img$PECAM1))+
  theme_minimal() + #去掉背景
  theme(panel.grid=element_blank(),aspect.ratio = 1)+ #去掉网格线
  theme(axis.text.x=element_blank(), #去掉y轴坐标
        axis.text.y=element_blank()) +#去掉x轴坐标
  #theme(plot.background = element_rect(fill = "black" ))+ 
  theme(panel.background = element_rect(fill = 'black', color = 'black'))+
  labs(x="",y = "",title = "")+theme(aspect.ratio = 1)+
  theme(legend.background = element_rect(fill = 'white', color = 'white'),
        legend.direction = "horizontal",#图例方向
        legend.position = "top"#图例位置
  )
#图4配色
ggplot(img, aes(x=imagecol, y=600-imagerow,color = PECAM1)) +
  geom_point(size = 1.5)+
  scale_color_gradientn(colors = c("#040103","#0D0B1D","#6F2663","#CA623C","#EAE69A","#F9F9B7"))+
  theme_minimal() + #去掉背景
  theme(panel.grid=element_blank(),aspect.ratio = 1)+ #去掉网格线
  theme(axis.text.x=element_blank(), #去掉y轴坐标
        axis.text.y=element_blank()) +#去掉x轴坐标
  #theme(plot.background = element_rect(fill = "black" ))+ 
  theme(panel.background = element_rect(fill = 'black', color = 'black'))+
  labs(x="",y = "",title = "")+theme(aspect.ratio = 1)+
  theme(legend.background = element_rect(fill = 'white', color = 'white'),
        legend.direction = "horizontal",#图例方向
        legend.position = "top"#图例位置
  )
ggsave("Normol_flow/PECAM1_expr.pdf",height = 6,width = 6)

#接下来是常规单细胞转录组分析流程：PCA降维、UMAP聚类；####
slide1<- RunPCA(slide1, assay = "SCT", verbose = FALSE)
slide1<- FindNeighbors(slide1, reduction = "pca", dims = 1:30)
slide1<- FindClusters(slide1, verbose = FALSE)
slide1<- RunUMAP(slide1, reduction = "pca", dims = 1:30)

#绘制UMAP分群图；
p9 <- DimPlot(slide1, reduction = "umap", label = TRUE)+theme(aspect.ratio = 1)
p9
ggsave("Normol_flow/slide1_umap.pdf",width = 6,height = 6)

#在切片图像上映射分群信息；
p10 <- SpatialDimPlot(slide1, label = TRUE, label.size = 3,pt.size.factor = 4.5)
p10
ggsave("Normol_flow/spatial_cluster.pdf",width = 6,height = 6)
p9+p10

####对所有cluster进行差异分析；
de_markers <- FindAllMarkers(slide1)
write.csv(de_markers,"Normol_flow/de_markers.csv")

########画dotplot图#######
genes_to_check = c(
  "CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC2", # T
  "GNLY", "NKG7", "KLRD1", "KLRB1", "PRF1","NCAM1", # NK
  "IGHG1", "IGKC", "JCHAIN", "IGHG3", "MZB1", # Plasma
  "CD79A", "CD37", "MS4A1", "BANK1", "TNFRSF13C",# B cell
  "CSF3R", "S100A8", "CXCL8", "MNDA", "S100A8", "S100A9", # Neutrophil
  "LYZ", "C1QB", "C1QA", "CD14", "FCGR3A", # Myeloid
  "MS4A2", "TPSB2", "HPGDS", "GATA2", "TPSAB1", "CPA3", "LTC4S", "RGS13", # Mast
  "DCN", "COL1A2", "COL1A1", "COL3A1", "RGS5", "TAGLN", "ACTA2", "MYL9", "TPM2", # Fibroblast
  "PECAM1", "VWF", "ENG", "PTPRB", "KDR", "PLVAP", # Endothelium
  "EPCAM","KRT19","KRT13", "CLDN4","FXYD3","S100P" # Epithelial
)
DotPlot(slide1, features = unique(genes_to_check),
        assay='SCT')  + coord_flip()+scale_color_gradient2(low="black",high="red",mid="white",midpoint = 0)
ggsave("Normol_flow/cluster_dotplot_genemarker.pdf",width = 10,height = 10)

######保存下数据
slide1<-slide1
save(slide1,file="slide1_normal_flow.Rdata")
load("slide1_normal_flow.Rdata")
data<-slide1@assays$SpaCET_TNKmajor
data1<-as.data.frame(data@counts)


table(slide1@meta.data$Mature_two)
# 设定Seurat对象的身份为"area"并绘制空间图
Idents(slide1) <- "Mature_two"
plot <- SpatialDimPlot(slide1, label = TRUE, label.size = 3, crop = TRUE, pt.size.factor = 4.5, stroke = 0.0)
plot
# 保存图像
ggsave("maturation/picture/cluster_粗分三群.pdf", plot, width = 5, height = 5)



