



#############harmony处理整个队列新流程
#######包
library(harmony)
library(devtools)
library(Seurat)
library(tidyverse)
library(patchwork)
#rm(list=ls())

#####加载数据
#load("BLCA_combine_merge.Rdata")
sce<-BLCA_combine13
rm(BLCA_combine13)

############################################################QC
#############step1   QC
#人和鼠的基因名字稍微不一样 
mito_genes=rownames(sce)[grep("^MT-", rownames(sce))] #把含有MT开头，也就是含有线粒体基因的行所在的基因名取出来
mito_genes #13个线粒体基因
sce=PercentageFeatureSet(sce, "^MT-", col.name = "percent_mito")
fivenum(sce@meta.data$percent_mito)#fivenum():返回五个数据：最小值、下四分位数、中位数、上四分位数、最大值。

#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA")

###nFea_nCount
VlnPlot(sce, group.by = "sample_id", features = feats, pt.size = 0.0001, ncol = 1) + 
  NoLegend()+ theme(axis.text.x=element_text(vjust=1,size=20,face = "bold"))
ggsave("nFea_nCount_before_QC.pdf",height = 10,width = 10)

####mito
feats <- c("percent_mito")#感觉线粒体基因比例有点高
VlnPlot(sce, group.by = "sample_id", features = feats, pt.size = 0.01, ncol = 1, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
ggsave("mito_beforeQC.pdf",height = 6,width = 10)

###nFeature_nCount_scatter
FeatureScatter(sce, "nCount_RNA", "nFeature_RNA", group.by = "sample_id", pt.size = 0.5)
ggsave("nFeature_nCount_scatter_beforeQC.pdf",height = 5,width = 7)

#根据上述指标，过滤低质量细胞/基因
#过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
##文章nature cancer肾癌的数据nFeature_RNA500-5000，nCount_RNA<50000,
selected_c <- WhichCells(sce, expression = nFeature_RNA > 500 & 
                           nFeature_RNA < 5000 & nCount_RNA > 1000 & nCount_RNA < 50000)###
selected_f <- rownames(sce)[Matrix::rowSums(sce@assays$RNA@counts > 0 ) > 10]
sce.filt <- subset(sce, features = selected_f, cells = selected_c)

#sce.filt <- subset(sce, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 1000 & nCount_RNA < 50000)
#x1<-sce@meta.data$nFeature_RNA
#class(x1)
#length(x1[x1>550])  

#过滤前/过滤后
dim(sce) 
dim(sce.filt) 

#过滤指标2:线粒体/核糖体基因比例(根据上面的violin图),
selected_mito <- WhichCells(sce.filt, expression = percent_mito < 20)#文章写了超过20%的去除
length(selected_mito)
sce.filt <- subset(sce.filt, cells = selected_mito)

dim(sce) #before
dim(sce.filt)#after
table(sce$sample_id)
table(sce.filt$sample_id) #看一下每个病人多少

#可视化过滤后的情况
feats <- c("nFeature_RNA", "nCount_RNA")
VlnPlot(sce.filt, group.by = "sample_id", features = feats, pt.size = 0.1, ncol = 1,raster=FALSE) + 
  NoLegend()
ggsave("nFeat_nCount_afterQC.pdf",height = 10,width = 10)

#指控后线粒体基因情况
feats <- c("percent_mito")
VlnPlot(sce.filt, group.by = "sample_id", features = feats, pt.size = 0.1, ncol = 1) + 
  NoLegend()
ggsave("mito_afterQC.pdf",height = 6,width = 10)
#保存数据
sce=sce.filt 
rm(sce.filt)
save(sce, file = "sce.filt.RData")


###############删除特定基因
####删除MT基因
mito_genes=rownames(sce)[grep("^MT-", rownames(sce))] #把含有MT开头，也就是含有线粒体基因的行所在的基因名取出来
mito_genes #13个线粒体基因
keep = c(!rownames(sce) %in% c(mito_genes))
scedeMT <- subset(x = sce,features =c(1:(dim(sce)[1]))[keep])
mito_genes=rownames(scedeMT)[grep("^MT-", rownames(scedeMT))] #把含有MT开头，也就是含有线粒体基因的行所在的基因名取出来
mito_genes #13个线粒体基因
###删除核糖体基因
ribo_genes=rownames(scedeMT)[grep("^Rp[sl]", rownames(scedeMT),ignore.case = T)]
ribo_genes
keep = c(!rownames(scedeMT) %in% c(ribo_genes))
scedeMT <- subset(x = scedeMT,features =c(1:(dim(scedeMT)[1]))[keep])
rm(sce)
BLCA_combine_sce.filt_degene<-scedeMT
rm(scedeMT)
save(BLCA_combine_sce.filt_degene,file="BLCA_combine_sce.filt_degene.Rdata")


########################################################step2
library(harmony)
library(devtools)
library(Seurat)
library(tidyverse)
library(patchwork)
#rm(list=ls())
load("BLCA_combine_sce.filt_degene.Rdata")

scRNA_harmony<-BLCA_combine_sce.filt_degene
rm(BLCA_combine_sce.filt_degene)
#########harmony整合多样本
DefaultAssay(scRNA_harmony)<-"RNA"
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
##RunHarmony的参数lambda值调小，整合力度变大，默认1，一般在0.5-2调整
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "sample_id")})
###运行结果保存在下面
#scRNA_harmony@reductions$harmony
####访问新的harmony
#scRNA_harmony_Embeddings<-Embeddings(scRNA_harmony,"harmony")


#################################选择PC
# Determine percent of variation associated with each PC
pct <- scRNA_harmony [["pca"]]@stdev / sum( scRNA_harmony [["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
ggsave('PC选择.pdf',width = 13,height = 5)
# Printing out the most variable genes driving PCs
print(x = scRNA_harmony [["pca"]],  dims = 1:25,  nfeatures = 5)




#降维聚类
table(scRNA_harmony$sample_id)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:15)
options(future.globals.maxSize = Inf)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 2)
#save(scRNA_harmony,file="scRNA_harmony_res2.Rdata")


#####看整合结果
######### Visualization the integration
table(scRNA_harmony@meta.data$seurat_clusters)
pp<-DimPlot(scRNA_harmony, reduction = "pca", group.by = "sample_id",raster=FALSE)

p0<-DimPlot(scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.2",raster=FALSE,label = T,label.size = 4)
ggsave('umap_by_cluster.pdf',width = 10,height = 10,p0)

p3<-DimPlot(scRNA_harmony, reduction = "umap", group.by = "sample_id",raster=FALSE)
ggsave('umap_by_sample_id.pdf',width = 10,height = 10,p3)

p4<-DimPlot(scRNA_harmony, reduction = "umap", group.by = "patient_id",raster=FALSE)
ggsave('umap_by_patient_id.pdf',width = 10,height = 10,p4)

ggsave('umap_pca_harmony.pdf',width = 20,height = 10,pp+p3)

p0<-DimPlot(scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.2",raster=FALSE)
ggsave('umap_all_id.pdf',width = 30,height = 10,p0+p3+p4)

#######差异基因
# find markers for every cluster compared to all remaining cells, report only the positive
harmony.markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(harmony.markers,"scRNA_harmony.markers.csv")




