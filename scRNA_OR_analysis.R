






###################
library("sscVis")
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")
library("dplyr")
library("tidyr")

table(scRNA_harmony$celltype_subtype)
table(scRNA_harmony$group_FL)
#加载函数及数据
#mouse_data <- readRDS("./mouse_data.rds")
source("all/OR/test_function.R")
source("all/OR/draw_analysis.R")
#数据分析
A <- do.tissueDist(cellInfo.tb = scRNA_harmony@meta.data,#这里的输入文件需要的是整理好的有分组和细胞类型的metadata文件
                   out.prefix = "./ORplot",#设定导出的文件名，自己设置
                   pdf.width = 5,#热图的宽设置
                   pdf.height = 8,#热图的高度设置
                   verbose=1,#设置为1文件以list形式存储
                   meta.cluster = 'celltype_subtype',#这里是细胞类型，也可以是seurat_clusters，名称没有要求，就是你细胞类型的列名
                   loc = 'group_FL',#这里就是分组，metadata中分组的列名，至于命名没有要求
                   z.hi=15)#热图legend最大值，根据实际情况自己设置

#查看并保存文件
A$OR.dist.mtx #做热图数据，OR值
A$p.dist.tb #p值
A$OR.dist.tb #OR值
A$count.dist.melt.ext.tb#组合表，adjust-pvalue等


###################################################################
#自己做图
color_palette <- c(
  "#f7fcfd",
  "#e0ecf4",
  "#bfd3e6",
  "#9ebcda",
  "#8c96c6",
  "#8c6bb1",
  "#88419d",
  "#810f7c",
  "#4d004b"
)
data <- A$count.dist.melt.ext.tb
write.csv(data, file = 'Bplasma/scRNA_harmony_data_celltype_subtype.csv')
library(ggplot2)
library(RColorBrewer)
data <- read.csv("Bplasma/scRNA_harmony_data_celltype_subtype.csv", header = T)
data$cid <- factor(data$cid, levels = c("nTLS", "FL"))

ggplot(data, aes(cid, rid)) + 
  geom_tile(aes(fill = OR), colour = "black", size = 0.6)+
  scale_fill_gradientn(name='OR',
                       colours= color_palette )+
  theme_minimal() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.y = element_text(size = 10,color = 'black')) + 
  scale_y_discrete(position = "right")
ggsave("Bplasma/scRNA_harmony_OR_热图_celltype_subtype.pdf",width =4,height = 4)

table(scRNA_harmony$patient_id)
# 加载必要的库
library(ggplot2)
library(RColorBrewer)
# 读取数据
data <- read.csv("Bplasma/scRNA_harmony_data_celltype_subtype.csv", header = T)
data$cid <- factor(data$cid, levels = c("nTLS", "FL"))
# 绘制热图并显示 OR 值
ggplot(data, aes(cid, rid)) + 
  geom_tile(aes(fill = OR), colour = "black", size = 0.6) +
  scale_fill_gradientn(name = 'OR',
                       colours = color_palette) +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = 'black')) + 
  scale_y_discrete(position = "right") +
  geom_text(aes(label = sprintf("%.2f", OR)), size = 3, color = "black") # 显示 OR 值

# 保存热图
ggsave("Bplasma/scRNA_harmony_OR_热图_celltype_subtype.pdf", width = 4, height = 4)
