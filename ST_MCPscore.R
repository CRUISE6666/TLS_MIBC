


######MCP打分################
lapply(c("MCP_score"), FUN = dir.create)
options(future.globals.maxSize = 100 * 1024^3)
##########MCP计算######
library(Seurat)
library(dplyr)
library(MCPcounter)
# 首先，为了适应MCPcounter，我们需要确保数据是在基因水平上进行注释的，且使用的是HUGO gene symbols。
# 提取Seurat对象的表达矩阵
expr_matrix <- slide1[["SCT"]]@data

# 使用MCPcounter估算每个细胞的MCP得分
probesets<-read.table("MCP_score/macro_score/probesets.txt",sep="\t",stringsAsFactors=FALSE,colClasses="character")
genes<-read.table("MCP_score/macro_score/genes.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
head(genes)
genes
###去隔壁代码加下基因ID
library(biomaRt)
# 选择你的mart和数据集。通常对于人类基因，你会使用ENSEMBL_MART_ENSEMBL和hsapiens_gene_ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 你的HUGO symbols向量
gene_symbols <- genes$`HUGO symbols`
# 查询Ensembl数据库来获取ENTREZID和ENSEMBL_ID
genes_info <- getBM(attributes = c('hgnc_symbol', 'entrezgene_id', 'ensembl_gene_id'),
                    filters = 'hgnc_symbol',
                    values = gene_symbols,
                    mart = mart)
# 由于有些基因可能没有对应的ENTREZ ID或ENSEMBL ID，结果可能会有NA
# 合并这些信息回到你的原始数据框，基于HUGO symbols
genes <- merge(genes, genes_info, by.x = 'HUGO symbols', by.y = 'hgnc_symbol', all.x = TRUE)
# 重新命名列以匹配你的要求
names(genes)[names(genes) == 'entrezgene_id'] <- 'ENTREZID'
names(genes)[names(genes) == 'ensembl_gene_id'] <- 'ENSEMBL_ID'
#定义分析过程中的特征类型，为了避免转换，选用symbol
featuresType <- c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID")[2]
#如果网络不佳proesets和genes也可以不提供
mcp_scores<- MCPcounter.estimate(expr_matrix,featuresType,
                                 probesets=probesets,
                                 genes=genes
)
#mcp_scores <- MCPcounter.estimate(expr_matrix, featuresType = "HUGO_symbols")
# 替换行名中的空格为下划线
rownames(mcp_scores) <- gsub(" ", "_", rownames(mcp_scores))
mcp_scores<-as.data.frame(mcp_scores)
# 将MCP得分加入到Seurat对象的元数据中
slide1[["MCP_macroscore"]] <- NULL
data<-mcp_scores
slide1[["MCP_macroscore"]] <- CreateAssayObject(data)
DefaultAssay(slide1) <- "MCP_macroscore"
celltypes <- rownames(slide1)
celltypes
####保存数据
saveRDS(mcp_scores, "MCP_score/macro_score/MCP_slide1.rds")
save(slide1,file="slide1_normal_flow.Rdata")

########对每个celltype画图##############
dir.create("MCP_score/macro_score/celltype", showWarnings = FALSE)
# 定义要绘制的基因列表
DefaultAssay(slide1) <- "MCP_macroscore"
genes <- celltypes
# 循环绘制多个基因的图并保存
for (gene in genes) {
  # 使用 SpatialFeaturePlot 绘制单个基因的图
  P <- SpatialFeaturePlot(slide1, features = gene)
  img <- P$data[, 1:3]
  plot <- ggplot(img, aes(x = imagecol, y = 600 - imagerow, color = .data[[gene]])) +
    geom_point(size = 1.5) +
    scale_color_gradientn(colors = c("#000000","#081d58", "#253494", "#225ea8", "#1d91c0", "#41b6c4",
                                     "#7fcdbb", "#c7e9b4", "#edf8b1", "#ffffd9")) +
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
  ggsave(paste0("MCP_score/macro_score/celltype/", gene, "_expr.pdf"), plot, height = 6, width = 6)
}



########对每个MCP_TLSscore画图##############
dir.create("MCP_score/celltype/red/less", showWarnings = FALSE)
# 定义要绘制的基因列表
DefaultAssay(slide1) <- "MCP_TLSscore"
celltypes <- rownames(slide1)
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
  ggsave(paste0("MCP_score/celltype/red/less/", gene, "_expr.pdf"), plot, height = 6, width = 6)
}

######针对TLS分组画小提琴图和箱型图#########################
library(ggplot2)
library(ggpubr)
###########提取打分以及分组
dataMCP<-slide1[["MCP_TLSscore"]]###修改提取的打分
data<-as.data.frame(dataMCP@counts)
data<-t(data)
data<-as.data.frame(data)
colnames(data) <- gsub("-", "_", colnames(data))
data$id<-row.names(data)
dataarea<-as.data.frame(slide1@meta.data)####可以只提取想要的分组，也可以都提取
dataarea$id<-row.names(dataarea)
scoredata<-inner_join(data,dataarea,by="id")
colnames(scoredata) <- gsub("-", "_", colnames(scoredata))
#####################
lapply(c("MCP_score/TLS_vlin","MCP_score/TLS_box"), FUN = dir.create)
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
  
  ggsave(paste0("MCP_score/TLS_box/", cell_type, "_box_合起来打分.pdf"), plot = p_box, width = 4, height = 4)
  
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
  
  ggsave(paste0("MCP_score/TLS_vlin/", cell_type, "_vlin_合起来打分.pdf"), plot = p_violin, width = 5, height = 4)
}

