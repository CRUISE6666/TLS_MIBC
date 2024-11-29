







rm(list=ls())
library(stringr)
library(pheatmap)
library(gplots)
library(grid)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 显示进程
display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
  
} 

# 计算组间统计差异
cross_subtype_compr <- function(expr = NULL,
                                subt = NULL,
                                subt.label = "Subtype",
                                two_sam_compr_method = "wilcox",
                                multi_sam_compr_method = "kruskal",
                                res.path = NULL) {
  
  if (!is.element(two_sam_compr_method, c("t.test", "wilcox"))) {stop("Two samples comparison should be t.test or wilcox!\n") }
  if (!is.element(multi_sam_compr_method, c("anova", "kruskal"))) {stop("multiple samples comparison should be kruskal or anova!\n") }
  
  subt.name <- unique(subt[,subt.label])
  n.subt <- length(subt.name)
  if(n.subt < 2) {stop("The number of subtype should be greater than 2!\n")}
  
  comprTab <- NULL
  
  # 两个亚型且为非参数检验
  if(n.subt == 2 & two_sam_compr_method == "wilcox") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp1 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[1]),,drop = F])])
      tmp2 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[2]),,drop = F])])
      wt <- wilcox.test(tmp1,tmp2)
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = wt$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 两个亚型且为参数检验
  if(n.subt == 2 & two_sam_compr_method == "t.test") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp1 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[1]),,drop = F])])
      tmp2 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[2]),,drop = F])])
      tt <- t.test(tmp1,tmp2)
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = tt$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 多个亚型且为非参数检验
  if(n.subt > 2 & multi_sam_compr_method == "kruskal") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp.list <- list()
      for (n in 1:n.subt) {
        tmp.list[[n]] <- data.frame(value = as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[n]),,drop = F])]),
                                    subt = subt.name[n],
                                    stringsAsFactors = F)
      }
      tmp <- do.call(rbind,tmp.list)
      kt <- kruskal.test(value ~ subt,data = tmp)
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = kt$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 多个亚型且为参数检验
  if(n.subt > 2 & multi_sam_compr_method == "anova") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp.list <- list()
      for (n in 1:n.subt) {
        tmp.list[[n]] <- data.frame(value = as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[n]),,drop = F])]),
                                    subt = subt.name[n],
                                    stringsAsFactors = F)
      }
      tmp <- do.call(rbind,tmp.list)
      at <- summary(aov(value ~ subt,data = tmp))
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = at[[1]][1,5],
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 调整p值
  comprTab$adjusted.p.value = p.adjust(comprTab$nominal.p.value,method = "BH")
  # 按p值排序
  #comprTab <- comprTab[order(comprTab$adjusted.p.value, decreasing = F),] 
  
  write.table(comprTab,file.path(res.path,"comprTab.txt"),sep = "\t",row.names = F,quote = F)
  return(comprTab)
}



##########表达矩阵，列名为样品名
data<-read.csv('RMYY_BLCA/RMYY_BLCA_MIBC_MCPcounter.csv')
row.names(data)<-data[,1]
data<-data[,-1]
mygene_data<-data
mygene_data<-as.data.frame(mygene_data)
#mygene_data<-as.data.frame(t(mygene_data))

## 读取分组信息，行名为样品名
clinical <- read.table("raw_data//clinical.txt", header = TRUE, sep = "\t")
Subtype <- clinical[,c(1,22,23)]
rownames(Subtype)<-Subtype[,1]
Subtype<-Subtype[,-1]
colnames(Subtype)[2]<-"Subtype"
colnames(Subtype)
# 查看各组名字和sample数量，画图时要用
table(Subtype$Subtype)

## 取表达矩阵和分组信息中共有的sample
com_sam <- intersect(colnames(mygene_data),rownames(Subtype))
mygene_data <- mygene_data[,com_sam]
Subtype <- Subtype[com_sam,,drop = F]
head(Subtype)













## 用前面的自定义函数计算组间统计差异
comprTab <- cross_subtype_compr(expr = mygene_data, # 或log2(mygene_data + 1)，如果用参数检验，请注意对数转化；若非参均可
                                subt = Subtype,
                                #two_sam_compr_method = "wilcox", # 两组"t.test", "wilcox"
                                multi_sam_compr_method = "kruskal", # 多组"anova", "kruskal"
                                res.path = ".")
# 用全部基因来画
n.show_top_gene <- nrow(mygene_data)
# 或者取top 20个基因来画
#n.show_top_gene <- 20 

# 按分组排序
subt.order <- Subtype[order(Subtype$Subtype),,drop = F]
indata <- mygene_data[comprTab$gene[1:n.show_top_gene],rownames(subt.order)]


# 开始画图
# 数据标准化和边界设置
plotdata <- t(scale(t(indata)))
plotdata[plotdata > 2] <- 2
plotdata[plotdata < -2] <- -2

# 调整行名
blank <- "    " # 行名和p值之间的间隔
p.value <- comprTab$adjusted.p.value[1:n.show_top_gene]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*",""))))
p.label <- formatC(p.value, # 将p值变成保留两位小数的科学计数法
                   format = "e",
                   digits = 2)

add.label <- str_pad(paste0(rownames(plotdata),sig.label), # 固定行名宽度并再右侧补齐" "
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")

annCol <- subt.order # 获得排序后的亚型注释信息，这里只有一个变量需要注释
colnames(annCol)[1] <- paste(str_pad(colnames(annCol)[1], # 注释列名补上"P-value"，宽度和刚才一致
                                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                                     side = "right"),
                             "P-value",
                             sep = blank)
annColors <- list(c("high"="lightblue", "low"="darkgreen")) # 如果有多个变量要注释颜色请补充c()
#annColors <- list(c("High"="lightblue", "Low"="darkgreen", "Medium"="pink")) # 如果有多个变量要注释颜色请补充c()
names(annColors) <- colnames(annCol)[1] # 如果有多个变量要注释颜色请补充每张list的name

# 绘制热图
pheatmap(cellheight = 15, cellwidth = 1,
         mat = plotdata, # 输入数据
         scale = "none", # 不标准化因为数据已经被标准化
         annotation_col = annCol, # 列注释信息
         annotation_colors = annColors, # 列注释对应的颜色
         cluster_cols = F, # 列不聚类
         cluster_rows = F, # 行不聚类
         show_colnames = F, # 不显示列名
         show_rownames = T, # 显示基因名
         #annotation_legend = F, # 不显示图例
         labels_row = paste(add.label, p.label, sep=blank), # 自定义样本名义blank作间隔
         fontfamily = "mono", # 关键，使用fixed font而不是proportional font
         filename = "heatmapPvalue_2.pdf")



library(pheatmap)
library(RColorBrewer)

# 选择一个颜色盘
color_palette <- brewer.pal(9, "RdBu")  # 这将产生一个蓝色系的颜色渐变

# 或者创建一个自定义的颜色渐变
# 例如，从蓝色到白色，再到红色
color_palette <- c("#2166ac", "#67a9cf", "#d1e5f0","#f7f7f7","#fddbc7",  "#ef8a62","#b2182b")

# 绘制热图
pheatmap(cellheight = 15, cellwidth = 7,
         mat = plotdata,  # 输入数据
         color = color_palette,  # 设置颜色渐变
         scale = "none",  # 不标准化因为数据已经被标准化
         annotation_col = annCol,  # 列注释信息
         annotation_colors = annColors,  # 列注释对应的颜色
         cluster_cols = FALSE,  # 列不聚类
         cluster_rows = FALSE,  # 行不聚类
         show_colnames = FALSE,  # 不显示列名
         show_rownames = TRUE,  # 显示基因名
         # annotation_legend = FALSE,  # 不显示图例（如果你想显示图例，删除这一行的注释）
         labels_row = paste(add.label, p.label, sep = " "),  # 自定义样本名，注意这里的空格不是'blank'而是实际的空格
         fontfamily = "mono",  # 关键，使用fixed font而不是proportional font
         filename = "heatmapPvalue_2.pdf"
)
