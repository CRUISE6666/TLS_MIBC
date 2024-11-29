

######MCP前期数据处理
probesets<-read.table("CXCL13_Tcell/probesets.txt",sep="\t",stringsAsFactors=FALSE,colClasses="character")
genes<-read.table("CXCL13_Tcell/genes.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
head(genes)
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


####################开始数据分析
#datasets <- c("IMvigor210", "IMvigor210", "IMvigor210", "IMvigor210", "IMvigor210", "IMvigor210", "IMvigor210", "IMvigor210", "IMvigor210")

lapply(c("CXCL13_Tcell/IMvigor210"), FUN = dir.create)
options(future.globals.maxSize = 100 * 1024^3)

lapply(c("CXCL13_Tcell/IMvigor210/MCP_best"), FUN = dir.create)
options(future.globals.maxSize = 100 * 1024^3)

########
#install.packages(c('devtools','curl')) 
library(curl)
library(devtools)
#install_github('ebecht/MCPcounter',ref='master', subdir='Source',force = TRUE)
library(MCPcounter)

######读取原始数据
load("data//IMvigor210//IMvigor210_MIBC_matrix.Rdata")

#save(dat,data="data//IMvigor210//IMvigor210_MIBC_matrix.Rdata")
dat<-dat
#save(dat,file="data//IMvigor210//IMvigor210_MIBC_matrix.Rdata")

test<-dat[1:10,1:10]
#dat<-aggregate(.~ID,mean,data=dat)

#row.names(dat)<-dat[,1]
#dat<-dat[,-1]
#test<-dat[1:10,1:10]


#如果网络不佳proesets和genes也可以不提供
results<- MCPcounter.estimate(dat,featuresType,
                              probesets=probesets,
                              genes=genes
)

######存储数据
write.csv(results,'CXCL13_Tcell/IMvigor210/IMvigor210_MIBC_MCPcounter.csv')
