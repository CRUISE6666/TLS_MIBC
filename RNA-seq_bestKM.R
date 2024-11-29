

###########################################best km计算
##加载必需的包
library(survival)
library(survminer)

####
#rm(list=ls())
#svdata <- read.csv("km-MCPcounter.csv",header = T,row.names = 1)
svdata<-MCP_KM
dim(svdata) #一共18个基因
head(svdata[1:3,1:6])

##如果你把所有基因都拿来作为输入，其中一些基因表达量没变化就不可能分开，需要先删掉
mean <- apply(svdata, 2, mean)
max <- apply(svdata, 2, max)
check <- max - mean
svdata <- svdata[,check!=0]

##对数据集的所有基因进行bestSeparation统计
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.25) #默认组内sample不能低于30%


head(res.cut)
class(res.cut)
##按照bestSeparation分高低表达
res.cat <- surv_categorize(res.cut)

##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  
  ##计算HR以及95%CI
  ##修改分组参照
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  
  #只画出p value<=0.05的基因，如果不想筛选，就删掉下面这行
  #if (p.val>0.05) next
  
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  #按照基因表达量从低到高排序，便于取出分界表达量
  svsort <- svdata[order(svdata[,i]),]
  
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #想要网格就运行这行
                      conf.int = T, #不画置信区间，想画置信区间就把F改成T
                      conf.int.style = "ribbon",#置信区间的类型，还可改为ribbon step
                      censor = F, #不显示观察值所在的位置
                      palette = c("#386cb0","#666666"), #线的颜色对应高、低
                      
                      legend.title = i,#基因名写在图例题目的位置
                      font.legend = 11,#图例的字体大小
                      #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                      
                      #在图例上标出高低分界点的表达量，和组内sample数量
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      
                      #在左下角标出pvalue、HR、95% CI
                      #太小的p value标为p < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  
  #如果想要一个图保存为一个pdf文件，就把下面这行前面的“#”删掉
  ggsave(paste0("CXCL13_Tcell/IMvigor210/", i, ".pdf"), width = 4, height = 4)
}
length(pl)

#######合并成一个图
res <- arrange_ggsurvplots(pl, 
                           print = T,
                           ncol = 3, nrow = 4)#每页纸画几列几行

#保存到pdf文件
ggsave("CXCL13_Tcell/IMvigor210/bestSurvPlot.pdf",res,width = 12,height = 16)


#########################################################把结果生成一个表格
# 统计作图所需的生存对象
my.surv <- Surv(res.cat$futime, res.cat$fustat)
# 初始化汇总表
results_df <- data.frame(Gene = character(), 
                         P_Value = numeric(), 
                         Hazard_Ratio = numeric(), 
                         CI_Lower = numeric(), 
                         CI_Upper = numeric(), 
                         stringsAsFactors = FALSE)

# 循环遍历每个基因
for (i in colnames(res.cat)[3:ncol(res.cat)]) {
  group <- res.cat[, i]
  fit <- survfit(my.surv ~ group)
  
  # 计算HR以及95% CI
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  
  # 将结果添加到汇总表
  results_df <- rbind(results_df, 
                      data.frame(Gene = i, 
                                 P_Value = ifelse(p.val < 0.001, NA, p.val), 
                                 Hazard_Ratio = HR, 
                                 CI_Lower = low95, 
                                 CI_Upper = up95))
}

# 查看汇总表
print(results_df)
# 保存汇总表为CSV文件
write.csv(results_df, "CXCL13_Tcell/IMvigor210/IMvigor210_KM_results.csv", row.names = FALSE)



##########################################################分组计算单因素和多因素cox#####
###先提取部分细胞进行分析
colnames(res.cat)
# Columns that end with "_all" plus "futime" and "fustat"
columns_all <- c("futime", "fustat", grep("_all$", colnames(res.cat), value = TRUE))
# Data frame with "_all" columns and "futime", "fustat"
df_all <- res.cat[, columns_all]
colnames(df_all)
# All column names except those in columns_all
remaining_columns <- setdiff(colnames(res.cat), columns_all)
# Data frame with the remaining columns plus "futime" and "fustat"
df_major <- res.cat[, c("futime", "fustat", remaining_columns)]
colnames(df_major)

#######分两步做
rt<-df_major
head(rt)
# Replace "high" with 1 and "low" with 0 in all columns except the first two
rt[, -c(1, 2)] <- lapply(rt[, -c(1, 2)], function(x) ifelse(x == "high", 1, ifelse(x == "low", 0, x)))
# View the first few rows of the modified data frame
head(rt)

#单因素独立预后分析
uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="CXCL13_Tcell/IMvigor210/CGA_BLCA_uniCox_major.txt",sep="\t",row.names=F,quote=F)

#多因素独立预后分析
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="CXCL13_Tcell/IMvigor210/IMvigor210_multiCox_major.txt",sep="\t",row.names=F,quote=F)


############绘制森林图函数############
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width = 6.3,height = 4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}
############绘制森林图函数############

bioForest(coxFile="CXCL13_Tcell/IMvigor210/IMvigor210_uniCox_major.txt",forestFile="CXCL13_Tcell/IMvigor210/IMvigor210_uniForest_major.pdf",forestCol="green")
bioForest(coxFile="CXCL13_Tcell/IMvigor210/IMvigor210_multiCox_major.txt",forestFile="CXCL13_Tcell/IMvigor210/IMvigor210_multiForest_major.pdf",forestCol="red")
