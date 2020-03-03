#挑选共有基因
#导入R包
library(dplyr)
library(org.Hvulgare1.eg.db)
library(readr)
library(clusterProfiler)
library(DOSE)
library(purrr)
library(tidyverse)
#设置字符串不变为因子
options(stringsAsFactors = F)
#setwd("/home/sda1/workspace/MetalRNA/10.GO_KEGG")
#读取文件夹中的重金属比较文件
list.files("/home/sda1/workspace/MetalRNA/13.DEG_analysis/")
FileName <- c(list.files(pattern = "*vs_CK$"))
#读取基因注释文件

f_anno <- read.table("anno_new.txt", 
                      sep = "\t", 
                      header = T, 
                      check.names = F,
                      fill = T, 
                      na.strings = "",
                      quote = "")

#基因文件取交集函数
co_gene <- function(f1,genelist,name){
  f1 <- dplyr::inner_join(f1, genelist,by = "GeneID")
  return (f1)
}
union_gene <- function(f1,genelist,name){
  f1 <- dplyr::full_join(f1, genelist,by = "GeneID")
  return (f1)
} 
#循环取集合
Name1 <- "U" #("I" or "U")
f_total <- data.frame(GeneID="")
取并集
for (i in 1:length(FileName)) {
  Name1 <- paste(Name1,FileName[i],sep = "_")
  DEG <- read.table(FileName[i], sep = "\t", header = TRUE, row.names = NULL, check.names = F, stringsAsFactors = F)
  print(Name1)
  colnames(DEG)[[1]] <- "GeneID"
  DEG_list_total <- as.data.frame(DEG)
  Name2 <- paste0(Name1,"_total")
  f_total <- union_gene(f_total,DEG_list_total, Name2)
  write.table(f_total$GeneID,file = paste0(Name2,"_GeneID"),sep = "\t",quote = F, row.names = F)
  f_geneid <- as.data.frame(f_total$GeneID)
  colnames(f_geneid)[1] <- "GeneID"
  f_write_anno_up <- dplyr::inner_join(f_total,f_anno, by = "GeneID")
  write.table(f_write_anno_up,file = paste0(Name2,"_anno"),sep = "\t",quote = F, row.names = F)
}



#取交集
Name1 <- "I" #("I" or "U")
f_total <-  read.table("U.metal.all.GeneID", 
                       sep = "\t", 
                       header = T, 
                       check.names = F,
                       fill = T, 
                       na.strings = "",
                       quote = "") 
#colnames(f_total)[1] <- "GeneID"
# for (i in 1:length(FileName)) {
#   Name1 <- paste(Name1,FileName[i],sep = "_")
#   DEG <- read.table(FileName[i], sep = "\t", header = TRUE, row.names = NULL, check.names = F, stringsAsFactors = F)
#   print(Name1)
#   colnames(DEG)[[1]] <- "GeneID" 
#   #全部基因取交集
#   DEG_list_total <- as.data.frame(DEG)
#   Name2 <- paste0(Name1,"_inter")
#   f_total <- co_gene(f_total,DEG_list_total, Name2)
#   write.table(f_total$GeneID,file = paste0(Name2,"_GeneID"),sep = "\t",quote = F, row.names = F)
#   f_geneid <- as.data.frame(f_total$GeneID)
#   colnames(f_geneid)[1] <- "GeneID"
#   f_write_anno_up <- dplyr::inner_join(f_anno, f_geneid, by = "GeneID")
#   write.table(f_write_anno_up,file = paste0(Name2,"_anno"),sep = "\t",quote = F, row.names = F)
# }

#GO富集分析函数
# Go_analysis <- function(gene_list, Name){
#   ego <- enrichGO(gene          = gene_list, #差异基因 vector
#                   keyType       = "GID",  #差异基因的 ID 类型，需要是 OrgDb 支持的
#                   OrgDb         = org.Hvulgare1.eg.db, #对应的OrgDb
#                   ont           = "ALL", #GO 分类名称，CC BP MF 
#                   pvalueCutoff  = 1, #Pvalue 阈值
#                   qvalueCutoff  = 1, #qvalue 阈值
#                   pAdjustMethod = "BH", #Pvalue 矫正方法
#                   readable      = FALSE) #TRUE 则展示SYMBOL，FALSE 则展示原来的ID
#   ego_results<-as.data.frame(ego)
#   write.table(ego_results, file = paste0(Name,"_ego_results.txt"), sep = "\t", quote = F)
#   pdf(paste0(Name,"_ego_barplot.pdf"))
#   p <- barplot(ego,showCategory = 20)
#   print(p)#这一步很关键
#   dev.off()
# }
#循环读取文件名称,然后交集匹配
#这里一定要注意的是DEG_list_up必须为数据框
Name1 <- "Inter"
####一个迭代函数
##FileName是向量
Iteration_name <- function(FileName){
  for (i in 1:length(FileName)) {
    Name1 <- paste(Name1,FileName[i],sep = "_")
    DEG <- read.table(FileName[i], sep = "\t", header = TRUE, row.names = NULL, check.names = F, stringsAsFactors = F)
    print(Name1)
    colnames(DEG)[[1]] <- "GeneID" 
    DEG_list_total <- as.data.frame(DEG)
    Name4 <- paste0(Name1,"_total")
    print(length(DEG_list_total$GeneID))
    f_total <- co_gene(f_total,DEG_list_total, Name4)
    f_num <- paste0(Name4,":",length(f_total$GeneID))
    write.table(f_num,file = "Metal_Inter_NUM",append = T,
                sep = "\t",quote = F, row.names = F)
    write.table(f_total,file = paste0(Name4,"_DegList"),sep = "\t",quote = F, row.names = F)
    f_write_anno_total <- dplyr::inner_join(f_total, f_anno, by = "GeneID")
    write.table(f_write_anno_total,file = paste0(Name4,"_anno"),sep = "\t", quote = F,row.names = F)
  }
}
#Iteration_name(FileName[1])

#Iteration_name(FileName)
#4和5这一步还么有解决已经解决，变量相加要加小括号
#又和python混了
#combn是排列组合
for (i in 1:5){
  test <- combn(FileName,3)
  for (j in 1:length(test[1,])){
    con_name <- test[,j]
    Iteration_name(con_name)
  }
}








