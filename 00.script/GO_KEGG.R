#(1)导入包及设置路径
#setwd("/home/sda1/workspace/MetalRNA/10.GO_KEGG")
#(2)导入org.Hvulgare1.eg.db,准备GO分析
library(org.Hvulgare1.eg.db)
library(readr)
library(clusterProfiler)
library(DOSE)
library(purrr)
library(tidyverse)
library(clusterProfiler)
library(dplyr)
#(3)设置字符串不为因子
options(stringsAsFactors = F)
#(4)引入参数,T的作用就是摈弃R自带的参数,使输入参数开始值为1
Args <- commandArgs(T)
f1 <- Args[1]
#(5)查看物种包里的基本信息
#columns(org.Hvulgare1.eg.db)
#keys(org.Hvulgare1.eg.db)
#AnnotationDbi::select(org.Hvulgare.eg.db, 
#                      keys = "HORVU0Hr1G000620",
#                      columns = c("GO","Pathway"))

#GO分析
# 导入需要进行富集分析的基因列表，并转换为向量
#正则表达式，是要好好看下啦
#正则获取目录
list.files(f1)
FileName <- c(list.files(pattern = "*vs_CK$"))

#GO_ANALYSIS函数
Go_analysis <- function(gene_list, Name){
  ego <- enrichGO(gene          = gene_list, #差异基因 vector
                  keyType       = "GID",  #差异基因的 ID 类型，需要是 OrgDb 支持的
                  OrgDb         = org.Hvulgare1.eg.db, #对应的OrgDb
                  ont           = "ALL", #GO 分类名称，CC BP MF 
                  pvalueCutoff  = 1, #Pvalue 阈值
                  qvalueCutoff  = 1, #qvalue 阈值
                  pAdjustMethod = "BH", #Pvalue 矫正方法
                  readable      = FALSE) #TRUE 则展示SYMBOL，FALSE 则展示原来的ID
  ego_results<-as.data.frame(ego)
  write.table(ego_results, file = paste0(Name,"_ego_results.txt"), sep = "\t", quote = F)
  pdf(paste0(Name,"_ego_barplot.pdf"))
  p <- barplot(ego,showCategory = 20)
  print(p)#这一步很关键
  dev.off()
}

#循环提取数据
for (i in 1:length(FileName)) {
  Name1 <- FileName[i]
  DEG <- read.table(FileName[i], sep = "\t", header = TRUE, row.names = NULL, check.names = F, stringsAsFactors = F)
  DEG_list_up <- as.character((DEG %>% 
  dplyr::filter(log2FoldChange >= 1 & padj <= 0.05))$row.names)
  Name_up <- paste0(Name1,"_up")
  Go_analysis(DEG_list_up, Name_up)
#此处注意数据框，是不能用在enrichGO里的,只能用向量
  DEG_list_down <- as.character((DEG %>% 
    dplyr::filter(log2FoldChange <= -1 & padj <= 0.05))$row.names)
  Name_down <- paste0(Name1,"_down")
  Go_analysis(DEG_list_down, Name_down)
  DEG_list_total <- as.character((DEG %>% 
    dplyr::filter(abs(log2FoldChange) >= 1 & padj <= 0.05))$row.names)
  Name_total <- paste0(Name1,"_total")
  Go_analysis(DEG_list_total, Name_total)
}
###ENRICHGO示例


#enrichGO多种图
#dotplot(ego)
#emapplot(ego)
#plotGOgraph(ego)
#KEGG分析

#pathway2gene <- AnnotationDbi::select(org.Hvulgare.eg.db, 
#                                      keys = keys(org.Hvulgare.eg.db), 
#                                      columns = c("Pathway","Ko")) %>%
#  na.omit() %>%
#  dplyr::select(Pathway, GID)

# 导入 Pathway 与名称对应关系

#load("kegg_info.RData")

#KEGG pathway 富集
#ekp <- enricher(gene_list, 
#                TERM2GENE = pathway2gene, 
#                TERM2NAME = pathway2name, 
#                pvalueCutoff = 1, 
#                qvalueCutoff = 1,
#                pAdjustMethod = "BH",
#                minGSSize = 1)

#ekp_results <- as.data.frame(ekp)

#barplot(ekp, showCategory=20, x = "GeneRatio")

#dotplot(ekp)

#emapplot(ekp)
#挑选共有基因
library("dplyr")
library(org.Hvulgare1.eg.db)
library(readr)
library(clusterProfiler)
library(DOSE)
library(purrr)
library(tidyverse)
library(clusterProfiler)
library(dplyr)
list.files("/home/sda1/workspace/MetalRNA/10.GO_KEGG")
FileName <- c(list.files(pattern = "*vs_CK$"))
f1 <- as.data.frame(read.table("all_gene.txt",sep = "\t",check.names = F, stringsAsFactors = F))
colnames(f1)[[1]] <- "GeneID" 
Go_analysis <- function(gene_list, Name){
  ego <- enrichGO(gene          = gene_list, #差异基因 vector
                  keyType       = "GID",  #差异基因的 ID 类型，需要是 OrgDb 支持的
                  OrgDb         = org.Hvulgare1.eg.db, #对应的OrgDb
                  ont           = "ALL", #GO 分类名称，CC BP MF 
                  pvalueCutoff  = 1, #Pvalue 阈值
                  qvalueCutoff  = 1, #qvalue 阈值
                  pAdjustMethod = "BH", #Pvalue 矫正方法
                  readable      = FALSE) #TRUE 则展示SYMBOL，FALSE 则展示原来的ID
  ego_results<-as.data.frame(ego)
  write.table(ego_results, file = paste0(Name,"_ego_results.txt"), sep = "\t", quote = F)
  pdf(paste0(Name,"_ego_barplot.pdf"))
  p <- barplot(ego,showCategory = 20)
  print(p)#这一步很关键
  dev.off()
}
co_gene <- function(f1,genelist,name){
  f1 <- dplyr::inner_join(f1, genelist,by = "GeneID")
  f2 <- as.character(f1$GeneID)
  Go_analysis(f2,name)
  return (f1)
}
Name1 <- "Union"
for (i in 1:length(FileName)) {
  Name1 <- paste(Name1,FileName[i],sep = "_")
  DEG <- read.table(FileName[i], sep = "\t", header = TRUE, row.names = NULL, check.names = F, stringsAsFactors = F)
  colnames(DEG)[[1]] <- "GeneID" 
  DEG_list_up <- as.character((DEG %>% 
                                 dplyr::filter(log2FoldChange >= 1 & padj <= 0.05))$row.names)
  Name2 <- paste0(Name1,"_up")
  f1 <- co_gene(f1,DEG_list_up, Name2)
  DEG_list_down <- as.character((DEG %>% 
                                   dplyr::filter(log2FoldChange <= -1 & padj <= 0.05))$row.names)
  Name3 <- paste0(Name1,"_down")
  f1 <- co_gene(f1,DEG_list_down, Name3)
  DEG_list_total <- as.data.frame(DEG %>% 
                    dplyr::filter(abs(log2FoldChange) >= 1 & padj <= 0.05))
  Name4 <- paste0(Name1,"_total")
  f1 <- co_gene(f1,DEG_list_total, Name4)
}

Iteration_name <- function(FileName) {
  for (i in 1:length(FileName)) {
    Name1 <- paste(Name1,FileName[i],sep = "_")
    DEG <- read.table(FileName[i], sep = "\t", header = TRUE, row.names = NULL, check.names = F, stringsAsFactors = F)
    colnames(DEG)[[1]] <- "GeneID" 
    DEG_list_up <- as.character((DEG %>% 
                                   dplyr::filter(log2FoldChange >= 1 & padj <= 0.05))$row.names)
    Name2 <- paste0(Name1,"_up")
    f1 <- co_gene(f1,DEG_list_up, Name2)
    DEG_list_down <- as.character((DEG %>% 
                                     dplyr::filter(log2FoldChange <= -1 & padj <= 0.05))$row.names)
    Name3 <- paste0(Name1,"_down")
    f1 <- co_gene(f1,DEG_list_down, Name3)
    DEG_list_total <- as.data.frame(DEG %>% 
                                      dplyr::filter(abs(log2FoldChange) >= 1 & padj <= 0.05))
    Name4 <- paste0(Name1,"_total")
    f1 <- co_gene(f1,DEG_list_total, Name4)
 }
}
Iteration_name(FileName)


