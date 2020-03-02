library(VennDiagram)
library(RColorBrewer)
library(grDevices)
#install.packages("VennDiagram")
# 设置工作路径
setwd("/Users/Davey/Desktop/VennDiagram/")
# 清除当前环境中的变量
rm(list=ls())
#输入数据
data = read.table("tcdb.txt",header = T,sep="\t", stringsAsFactors = F, check.names = F)
head(data)
#绘制五个
venn.plot <- venn.diagram(
  list(Cd_CK=data$Cd,Cu_CK=data$Cu,Fe_CK=data$Fe,Mn_CK=data$Mn,Zn_CK=data$Zn),fill=c(brewer.pal(7,"Set1")[1:5]),
  filename = "TCDB_venn.tiff",
  #lty = "dotted",
  lwd = 1,
  col = "transparent",  #"transparent",
  #fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  alpha = 0.60,
  cat.col = c(fill=c(brewer.pal(7,"Set1")[1:5])),
  cat.cex = 1.2,
  cat.fontface = "bold",
  margin = 0.10,
  cex = 1.2,
  cat.dist = c(0.2,0.25,0.2,0.2,0.25)#此顺序是按照输入的向量顺序
  #cat.pos = 0.07
)
#若想生成pdf,按照如下命令行
pdf(file="TCDB_venn.pdf")
grid.draw(venn.plot)
dev.off()
