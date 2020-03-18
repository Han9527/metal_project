#导入包
library('ggplot2')
library('dplyr')
library('reshape2')
library("Rmisc")
library("lubridate")
library('scales')
library('RColorBrewer')
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(agricolae)
#install.packages("agricolae")
#导入数据
T_data <- read.table("./test_box.txt", 
                     header = T, 
                     stringsAsFactors = F, 
                     check.names = F)


#计算三个重复的标准差和标准误
T_se <- summarySE(data = T_data, 
                  measurevar = "value", 
                  groupvars = "treat")
head(T_se)
#计算多重比较
#(1)方差分析
tmp <- aov(value ~ treat, T_data)
#(2)求多重比较结果
T_mul <- LSD.test(tmp, 'treat', p.adj = 'bonferroni')
print(T_mul$groups)
#ggplot2绘图
#全局设置科研模板

theme_set(theme_bw())
#palette<-c(brewer.pal(7,"Set1")[c(1,2,3,4,5,6)])
#柱状图的绘制
# compaired <- list(c("CK", "-Fe"), 
#                   c("CK","-Cu"), 
#                   c("CK","-Zn"),
#                   c("CK","-Mn"),
#                   c("CK","+Cd"))
#stat_compare_means(comparisons = compaired,method = "wilcox.test") +

T_se$treat <- factor(T_se$treat,
                       levels = c("CK","+Cd","-Mn","-Fe","-Zn","-Cu"))
#label <- c("a", "b", "a", "a", "a", "a")
p2 <- ggplot(T_se,aes(x = treat, y = value, colours =  treat)) +
  geom_bar(stat = "identity",
           width = 0.25,
           position = position_dodge(),
           aes(fill = treat)) +
  geom_errorbar(aes(ymin=value-T_se$se, ymax=value+T_se$se), 
                width=.10,
                position = position_dodge(), ) +
  scale_fill_brewer(palette = "Set1") +
  labs(x="") +
  labs(y="Metal Content") +
  #theme_classic(base_size = 16, base_family = 'momo') +
  #网格线
  theme(panel.grid.major =element_blank(),
        panel.grid.minor =element_blank(),
        #边框线
        panel.background = element_rect(colour = "black", size = 2),
        text=element_text(family="serif",colour = "black"),
        axis.text=element_text(color="black")) +
        coord_cartesian(ylim = c(0,7)) +
        ##这个可以去掉与X轴间隙
        scale_x_discrete(expand=c(0.1, 0)) +
        scale_y_continuous(expand = c(0.02,0)) + 
        ##这个可以去掉与Y轴间隙
        #scale_x_continuous(limits=c(0.6,2.4),breaks=c(1,2),
                     #labels=c("农气观测一级站","农气观测二级站"))
        #scale_x_discrete(expand = c(0,0)) +
        #scale_x_discrete(limits = c(0.4,3.1), breaks = c(0.5,1,1.5,2,2.5,3.0), 
                           #labels = c("CK","+Cd","-Mn","-Fe","-Zn","-Cu")
                           #) +
        #修改坐标轴大小
        theme(axis.text = element_text(face="bold", color="black", size=12),
              axis.title = element_text(face="bold", color="black", size=12),
              legend.text = element_text(face="bold", color="black", size=12),
              legend.title = element_text(face="bold", color="black", size=12))
  

ggsave(p2,filename = "HvKFB39.tiff",width = 4, height = 3.19, units = "in", dpi = 600)
ggsave(p2,filename = 'HvKFB39.pdf', width = 4, height = 3.19, units = "in")




#函数成品完成

metal_fun <- function(metal_con) {
  T_data <- read.table(metal_con, 
                       header = T, 
                       stringsAsFactors = F, 
                       check.names = F)
  
  #计算三个重复的标准差和标准误
  T_se <- summarySE(data = T_data, 
                    measurevar = "value", 
                    groupvars = "treat")
  #计算多重比较
  #(1)方差分析
  tmp <- aov(value ~ treat, T_data)
  #(2)求多重比较结果
  T_mul <- LSD.test(tmp, 'treat', p.adj = 'bonferroni')
  write.table(T_mul$groups, file = paste0(metal_con,"multi"),
              sep = "\t", row.names = F, quote = F)
  T_se$treat <- factor(T_se$treat,
                       levels = c("CK","+Cd","-Mn","-Fe","-Zn","-Cu"))
  #ggplot2绘图
  #全局设置科研模板
  p1 <- ggplot(T_se,aes(x = treat, y = value, colours =  treat)) +
    geom_bar(stat = "identity", 
             position = position_dodge(),
             aes(fill = treat)) +
    geom_errorbar(aes(ymin=value-T_se$se, ymax=value+T_se$se), 
                  width=.15,
                  position = position_dodge(0.9)) +
    scale_fill_brewer(palette = "Set1") +
    labs(x="") +
    labs(y="Metal Content") +
    #theme_classic(base_size = 16, base_family = 'momo') +
    #网格线
    theme(panel.grid.major =element_blank(),
          panel.grid.minor =element_blank(),
          #边框线
          panel.background = element_rect(colour = "black", size = 2),
          text=element_text(family="serif",colour = "black"),
          axis.text=element_text(color="black")) +
    coord_cartesian(expand = T) +
    ##这个可以去掉与X轴间隙
    scale_y_continuous(expand = c(0.02,0.02)) + 
    ##这个可以去掉与Y轴间隙
    scale_x_discrete(expand = c(0.09,0.09)) +
    #修改坐标轴大小
    theme(axis.text = element_text(face="bold", color="black", size=12),
          axis.title = element_text(face="bold", color="black", size=14),
          legend.text = element_text(face="bold", color="black", size=14),
          legend.title = element_text(face="bold", color="black", size=14))

 ggsave(p1,filename = paste0(metal_con,".png"),dpi = 600)
 ggsave(p1,filename = paste0(metal_con,'.pdf'))
}

metal_fun("test_box.txt")
metal_fun("Cd_con.txt")
metal_fun("Fe_con.txt")
metal_fun("Mn_con.txt")
metal_fun("Zn_con.txt")
metal_fun("Cu_con.txt")
#(2)宅柱+小图
metal_fun1 <- function(metal_con) {
  T_data <- read.table(metal_con, 
                       header = T, 
                       stringsAsFactors = F, 
                       check.names = F)
  theme_set(theme_bw())
  #计算三个重复的标准差和标准误
  T_se <- summarySE(data = T_data, 
                    measurevar = "value", 
                    groupvars = "treat")
  #计算多重比较
  #(1)方差分析
  tmp <- aov(value ~ treat, T_data)
  #(2)求多重比较结果
  T_mul <- LSD.test(tmp, 'treat', p.adj = 'bonferroni')
  write.table(T_mul$groups, file = paste0(metal_con,"multi"),
              sep = "\t", row.names = F, quote = F)
  T_se$treat <- factor(T_se$treat,
                       levels = c("CK","+Cd","-Mn","-Fe","-Zn","-Cu"))
  #ggplot2绘图
  #order(T_se$value, decreasing = T)
  #max(T_se$value)
  #全局设置科研模板
  p2 <- ggplot(T_se,aes(x = treat, y = value, colours =  treat)) +
    geom_bar(stat = "identity",
             width = 0.25,
             position = position_dodge(),
             aes(fill = treat)) +
    geom_errorbar(aes(ymin=value-T_se$se, ymax=value+T_se$se), 
                  width=.10,
                  position = position_dodge(), ) +
    scale_fill_brewer(palette = "Set1") +
    labs(x="") +
    labs(y="Metal Content") +
    #theme_classic(base_size = 16, base_family = 'momo') +
    #网格线
    theme(panel.grid.major =element_blank(),
          panel.grid.minor =element_blank(),
          #边框线
          panel.background = element_rect(colour = "black", size = 2),
          text=element_text(family="serif",colour = "black"),
          axis.text=element_text(color="black")) +
    coord_cartesian(ylim = c(0,max(T_se$value)*1.2)) +
    ##这个可以去掉与X轴间隙
    scale_x_discrete(expand=c(0.1, 0)) +
    scale_y_continuous(expand = c(0.02,0)) + 
    ##这个可以去掉与Y轴间隙
    #scale_x_continuous(limits=c(0.6,2.4),breaks=c(1,2),
    #labels=c("农气观测一级站","农气观测二级站"))
    #scale_x_discrete(expand = c(0,0)) +
    #scale_x_discrete(limits = c(0.4,3.1), breaks = c(0.5,1,1.5,2,2.5,3.0), 
    #labels = c("CK","+Cd","-Mn","-Fe","-Zn","-Cu")
    #) +
    #修改坐标轴大小
    theme(axis.text = element_text(face="bold", color="black", size=12),
          axis.title = element_text(face="bold", color="black", size=12),
          legend.text = element_text(face="bold", color="black", size=12),
          legend.title = element_text(face="bold", color="black", size=12))
  
  
  ggsave(p2,filename = paste0(metal_con,".tiff"),width = 4, height = 3.19, units = "in", dpi = 600)
  ggsave(p2,filename = paste0(metal_con,'.pdf'), width = 4, height = 3.19, units = "in")

}


metal_fun1("test_box.txt")
metal_fun1("Cd_con.txt")
metal_fun1("Fe_con.txt")
metal_fun1("Mn_con.txt")
metal_fun1("Zn_con.txt")
metal_fun1("Cu_con.txt")