
rm(list = ls())
setwd('')
if (! dir.exists("./02_ssGSEA/")){
  dir.create("./02_ssGSEA")
}
setwd("./02_ssGSEA")

library(GSVA)

gene_set <- read.csv('../00_rawdata/mtUPR_RGs.csv') 
head(gene_set)
gene_set <- data.frame(symbols = gene_set$symbol, pathway =rep('mtUPR_RGs',length(gene_set$symbol.)))  

fpkm_log2 <- read.csv('../00_rawdata/dat(GSE58294).csv',row.names = 1,check.names = F)

dat.final2 <- as.matrix(fpkm_log2)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])

# ssgsea_score = gsva(dat.final2, gene_list, 
#                     method = "ssgsea", 
#                     ssgsea.norm = TRUE, 
#                     verbose = TRUE)
ssgsea_par <- ssgseaParam(
  exprData = dat.final2,               # 表达矩阵
  geneSets = gene_list,                # 基因集
  normalize = TRUE                    # 是否进行规范化
)

# 计算ssGSEA分数
ssgsea_score <- gsva(ssgsea_par)
                 
write.csv(ssgsea_score,
          file = "02.ssgsea_result.csv",
          quote = F)

save.image('ssgsea.RData')

rm(list = ls())
setwd('')
if (! dir.exists("./02_ssGSEA/")){
  dir.create("./02_ssGSEA")
}
setwd("./02_ssGSEA")

load('./ssgsea.RData')

group <- read.csv('../00_rawdata/group(GSE58294).csv')

colnames(group) <- c('sample', 'Group')
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install('dplyr')
library(dplyr)
group <- group %>%
  dplyr::mutate(Group = if_else(Group == "Cardioembolic", "IS", Group))

ssgsea_result <- ssgsea_score

ssgsea_result <- data.frame(t(ssgsea_result))
ssgsea_result$sample <- rownames(ssgsea_result)
identical(ssgsea_result$sample, group$sample)

ssgsea_result <- merge(ssgsea_result, group, by.x = 'sample', by.y = 'sample')

library(tidyr)
ssgsea_result<-gather(ssgsea_result,
                      key = pathway,
                      value = score,
                      -c("sample",'Group'))
plot.dat <- ssgsea_result

library(ggpubr)
library(ggplot2)
table(plot.dat$Group
      )
plot.dat$Group <- factor(plot.dat$Group, levels = c('Control','IS'))

exp_plot <- ggplot(plot.dat,aes(x = Group, y = score, color = Group)) +
  geom_violin(trim=F,color='black',aes(fill=Group)) + 
  
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA)+ 
  scale_fill_manual(values= c("#69bcce","#dc4f3e"), name = "Group")+
  labs(title="", x="", y = "score",size=20) +
  scale_colour_manual(values = c("#ff9a9b","#d2e099"))+
  
  stat_compare_means(comparisons = list(c('Control','IS')),
                     method = "wilcox.test",
                     
                     label = 'p.signif',
                     
                     label.y = 10, y.position = 14
                    
                    )+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(family='Times'))

exp_plot
ggsave(filename = '02.score_wilcox.pdf',exp_plot,w=5,h=5)
ggsave(filename = '02.score_wilcox.png',exp_plot,w=5,h=5,dpi = 600)

