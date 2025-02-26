rm(list = ls())
setwd('')
if (! dir.exists("./06_expression")){
  dir.create("./06_expression")
}
setwd("./06_expression")

library(tidyverse)
library(dplyr)
gene.all <- read.csv('../05_LASSO_RF_/03.Intersection_gene.csv') %>% as.data.frame()
class(gene.all)
colnames(gene.all) <- 'symbol'
gene.all <- data.frame(symbol=unique(gene.all$symbol))

dat<-read.csv('../00_rawdata/dat(GSE16561).csv',row.names = 1,check.names = F) #%>% lance::lc.tableToNum()
group <- read.csv('../00_rawdata/group(GSE16561).csv')

colnames(group) <- c('sample', 'group')

table(group$group)

hub_exp<-dat[gene.all$symbol,group$sample]
hub_exp2<-hub_exp

hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

Control.sample <- group$sample[which(group$group=='Control')]

hub_exp2$Group<-ifelse(hub_exp2$sample%in%Control.sample,'Control','IS')
hub_exp2<-na.omit(hub_exp2)

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)

stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')

stat.test$p1<-ifelse(stat.test$p<0.0001,"****",
                     ifelse(stat.test$p<0.001,"***",
                            ifelse(stat.test$p<0.01,"**",
                                   ifelse(stat.test$p<0.05,"*",'ns'))))

library(ggplot2)
library(ggpubr)
hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Control','IS'))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  
  
  stat_boxplot(geom="errorbar", 
               width=0.2,
               color="black",
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.8,
               aes(fill=Group),
               color="black",
               position=position_dodge(0.9)
               
  )+ 
  
  
  
  
  

scale_fill_manual(values= c("#69bcce","#fed53b"), name = "Group")+
  labs(title="Expression in GSE16561", x="", y = "expression",size=20) +
  scale_colour_manual(values = c("#69bcce","#fed53b"))+
  
  
  
  stat_pvalue_manual(stat.test,
                     y.position = c(10,11,16,15),
                     size = 3.2,
                     family = "Times",
                     label = "p1",  
                     
                     face = "bold")+
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
        text=element_text(family='Times'))+
  facet_wrap(~Symbol,scales = "free",nrow = 1,strip.position = "top") 
exp_plot
ggsave(filename = '01.expression1.pdf',exp_plot,w=8,h=5)
ggsave(filename = '01.expression1.png',exp_plot,w=8,h=5,dpi = 600)

write.csv(gene.all,'03.hubgene.csv',row.names = F)

rm(list = ls())
# setwd('')
# if (! dir.exists("./06_expression")){
#   dir.create("./06_expression")
# }
# setwd("./06_expression")

library(tidyverse)
library(dplyr)
gene.all <- read.csv('../05_LASSO_RF_/03.Intersection_gene.csv') %>% as.data.frame()
class(gene.all)
colnames(gene.all) <- 'symbol'
gene.all <- data.frame(symbol=unique(gene.all$symbol))

dat<-read.csv('../00_rawdata/dat(GSE58294).csv',row.names = 1,check.names = F)# %>% lance::lc.tableToNum()
group <- read.csv('../00_rawdata/group(GSE58294).csv')

colnames(group) <- c('sample', 'group')

table(group$group)

hub_exp<-dat[gene.all$symbol,group$sample]
hub_exp2<-hub_exp

hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

Control.sample <- group$sample[which(group$group=='Control')]

hub_exp2$Group<-ifelse(hub_exp2$sample%in%Control.sample,'Control','IS')

table(hub_exp2$Group)

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')

stat.test$p1<-ifelse(stat.test$p<0.0001,"****",
                     ifelse(stat.test$p<0.001,"***",
                            ifelse(stat.test$p<0.01,"**",
                                   ifelse(stat.test$p<0.05,"*",'ns'))))

library(ggplot2)
library(ggpubr)
hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Control','IS'))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  
  
  stat_boxplot(geom="errorbar", 
               width=0.2,
               color="black",
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.8,
               aes(fill=Group),
               color="black",
               position=position_dodge(0.9)
               
  )+ 
  
  
  
  
  
  
  
  
  
  scale_fill_manual(values= c("#69bcce","#fed53b"), name = "Group")+
  labs(title="Expression in GSE58294", x="", y = "expression",size=20) +
  scale_colour_manual(values = c("#69bcce","#fed53b"))+
  
  
  
  stat_pvalue_manual(stat.test,
                     y.position = c(3,5,6,5),
                     size = 3.2,
                     family = "Times",
                     label = "p1",  
                     
                     face = "bold")+
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
        text=element_text(family='Times'))+
  facet_wrap(~Symbol,scales = "free",nrow = 1,strip.position = "top") 
exp_plot
ggsave(filename = '02.expression.pdf',exp_plot,w=8,h=5)
ggsave(filename = '02.expression.png',exp_plot,w=8,h=5,dpi = 600)

