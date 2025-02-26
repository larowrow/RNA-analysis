rm(list = ls())

setwd('')
if (!dir.exists("./07_ROC")) {
  dir.create("./07_ROC")
}
setwd("./07_ROC")
# library(lance)

library(tidyverse)

hubgene <- read.csv('../06_expression/03.hubgene.csv')
dat<-read.csv('../00_rawdata/dat(GSE16561).csv',row.names = 1,check.names = F)

group<-read.csv('../00_rawdata/group(GSE16561).csv')
colnames(group) <- c('sample', 'type')
table(group$type)

group$type = factor(group$type, levels = c("Control", "Stroke"))
levels(group$type)[levels(group$type) == "Stroke"] <- "IS"

hub_exp<-dat[hubgene$symbol,group$sample]
library(ROCR)
library(ggplot2)
hub_exp2<-t(hub_exp)
project <- 'ROC'

identical(rownames(hub_exp2), group$sample)

res <- data.frame()
library(pROC)
setwd('')
for (i in c(1:length(hubgene$symbol))) {
  roc<-roc(group$type,hub_exp2[,i],levels=c("Control", "IS"))
  print(paste0(colnames(hub_exp2)[i],':',roc$auc))
  roc$auc <- round(roc$auc,3)
  png(paste0('07_ROC/01.',colnames(hub_exp2)[i],'_',project,".png"),width = 4,height = 4,units = 'in',res = 600, family='Times')
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.3f',
       
       
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       
       main=paste0(colnames(hub_exp2)[i]),
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  pdf(paste0('07_ROC/01.',colnames(hub_exp2)[i],'_',project,".pdf"),width = 4,height = 4, family='Times')
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.3f',
       
       
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       
       main=paste0(colnames(hub_exp2)[i]),
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  roc <- data.frame(symbol=colnames(hub_exp2)[i], auc=roc$auc)
  res <- rbind(res,roc)
}

hubgene1 <-hubgene[-1, ]
hubgene1<-as.data.frame(hubgene1)
colnames(hubgene1)<-"symbol"
write.csv(hubgene1,file = "0.1_keygene.csv",row.names=FALSE)

rm(list = ls())

setwd('')
if (!dir.exists("./07_ROC")) {
  dir.create("./07_ROC")
}
setwd("./07_ROC")
# library(lance)

library(tidyverse)
hubgene <- read.csv('../06_expression/03.hubgene.csv')
dat<-read.csv('../00_rawdata/dat(GSE58294).csv',row.names = 1,check.names = F)

group<-read.csv('../00_rawdata/group(GSE58294).csv')
colnames(group) <- c('sample', 'type')
table(group$type)

group$type = factor(group$type, levels = c("Control", "Cardioembolic"))
levels(group$type)[levels(group$type) == "Cardioembolic"] <- "IS"

hub_exp<-dat[hubgene$symbol,group$sample]
library(ROCR)
library(ggplot2)
hub_exp2<-t(hub_exp)
project <- 'ROC'

identical(rownames(hub_exp2), group$sample)

res <- data.frame()
library(pROC)
setwd('')
for (i in c(1:length(hubgene$symbol))) {
  roc<-roc(group$type,hub_exp2[,i],levels=c("Control", "IS"))
  print(paste0(colnames(hub_exp2)[i],':',roc$auc))
  roc$auc <- round(roc$auc,3)
  png(paste0('07_ROC/02.',colnames(hub_exp2)[i],'_',project,".png"),width = 4,height = 4,units = 'in',res = 600, family='Times')
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.3f',
       
       
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       
       main=paste0(colnames(hub_exp2)[i]),
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  pdf(paste0('07_ROC/02.',colnames(hub_exp2)[i],'_',project,".pdf"),width = 4,height = 4, family='Times')
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.3f',
       
       
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       
       main=paste0(colnames(hub_exp2)[i]),
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  roc <- data.frame(symbol=colnames(hub_exp2)[i], auc=roc$auc)
  res <- rbind(res,roc)
}

