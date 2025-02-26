rm(list = ls())
setwd('')
if (! dir.exists("./08_nomogram")){
  dir.create("./08_nomogram")
}
setwd("./08_nomogram")
# library(lance)

library(tidyverse)

hubgene <- read.csv('../07_ROC/0.1_keygene.csv')
dat<-read.csv('../00_rawdata/dat(GSE58294).csv',row.names = 1,check.names = F)

group<-read.csv('../00_rawdata/group(GSE58294).csv')
colnames(group) <- c('sample', 'type')
group$type = factor(group$type, levels = c("Control", "Cardioembolic"))
levels(group$type)[levels(group$type) == "Cardioembolic"] <- "IS"
table(group$type)
group$y<-ifelse(group$type=='IS',1,0)
group$type = factor(group$type, levels = c("Control", "IS"))
hub_exp<-dat[hubgene$symbol,group$sample]
hub_exp2<-t(hub_exp)
d <- merge(hub_exp2, group, by.x = "row.names", by.y = "sample")
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install('rms')
library(rms)
ddist <- datadist(d)
options(datadist='ddist')

f <- lrm(y~ MCEMP1+CACNA1E+CLEC4D,
         
         data = d, x = TRUE, y = TRUE)
f
summary(f)

nom.cox <- nomogram(f,  
                    fun=plogis, 
                    funlabel="Risk of IS",
                    lp=F , 
                    fun.at=c(0.1,seq(0.3,0.7,by=0.2),0.9)  
)
plot(nom.cox, cex.axis  = 1.2, cex.var = 1.7)

png(filename = "01.nomogram_line_points.png", height = 5, width = 10,units='in',res=600, family = "Times")
par(family = "Times",font=2)
plot(nom.cox, cex.axis  = 0.8 ,cex.var = 1.5)
dev.off()
pdf(file = "01.nomogram_line_points.pdf", height = 5, width = 10, family = "Times")
par(family = "Times",font=2)
plot(nom.cox, cex.axis  = 0.8, cex.var = 1.5)
dev.off()

library(ResourceSelection)

hl1 <- ResourceSelection::hoslem.test(d$y,predict(f,d),g=10)
hl1

hl2 <- stats::resid(f,"gof")
hl2
pval <- signif(hl2[5], 3)
pval

set.seed(3)
cal1 <- calibrate(f, cmethod='KM', method='boot', B=1000)

pdf("02.calibrate.pdf",width=9,height=6, family = "Times")
par(mar = c(6,5,2,2))
plot(cal1, lwd=2, lty=1, 
     cex.lab=1.4, cex.axis=1.4, cex.main=1.8, cex.sub=1.4, 
     xlim=c(0, 1), ylim= c(0, 1), 
     xlab="Nomogram-Predicted Probability", 
     ylab="Actual Probability", 
     col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
     legend=FALSE)
lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
abline(0, 1, lty=3, lwd=2) 
legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"), cex=1.2,
       lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")

text(x = 0.34, y = 0.8, as.expression(bquote(italic('p')==.(pval))),cex=1.0)
dev.off()
png("02.calibrate.png", width=9,height=6, units = "in", res = 600, family = "Times")
par(mar = c(6,5,2,2))
plot(cal1, lwd=2, lty=1, 
     cex.lab=1.4, cex.axis=1.4, cex.main=1.8, cex.sub=1.4, 
     xlim=c(0, 1), ylim= c(0, 1), 
     xlab="Nomogram-Predicted Probability", 
     ylab="Actual Probability", 
     col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
     legend=FALSE)
lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
abline(0, 1, lty=3, lwd=2) 
legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"), cex=1.2,
       lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")

text(x = 0.34, y = 0.8, as.expression(bquote(italic('p')==.(pval))),cex=1.0)
dev.off()

predicted<-predict(f,newdata = d)

library(PRROC)
roc_curve = roc.curve(predicted, weights.class0 =  d$type == "IS", curve = T)

pdf(file='03.roc.pdf', width = 8, height = 8, family='Times')
par(pin = c(4,4), mar = c(6,6,6,1))
plot(roc_curve, auc.main = T, legend = F, color = 'darkblue', xlab = "1-Specificity", asp = 1 ,cex.axis=1.8,  
     cex.lab=2.0,   
     cex.main=2.0,   
     main='',
     font.lab = 2,
     font.main = 2,
     font.sub =2)
abline(0,1); dev.off()
png(file='03.roc.png', width = 8, height = 8, res = 600, units = "in", bg = "white",family='Times')
par(pin = c(4,4), mar = c(6,6,6,1))
plot(roc_curve, auc.main = T, legend = F, color = 'darkblue', xlab = "1-Specificity", asp = 1 ,cex.axis=1.8,  
     cex.lab=2.0,   
     cex.main=2.0,   
     main='',
     font.lab = 2,
     font.main = 2,
     font.sub =2)
abline(0,1); dev.off()

 
 
library(ggDCA)
data  <-ggDCA:: dca(f,
                     model.names =c('Cardioembolic'))

 # options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
 # if (!require("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")
 # BiocManager::install('rmda')
library(rmda)
d$predicted<-predicted
Nomogram<-decision_curve(y~ predicted,data= d,
                         family = binomial(link ='logit'),
                         thresholds= seq(0,0.8, by = 0.01),
                         confidence.intervals = 0.95,
                         study.design = 'case-control',
                         population.prevalence = 0.3)

List<- list(Nomogram)
pdf('04.dca1.pdf',w=10,h=9,family='Times')
par(pin = c(4,4), mar = c(6,6,6,4),font.lab=2,cex.lab=2,cex.axis=2)
plot_decision_curve(List,
                    curve.names=c('Nomogram'),
                    cost.benefit.axis =FALSE,
                    
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()

png('04.dca1.png',w=10,h=9,family='Times',units='in',res=600)
par(pin = c(4,4), mar = c(6,6,6,4),font.lab=2,cex.lab=2,cex.axis=2)
plot_decision_curve(List,
                    curve.names=c('Nomogram'),
                    cost.benefit.axis =FALSE,
                    
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()

