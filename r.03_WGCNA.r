rm(list = ls())
library(dplyr)
library(magrittr)
# library(lance)
rm(list = ls())
setwd('')
if (! dir.exists('./03_WGCNA')){
  dir.create('./03_WGCNA')
}
setwd('./03_WGCNA')
dat_expr <- read.csv('../00_rawdata/dat(GSE58294).csv',row.names = 1,check.names = F) #%>% lance::lc.tableToNum()

group <- read.csv('../02_ssGSEA/02.ssgsea_result.csv',check.names = F,row.names = 1)

group<- as.data.frame(t(group))

str(group)

eset<-dat_expr

datExprOri=as.data.frame(t(eset))

nGenes0 = ncol(datExprOri)
nSamples0 = nrow(datExprOri)
nGenes0 ;nSamples0

library(WGCNA)
gsg = goodSamplesGenes(datExprOri, verbose = 3)  
gsg$allOK   

nGenes1<-ncol(datExprOri)      
nSamples1 <- nrow(datExprOri) 

nGenes1 ;nSamples1

tree=hclust(dist(datExprOri),method ='complete')   

pdf(file='01.sampleClustering.pdf',w=10,h=7)

par(pin = c(4, 4), mar = c(6, 6, 6, 1), family = "serif")

plot(tree,xlab="", sub="", main="Sample Clustering",
     labels=F,
     cex=1.0,  
     font=2,
     cex.axis=1.6,  
     cex.lab=1.8,   
     cex.main=1.8,   
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
abline(h =400, col = 'red') 
dev.off()

png(file='01.sampleClustering.png',w=10,h=7,units='in',res=600,bg='white')
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plot(tree,xlab="", sub="", main="Sample Clustering",
     labels=F,
     cex=1.0,  
     font=2,
     cex.axis=1.6,  
     cex.lab=1.8,   
     cex.main=1.8,   
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)

abline(h =400 ,col = 'red')
dev.off()

clust = cutreeStatic(tree, cutHeight = 400, minSize = 10)

table(clust)  

datExpr = datExprOri[clust == 1, ]
nGenes = ncol(datExpr)      
nSample =nrow(datExpr) 
nGenes ;nSample

SampleName<-rownames(datExpr)

condition <- group

datTraits<-condition
identical(rownames(datTraits), rownames(datExpr))
datExpr[1:4,1:4]

traitColors = numbers2colors(datTraits, signed = FALSE)

tree2=hclust(dist(datExpr),method ='complete')   

pdf(file='02.sampleClustering2.pdf',w=10,h=7)
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plotDendroAndColors(tree2,
                    traitColors,
                    dendroLabels = FALSE,
                    hang = 0.03,
                    cex.dendroLabels = 1,
                    
                    font=2,
                    guideHang = 0.05,
                    cex.axis=1.4,  
                    cex.lab=1.6,   
                    cex.main=1.6,   
                    font.axis = 2,
                    font.lab = 2,
                    font.main = 2,
                    font.sub =2,
                    cex.colorLabels=1,
                    groupLabels = colnames(datTraits),
                    main = "Sample Clustering and trait heatmap")
dev.off()

png(file='02.sampleClustering2.png',w=10,h=7,units='in',res=600,bg='white')
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plotDendroAndColors(tree2,
                    traitColors,
                    dendroLabels = FALSE,
                    hang = 0.03,
                    cex.dendroLabels = 1,
                    
                    font=2,
                    guideHang = 0.05,
                    cex.axis=1.4,  
                    cex.lab=1.6,   
                    cex.main=1.6,   
                    font.axis = 2,
                    font.lab = 2,
                    font.main = 2,
                    font.sub =2,
                    cex.colorLabels=1,
                    groupLabels = colnames(datTraits),
                    main = "Sample Clustering and trait heatmap")
dev.off()

powers <- c(seq(1, 10, by=1), seq(12, 20, by=2)) 

NUM <- 0.85  
enableWGCNAThreads()

sft = pickSoftThreshold(datExpr, RsquaredCut = NUM, powerVector = powers, verbose = 5)
sft$powerEstimate 

pdf('03.softThreshold.pdf',w=12,h=8)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif')
cex1=1.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     cex.axis=1.6,  
     cex.lab=1.8,   
     cex.main=1.8,   
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2,
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=NUM,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"),
     cex.axis=1.6,  
     cex.lab=1.8,   
     cex.main=1.8,   
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

png('03.softThreshold.png',w=12,h=8,units='in',res=600,bg='white')
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif')
cex1=1.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     cex.axis=1.6,  
     cex.lab=1.8,   
     cex.main=1.8,   
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2,
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=NUM,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"),
     cex.axis=1.6,  
     cex.lab=1.8,   
     cex.main=1.8,   
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

cor <- WGCNA::cor

net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,  
  minModuleSize =100,                    
  deepSplit = 4,                         
  mergeCutHeight = 0.25,                 
  numericLabels = TRUE,                  
  networkType  = "signed",               
  maxBlockSize = ncol(datExpr),
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "FPKM-TOM",
  loadTOMs = TRUE,
  verbose = 3
)

save.image('Modules.RData')
load('./Modules.RData')

cor<-stats::cor
table(net$colors)     

moduleColors = labels2colors(net$colors)
moduleColors

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); 

useMEs=subset(MEs, select = -c(MEgrey))  

mergedColors = labels2colors(net$colors)
mergedColors[net$blockGenes[[1]]]
pdf("04.wgcna.dendroColors.pdf",height = 7,width = 9)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif')
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    cex.axis=1.6,  
                    cex.lab=1.8,   
                    cex.main=1.8,   
                    font.axis = 2,
                    font.lab = 2,
                    font.main = 2,
                    cex.colorLabels=1,
                    font.sub =2)

dev.off()

png("04.wgcna.dendroColors.png",height = 7,width = 9,units='in',res=600)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif')
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    cex.axis=1.6,  
                    cex.lab=1.8,   
                    cex.main=1.8,   
                    font.axis = 2,
                    font.lab = 2,
                    font.main = 2,
                    cex.colorLabels=1,
                    font.sub =2)

dev.off()

load('FPKM-TOM-block.1.RData')
TOM=as.matrix(TOM)
save.image(file="WGCNA.RData")
load('WGCNA.RData')

moduleTraitCor =cor(useMEs,datTraits, use = "p") 
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, length(SampleName))

textMatrix =  paste(signif(moduleTraitCor, 3), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
class(moduleTraitCor)
colnames(moduleTraitCor) <- 'mtUPR-RGs'

pdf("05.wgcna.Module-trait.heatmap.pdf", width = 8, height =10)
par(mar = c(3, 12, 3, 1.5),family='serif')  
labeledHeatmap(Matrix = moduleTraitCor,   
               xLabels = 'mtUPR-RGs', 
               yLabels = names(useMEs),
               ySymbols = names(useMEs),
               colorLabels = FALSE,
               colors =  blueWhiteRed(50),
               textMatrix = textMatrix, 
               setStdMargins = FALSE,
               xLabelsAngle = 0,
               xLabelsAdj = 0.5,
               cex.text =1.1,cex.lab.y = 1.5,cex.legendLabel =1.5,cex.lab.x=1.5,
               zlim = c(-1,1),
               main = "Module-trait relationships",
               cex.lab = 1.2,
               font.lab.x = 2, font.lab.y = 2 )

dev.off()

png("05.wgcna.Module-trait.heatmap.png", width = 8, height =10,unit='in',res=600,bg='white')
par(mar = c(3, 12, 3, 1.5),family='serif')  
labeledHeatmap(Matrix = moduleTraitCor,   
               xLabels = 'mtUPR-RGs', 
               yLabels = names(useMEs),
               ySymbols = names(useMEs),
               colorLabels = FALSE,
               colors =  blueWhiteRed(50),
               textMatrix = textMatrix, 
               setStdMargins = FALSE,
               xLabelsAngle = 0,
               xLabelsAdj = 0.5,
               cex.text =1.1,cex.lab.y = 1.5,cex.legendLabel =1.5,cex.lab.x=1.5,
               zlim = c(-1,1),
               main = "Module-trait relationships",
               cex.lab = 1.2,
               font.lab.x = 2, font.lab.y = 2 )
dev.off()

load('WGCNA1.RData')

kME=signedKME(datExpr, useMEs, outputColumnName = "kME", corFnc = "cor", corOptions = "use = 'p'")
geneTraitSignificance  = as.data.frame(cor(datExpr, datTraits, use = "p"))
write.table(geneTraitSignificance,"GS.txt",quote = F,sep = "\t",row.names = T)
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 70))
names(geneTraitSignificance) = paste("GS.", colnames(datTraits), sep="")
names(GSPvalue) = paste("GS.", colnames(sample), sep="")
modNames = substring(names(MEs), 3)

modNames1<-c("turquoise")

for (module in modNames1){
  if(module== "grey"){ next }
  column = match(module, modNames); 
  pdf(paste("06.GS_MM.", module, ".pdf", sep=""),height = 6,width = 7,family='Times')
  
  moduleGenes = moduleColors==module;
  par(mfrow = c(1,1))
  verboseScatterplot(abs(kME[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance",
                     main = paste("Module membership vs gene significance

"),cex.main = 1.2, cex.lab = 1.2, pch=19,cex.axis = 1.2, col = module,font.lab=2)
  abline(v=0.4,lwd=2,col="red")
  abline(h=0.45,lwd=2,col="red")
  dev.off()
}

for (module in modNames1){
  if(module== "grey"){ next }
  column = match(module, modNames); 
  png(paste("06.GS_MM.", module, ".png", sep=""),height = 6,width = 7,family='Times',units='in',res=600)
  
  moduleGenes = moduleColors==module;
  par(mfrow = c(1,1))
  verboseScatterplot(abs(kME[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance",
                     main = paste("Module membership vs gene significance

"),cex.main = 1.2, cex.lab = 1.2, pch=19,cex.axis = 1.2, col = module,font.lab=2)
  abline(v=0.4,lwd=2,col="red")
  abline(h=0.45,lwd=2,col="red")
  dev.off()
}

column = match("turquoise", modNames);
moduleGenes = moduleColors=="turquoise";
moduleColors_plot <- data.frame(moduleColors)
write.csv(moduleColors_plot,file = "07.module_gene_colors.csv",quote = F,row.names = T)
write.csv(kME,file = "07.MM.csv",quote = F,row.names = T)

mid_keep = abs(kME[moduleGenes, column])>0.4 & abs(geneTraitSignificance[moduleGenes, 1]) >0.45
mid_kme<- kME[moduleGenes, ]
mid_gene<-rownames(mid_kme[mid_keep,])

write.csv(rownames(mid_kme),file = "08.turquoise_gene.all.csv")
write.csv(mid_gene,file = "08.turquoise_gene.sig.csv")

