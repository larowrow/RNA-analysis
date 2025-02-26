rm(list = ls())
setwd('')
if (! dir.exists("./09_GSEA")){
  dir.create("./09_GSEA")
}

setwd("./09_GSEA")
library(magrittr)
expr <- read.csv(file = "../00_rawdata/dat(GSE58294).csv", check.names = F, row.names = 1) #%>% lance::lc.tableToNum()
gene <- data.frame(gene = c( "MCEMP1","CACNA1E","CLEC4D"))
gene <- gene$gene
group <- read.csv('../00_rawdata/group(GSE58294).csv')
colnames(group) <- c('sample', 'type')
group$type = factor(group$type, levels = c("Control", "Cardioembolic"))
levels(group$type)[levels(group$type) == "Cardioembolic"] <- "IS"
gene_expr <- expr[gene,group$sample] %>% t %>% as.data.frame
allgene_expr <- expr[,group$sample] %>% t %>% as.data.frame
identical(rownames(gene_expr),rownames(allgene_expr))
cor_res <- cor(gene_expr, allgene_expr,method = "spearman") %>% t %>% as.data.frame 

# BiocManager::install('clusterProfiler')

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

j <- 0
for (num in 1:length(gene)){
  i <- gene[num]
  
  print(i)
  j <- j+1
  dat <- data.frame(row.names = rownames(cor_res),
                    cor=cor_res[,i])
  data_sort <- dat %>% dplyr::arrange(desc(cor)) %>% na.omit()
  
  geneList <- data_sort$cor
  names(geneList) <- rownames(data_sort) 
  
  

kegg_set<- read.gmt("c2.cp.kegg.v2023.1.Hs.symbols.gmt")
  
kegg_list <- kegg_set %>% split(.$term) %>% lapply( "[[", 2)
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 1)
kegg_result <- kegg_gsea@result
kegg_result <- kegg_result[kegg_result$p.adjust < 0.05,]
kegg_result <- kegg_result[order(abs(kegg_result$NES), decreasing = T),]
kegg_result.write <- kegg_result[,-c(8,9,10)]
print(dim(kegg_result.write))
write.csv(kegg_result.write,paste0('0',j,'.',i,'.GSEA(kegg)_result.csv'),row.names = F)
  
  
terms <- kegg_result$ID[1:3]
library(GseaVis)
  
p = gseaNb(object = kegg_gsea,
             geneSetID = terms,
             curveCol = (rainbow(10)),
             subPlot = 3,
             termWidth = 50,
             legend.position = 'top',
             
             pvalX = 0.05,pvalY = 0.05)
  
  fn1 = paste0("GSE58294", j, ".", i, ".png")
  fn2 = paste0("GSE58294", j, ".", i,".pdf")
  png(fn1,w=10,h=6,units = "in",res = 600, family='Times')
  print(p)
  dev.off()
  pdf(fn2,w=10,h=6,family='Times')
  print(p)
  dev.off()
  
}

CACNA1E<-read_csv("02.CACNA1E.GSEA(kegg)_result.csv")
CACNA1E_1<-data.frame(CACNA1E = CACNA1E$ID)

mce<-read_csv("01.MCEMP1.GSEA(kegg)_result.csv")
mce_1<-data.frame(mce = mce$ID)

cle<-read_csv("03.CLEC4D.GSEA(kegg)_result.csv")
cle_1<-data.frame(cle = cle$ID)                  

intersection_df <- intersect(intersect(CACNA1E_1, mce_1), cle_1)  
print(intersection_df)
print(getwd())  
