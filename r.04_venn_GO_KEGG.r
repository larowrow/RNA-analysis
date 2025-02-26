
rm(list = ls())
setwd('')
if (! dir.exists("./04_venn_GO_KEGG")){
  dir.create("./04_venn_GO_KEGG/")
}
setwd("./04_venn_GO_KEGG/")

DEG <- read.csv("../01_DEGs/01.DEG_sig.csv")
DEG <- DEG$X

wgcna_turquoise <- read.csv("../03_WGCNA./08.turquoise_gene.sig.csv") 
wgcna_turquoise <- wgcna_turquoise$x

venn_gene <- intersect(DEG, wgcna_turquoise)
write.csv(venn_gene,file = '01.venn_gene.csv',row.names = F)

library(grid)
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install('ggvenn')
library(ggvenn)
library(ggplot2)
mydata<-list('WGCNA'=wgcna_turquoise,
             'DEGs'=DEG)
pdf('01.venn.pdf',w=7,h=7,family='Times')
ggvenn(mydata,c('WGCNA','DEGs'),
       fill_color = c("#49afca","#fed53b"),
       show_percentage = T,
       stroke_alpha = 0.5,
       stroke_size = 1,
       text_size = 5,
       stroke_color="white",
       stroke_linetype="dashed",
       set_name_color='black',
       text_color = 'black')
dev.off()

png('01.venn.png',w=7,h=7,units = 'in',res = 600,family='Times')
ggvenn(mydata,c('WGCNA','DEGs'),
       fill_color = c("#49afca","#fed53b"),
       show_percentage = T,
       stroke_alpha = 0.5,
       stroke_size = 1,
       text_size = 5,
       stroke_color="white",
       stroke_linetype="dashed",
       set_name_color='black',
       text_color = 'black')
dev.off()

rm(list = ls())
library(dplyr)
library(magrittr)
# library(lance)
# rm(list = ls())
# setwd('')
# if (! dir.exists('./04_venn_GO_KEGG')){
#   dir.create('./04_venn_GO_KEGG')
# }
# setwd('./04_venn_GO_KEGG')

library(clusterProfiler)
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install('clusterProfiler')
# install.packages("igraph")
# BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)

gene <- read.csv('01.venn_gene.csv')

gene_transform <- bitr(gene$x, 
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")

ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE)

go_result <- data.frame(ego)

go_result <- go_result[go_result$pvalue<0.05,]

go_result <- go_result[order(go_result$pvalue,decreasing = F),]
dim(go_result)

table(go_result$ONTOLOGY)

go_result.write <- go_result[,-c(4,5,8)]
write.csv(go_result.write,file = "02.DEGs_WGCNA_GO.csv",quote = T,row.names = F)

display_number = c(5, 5, 5)  

BP <- go_result[which(go_result$ONTOLOGY=='BP'),]
BP <- BP[order(BP$pvalue, decreasing = F),,drop=F] 

CC <- go_result[which(go_result$ONTOLOGY=='CC'),]
CC <- CC[order(CC$pvalue, decreasing = F),,drop=F] 

MF <- go_result[which(go_result$ONTOLOGY=='MF'),]
MF <- MF[order(MF$pvalue, decreasing = F),,drop=F] 

display_number = c(5, 5, 5)
go_result_BP = as.data.frame(BP)[1:display_number[1], ]
go_result_CC = as.data.frame(CC)[1:display_number[2], ]
go_result_MF = as.data.frame(MF)[1:display_number[3], ]

go_enrich = data.frame(
  ID=c(go_result_BP$ID, go_result_CC$ID, go_result_MF$ID),  
  Description=c(go_result_BP$Description,go_result_CC$Description,go_result_MF$Description),
  Count=c(go_result_BP$Count, go_result_CC$Count, go_result_MF$Count), 
  type=factor(c(rep("Biological Process", display_number[1]), 
                rep("Cellular Component", display_number[2]),
                rep("Molecular Function", display_number[3])),
              levels=c("Biological Process", "Cellular Component","Molecular Function" )))

go_enrich$type_order = factor(go_enrich$Description,levels=go_enrich$Description,ordered = T)
head(go_enrich)
p <- ggplot(go_enrich,
            aes(x=type_order,y=Count, fill=type)) +  
  geom_bar(stat="identity", width=0.8) +  
  scale_fill_manual(values = c("cadetblue1", "#FF6666", "#b2e7cb") ) + 
  xlab("") + 
  ylab("Counts") +  
  labs(title = "GO Terms Enrich")+ 
  theme_bw() +
  theme(text=element_text(family="Times",size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")+
  theme(axis.text.x=element_text(family="Times",face = "bold", color="black",angle = 80,vjust = 1, hjust = 1 )) 
p

ggsave(filename = '02.DEGs_WGCNA_RRGs_GO.pdf',w=10,h=10)
ggsave(filename = '02.DEGs_WGCNA_RRGs_GO.png',w=10,h=10,units='in', dpi = 600)
dev.off()

kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)

kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kk_result <- data.frame(kk)
kk_result <- kk_result[kk_result$pvalue<0.05,]
kk_result <- kk_result[order(kk_result$pvalue,decreasing = F),]
dim(kk_result)

kk_result.write <- kk_result[,-c(3,4,7)]
write.csv(kk_result.write,file = "DEGs_WGCNA_RRGs_KEGG.csv",quote = T,row.names = F)

colnames(kk_result)
kk_result <- head(kk_result,12)
colnames(kk_result)
kk_result <- head(kk_result,14)
p <- ggplot(kk_result,aes(x=-log10(pvalue),y=Description))+
  xlim(1,2.0)+
  geom_point(aes(size=Count,color=Description),alpha=2)+
  
  
  
  
  
  scale_size(range = c(5,10))+
  theme_bw()+
  theme(legend.position = c('none'))+
  ggrepel::geom_text_repel(aes(label=Description),size=3,segment.color='black',show.legend = T)+
  theme(axis.text.y = element_blank())
p
pdf('03.DEGs_WGCNA_RRGs_KEGG.pdf',w=5,h=4,family='Times')
p
dev.off()
png('03.DEGs_WGCNA_RRGs_KEGG.png',w=5,h=4,family='Times',units='in',res=600)
p
dev.off()

