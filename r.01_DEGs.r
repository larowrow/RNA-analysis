rm(list = ls())
setwd('')
if (! dir.exists('./01_DEGs')){
  dir.create('./01_DEGs')
}
setwd('./01_DEGs')

library(magrittr)
library(stringr)
# library(lance)
BiocManager::install('limma')
library(limma)

dat_expr <- read.csv('../00_rawdata/dat(GSE58294).csv',row.names = 1)

group <- read.csv('../00_rawdata/group(GSE58294).csv')
colnames(group) <- c('sample','type')
dat_expr = dat_expr[,group$sample]
table(group$type)

group$type = factor(group$type, levels = c("Control", "Cardioembolic"))
levels(group$type)[levels(group$type) == "Cardioembolic"] <- "IS"
design.mat = cbind(control = ifelse(group$type == "Control", 1, 0), 
                   Cardioembolic = ifelse(group$type == "IS", 1, 0))

print(design.mat)

contrast.mat = makeContrasts(contrasts="Cardioembolic-control", levels=design.mat) 
contrast.mat

fit = lmFit(dat_expr, design.mat)

fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")

DEG=na.omit(fit)
logFC_cutoff <- 1.5
colnames(DEG)
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val<0.05& abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  
)
sig_diff <- subset(DEG,
                     DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)
dim(sig_diff)

summary(sig_diff$change)

DEG_write <- DEG[,-c(2,3,6)]
write.csv(DEG_write,file = "01.DEG_all.csv")
sig_diff_write <- sig_diff[,-c(2,3,6)]
write.csv(sig_diff_write, file = "01.DEG_sig.csv")

dat_rep<-DEG[rownames(DEG)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],5),
                              head(sig_diff[order(sig_diff$logFC,decreasing = F),],5))),]

library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)

volcano_plot<- ggplot(data = DEG,
                      aes(x = logFC,
                          y = -log10(adj.P.Val),
                          color =change)) +
  scale_color_manual(values = c("blue", "darkgray","red")) +
  
  
  
  geom_point(size = 2.4, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff),
             lty = 4,
             col = "gray40",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "gray40",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log2(Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
pdf('01.volcano.pdf',  w=8,h=6,family='Times')
volcano_plot
dev.off()
png('01.volcano.png',  w=8,h=6,family='Times', units='in',res=600)
volcano_plot
dev.off()

library(ggnewscale)
library(reshape2)
dat_rep<-dat_rep[order(dat_rep$change),]

str(group)
colnames(group) <- c('sample', 'type')
group <- group[order(group$type,decreasing = F),]
rt <- dat_expr
diff <- rt[rownames(dat_rep),group$sample] %>% data.frame()
diff <- data.frame(t(scale(t(diff))))
diff[diff < (-2)] <- (-2)
diff[diff > 2] <- 2

annotation_col <- data.frame(group=group$type)
rownames(annotation_col) <- colnames(diff)
annotation_col$group <- factor(annotation_col$group,levels = c("IS","Control"))

annotation_row <- data.frame(dat_rep$change)
rownames(annotation_row)=rownames(dat_rep)
colnames(annotation_row) <- " "

diff$gene <- rownames(diff)
diff <- melt(diff)
str(diff)

res <- diff %>% dplyr::filter(variable  == colnames(rt)[1])

{
  res$ang<-NA
  res$ang[1]<- -30
  for (i in (2:nrow(annotation_row))){
    res$ang[i]<-  res$ang[i-1]-c(360/24)
  }
  
  
  res$hjust <- 0
  res$hjust[which(res$ang < -90)] <- 1
  res$ang[which(res$ang < -90)] <- (180+res$ang)[which(res$ang < -90)]
  
  
  
  
  
  
  
  
  
  
  
  diff$var <- rep(1:nrow(annotation_row),nrow(annotation_col))
  range(diff$value)
  median(diff$value)
  
  
  
  
  
  
  
  
}

p <-  ggplot() +
  geom_bar(data = annotation_col,stat = 'identity',
           aes(x = 0.25,y = 1,fill = group),
           width = 0.5,
           color = NA) +
  scale_fill_manual(name = 'group',
                    values = c(Control="#FACF94FF",IS ="#D6E7A3"))+
  new_scale("fill") +
  geom_tile(data = diff[which(diff$variable == group$sample[1]),],
            aes_string(x = 1:nrow(annotation_row),y = 0.5,fill = diff[which(diff$variable == group$sample[1]),]$value),
            color = 'white')

for (i in 2:nrow(annotation_col)) {
  
  
  p <- p+geom_tile(data = diff[which(diff$variable == group$sample[i]),],
                   aes_string(x = 1:nrow(annotation_row),y = i-0.5,fill = diff[which(diff$variable == group$sample[i]),]$value),
                   color = 'white')
}

p1<- p+
  
  scale_fill_gradient2(name=" ",
                       midpoint = median(diff$value),
                       low = 'blue',
                       mid = "white",
                       high = 'red') +
  coord_polar(theta = 'x') +
  theme_void() +
  xlim(-4,nrow(annotation_row)+1) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.margin = margin(10,10, 10, 10))+
  geom_text(data = res,
            aes(x = as.numeric(rownames(res)),
                y = nrow(annotation_col)+1,
                label = gene, angle = ang, hjust = hjust),
            size = 3)+
  ylim(-8,nrow(annotation_col)+1)
p1
dev.off()

ggsave('04.heatmap.pdf', p1, width = 13, height = 7)
ggsave('04.heatmap.png', p1, width = 13, height = 7)

