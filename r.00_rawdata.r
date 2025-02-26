
rm(list = ls())
setwd('')
if (! dir.exists('./00_rawdata')) {
  dir.create('./00_rawdata')
}
setwd("./00_rawdata")
getwd()
library(GEOquery)
library(Biobase)
library(tidyverse)
gset1<-getGEO('GSE58294',
              destdir = '.',
              GSEMatrix = T,
              getGPL = F)
expr1<-as.data.frame(exprs(gset1[[1]]))
gset1
gpl1<-getGEO("GPL570",destdir = '.')

gpl1<-Table(gpl1)
colnames(gpl1)

colnames(gpl1) <- gsub(" ", "_", colnames(gpl1))
colnames(gpl1)
gpl1$`Gene_Symbol`
colnames(gpl1)
str(gpl1)
head(gpl1)

library(tidyverse)
probe2symobl_1<-gpl1 %>%
  dplyr::select('ID','Gene_Symbol')%>%filter('Gene_Symbol'!='')%>%
 separate('Gene_Symbol',c("Gene_Symbol","drop"),sep = ' /// ')%>%
 dplyr::select(-drop)

probe2symobl_1=probe2symobl_1[probe2symobl_1$'Gene_Symbol'!='',] %>% na.omit()
probe2symobl_1 <- probe2symobl_1[!probe2symobl_1$'Gene_Symbol' == "previous version conserved probe", ]
dat_expr_1<-expr1
dat_expr_1$ID<-rownames(dat_expr_1)
dat_expr_1$ID<-as.character(dat_expr_1$ID)
probe2symobl_1$ID<-as.character(probe2symobl_1$ID)

dat_expr_1<-dat_expr_1 %>%
  inner_join(probe2symobl_1,by='ID')%>%
  dplyr::select(-ID)%>%     
  dplyr::select('Gene_Symbol',everything())%>%     
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(Gene_Symbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])   

# pcg <- read.delim2('/data/nas3/luchensi/data/PCG(1).xls(v22)')
pcg <- read.delim2('PCG(1).xls(v22)')
dat_expr_1 <- dat_expr_1[rownames(dat_expr_1)%in%pcg$gene_name,] %>% na.omit()

dim(dat_expr_1)
pd1<-pData(gset1[[1]])
colnames(pd1)
group1 <- data_frame(sample = rownames(pd1), type = pd1$title)
group1$group<- str_split_fixed(group1$type,"[_]",3)[,2]
group1<- group1[,-2]
  

identical(group1$sample, colnames(dat_expr_1))
write.csv(dat_expr_1, file = 'dat(GSE58294).csv',row.names = T)
write.csv(group1, file = 'group(GSE58294).csv',row.names = F)

rm(list = ls())

library(GEOquery)
library(Biobase)
library(tidyverse)
gset1<-getGEO('GSE16561',
              destdir = '.',
              GSEMatrix = T,
              getGPL = F)
expr1<-as.data.frame(exprs(gset1[[1]]))
gset1

dat_expr_1<-read_table("GSE16561_RAW.txt")
colnames(dat_expr_1)[1] <- "ID"

gpl1<-getGEO("GPL6883",destdir = '.')

gpl1<-Table(gpl1)
colnames(gpl1)

gpl1$`Symbol`
colnames(gpl1)
str(gpl1)
head(gpl1)

library(tidyverse)

probe2symobl_1<-gpl1 %>%
  dplyr::select('ID','Symbol')%>%filter('Symbol'!='')%>%
  separate('Symbol',c("Symbol","drop"),sep = ' /// ')%>%
  dplyr::select(-drop)

probe2symobl_1=probe2symobl_1[probe2symobl_1$'Symbol'!='',] %>% na.omit()
probe2symobl_1 <- probe2symobl_1[!probe2symobl_1$'Symbol' == "previous version conserved probe", ]

dat_expr_1$ID<-as.character(dat_expr_1$ID)

probe2symobl_1$ID<-as.character(probe2symobl_1$ID)

probe2symobl_1<-dat_expr_1 %>%
  inner_join(probe2symobl_1,by='ID')%>%
  dplyr::select(-ID)%>%     
  dplyr::select('Symbol',everything())%>%     
  mutate(rowMean=rowMeans(.[-1]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(Symbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])   

pcg <- read.delim2('PCG(1).xls(v22)')
probe2symobl_1 <- probe2symobl_1[rownames(probe2symobl_1)%in%pcg$gene_name,] %>% na.omit()
probe2symobl_1<- log2(probe2symobl_1+1)

a <- gset1[[1]]
pd <- pData(a)
group = data.frame(sample = pd$title,group = pd$title)
group$group <- str_split_fixed(group$group,"[_]",2)[,2]
identical(group$sample, colnames(probe2symobl_1))
write.csv(probe2symobl_1, file = 'dat(GSE16561).csv',row.names = T)
write.csv(group, file = 'group(GSE16561).csv',row.names = F)

