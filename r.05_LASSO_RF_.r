
rm(list = ls())

setwd('')
if (!dir.exists("./05_LASSO_RF_")) {
  dir.create("./05_LASSO_RF_")
}
setwd("./05_LASSO_RF_")

library(glmnet)
library(dplyr)
library(ggplot2)

intersect <- read.csv('../04_venn_GO_KEGG/01.venn_gene.csv')
names(intersect) <- c("gene_name")

dat_expr<-read.csv("../00_rawdata/dat(GSE58294).csv", row.names = 1,check.names = F) # %>% lance::lc.tableToNum()
group<-read.csv('../00_rawdata/group(GSE58294).csv')
colnames(group) <- c('sample', 'type')
table(group$type)
group$type <- ifelse(group$type == 'Control', 'Control', 'IS')
dat <- dat_expr[intersect$gene_name, group$sample] %>% t() %>% as.data.frame()
dat <- merge(dat, group, by.x = 'row.names', by.y = 'sample') %>% tibble::column_to_rownames(var = 'Row.names')

dat$type <- factor(dat$type, levels = c('Control', 'IS'))

set.seed(3)

cvfit <- cv.glmnet(as.matrix(dat[-ncol(dat)]),
                   dat$type,
                   nfold = 5,
                   family = "binomial")

coef.min <- coef(cvfit, s = "lambda.min")
active.min <- which(coef.min@i != 0)
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i + 1]

cv_model  <- cvfit
best_lambda <- cv_model$lambda.min
lambda_lse <- cv_model$lambda.1se

lambda_values <- cv_model$lambda
coef_list <- lapply(lambda_values, function(lambda) coef(cv_model, s = lambda))
coef_matrices <- lapply(coef_list, function(coef) as.matrix(coef))
coef_df_list <- lapply(coef_matrices, function(matrix) as.data.frame(t(matrix)))
coef_df <- do.call(rbind, coef_df_list)
coef_df$lambda <- rep(lambda_values, each = nrow(coef_df_list[[1]]))
coef_df <- coef_df[, -1]

colnames(coef_df) <- c(colnames(coef_df)[-ncol(coef_df)], "lambda")
library(reshape2)
coef_long <- melt(coef_df, id.vars = "lambda", variable.name = "Variable", value.name = "Coefficient")
coef_long$log_lambda <- log(coef_long$lambda)

library(ggplot2)

table(coef_long$Variable)
p1 <- ggplot(coef_long, aes(x = log_lambda, y = Coefficient, color = Variable)) + 
  geom_vline(xintercept = log(best_lambda), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_vline(xintercept = log(lambda_lse), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_line(size = 1) +
  xlab("Log Lambda") + 
  ylab("Coefficients") + 
  theme_bw(base_rect_size = 2) +
  
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 15, color = 'black'), 
        axis.text = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'), 
        legend.position = 'right') +
  annotate('text', x = log(best_lambda),
           y = 2.5, 
           label = paste0('log(lambda.min)\n', round(log(best_lambda), 4)), 
           color = 'black', size = 5, family = 'serif') +
  annotate('text', x = log(lambda_lse),
           y = 1.7, 
           label = paste0('log(lambda.lse)\n', round(log(lambda_lse), 4)), 
           color = 'black', size = 5, family = 'serif') +
  guides(col = guide_legend(ncol = 4))
p1
ggsave('01.lasso_coefficients_venalty.pdf',p1,w=11,h=7)
ggsave('01.lasso_coefficients_venalty.png',p1,w=11,h=7)

xx <- data.frame(
  lambda = cvfit$lambda,
  cvm = cvfit$cvm,
  cvsd = cvfit$cvsd,
  cvup = cvfit$cvup,
  cvlo = cvfit$cvlo,
  nozezo = cvfit$nzero
)
xx$ll <- log(xx$lambda)
xx$NZERO <- paste0(xx$nozezo, ' vars')
dev.off()

p2 <- ggplot(xx, aes(ll, cvm, color = NZERO)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.05, size = 1) +
  geom_vline(xintercept = log(best_lambda), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_vline(xintercept = log(lambda_lse), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_point(size = 2) +
  xlab("Log Lambda") +
  ylab('Partial Likelihood Deviance') +
  theme_bw(base_rect_size = 1.5) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 15, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = 'black'),
    legend.position = 'bottom'
  ) +
  annotate('text', x = log(best_lambda), y = 0.6, 
           label = paste0('log(lambda.min)\n', round(log(best_lambda), 4)), 
           color = 'black', size = 4, family = 'serif') +
  annotate('text', x = log(lambda_lse), y =0.8, 
           label = paste0('log(lambda.1se)\n', round(log(lambda_lse), 4)), 
           color = 'black', size = 4, family = 'serif') +
  guides(col = guide_legend(ncol = 9))
p2
ggsave('02.lasso_model.pdf',p2,w=11,h=7)
ggsave('02.lasso_model.png',p2,w=11,h=7)

lassogenes <- as.data.frame(as.data.frame(lasso_geneids)[2:16,])

colnames(lassogenes) <- 'gene'
write.csv(lassogenes,file = '03.lasso_gene.csv',row.names = F,quote = F)

rm(list = ls())

# setwd('')
# if (!dir.exists("./05_LASSO_RF_")) {
#   dir.create("./05_LASSO_RF_")
# }
# setwd("./05_LASSO_RF_/")

library(tidyverse)
library(randomForest)
dat_expr <- read.csv('../00_rawdata/dat(GSE58294).csv',row.names = 1)
group <- read.csv('../00_rawdata/group(GSE58294).csv')
colnames(group) <- c('sample', 'type')
table(group$type)
group$type <- ifelse(group$type == 'Control', 'Control', 'IS')
gene <- read.csv('../04_venn_GO_KEGG/01.venn_gene.csv')

dat <- data.frame(t(dat_expr[rownames(dat_expr) %in% gene$x,]))
dat$sample <- rownames(dat)
dat <- merge(dat, group, by = 'sample')
rownames(dat) <- dat$sample 
dat <- dat %>% dplyr::select(-sample)
dat$type<-factor(dat$type,levels = c('Control','IS'))
str(dat)
set.seed(1)
cv_error <- c()
ntree_range <- seq(5, 100, by = 2)
for (ntree in ntree_range) {
  rf_model <- randomForest(type ~ ., data = dat, ntree = ntree, mtry=5,
                           importance = TRUE)
  oob_error <- rf_model$err.rate[ntree, "OOB"]
  cv_error <- c(cv_error, oob_error)
}

cv_results <- data.frame(ntree = ntree_range, error = cv_error)

p <- ggplot(cv_results, aes(x = ntree, y = error)) +
  geom_line(color = "#1f77b4", size = 1) +  
  geom_point(color = "#d62728", size = 3, shape = 21, fill = "white") +  
  labs(title = "Cross-Validation Error vs. Number of Trees",
       x = "Number of Trees",
       y = "OOB Error Rate") +
  geom_vline(xintercept = 83, size = 0.8, color = '#EB4B17', linetype = 2) +
  theme_bw() +
  annotate("text", x = 10, y = 0.005, label = "ntree = 83 ", colour = "#EB4B17",
           size = 3.5, fontface = "bold") +
  theme(
    text = element_text(size = 14),  
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.title.x = element_text(size = 12, face = "bold"),  
    axis.title.y = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),  
    panel.grid.major = element_line(color = "#e0e0e0", linetype = "dashed"),  
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "#f9f9f9")  
  )
ggsave('04_RF_select_ntree.png', p, w = 12, h = 6, dpi = 600)
ggsave('04_RF_select_ntree.pdf', p, w = 12, h = 6)

optimal_ntree <- cv_results$ntree[which.min(cv_results$error)]
cat("Optimal number of trees:", optimal_ntree, "\n")
final_rf_model <- randomForest(type ~ ., data = dat, ntree = optimal_ntree, importance = TRUE)
print(final_rf_model)

importance_scores <- importance(final_rf_model)
importance_df <- data.frame(Feature = rownames(importance_scores),
                            Importance = importance_scores[, "MeanDecreaseGini"])
importance_df[order(importance_df$Importance),]

write.csv(importance_df,file = '02.RF_importance.csv',row.names = F,quote = F)

importance_df_selected <- importance_df[order(-importance_df$Importance), ][c(1:8),]
write.csv(importance_df_selected[c(1:8),],file = '02.RF_Gene.csv',row.names = F,quote = F)

library(ggplot2)
p1 <- ggplot(importance_df_selected, aes(x = reorder(Feature, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_gradient(low = "#92C5DE", high = "#DC0000B2") +  
  labs(x = "Gene", y = "Importance", title = "Feature Importance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),  
    axis.title.x = element_text(size = 14, face = "bold"),  
    axis.title.y = element_text(size = 14, face = "bold"),  
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    panel.grid.major = element_line(color = "#e0e0e0", linetype = "dashed"),  
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "#f9f9f9")  
  )
p1
ggsave('02_RF_importance.png', p1, w = 11, h = 6, dpi = 600)
ggsave('02_RF_importance.pdf', p1, w = 11, h = 6)

# setwd('05_LASSO_RF_')
LASSO <- read.csv('03.lasso_gene.csv')
LASSO<-na.omit(LASSO)
LASSO<-LASSO$gene

randomForest<- read.csv('02.RF_Gene.csv')
randomForest<-randomForest$Feature
DEERG<- intersect(LASSO,randomForest)

write.table(DEERG,file = '03.Intersection_gene.csv',row.names = F,quote = F)

library(grid)
library(ggvenn)
library(ggplot2)
mydata<-list('LASSO'=LASSO,
             'RF'=randomForest  )
pdf('03.Intersection_gene.pdf',w=7,h=7,family='Times')
ggvenn(mydata,c('LASSO', 'RF'),
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
png('03.Intersection_gene.png',w=7,h=7, units = 'in',res = 600,family='Times')
ggvenn(mydata,c('LASSO', 'RF'),
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

