# ���ù���·��
work.path <- "i:/genomicdata/External/Xlu/MLClassifcation"; setwd(work.path) 

# ��������·��
code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")

# �粻������Щ·���򴴽�·��
if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

# BiocManager::install("mixOmics")
# BiocManager::install("survcomp")
# devtools::install_github("binderh/CoxBoost")
# install.packages("randomForestSRC")
# install.packages("snowfall")

# ������Ҫʹ�õ�R��
library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)

# ����ģ��ѵ���Լ�ģ�������Ľű�
source(file.path(code.path, "ML.R"))

## Training Cohort ---------------------------------------------------------
# ѵ��������������Ϊ���򣨸���Ȥ�Ļ��򼯣�����Ϊ�����ı�����󣨻���������Լ�������ͬ���ͣ���ͬΪSYMBOL��ENSEMBL�ȣ�
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_expr <- Train_expr[rowSums(Train_expr > 0) > ncol(Train_expr) * 0.1, ] # �޳������ޱ���Ļ������⽨ģ���̱���
# ѵ����������������Ϊ��������Ϊ�����Ϣ�����ݿ�
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_surv <- Train_surv[Train_surv$OS.time > 0, c("OS", "OS.time")] # ��ȡOS����0������
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -------------------------------------------------------
# ���Լ�����������Ϊ���򣨸���Ȥ�Ļ��򼯣�����Ϊ�����ı�����󣨻�������ѵ����������ͬ���ͣ���ͬΪSYMBOL��ENSEMBL�ȣ�
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
# ���Լ�������������Ϊ��������Ϊ�����Ϣ�����ݿ�
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_surv <- Test_surv[Test_surv$OS.time > 0, c("Cohort","OS", "OS.time")] # ��ȡOS.time����0������
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

# ��ȡ��ͬ����
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # ����ģ�͵ı�������Ϊ��������Ϊ����
Test_expr <- t(Test_expr[comgene,]) # ����ģ�͵ı�������Ϊ��������Ϊ����

# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
# �˴���¼��Ҫ���е�ģ�ͣ���ʽΪ���㷨1����[�㷨����]+�㷨2����[�㷨����]
# Ŀǰ����StepCox��RunEnet֧�������㷨����
methods <- read.xlsx(file.path(code.path, "methods.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)
head(methods)

## Train the model --------------------------------------------------------

classVar = "OS"
Train_surv[[classVar]] <- ifelse(test = Train_surv[[classVar]] == "dead", yes = 1, no = 0)
Test_surv[[classVar]] <- ifelse(test = Test_surv[[classVar]] == "dead", yes = 1, no = 0)


## Pre-training 
Variable = colnames(Train_expr)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))
preTrain.method

preTrain.var <- list()
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, # ����ѧϰ����
                                 Train_set = Train_expr, # ѵ������Ǳ��Ԥ���ֵ�ı���
                                 Train_label = Train_surv, # ѵ������������
                                 mode = "Variable",       # ����ģʽ��Variable(ɸѡ����)��Model(��ȡģ��)
                                 classVar = classVar) # ����ѵ����������������������Train_surv��
}
preTrain.var[["simple"]] <- colnames(Train_expr)

## Model training
model <- list()
set.seed(seed = 123)
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method # �����㷨����
  method <- strsplit(method, "\\+")[[1]] # �������㷨����
  
  if (length(method) == 1) method <- c("simple", method)
 
  Variable = preTrain.var[[method[1]]]
  Train_set = Train_expr[, Variable]
  Train_label = Train_surv
  model[[method_name]] <- RunML(method = method[2], 
                                Train_set = Train_set, 
                                Train_label = Train_surv,
                                mode = "Model",
                                classVar = classVar)
}
saveRDS(model, file.path(res.path, "model.rds"))

## Evaluate the model -----------------------------------------------------

# ��ȡ�ѱ����ģ���б�
model <- readRDS(file.path(res.path, "model.rds"))
methodsValid <- names(model)

# ���ݸ�������������������������
# linear Predictor
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                       new_data = Test_expr)
}

# ��ȡ��ɸѡ�ı������б���ʽ��
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- model[[method]]$subFeature
}

# ��ȡ��ɸѡ�ı��������ݿ��ʽ��
fea_df <- lapply(model, function(fit){
  data.frame(fit$subFeature)
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file.path(res.path, "fea_df.txt"), # ���У������㷨�Լ��㷨��ɸѡ���ı���
            sep = "\t", row.names = F, col.names = T, quote = F)

# �Ը�ģ�ͼ���C-index
AUC_list <- list()
for (method in methods){
  AUC_list[[method]] <- RunEval(fit = model[[method]], # Ԥ��ģ��
                                Test_set = Test_expr, # ���Լ�Ԥ�������Ӧ������ѵ���������еı���������ᱨ��
                                Test_label = Test_surv, # ѵ�����������ݣ�Ӧ������ѵ���������еı���������ᱨ��
                                Train_set = Train_expr, # ����Ҫͬʱ����ѵ�����������ѵ���������ף�������NULL
                                Train_label = Train_surv, # ����Ҫͬʱ����ѵ�����������ѵ�����������ݣ�������NULL
                                Train_name = "TCGA", # ����Ҫͬʱ����ѵ�������ɸ���ѵ�����ı�ǩ�����򰴡�Training������
                                cohortVar = "Cohort", # ��Ҫ������ָ�����еı��������б��������ָ��[Ĭ��Ϊ��Cohort��]������ᱨ��
                                classVar = classVar) # ��������������״̬�����������Test_surv��
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(AUC_mat, file.path(res.path, "AUC_mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean)           # ����ÿ���㷨�����ж�����ƽ��C-index
avg_AUC <- sort(avg_AUC, decreasing = T)     # �Ը��㷨C-index�ɸߵ�������
AUC_mat <- AUC_mat[names(avg_AUC), ]      # ��C-index��������

fea_sel <- fea_list[[rownames(AUC_mat)[1]]] # ����ģ�ͣ����Լ�Cָ����ֵ�����ɸѡ������

avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3)) # ������λС��


CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired") # ���ö�����ɫ
names(CohortCol) <- colnames(AUC_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(AUC_mat, # ������
                    avg_AUC, # �����״ͼ
                    CohortCol, "steelblue", # �б�ǩ��ɫ���Ҳ���״ͼ��ɫ
                    cellwidth = cellwidth, cellheight = cellheight, # ��ͼÿ��ɫ��ĳߴ�
                    cluster_columns = F, cluster_rows = F) # �Ƿ�����н��о���

pdf(file.path(fig.path, "AUC.pdf"), width = cellwidth * ncol(AUC_mat) + 3, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())