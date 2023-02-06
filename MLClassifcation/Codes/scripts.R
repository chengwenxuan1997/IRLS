# 设置工作路径
work.path <- "i:/genomicdata/External/Xlu/MLClassifcation"; setwd(work.path) 

# 设置其他路径
code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")

# 如不存在这些路径则创建路径
if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

# BiocManager::install("mixOmics")
# BiocManager::install("survcomp")
# devtools::install_github("binderh/CoxBoost")
# install.packages("randomForestSRC")
# install.packages("snowfall")

# 加载需要使用的R包
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

# 加载模型训练以及模型评估的脚本
source(file.path(code.path, "ML.R"))

## Training Cohort ---------------------------------------------------------
# 训练集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与测试集保持相同类型，如同为SYMBOL或ENSEMBL等）
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_expr <- Train_expr[rowSums(Train_expr > 0) > ncol(Train_expr) * 0.1, ] # 剔除大量无表达的基因，以免建模过程报错
# 训练集生存数据是行为样本，列为结局信息的数据框
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_surv <- Train_surv[Train_surv$OS.time > 0, c("OS", "OS.time")] # 提取OS大于0的样本
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -------------------------------------------------------
# 测试集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与训练集保持相同类型，如同为SYMBOL或ENSEMBL等）
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
# 测试集生存数据是行为样本，列为结局信息的数据框
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_surv <- Test_surv[Test_surv$OS.time > 0, c("Cohort","OS", "OS.time")] # 提取OS.time大于0的样本
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

# 提取相同基因
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因

# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
# 此处记录需要运行的模型，格式为：算法1名称[算法参数]+算法2名称[算法参数]
# 目前仅有StepCox和RunEnet支持输入算法参数
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
  preTrain.var[[method]] = RunML(method = method, # 机器学习方法
                                 Train_set = Train_expr, # 训练集有潜在预测价值的变量
                                 Train_label = Train_surv, # 训练集生存数据
                                 mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                 classVar = classVar) # 用于训练的生存变量，必须出现在Train_surv中
}
preTrain.var[["simple"]] <- colnames(Train_expr)

## Model training
model <- list()
set.seed(seed = 123)
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method # 本轮算法名称
  method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称
  
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

# 读取已保存的模型列表
model <- readRDS(file.path(res.path, "model.rds"))
methodsValid <- names(model)

# 根据给定表达量计算样本风险评分
# linear Predictor
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                       new_data = Test_expr)
}

# 提取所筛选的变量（列表格式）
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- model[[method]]$subFeature
}

# 提取所筛选的变量（数据框格式）
fea_df <- lapply(model, function(fit){
  data.frame(fit$subFeature)
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file.path(res.path, "fea_df.txt"), # 两列，包含算法以及算法所筛选出的变量
            sep = "\t", row.names = F, col.names = T, quote = F)

# 对各模型计算C-index
AUC_list <- list()
for (method in methods){
  AUC_list[[method]] <- RunEval(fit = model[[method]], # 预后模型
                                Test_set = Test_expr, # 测试集预后变量，应当包含训练集中所有的变量，否则会报错
                                Test_label = Test_surv, # 训练集生存数据，应当包含训练集中所有的变量，否则会报错
                                Train_set = Train_expr, # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                Train_label = Train_surv, # 若需要同时评估训练集，则给出训练集生存数据，否则置NULL
                                Train_name = "TCGA", # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                cohortVar = "Cohort", # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                classVar = classVar) # 用于评估的生存状态，必须出现在Test_surv中
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(AUC_mat, file.path(res.path, "AUC_mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean)           # 计算每种算法在所有队列中平均C-index
avg_AUC <- sort(avg_AUC, decreasing = T)     # 对各算法C-index由高到低排序
AUC_mat <- AUC_mat[names(avg_AUC), ]      # 对C-index矩阵排序

fea_sel <- fea_list[[rownames(AUC_mat)[1]]] # 最优模型（测试集C指数均值最大）所筛选的特征

avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3)) # 保留三位小数


CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired") # 设置队列颜色
names(CohortCol) <- colnames(AUC_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(AUC_mat, # 主矩阵
                    avg_AUC, # 侧边柱状图
                    CohortCol, "steelblue", # 列标签颜色，右侧柱状图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path(fig.path, "AUC.pdf"), width = cellwidth * ncol(AUC_mat) + 3, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())
