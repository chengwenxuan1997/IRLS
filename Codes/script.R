# work.path <- "D:/Other/External/IRLS";setwd(work.path)
work.path <- "j:/Other/External/IRLS";setwd(work.path)
code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")
pkg.path <- file.path(work.path, "Packages")
ann.path <- file.path(work.path, "Annotations")

if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(pkg.path)) dir.create(pkg.path)
if (!dir.exists(ann.path)) dir.create(ann.path)

# BiocManager::install("mixOmics")
# BiocManager::install("survcomp")
# devtools::install_github("binderh/CoxBoost")
# devtools::install_github('seandavi/GEOquery')
# install.packages("randomForestSRC")
# install.packages("snowfall")
# devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")


library(GEOquery)
library(openxlsx)
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(rtracklayer)
library(seqinr)
library(plyr)
library(AnnotationDbi)
library(sva)
library(survival)
library(randomForestSRC)
library(SummarizedExperiment)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)

source(file.path(code.path, "ML.R"))


# Data Preparation --------------------------------------------------------

org.db <- org.Hs.eg.db
#
protein_pseudo <- c("IG_C_pseudogene", "IG_V_pseudogene", "polymorphic_pseudogene",
                    "processed_pseudogene", "protein_coding", "pseudogene", "TR_V_pseudogene",
                    "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene",
                    "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene",
                    "unitary_pseudogene", "unprocessed_pseudogene")

## annotation --------------------------------------------------------
# 获取TCGA数据的基因注释，CRC队列的注释为Data Release 32~25，且官网声称现在使用的是Dr 36
# 所以此处使用Dr36的gtf，可从https://www.gencodegenes.org/human/releases.html进行下载
# gtf文件存放在Annotations文件夹下
# 运行一次即可
# 读取gtf --> 根据ensembl号进行去重 --> 仅保留重要列 -->输出为txt
# 结果存放在Annotations/TCGA.anno.txt中
# TCGA.anno <- import(file.path(ann.path, "gencode.v36.annotation.gtf"))
# TCGA.anno <- as.data.frame(TCGA.anno)
# TCGA.anno <- TCGA.anno[!duplicated(TCGA.anno$gene_id),
#                        c("gene_id", "gene_type", "gene_name")]
# TCGA.anno$gene_id <- do.call(rbind, strsplit(TCGA.anno$gene_id, "\\."))[, 1]
# write.table(TCGA.anno, file.path(ann.path, "TCGA.anno.txt"),
#             sep = "\t", col.names = T, row.names = F, quote = F)
#

# 对微阵列数据进行重注释，使用与gencode.v36.annotation.gtf配对的gencode.v36.transcripts.fa文件
# fa文件的下载地址和存放地址均与gtf文件一致
# 探针fasta文件可从GEO数据库或者生产商官网上获得
# 本芯片fasta文件从https://www.thermofisher.cn/cn/zh/home.html获得，放在Annotations文件夹下
# 运行一次即可
# system(
#   paste(file.path(pkg.path, "seqmap-1.0.12-windows.exe"), 0,  # seqmap软件位置
#         file.path(ann.path, "HG-U133_Plus_2.probe_fasta"),    # 探针fasta文件位置
#         file.path(ann.path, "gencode.v36.transcripts.fa"),    # gencode文件位置
#         file.path(ann.path, "probe.reannotation"),            # 重注释结果
#         "/output_all_matches")
# )

# 对seqmap的结果进行处理
# 提取比对结果的探针名称(probe_id的第三列)，对应基因信息(trans_id)
# 结果存放于Annotations/probe.anno.txt中
# probe.anno <- read.delim(file.path(ann.path, "probe.reannotation"), sep = "\t")
# probe.anno$probe <- do.call(rbind, strsplit(probe.anno$probe_id, ":|;"))[, 3]
# probe.anno <- unique(probe.anno[, c("probe", "trans_id")])
# tmp <- do.call(rbind, strsplit(x = probe.anno$trans_id, split = "\\|"))
# probe.anno$ENSEMBL <- tmp[, 2]
# probe.anno$transcript <- tmp[, 5]
# probe.anno$GeneSymbol <- tmp[, 6]
# probe.anno$EntrezID <- tmp[, 7]
# probe.anno$BioType <- tmp[, 8]
# probe.anno <- probe.anno[, c("probe", "ENSEMBL", "transcript", "GeneSymbol", "EntrezID", "BioType")]
# probe.anno$ENSEMBL <- do.call(rbind, strsplit(probe.anno$ENSEMBL, "\\."))[, 1]
# write.table(probe.anno, file.path(ann.path, "probe.anno.txt"),
#             sep = "\t", col.names = T, row.names = F, quote = F)

# 筛选出TCGA队列中属于lncRNA的基因
# 结果存放于Annotatins/TCGA.lnc.anno.txt
# TCGA.anno <- read.table(file.path(ann.path, "TCGA.anno.txt"), header = T)
# table(TCGA.anno$gene_type)[c("protein_coding", "lncRNA")]
# TCGA.lnc.anno <- TCGA.anno[TCGA.anno$gene_type == "lncRNA", ]
# write.table(TCGA.lnc.anno, file.path(ann.path, "TCGA.lnc.anno.txt"),
#             sep = "\t", col.names = T, row.names = F, quote = F)

# 筛选出微阵列队列中属于lncRNA的基因
# 需要删除被比对为编码基因或同时比对到多个lncRNA的探针
# 结果存放于probe.lnc.anno.txt
# probe.anno <- read.table(file.path(ann.path, "probe.anno.txt"), header = T)
# table(probe.anno$BioType)[c("protein_coding", "lncRNA")]
# protein_probe <- probe.anno$probe[probe.anno$BioType %in% protein_pseudo] #被比对为编码基因的探针
# probe.anno <- probe.anno[!probe.anno$probe %in% protein_probe, ]
# probe.anno <- probe.anno[probe.anno$BioType == "lncRNA", ]
# duplicate_probe <- probe.anno$probe[duplicated(probe.anno$probe)] # 被比对到多个lncRNA的探针
# probe.lnc.anno <- probe.anno[!probe.anno$probe %in% duplicate_probe, ]
# write.table(probe.lnc.anno, file.path(ann.path, "probe.lnc.anno.txt"),
#             sep = "\t", col.names = T, row.names = F, quote = F)


## Signature ---------------------------------------------------------------

# 以原文使用的43个lncRNA为例
sig <- read.table(file.path(res.path, "signature.txt"))
sig

TCGA.lnc.anno <- read.table(file.path(ann.path, "TCGA.lnc.anno.txt"), header = T)
probe.lnc.anno <- read.table(file.path(ann.path, "probe.lnc.anno.txt"), header = T)
sig <- sig$V2
# 如果不保留比对到多个lncRNA的探针，仅有18个可以匹配（保留的话有33个）
sig <- intersect(sig, intersect(TCGA.lnc.anno$gene_id, probe.lnc.anno$ENSEMBL))

## Training Cohort ---------------------------------------------------------

# # 使用TCGAbiolinks下载TCGA队列转录组数据
# project <- paste0("TCGA-", c("COAD", "READ"))
# query <- GDCquery(
#   project = project,
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification",
#   workflow.type = "STAR - Counts"
# )
# # GDCdownload(query, method = "api",
# #             directory = file.path(data.path, "GDCdata"))
# # TCGA.rnaseq <- GDCprepare(query, directory = file.path(data.path, "GDCdata"))
# # saveRDS(TCGA.rnaseq, file.path(data.path, "TCGA.rnaseq.rds"))
#
# TCGA.rnaseq <- readRDS(file.path(data.path, "TCGA.rnaseq.rds"))
# TCGA.expr = assay(TCGA.rnaseq, i = "tpm_unstrand")
# rownames(TCGA.expr) <- do.call(rbind, strsplit(rownames(TCGA.expr), "\\."))[, 1]
# TCGA.expr = TCGA.expr[, TCGA.rnaseq$is_ffpe == F]
# TCGA.expr = TCGA.expr[, substr(colnames(TCGA.expr), 14, 16) == "01A"]
#
# # 对TCGA样本按患者(patient)、样本(sample)和plate进行排序
# # 仅选取最新plate（编号最大）的01A样本
# TCGA.sample <- data.frame("name" = colnames(TCGA.expr))
# TCGA.sample$patient <- substr(colnames(TCGA.expr), 1, 12)
# TCGA.sample$sample <- substr(colnames(TCGA.expr), 14, 16)
# TCGA.sample$plate <- substr(colnames(TCGA.expr), 22, 25)
# TCGA.sample <- arrange(TCGA.sample, TCGA.sample$patient, TCGA.sample$sample, desc(TCGA.sample$plate))
# TCGA.sample <- TCGA.sample[!duplicated(TCGA.sample$patient), ]
# TCGA.sample <- TCGA.sample[TCGA.sample$sample == "01A", ]
# # View(TCGA.sample[TCGA.sample$patient %in% a, ])
# TCGA.expr <- TCGA.expr[, colnames(TCGA.expr) %in% TCGA.sample$name]
#
#
#
# # 生存数据
# # https://xenabrowser.net/datapages/
# # TCGA Colon and Rectal Cancer (COADREAD)
# # Curated survival data
# TCGA.clin <- read.table(file.path(data.path, "TCGA-CRC_Clinic.txt"), header = T, sep = "\t")
# TCGA.clin$sample <- NULL
# TCGA.clin <- TCGA.clin[!duplicated(TCGA.clin$X_PATIENT), ]
# TCGA.clin <- TCGA.clin[!is.na(TCGA.clin$OS.time), ]
# sum(substr(colnames(TCGA.expr), 1, 12) %in% TCGA.clin$X_PATIENT)

# # 对表达数据和生存数据取交集
# TCGA.expr <- TCGA.expr[, substr(colnames(TCGA.expr), 1, 12) %in% TCGA.clin$X_PATIENT]
# TCGA.clin <- TCGA.clin[match(substr(colnames(TCGA.expr), 1, 12), TCGA.clin$X_PATIENT), ]
# rownames(TCGA.clin) <- TCGA.clin$X_PATIENT
# write.table(TCGA.clin, file.path(res.path, "TCGA.clin.txt"),
#             sep = "\t", row.names = F, col.names = T, quote = F)
#
# # 基于logTPM对COAD和READ的数据去批次
# # 结果存放在Results/combat.TCGA.expr.rds
# batch <- setNames(object = TCGA.rnaseq$project_id, nm = TCGA.rnaseq$barcode)
# batch <- batch[colnames(TCGA.expr)]
# table(batch)
# TCGA.expr <- log(TCGA.expr+1)
# combat.TCGA.expr <- ComBat(dat = TCGA.expr, batch = batch)
# saveRDS(combat.TCGA.expr, file.path(res.path, "combat.TCGA.expr.rds"))

# 制作训练集
combat.TCGA.expr <- readRDS(file.path(res.path, "combat.TCGA.expr.rds"))
TCGA.clin <- read.table(file.path(res.path, "TCGA.clin.txt"), header = T, sep = "\t")
Train_clin <- TCGA.clin[TCGA.clin$OS.time>0, c("OS", "OS.time")]
Train_expr <- t(combat.TCGA.expr[sig, TCGA.clin$OS.time>0])

## Validation Cohort -------------------------------------------------------

# Sys.setlocale( 'LC_ALL','C' )
# meta_expr = meta_clin = batch = list()
# for (project in c("GSE17536", "GSE17537", "GSE29621", "GSE38832", "GSE39582", "GSE72970")){
#   # project = "GSE72970"
#   cat(project, "\n")
#   filename <- file.path(data.path, paste0(project, "_series_matrix.txt.gz"))
#   geo <- getGEO(filename = filename)
#   clin <- as.data.frame(geo@phenoData@data)
#
#   meta_expr[[project]] <- geo@assayData$exprs
#   meta_clin[[project]] <- clin
#   batch[[project]] <- rep(project, ncol(geo))
# }
#
# batch <- unlist(batch)
# meta_expr <- do.call(cbind, meta_expr)
# combat.meta_expr <- ComBat(dat = meta_expr, batch = batch)
# saveRDS(combat.meta_expr, file.path(res.path, "combat.meta_expr.rds"))
#
# # 手动提取出OS和OS.time
#
# ## GSE17536
# {
#   clin <- meta_clin$GSE17536
#   colnames(clin)
#   clin <- clin[, c("overall survival follow-up time:ch1", "overall_event (death from any cause):ch1")]
#   colnames(clin) <- c("OS.time", "OS")
#   table(clin$OS)
#   clin$OS <- ifelse(test = clin$OS == "death", yes = 1, no = 0)
#   clin$OS.time <- as.numeric(clin$OS.time)
#   meta_clin$GSE17536 <- clin
# }
#
# ## GSE17537
# {
#   clin <- meta_clin$GSE17537
#   colnames(clin)
#   clin <- clin[, c("overall survival follow-up time:ch1", "overall_event (death from any cause):ch1")]
#   colnames(clin) <- c("OS.time", "OS")
#   table(clin$OS)
#   clin$OS <- ifelse(test = clin$OS == "death", yes = 1, no = 0)
#   clin$OS.time <- as.numeric(clin$OS.time)
#   meta_clin$GSE17537 <- clin
# }
#
# ## GSE29621
# {
#   clin <- meta_clin$GSE29621
#   colnames(clin)
#   clin <- clin[, c("overall survival (os):ch1", "os event:ch1")]
#   colnames(clin) <- c("OS.time", "OS")
#   table(clin$OS)
#   clin$OS <- ifelse(test = clin$OS == "dead", yes = 1, no = 0)
#   clin$OS.time <- as.numeric(clin$OS.time)
#   meta_clin$GSE29621 <- clin
# }
#
# ## GSE38832 无OS，使用DSS代替
# {
#   clin <- meta_clin$GSE38832
#   colnames(clin)
#   clin <- clin[, c("dss_time (disease specific survival time, months):ch1", "dss_event (disease specific survival):ch1")]
#   colnames(clin) <- c("OS.time", "OS")
#   table(clin$OS)
#   clin$OS <- ifelse(test = clin$OS == "1 (death from cancer)", yes = 1, no = 0)
#   clin$OS.time <- as.numeric(clin$OS.time)
#   meta_clin$GSE38832 <- clin
# }
#
# ## GSE39582
# {
#   clin <- meta_clin$GSE39582
#   colnames(clin)
#   clin <- clin[, c("os.delay (months):ch1", "os.event:ch1")]
#   colnames(clin) <- c("OS.time", "OS")
#   table(clin$OS)
#   clin <- clin[clin$OS != "N/A", ]
#   clin$OS <- ifelse(test = clin$OS == "1", yes = 1, no = 0)
#   clin$OS.time <- as.numeric(clin$OS.time)
#   meta_clin$GSE39582 <- clin
# }
#
# ## GSE72970
# {
#   clin <- meta_clin$GSE72970
#   colnames(clin)
#   clin <- clin[, c("os:ch1", "os censored:ch1")]
#   colnames(clin) <- c("OS.time", "OS")
#   table(clin$OS)
#   clin$OS <- ifelse(test = clin$OS == "1", yes = 1, no = 0)
#   clin$OS.time <- as.numeric(clin$OS.time)
#   meta_clin$GSE72970 <- clin
# }
#
# meta_clin <- do.call(rbind, meta_clin)
# meta_clin$GEO <- do.call(rbind, strsplit(rownames(meta_clin), "\\."))[, 1]
# meta_clin$Sample <- do.call(rbind, strsplit(rownames(meta_clin), "\\."))[, 2]
# meta_clin <- meta_clin[, c("Sample", "GEO", "OS", "OS.time")]
# write.table(meta_clin, file.path(res.path, "meta_clin.txt"),
#             sep = "\t", row.names = F, col.names = T, quote = F)

## 制作测试集
combat.meta_expr <- readRDS(file.path(res.path, "combat.meta_expr.rds"))
meta_clin <- read.table(file.path(res.path, "meta_clin.txt"), header = T, sep = "\t")
Test_expr <- t(combat.meta_expr[probe.lnc.anno$probe[match(sig, probe.lnc.anno$ENSEMBL)],
                                meta_clin$Sample])
colnames(Test_expr) <- sig
Test_clin <- meta_clin


# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
# 此处记录需要运行的模型，格式为：算法1名称[算法参数]+算法2名称[算法参数]
# 目前仅有StepCox和RunEnet支持输入算法参数
methods <- read.xlsx(file.path(data.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)
head(methods)

## Train the model --------------------------------------------------------
model <- list()
set.seed(seed = 777)
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method # 本轮算法名称
  method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称

  Variable = colnames(Train_expr) # 最后用于构建模型的变量
  for (i in 1:length(method)){
    if (i < length(method)){
      selected.var <- RunML(method = method[i], # 机器学习方法
                            Train_expr = Train_expr, # 训练集有潜在预测价值的变量
                            Train_clin = Train_clin, # 训练集生存数据
                            mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                            timeVar = "OS.time", statusVar = "OS") # 用于训练的生存变量，应当出现在Train_clin中
      if (length(selected.var)>5) Variable <- intersect(Variable, selected.var)
    }else{
      model[[method_name]] <- RunML(method = method[i],
                                    Train_expr = Train_expr[, Variable],
                                    Train_clin = Train_clin,
                                    mode = "Model",
                                    timeVar = "OS.time", statusVar = "OS")
    }
  }
}
saveRDS(model, file.path(res.path, "model.rds"))

## Evaluate the model -----------------------------------------------------

model <- readRDS(file.path(res.path, "model.rds"))
# 对各模型计算C-index
Cindexlist <- list()
for (method in methods){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], # 预后模型
                                  Test_expr = Test_expr, # 测试集预后变量，应当包含训练集中所有的变量，否则会报错
                                  Test_clin = Test_clin) # 训练集生存数据，应当包含训练集中所有的变量，否则会报错
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "Cindex_mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

Cindex_mat <- read.table(file.path(res.path, "Cindex_mat.txt"),
                         sep = "\t", row.names = 1)
avg_Cindex <- apply(Cindex_mat, 1, mean)           # 计算每种算法在所有队列中平均C-index
avg_Cindex <- sort(avg_Cindex, decreasing = T)     # 对各算法C-index由高到低排序
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]      # 对C-index矩阵排序

avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3))
row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8,
                                          gp = gpar(fill = "#1691CD", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          width = unit(2, "cm")))

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Set1")
names(CohortCol) <- colnames(Cindex_mat)
col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = CohortCol))

cellwidth = 1; cellheight = 0.5
pdf(file.path(fig.path, "Cindex.pdf"),
    width = cellwidth * ncol(Cindex_mat) + 3,
    height = cellheight * nrow(Cindex_mat) * 0.45)
Heatmap(as.matrix(Cindex_mat), name = "C-index",
        right_annotation = row_ha, top_annotation = col_ha,
        col = c("#1CB8B2", "#FFFFFF", "#EEB849"),
        cluster_columns = F, cluster_rows = F,
        rect_gp = gpar(col = "#CCCCCA", lwd = 1),
        show_column_names = F, show_row_names = T,
        width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
        height = unit(cellheight * nrow(Cindex_mat), "cm"),
        column_split = colnames(Cindex_mat), column_title=NULL,
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                    x, y, gp = gpar(fontsize = 10))
        })
dev.off()

fit = model[["survivalSVM"]]
RS <- switch(
  EXPR = class(fit)[1],
  "coxnet"      = coef(fit)[, 1],
  "coxph"       = coef(fit),
  "survivalsvm" = setNames(object = as.vector(t(fit$model.fit$SV) %*% fit$model.fit$Beta), 
                           nm = fit$var.names),
  "CoxBoost"    = coef(fit),
  "superpc"     = fit$feature.scores,
  "plsRcoxmodel" = setNames(object = fit$Coeffs[, 1],
                            nm = rownames(fit$Coeffs)),
  "rfsrc"        = setNames(object = var.select(fit, verbose = F)$varselect$vimp,
                            nm = rownames(var.select(fit, verbose = F)$varselect)),
  "gbm"          = setNames(object = summary.gbm(fit)$rel.inf,
                            nm = summary.gbm(fit)$var)
)

Predict.out <- Test_surv
Predict.out$RS <- as.vector(RS)
Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
f <- as.formula(paste0("Surv(", timeVar,",",statusVar,")~RS"))
unlist(lapply(Predict.out, function(data){
  unname(summary(coxph(formula = f,
                       data = data))$concordance["C"])
}))