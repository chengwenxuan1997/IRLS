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
# ??????TCGA????????????????????????CRC??????????????????Data Release 32~25????????????????????????????????????Dr 36
# ??????????????????Dr36???gtf?????????https://www.gencodegenes.org/human/releases.html????????????
# gtf???????????????Annotations????????????
# ??????????????????
# ??????gtf --> ??????ensembl??????????????? --> ?????????????????? -->?????????txt
# ???????????????Annotations/TCGA.anno.txt???
# TCGA.anno <- import(file.path(ann.path, "gencode.v36.annotation.gtf"))
# TCGA.anno <- as.data.frame(TCGA.anno)
# TCGA.anno <- TCGA.anno[!duplicated(TCGA.anno$gene_id),
#                        c("gene_id", "gene_type", "gene_name")]
# TCGA.anno$gene_id <- do.call(rbind, strsplit(TCGA.anno$gene_id, "\\."))[, 1]
# write.table(TCGA.anno, file.path(ann.path, "TCGA.anno.txt"),
#             sep = "\t", col.names = T, row.names = F, quote = F)
#

# ?????????????????????????????????????????????gencode.v36.annotation.gtf?????????gencode.v36.transcripts.fa??????
# fa??????????????????????????????????????????gtf????????????
# ??????fasta????????????GEO???????????????????????????????????????
# ?????????fasta?????????https://www.thermofisher.cn/cn/zh/home.html???????????????Annotations????????????
# ??????????????????
# system(
#   paste(file.path(pkg.path, "seqmap-1.0.12-windows.exe"), 0,  # seqmap????????????
#         file.path(ann.path, "HG-U133_Plus_2.probe_fasta"),    # ??????fasta????????????
#         file.path(ann.path, "gencode.v36.transcripts.fa"),    # gencode????????????
#         file.path(ann.path, "probe.reannotation"),            # ???????????????
#         "/output_all_matches")
# )

# ???seqmap?????????????????????
# ?????????????????????????????????(probe_id????????????)?????????????????????(trans_id)
# ???????????????Annotations/probe.anno.txt???
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

# ?????????TCGA???????????????lncRNA?????????
# ???????????????Annotatins/TCGA.lnc.anno.txt
# TCGA.anno <- read.table(file.path(ann.path, "TCGA.anno.txt"), header = T)
# table(TCGA.anno$gene_type)[c("protein_coding", "lncRNA")]
# TCGA.lnc.anno <- TCGA.anno[TCGA.anno$gene_type == "lncRNA", ]
# write.table(TCGA.lnc.anno, file.path(ann.path, "TCGA.lnc.anno.txt"),
#             sep = "\t", col.names = T, row.names = F, quote = F)

# ?????????????????????????????????lncRNA?????????
# ????????????????????????????????????????????????????????????lncRNA?????????
# ???????????????probe.lnc.anno.txt
# probe.anno <- read.table(file.path(ann.path, "probe.anno.txt"), header = T)
# table(probe.anno$BioType)[c("protein_coding", "lncRNA")]
# protein_probe <- probe.anno$probe[probe.anno$BioType %in% protein_pseudo] #?????????????????????????????????
# probe.anno <- probe.anno[!probe.anno$probe %in% protein_probe, ]
# probe.anno <- probe.anno[probe.anno$BioType == "lncRNA", ]
# duplicate_probe <- probe.anno$probe[duplicated(probe.anno$probe)] # ??????????????????lncRNA?????????
# probe.lnc.anno <- probe.anno[!probe.anno$probe %in% duplicate_probe, ]
# write.table(probe.lnc.anno, file.path(ann.path, "probe.lnc.anno.txt"),
#             sep = "\t", col.names = T, row.names = F, quote = F)


## Signature ---------------------------------------------------------------

# ??????????????????43???lncRNA??????
sig <- read.table(file.path(res.path, "signature.txt"))
sig

TCGA.lnc.anno <- read.table(file.path(ann.path, "TCGA.lnc.anno.txt"), header = T)
probe.lnc.anno <- read.table(file.path(ann.path, "probe.lnc.anno.txt"), header = T)
sig <- sig$V2
# ??????????????????????????????lncRNA??????????????????18?????????????????????????????????33??????
sig <- intersect(sig, intersect(TCGA.lnc.anno$gene_id, probe.lnc.anno$ENSEMBL))

## Training Cohort ---------------------------------------------------------

# # ??????TCGAbiolinks??????TCGA?????????????????????
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
# # ???TCGA???????????????(patient)?????????(sample)???plate????????????
# # ???????????????plate?????????????????????01A??????
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
# # ????????????
# # https://xenabrowser.net/datapages/
# # TCGA Colon and Rectal Cancer (COADREAD)
# # Curated survival data
# TCGA.clin <- read.table(file.path(data.path, "TCGA-CRC_Clinic.txt"), header = T, sep = "\t")
# TCGA.clin$sample <- NULL
# TCGA.clin <- TCGA.clin[!duplicated(TCGA.clin$X_PATIENT), ]
# TCGA.clin <- TCGA.clin[!is.na(TCGA.clin$OS.time), ]
# sum(substr(colnames(TCGA.expr), 1, 12) %in% TCGA.clin$X_PATIENT)

# # ???????????????????????????????????????
# TCGA.expr <- TCGA.expr[, substr(colnames(TCGA.expr), 1, 12) %in% TCGA.clin$X_PATIENT]
# TCGA.clin <- TCGA.clin[match(substr(colnames(TCGA.expr), 1, 12), TCGA.clin$X_PATIENT), ]
# rownames(TCGA.clin) <- TCGA.clin$X_PATIENT
# write.table(TCGA.clin, file.path(res.path, "TCGA.clin.txt"),
#             sep = "\t", row.names = F, col.names = T, quote = F)
#
# # ??????logTPM???COAD???READ??????????????????
# # ???????????????Results/combat.TCGA.expr.rds
# batch <- setNames(object = TCGA.rnaseq$project_id, nm = TCGA.rnaseq$barcode)
# batch <- batch[colnames(TCGA.expr)]
# table(batch)
# TCGA.expr <- log(TCGA.expr+1)
# combat.TCGA.expr <- ComBat(dat = TCGA.expr, batch = batch)
# saveRDS(combat.TCGA.expr, file.path(res.path, "combat.TCGA.expr.rds"))

# ???????????????
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
# # ???????????????OS???OS.time
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
# ## GSE38832 ???OS?????????DSS??????
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

## ???????????????
combat.meta_expr <- readRDS(file.path(res.path, "combat.meta_expr.rds"))
meta_clin <- read.table(file.path(res.path, "meta_clin.txt"), header = T, sep = "\t")
Test_expr <- t(combat.meta_expr[probe.lnc.anno$probe[match(sig, probe.lnc.anno$ENSEMBL)],
                                meta_clin$Sample])
colnames(Test_expr) <- sig
Test_clin <- meta_clin


# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
# ??????????????????????????????????????????????????????1??????[????????????]+??????2??????[????????????]
# ????????????StepCox???RunEnet????????????????????????
methods <- read.xlsx(file.path(data.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)
head(methods)

## Train the model --------------------------------------------------------
model <- list()
set.seed(seed = 777)
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method # ??????????????????
  method <- strsplit(method, "\\+")[[1]] # ?????????????????????

  Variable = colnames(Train_expr) # ?????????????????????????????????
  for (i in 1:length(method)){
    if (i < length(method)){
      selected.var <- RunML(method = method[i], # ??????????????????
                            Train_expr = Train_expr, # ???????????????????????????????????????
                            Train_clin = Train_clin, # ?????????????????????
                            mode = "Variable",       # ???????????????Variable(????????????)???Model(????????????)
                            timeVar = "OS.time", statusVar = "OS") # ?????????????????????????????????????????????Train_clin???
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
# ??????????????????C-index
Cindexlist <- list()
for (method in methods){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], # ????????????
                                  Test_expr = Test_expr, # ?????????????????????????????????????????????????????????????????????????????????
                                  Test_clin = Test_clin) # ?????????????????????????????????????????????????????????????????????????????????
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "Cindex_mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

Cindex_mat <- read.table(file.path(res.path, "Cindex_mat.txt"),
                         sep = "\t", row.names = 1)
avg_Cindex <- apply(Cindex_mat, 1, mean)           # ??????????????????????????????????????????C-index
avg_Cindex <- sort(avg_Cindex, decreasing = T)     # ????????????C-index??????????????????
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]      # ???C-index????????????

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