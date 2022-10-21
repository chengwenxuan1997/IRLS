# 调用各机器学习算法的桥接函数
RunML <- function(method, Train_expr, Train_clin, mode = "Model", timeVar = "OS.time", statusVar = "OS", ...){
  # Variable模式下返回筛选后的变量，Model模式下返回预测模型
  # 以"Enet [alpha=0.4]"为例
  method = gsub(" ", "", method) # 去除参数中的空格，得到Enet [alpha=0.4]
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)  #获取算法名称，即Enet
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method) #获取算法参数，即alpha=0.4

  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "StepCox" = list("direction" = method_param),
    NULL
  )
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_expr), " Variables")

  args = list("Train_expr" = Train_expr,
              "Train_clin" = Train_clin,
              "mode" = mode,
              "timeVar" = timeVar, "statusVar" = statusVar)
  args = c(args, method_param)

  obj <- do.call(what = paste0("Run", method_name),
                 args = args)   #调用机器学习算法

  if(mode == "Variable"){
    message(length(obj), " Variables retained;\n")
  }else{message("\n")}
  return(obj)
}

RunEnet <- function(Train_expr, Train_clin, mode, timeVar, statusVar, alpha){
  cv.fit = cv.glmnet(x = Train_expr,
                     y = Surv(Train_clin[[timeVar]], Train_clin[[statusVar]]),
                     family = "cox", alpha = alpha, nfolds = 10)
  fit = glmnet(x = Train_expr,
               y = Surv(Train_clin[[timeVar]], Train_clin[[statusVar]]),
               family = "cox", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  co = coef(fit)[, 1]
  if (mode == "Variable") return(names(co)[co!=0])
}

RunLasso <- function(Train_expr, Train_clin, mode, timeVar, statusVar){
  RunEnet(Train_expr, Train_clin, mode, timeVar, statusVar, alpha = 1)
}

RunRidge <- function(Train_expr, Train_clin, mode, timeVar, statusVar){
  RunEnet(Train_expr, Train_clin, mode, timeVar, statusVar, alpha = 0)
}

RunStepCox <- function(Train_expr, Train_clin, mode, timeVar, statusVar, direction){
  fit <- step(coxph(formula = Surv(Train_clin[[timeVar]], Train_clin[[statusVar]]) ~ .,
                    data = as.data.frame(Train_expr)),
              direction = direction, trace = 0)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(names(coef(fit)))
}

RunsurvivalSVM <- function(Train_expr, Train_clin, mode, timeVar, statusVar){
  fit = survivalsvm(formula = Surv(Train_clin[[timeVar]], Train_clin[[statusVar]]) ~ .,
                    data= as.data.frame(Train_expr),
                    gamma.mu = 1)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(fit["var.names"])
}

RunCoxBoost <- function(Train_expr, Train_clin, mode, timeVar, statusVar){
  pen <- optimCoxBoostPenalty(time = Train_clin[[timeVar]],
                              status = Train_clin[[statusVar]],
                              x = Train_expr,
                              trace = F, start.penalty=500, parallel = F)
  cv.res <- cv.CoxBoost(time = Train_clin[[timeVar]],
                        status = Train_clin[[statusVar]],
                        x = Train_expr,
                        maxstepno=500, K=10, type="verweij", penalty=pen$penalty)
  fit <- CoxBoost(time = Train_clin[[timeVar]],
                  status = Train_clin[[statusVar]],
                  x = Train_expr,
                  stepno = cv.res$optimal.step,
                  penalty = pen$penalty)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(names(coef(fit)[abs(coef(fit))>0]))
}

RunSuperPC <- function(Train_expr, Train_clin, mode, timeVar, statusVar){
  data <- list(x = t(Train_expr),
               y = Train_clin[[timeVar]],
               censoring.status = Train_clin[[statusVar]],
               featurenames = colnames(Train_expr))
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  cv.fit <- suppressWarnings(superpc.cv(fit, data,
                                        n.threshold = 20,#default
                                        n.fold = 10,
                                        n.components = 3,
                                        min.features = 5,
                                        max.features = nrow(data$x),
                                        compute.fullcv = TRUE,
                                        compute.preval = TRUE))
  fit$threshold <- cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])]
  fit$data <- data
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(names(fit$feature.scores)[abs(fit$feature.scores)>0.5])
}

RunplsRcox <- function(Train_expr, Train_clin, mode, timeVar, statusVar){
  data <- list(x = Train_expr,
               time = Train_clin[[timeVar]],
               status = Train_clin[[statusVar]])
  cv.plsRcox.res = cv.plsRcox(data = data,
                              nt=10, verbose = FALSE)
  fit <- plsRcox(Xplan = data$x,
                 time = data$time,
                 event = data$status,
                 nt = as.numeric(cv.plsRcox.res[5]),
                 verbose = F)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(rownames(fit$Coeffs))
}

RunRSF <- function(Train_expr, Train_clin, mode, timeVar, statusVar){
  rf_nodesize = 2
  fit <- rfsrc(formula = formula(paste0("Surv(", timeVar, ", ", statusVar, ")", "~.")),
               data = cbind(Train_expr, Train_clin),
               ntree = 1000, nodesize = rf_nodesize,##该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(var.select(fit, verbose = F)$topvars)
}

RunGBM <- function(Train_expr, Train_clin, mode, timeVar, statusVar){
  fit <- gbm(formula = Surv(Train_clin[[timeVar]], Train_clin[[statusVar]]) ~ .,
             data = as.data.frame(Train_expr),
             distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  best <- which.min(fit$cv.error)
  fit <- gbm(formula = Surv(Train_clin[[timeVar]], Train_clin[[statusVar]]) ~ .,
             data = as.data.frame(Train_expr),
             distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001, n.cores = 8)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(rownames(summary.gbm(fit))[summary.gbm(fit)$rel.inf>0])
}

RunEval <- function(fit, Test_expr, Test_clin){
  new_data <- Test_expr[, fit$subFeature]
  RS <- switch(
    EXPR = class(fit)[1],
    "coxnet"   = predict(fit, as.matrix(new_data)),
    "coxph"       = predict(fit, as.data.frame(new_data)),
    "survivalsvm" = predict(fit, as.data.frame(new_data))$predicted,
    "CoxBoost"    = predict(fit, as.data.frame(new_data)),
    "superpc"     = superpc.predict(object = fit,
                                    data = fit$data,
                                    newdata = list(x = t(as.matrix(new_data))),
                                    threshold = fit$threshold,
                                    n.components = 1)$v.pred,
    "plsRcoxmodel" = predict(fit, as.data.frame(new_data)),
    "rfsrc"       = predict(fit, as.data.frame(new_data))$predicted,
    "gbm"         = predict(fit, as.data.frame(new_data))
  )

  Predict.out <- Test_clin
  Predict.out$RS <- as.vector(RS)
  Predict.out <- split(x = Predict.out, f = Predict.out$GEO)
  unlist(lapply(Predict.out, function(data){
    unname(summary(coxph(formula = Surv(OS.time,OS)~RS,
                         data = data))$concordance["C"])
  }))
}
