#!/usr/bin/env Rscript
# R version 4.0.0

## Classifier from a Polygenic risk score (PRS) to predict prostate cancer.##

libraries <- c("caret", "MLeval", "mlbench", 
               "ROCR", "pROC", "unix")

for (lib in libraries) {
  if (require(package = lib, character.only = TRUE)) {
    successful <- "Successful"
  } else {
    installing <- "Installing"
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkgs = lib, suppressUpdates = T)
    library(lib, character.only = TRUE )
  }
}

setwd("/healthyvscancer/lassosum")
set.seed(2023)

#### Healthy vs Cancer #####
lass           <- read.table("lass_selection.txt", header = T)
colnames(lass) <- c("FID", "IID", "pheno", "PRS", "PC1", "PC2", "Condition")
lass$IID       <- gsub(pattern = "(.+\\/)|(_.+)", x = lass$IID, replacement = "")
lass$pheno     <- as.factor(lass$pheno)
lass$condition <- ifelse(lass$pheno == 2, "Cancer", "Healthy")
lass$OR        <- exp(lass$PRS)
PRS_data       <- lass


# Split the data into Training: 80%, and Validation: 20%.
PRS_data$pheno <- as.numeric(PRS_data$pheno)
PRS_data$pheno <- ifelse(PRS_data$pheno == 2, 1, 0)
train_knn      <- createDataPartition(PRS_data[,"pheno"], times = 1, p = 0.8,list = F) 
data.trn       <- PRS_data[train_knn,] 
data.trn$pheno <- as.factor(data.trn$pheno)
data.tst       <- PRS_data[-train_knn,] 
data.tst$pheno <- as.factor(data.tst$pheno)

# 10 fold Cross-validation step with multi class summary to improve the model"s accuracy.
ctrl  <- trainControl(method          = "cv", 
                      number          = 10, 
                      savePredictions = T, 
                      summaryFunction = twoClassSummary, 
                      classProbs      = T) 

# Set all methods available in the caret package.
mlmethods <- c("glm", "lda", "qda", 
               "knn", "rpart", "rf", 
               "gbm", "nnet", "svmLinear", 
               "svmRadial", "xgbLinear")


# Train the model with the selected methods.
# Square root of the number of observations (nrow(PRS_data)) determined the number of Ks to use.
fit_model <- list()
pred      <- list()

for(i in 1:length(mlmethods)){
  tryCatch({
    fit_model[[i]] <- train(as.factor(Condition) ~ PRS, 
                            data       = data.trn, 
                            method     = mlmethods[i], 
                            trControl  = ctrl, 
                            preProcess = c("center","scale"), 
                            tuneLength = 20, 
                            metric     = "Accuracy")
    pred[[i]]      <- confusionMatrix(predict(fit_model[[i]], 
                                              newdata = data.tst), 
                                      as.factor(data.tst$Condition)
                                      )
    print(paste0("Prediction done for model ", mlmethods[i]))
  }, error = function(e) {
    pred[[i]] <- e$message
    print(paste0("Error for model ", mlmethods[i]))
  })
  gc()
}


# Compared the SP and SE between the confusion matrix the roc curves, for the models that worked).
comparisons <- data.frame(Models  = c("glm", "lda", "qda", "knn", "rpart", "rf", "gbm", "nnet"),
                          se_test = rep(NA, length(mlmethods)), 
                          se_roc  = rep(NA, length(mlmethods)), 
                          sp_test = rep(NA, length(mlmethods)), 
                          sp_roc  = rep(NA, length(mlmethods)))


for(id in 1:length(mlmethods)){
  comparisons$se_test[id] <- pred[[id]]$byClass[1]
  comparisons$sp_test[id] <- pred[[id]]$byClass[2]
  comparisons$se_roc[id] <- roc_lda[[id]]$sensitivities[2]
  comparisons$sp_roc[id] <- roc_lda[[id]]$specificities[2]
}

# Calculate the Youden index for the validation data.
optimal.cutpoint.Youdenprs <- optimal.cutpoints(X           = "PRS", 
                                                status      = "pheno", 
                                                tag.healthy = 0, 
                                                methods     = "Youden", 
                                                data        = data.tst, 
                                                pop.prev    = NULL, 
                                                control     = control.cutpoints(), 
                                                ci.fit.     = FALSE, 
                                                conf.level  = 0.95, 
                                                trace       = FALSE)


# Plot the roc curves with several thresholds using probabilities.
roc_list  <- list()
pred_list <- list()
perf_list <- list()
auc_list  <- list()

for(id in 1:length(mlmethods)){
  roc_list[[id]]  <- predict(fit_model[[id]], newdata = data.tst, type = "prob")
  pred_list[[id]] <- prediction(roc_list[[id]][,1], 
                                data.tst$Condition, 
                                label.ordering = c("Healthy", "Cancer")
                                )
  perf_list[[id]] <- performance(pred_list[[id]], measure = "tpr", x.measure = "fpr")
  auc_list[[id]]  <- performance(pred_list[[id]], measure = "auc")@y.values
}

plot(perf_list[[1]], 
     main = "Healthy vs cancer", 
     col  = "blue", 
     lwd  = 2)
plot(perf_list[[2]], 
     col = "red", 
     lwd = 2, 
     add = T)
plot(perf_list[[3]], 
     col = "black", 
     lwd = 2, 
     add = T)
plot(perf_list[[4]], 
     col = "magenta", 
     lwd = 2, 
     add = T)
plot(perf_list[[5]], 
     col = "green", 
     lwd = 2, 
     add = T)
plot(perf_list[[6]], 
     col = "brown", 
     lwd = 2, 
     add = T)
abline(a   = 0, 
       b   = 1, 
       lty = 2, 
       lwd = 2, 
       col = "gray")
legend("bottomright", legend = paste(mlmethods,"-", 
                                     paste0("AUC = ", 
                                            round(unlist(auc_list), 3))
                                     ), 
       col = c("blue", "red", "black", "magenta", "green", "brown"), 
       lty = 1)

# Plot The best model.
plot(perf_list[[5]], 
     main         = "ROC Curve from rpart", 
     col          = "green", 
     lwd          = 2,
     print.auc    = T, 
     auc.polygon  = T, 
     ci           = T)
abline(a   = 0, 
       b   = 1, 
       lty = 2, 
       lwd = 2, 
       col = "gray")
legend("bottomright", paste("AUC = ", 
                            round(as.numeric(auc_list[[5]]), 3)), 
       bty = "n")


# Compare the model's preditions.
compare_models <- resamples(pred)
dotplot(compare_models) 

###### NEURAL NETWORK ######
#Neural network was run with different hidden layers.
n <- neuralnet(pheno ~ PRS, 
               data          = data.trn, 
               hidden        = 1,
               algorithm     = "rprop+",
               rep           = 5,
               err.fct       = "ce", 
               linear.output = F,
               stepmax       = 150000)
# Select the best rep value based on the lowest error rate.
for (id in 1:10) {
  plot(n, rep = id)
}


# Compute the predicted values for Healthy (0) and cancer (1) on the validation data.
out.    <- predict(n, data.tst, rep = 5)
out[,1] <- ifelse(out[,1] > 0.5, 0, 1)
tab     <- table(out[,1], as.factor(data.tst$pheno))
cm_nn   <- confusionMatrix(tab)

