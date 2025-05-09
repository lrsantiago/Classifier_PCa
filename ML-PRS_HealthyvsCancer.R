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

set.seed(119)

# Split the data into Training: 70%, and Validation: 30%.
train_knn      <- createDataPartition(PRS_data[,"pheno"], times = 1, p = 0.7,list = F)
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
               "svmRadial")


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


# Predict aggressive cancer on the test data and calculate the AUC for each model.
 pred_roc <- list()
 roc_lda <- list()
 for(id in 1:length(mlmethods)){
   pred_roc[[id]] <- predict(fit_model[[id]], newdata = data.tst, type = "raw")
   pred_roc[[id]] <- ifelse(pred_roc[[id]] == "Healthy", 0, 1)
   roc_lda[[id]] <- roc(as.factor(data.tst$Condition),
                        pred_roc[[id]],
                        levels = c("Healthy", "Cancer"),
                        auc = T,
                        ci = T,
                        direction = "<")
}


 # Get the sensitivity, specificity, accuracy, kappa, and AUC for the models.
comparisons <- data.frame(Models    = mlmethods,
                           se       = rep(NA, length(mlmethods)),
                           sp       = rep(NA, length(mlmethods)),
                           accuracy = rep(NA, length(mlmethods)),
                           kappa.   = rep(NA, length(mlmethods)),
                           auc      = rep(NA, length(mlmethods)))

for(id in 1:length(mlmethods)){
  comparisons$se[id]       <- pred[[id]]$byClass[1]
  comparisons$sp[id]       <- pred[[id]]$byClass[2]
  comparisons$accuracy[id] <- pred[[id]]$overall[1]
  comparisons$auc[id]      <- roc_lda[[1]]$auc
  comparisons$kappa[id].   <- pred[[id]]$overall[2]
}
