#!/usr/bin/env Rscript
# R version 4.0.0
## Classifier using PSA and age to predict prostate cancer aggressiveness.

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

# Prepare the data.
data_all               <- read.csv("TAPG_TURP_DATA.csv")
data_all$pheno         <- as.factor(ifelse(data_all$pheno == 1, 0, 1))
data_all               <- subset(data_all, !is.na(CAPRA))
surv_data              <- data_all[, c("patientid", "Status_2011", "Age", 
                                       "survival2011", "PSA", "LABNO", 
                                       "Gleason", "pheno", "CAPRA")]
surv_data$Age          <- as.integer(surv_data$Age)
surv_data$survival2011 <- as.integer(surv_data$survival2011)

# Split the data into Training: 70%, and Validation: 30%.
surv_data$pheno <- as.numeric(surv_data$pheno)
surv_data$pheno <- ifelse(surv_data$pheno == 2, 1, 0)

set.seed(110)
train            <- createDataPartition(surv_data$pheno, times = 1, p = 0.7, list = F)
train.data       <- surv_data[train,]
train.data$pheno <- as.factor(train.data$pheno)
val.data         <- surv_data[-train,]
val.data$pheno   <- as.factor(val.data$pheno)

# 10 fold cross-validation with multi class summary to improve the model's accuracy.
ctrl  <- trainControl(method          = "cv", 
                      number          = 10, 
                      savePredictions = T, 
                      summaryFunction = twoClassSummary, 
                      classProbs      = T) 

# Select six methods available in the caret package
mlmethods <- c("glm", "lda", "qda", "knn", "rpart", "rf")

# Train the model with the selected methods.
# Based on the  square root of the number of observations (nrow(data_all)), ks was set to 50.
fit_model <- list()
pred      <- list()


for(i in 1:length(mlmethods)){
  tryCatch({
    fit_model[[i]] <- train(as.factor(Status_2011) ~ PSA + age, 
                            data       = train.data, 
                            method.    = mlmethods[i], 
                            trCont.rol = ctrl, 
                            preProcess = c("center","scale"), 
                            tuneLength = 50, 
                            metric.    = "Spec")
    pred[[i]]      <- confusionMatrix(predict(fit_model[[i]], 
                                              newdata = val.data), 
                                      as.factor(val.data$Status_2011)
                                      )
    print(paste0("Prediction done for model", mlmethods[i]))
  }, error = function(e) {
    pred[[i]] <- e$message
    print(paste0("Error for model ", mlmethods[i]))
  })
  gc()
}

# Check the performance of the models.
results <- resamples(fit_model)
dotplot(results) 

# Plot the roc curve with several thresholds using probabilities.
roc_list  <- list()
pred_list <- list()
perf_list <- list()
auc_list  <- list()

for(id in 1:length(mlmethods)){
  roc_list[[id]]  <- predict(fit_model[[id]], 
                             newdata = val.data, 
                             type    = "prob")
  pred_list[[id]] <- prediction(roc_list[[id]][,1], 
                                val.data$Status_2011, 
                                label.ordering = c("DEAD", "ALIVE")
                                )
  perf_list[[id]] <- performance(pred_list[[id]], 
                                 measure   = "tpr", 
                                 x.measure = "fpr")
  auc_list[[id]]  <- performance(pred_list[[id]], 
                                 measure = "auc")@y.values
}

plot(perf_list[[1]], main = "ROC Curves", col = "blue", lwd = 2)
plot(perf_list[[2]], col = "red", lwd = 2, add = T)
plot(perf_list[[3]], col = "orange", lwd = 2, add = T)
plot(perf_list[[4]], col = "black", lwd = 2, add = T)
plot(perf_list[[5]], col = "green", lwd = 2, add = T)
plot(perf_list[[6]], col = "brown", lwd = 2, add = T)
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "gray")
legend("bottomright", legend = paste(mlmethods,"-", 
                                     paste0("AUC = ", 
                                            round(unlist(auc_list), 3))), 
       col = c("blue", "red", "orange", "black", "green", "brown"), 
       lty = 1)

# qda had the best performance.