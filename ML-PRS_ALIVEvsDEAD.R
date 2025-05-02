#!/usr/bin/env Rscript
# R version 4.0.0
####### Classifier from a Polygenic risk score (PRS) to predict prostate cancer aggressiveness.#######

libraries <- c('caret', 'MLeval', 'mlbench', 
               'ROCR', 'pROC', 'unix')

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

setwd("controlvscases/lassosum")
set.seed(2023)

##### ALIVE vs DEAD #####
lasscc           <- read.table("lass_selection.txt", header = T)
lasscc$IID       <- gsub(pattern = '(.+\\/)|(_.+)', x = lasscc$IID, replacement = '')
lasscc$pheno     <- as.factor(lasscc$pheno)
lasscc$condition <- ifelse(lasscc$pheno == 1, 'Controls', 'Cases')
lasscc$OR        <- exp(lasscc$PRS)
PRS_datac        <- lasscc

# Split the data into Training: 70%, and Validation: 30%.
PRS_datac$pheno     <- as.numeric(PRS_datac$pheno)
PRS_datac$pheno     <- ifelse(PRS_datac$pheno == 2, 1, 0)
PRS_datac$condition <- ifelse(PRS_datac$condition == 'Controls', 'ALIVE', 'DEAD') 

trainc              <- createDataPartition(PRS_datac$pheno, times = 1, p = 0.7, list = F)
trainc.data         <- PRS_datac[trainc,]
trainc.data$pheno   <- as.factor(trainc.data$pheno)
valc.data           <- PRS_datac[-trainc,]
valc.data$pheno     <- as.factor(valc.data$pheno)

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
fit_modelc <- list()
predc      <- list()


for(i in 1:length(mlmethods)){
  tryCatch({
    fit_modelc[[i]] <- train(as.factor(condition) ~ PRS+PSA+Gleason, 
                             data.      = trainc.data, 
                             method.    = mlmethods[i], 
                             trControl  = ctrl, 
                             preProcess = c("center","scale"),
                             tuneLength = 20, 
                             metric.    = 'Accuracy')
    predc[[i]]      <- confusionMatrix(predict(fit_modelc[[i]], 
                                               newdata = valc.data),
                                       as.factor(valc.data$condition)
                                       )
    print(paste0("Prediction done for model ", mlmethods[i]))
  }, error = function(e) {
    predc[[i]] <- e$message
    print(paste0("Error for model ", mlmethods[i]))
  })
  gc()
}

# Plot the roc curve with several thresholds using probabilities for the models that worked 

mlmethods  <- c("glm", "lda", "qda", "knn", "rpart", "rf")
roc_listc  <- list()
predc_list <- list()
perfc_list <- list()
auc_listc  <- list()

for(id in 1:length(mlmethods)){
  roc_listc[[id]]  <- predict(fit_modelc[[id]], 
                              newdata = valc.data, 
                              type    = "prob")
  predc_list[[id]] <- prediction(roc_listc[[id]][,2], 
                                 valc.data$condition, 
                                 label.ordering = c("ALIVE", "DEAD")
                                 )
  perfc_list[[id]] <- performance(predc_list[[id]], 
                                  measure   = "tpr", 
                                  x.measure = "fpr")
  auc_listc[[id]]  <- performance(predc_list[[id]], 
                                  measure = "auc")@y.values
}

plot(perfc_list[[1]], main = "Indolent vs Aggressive", col = "blue", lwd = 2)
plot(perfc_list[[2]], col = "red", lwd = 2, add = T)
plot(perfc_list[[3]], col = "black", lwd = 2, add = T)
plot(perfc_list[[4]], col = "magenta", lwd = 2, add = T)
plot(perfc_list[[5]], col = "green", lwd = 2, add = T)
plot(perfc_list[[6]], col = "brown", lwd = 2, add = T)
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "gray")
legend("bottomright", legend = paste(mlmethods,"-", paste0("AUC = ", 
                                                           round(unlist(auc_listc), 3))
                                     ), 
       col = c("blue", "red", "black", "magenta", "green", "brown"), 
       lty = 1)


## Neural network ##

# Neural network method with 70% training and 30% validation split.
nc <- neuralnet(pheno ~ PRS+PSA+age, 
                 data          = trainc.data, 
                 hidden        = 4,
                 algorithm     = "rprop+",
                 rep           = 3,
                 err.fct       = "ce", 
                 linear.output = F,
                 stepmax       = 1500000)

# Select the best rep value based on the lowest error rate.
for (id in 1:10) {
  plot(n, rep = id)
}

# Compute the predicted values for ALIVE (0) and DEAD (1) on the validation data.
outc     <- predict(nc, valc.data, rep = 3).
outc[,1] <- ifelse(outc[,1] > 0.5, 0, 1)
tabc     <- table(outc[,1], valc.data$pheno)
cm_nn    <- confusionMatrix(tabc)
