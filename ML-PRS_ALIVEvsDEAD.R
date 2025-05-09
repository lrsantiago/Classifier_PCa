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
trainc              <- createDataPartition(PRS_datac$pheno,
                                           times = 1,
                                           p     = 0.7,
                                           list  = F)
trainc.data         <- PRS_datac[trainc,]
trainc.data$pheno   <- as.factor(trainc.data$pheno)
valc.data           <- PRS_datac[-trainc,]
valc.data$pheno     <- as.factor(valc.data$pheno)

# 10 fold Cross-validation with multi class summary to improve the model"s accuracy.
ctrl  <- trainControl(method          = "cv",
                      number          = 10,
                      savePredictions = T,
                      summaryFunction = twoClassSummary,
                      classProbs      = T)

# Set the methods available in the caret package.
mlmethods <- c("glm", "lda", "qda", "knn",
               "rpart", "rf", "gbm", "nnet",
               "svmLinear", "svmRadial")

# Train the model with the selected methods.
fit_modelc <- list()
predc      <- list()


for(i in 1:length(mlmethods)){
  tryCatch({
    fit_modelc[[i]] <- train(as.factor(condition) ~ PRS+PSA+age,
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

# Predict aggressive cancer on the test data and calculate the AUC for each model.
predc_roc <- list()
roc_lda  <- list()

for(id in 1:length(mlmethods)){
  predc_roc[[id]] <- predict(fit_modelc[[id]], newdata = valc.data, type = "raw")
  predc_roc[[id]] <- ifelse(predc_roc[[id]] == "Controls", 0, 1)
  roc_lda[[id]]  <- roc(as.factor(valc.data$condition),
                       predc_roc[[id]],
                       levels    = c("Controls", "Cases"),
                       auc       = T,
                       ci        = T,
                       direction = "<")
}


# Get the sensitivity, specificity, accuracy, kappa, and AUC for the models.
comparisons <- data.frame(Models   = mlmethods,
                          se       = rep(NA, length(mlmethods)),
                          sp       = rep(NA, length(mlmethods)),
                          accuracy = rep(NA, length(mlmethods)),
                          kappa    = rep(NA, length(mlmethods)),
                          auc      = rep(NA, length(mlmethods)))

for(id in 1:length(mlmethods)){
  comparisons$se[id]       <- predc[[id]]$byClass[1]
  comparisons$sp[id]       <- predc[[id]]$byClass[2]
  comparisons$accuracy[id] <- predc[[id]]$overall[1]
  comparisons$auc[id]      <- roc_lda[[1]]$auc
  comparisons$kappa[id]    <- predc[[id]]$overall[2]
}


set.seed(120)
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
