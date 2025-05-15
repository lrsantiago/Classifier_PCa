#!/usr/bin/env Rscript
# R version 4.1.0
## Classifier using PSA and age to predict prostate cancer aggressiveness.

libraries <- c("caret", "MLeval", "mlbench",
               "ROCR", "pROC", "unix", "neuralnet")

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

setwd("controlvscases")

# Prepare the data.
data_all             <- read.csv("TAPG_TURP_DATA.csv")
data_all             <- data_all[-1]
data_all$pheno       <- as.factor(ifelse(data_all$Status_2011 == 'ALIVE', 0, 1))
data_all             <- subset(data_all, !is.na(CAPRA))
surv_data            <- data_all[, c("patientid", "Status_2011", "Age",
                                     "survival2011", "PSA", "LABNO",
                                     "Gleason", "pheno", "CAPRA")]

surv_data$Age          <- as.integer(surv_data$Age)
surv_data$survival2011 <- as.integer(surv_data$survival2011)
surv_data$pheno        <- as.numeric(surv_data$pheno)
surv_data$pheno        <- ifelse(surv_data$pheno == 2, 1, 0)

# Split the data into Training: 70%, and Validation: 30%.
set.seed(333)
train            <- createDataPartition(surv_data$pheno, times = 1, p = 0.7, list = F)
train.data       <- surv_data[train,]
train.data$pheno <- as.factor(train.data$pheno)
val.data         <- surv_data[-train,]
val.data$pheno   <- as.factor(val.data$pheno)


# Neural network method with 70% training and 30% validation split.
nc <- list()
for(i in 1:3){
  nc[[i]] <- tryCatch({neuralnet(pheno ~ PSA+Age,
                                 data          = train.data,
                                 hidden        = i,
                                 algorithm     = "rprop+",
                                 rep           = 5,
                                 err.fct       = "ce",
                                 act.fct       = "logistic",
                                 linear.output = F,
                                 stepmax       = 1500000)},
                      error = function(e) {
                        message("length zero: ", e$message)
                        return(NULL)})
}


# Select the best rep value based on the lowest error rate.
for(i in 1:5) {
  png(file = sprintf("hl%s.png", i))
  plot(nc[[1]], rep = i)
  dev.off()
}

reps <- c(2,2,1)

# Compute the predicted values for ALIVE (0) and DEAD (1) on the validation data.
outc  <- list()
tabc  <- list()
cm_nn <- list()
auc   <- list()

for(id in 1:3){
  outc[[id]]     <- predict(nc[[id]], val.data, rep = reps[id])
  outc[[id]][,1] <- ifelse(outc[[id]][,1] > 0.5, 0, 1)
  u              <- union(outc[[id]][,1], val.data$pheno)
  tabc[[id]]     <- table(factor(outc[[id]][,1], u), factor(val.data$pheno, u))
  cm_nn[[id]]    <- confusionMatrix(tabc[[id]])
  auc[[id]]      <- roc(as.factor(val.data$Status_2011),
                        outc[[id]][,1],
                        levels    = c("ALIVE", "DEAD"),
                        auc       = T,
                        ci        = T, direction = "<")
}


# Get the sensitivity, specificity, accuracy, kappa, and AUC for the model.
comparisons_nc <- data.frame(Models = paste0("ANN-", 1:3, " hidden layers"),
                             se        = rep(NA, length(reps)),
                             sp        = rep(NA, length(reps)),
                             accuracy  = rep(NA, length(reps)),
                             kappa     = rep(NA, length(reps)),
                             auc       = rep(NA, length(reps)))

for(id in 1:3){
  comparisons_nc$se[id]       <- cm_nn[[id]]$byClass[1]
  comparisons_nc$sp[id]       <- cm_nn[[id]]$byClass[2]
  comparisons_nc$accuracy[id] <- cm_nn[[id]]$overall[1]
  comparisons_nc$auc[id]      <- auc[[id]]$auc
  comparisons_nc$kappa[id]    <- cm_nn[[id]]$overall[2]
}


# Neural network method with 70% training and 30% validation split with PSA + Gleason.
nc2 <- list()
for(i in 1:3){
  nc2[[i]] <- tryCatch({neuralnet(pheno ~ PSA+Gleason,
                                  data          = train.data,
                                  hidden        = i,
                                  algorithm     = "rprop+",
                                  rep           = 5,
                                  err.fct       = "ce",
                                  act.fct       = "logistic",
                                  linear.output = F,
                                  stepmax       = 1500000)},
                       error = func2tion(e) {
                         message("length zero: ", e$message)
                         return(NULL)})
}


# Select the best rep value based on the lowest error rate.
for(i in 1:5) {
  png(file = sprintf("hl%s.png", i))
  plot(nc2[[1]], rep = i)
  dev.off()
}

reps <- c(2,2,1)

# Compute the predicted values for ALIVE (0) and DEAD (1) on the validation data.
outc2  <- list()
tabc2  <- list()
c2m_nn <- list()
auc2   <- list()

for(id in 1:3){
  outc2[[id]]     <- predict(nc2[[id]], val.data, rep = reps[id])
  outc2[[id]][,1] <- ifelse(outc2[[id]][,1] > 0.5, 0, 1)
  u              <- union(outc2[[id]][,1], val.data$pheno)
  tabc2[[id]]     <- table(factor(outc2[[id]][,1], u), factor(val.data$pheno, u))
  c2m_nn[[id]]    <- confusionMatrix(tabc2[[id]])
  auc2[[id]]      <- roc(as.factor(val.data$Status_2011),
                         outc2[[id]][,1],
                         levels    = c("ALIVE", "DEAD"),
                         auc       = T,
                         ci        = T, direction = "<")
}


# Get the sensitivity, specificity, accuracy, kappa, and AUC for the model.
comparisons_nc2 <- data.frame(Models = paste0("ANN-", 1:3, " hidden layers"),
                              se        = rep(NA, length(reps)),
                              sp        = rep(NA, length(reps)),
                              accuracy  = rep(NA, length(reps)),
                              kappa     = rep(NA, length(reps)),
                              auc       = rep(NA, length(reps)))

for(id in 1:3){
  comparisons_nc2$se[id]       <- c2m_nn[[id]]$byClass[1]
  comparisons_nc2$sp[id]       <- c2m_nn[[id]]$byClass[2]
  comparisons_nc2$accuracy[id] <- c2m_nn[[id]]$overall[1]
  comparisons_nc2$auc[id]      <- auc2[[id]]$auc
  comparisons_nc2$kappa[id]    <- c2m_nn[[id]]$overall[2]
}


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

mlmethods <- c("glm", "lda", "qda", "knn", "rpart", "rf",
               "gbm", "nnet", "svmLinear", "svmRadial")

# Train the model with the selected methods using PSA + Age.
fit_model <- list()
pred      <- list()


for(i in 1:length(mlmethods)){
  tryCatch({
    fit_model[[i]] <- train(as.factor(Status_2011) ~ PSA + Age,
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

# Predict aggressive cancer on the test data and calculate the AUC for each model.
predcg_roc <- list()
rocg_lda  <- list()

for(id in 1:length(mlmethods)){
  predcg_roc[[id]] <- predict(fit_model[[id]], newdata = val.data, type = "raw")
  predcg_roc[[id]] <- ifelse(predcg_roc[[id]] == "ALIVE", 0, 1)
  rocg_lda[[id]]  <- roc(as.factor(val.data$Status_2011),
                         predcg_roc[[id]],
                         levels    = c("ALIVE", "DEAD"),
                         auc       = T,
                         ci        = T,
                         direction = "<")
}


# Get the sensitivity, specificity, accuracy, kappa, and AUC for the models.
comparisons_cg <- data.frame(Models   = mlmethods,
                             accuracy = rep(NA, length(mlmethods)),
                             se       = rep(NA, length(mlmethods)),
                             sp       = rep(NA, length(mlmethods)),
                             kappa    = rep(NA, length(mlmethods)),
                             auc      = rep(NA, length(mlmethods)))

for(id in 1:length(mlmethods)){
  comparisons_cg$se[id]       <- pred[[id]]$byClass[1]
  comparisons_cg$sp[id]       <- pred[[id]]$byClass[2]
  comparisons_cg$accuracy[id] <- pred[[id]]$overall[1]
  comparisons_cg$auc[id]      <- rocg_lda[[id]]$auc
  comparisons_cg$kappa[id]    <- pred[[id]]$overall[2]
}


# Train the model with the selected methods using PSA + Gleason.
fit_model2 <- list()
pred2      <- list()


for(i in 1:length(mlmethods)){
  tryCatch({
    fit_model2[[i]] <- train(as.factor(Status_2011) ~ PSA + Gleason,
                             data       = train.data,
                             method.    = mlmethods[i],
                             trCont.rol = ctrl,
                             preProcess = c("center","scale"),
                             tuneLength = 50,
                             metric.    = "Spec")
    pred2[[i]]      <- confusionMatrix(predict(fit_model2[[i]],
                                               newdata = val.data),
                                       as.factor(val.data$Status_2011)
    )
    print(paste0("Prediction done for model", mlmethods[i]))
  }, error = function(e) {
    pred2[[i]] <- e$message
    print(paste0("Error for model ", mlmethods[i]))
  })
  gc()
}

# Check the performance of the models.
results2 <- resamples(fit_model2)
dotplot(results2)

# Predict aggressive cancer on the test data and calculate the AUC for each model.
predcg2_roc <- list()
rocg2_lda   <- list()

for(id in 1:length(mlmethods)){
  predcg2_roc[[id]] <- predict(fit_model2[[id]], newdata = val.data, type = "raw")
  predcg2_roc[[id]] <- ifelse(predcg2_roc[[id]] == "ALIVE", 0, 1)
  rocg2_lda[[id]]   <- roc(as.factor(val.data$Status_2011),
                           predcg2_roc[[id]],
                           levels    = c("ALIVE", "DEAD"),
                           auc       = T,
                           ci        = T,
                           direction = "<")
}


# Get the sensitivity, specificity, accuracy, kappa, and AUC for the models.
comparisons_cg2 <- data.frame(Models   = mlmethods,
                              accuracy = rep(NA, length(mlmethods)),
                              se       = rep(NA, length(mlmethods)),
                              sp       = rep(NA, length(mlmethods)),
                              kappa    = rep(NA, length(mlmethods)),
                              auc      = rep(NA, length(mlmethods)))

for(id in 1:length(mlmethods)){
  comparisons_cg2$se[id]       <- pred2[[id]]$byClass[1]
  comparisons_cg2$sp[id]       <- pred2[[id]]$byClass[2]
  comparisons_cg2$accuracy[id] <- pred2[[id]]$overall[1]
  comparisons_cg2$auc[id]      <- rocg2_lda[[id]]$auc
  comparisons_cg2$kappa[id]    <- pred2[[id]]$overall[2]
}

q(save = "no")
