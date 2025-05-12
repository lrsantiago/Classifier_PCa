#!/usr/bin/env Rscript
# R version 4.1.0
####### Classifier from a Polygenic risk score (PRS) to predict prostate cancer aggressiveness.#######

libraries <- c("caret", "MLeval", "mlbench",
               "ROCR", "pROC", "unix", "party",
               "neuralnet")

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

##### ALIVE vs DEAD #####
# Prepare the data set.
lasscc              <- read.table("lass_selection.txt", header = T)
lasscc$IID          <- gsub(pattern = "(.+\\/)|(_.+)", x = lasscc$IID, replacement = "")
lasscc$pheno        <- as.factor(lasscc$pheno)
lasscc$condition    <- ifelse(lasscc$pheno == 1, "Controls", "Cases")
lasscc$OR           <- exp(lasscc$PRS)
PRS_datac           <- lasscc
PRS_datac$pheno     <- as.numeric(PRS_datac$pheno)
PRS_datac$pheno     <- ifelse(PRS_datac$pheno == 2, 1, 0)
PRS_datac$condition <- ifelse(PRS_datac$condition == "Controls", "ALIVE", "DEAD")
PRS_datac$CAPRA     <- master$CAPRA[match(PRS_datac$IID, master$LABNO)]

set.seed(2023)
# Split the data into Training: 70%, and Validation: 30%.
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

# Train the model with the selected methods with PRS + PSA + age.
fit_modelc <- list()
predc      <- list()


for(i in 1:length(mlmethods)){
  tryCatch({
    fit_modelc[[i]] <- train(as.factor(condition) ~ PRS+PSA+age,
                             data       = trainc.data,
                             method     = mlmethods[i],
                             trControl  = ctrl,
                             preProcess = c("center","scale"),
                             tuneLength = 20,
                             metric     = "Accuracy")
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
  predc_roc[[id]] <- ifelse(predc_roc[[id]] == "ALIVE", 0, 1)
  roc_lda[[id]]  <- roc(as.factor(valc.data$condition),
                       predc_roc[[id]],
                       levels    = c("ALIVE", "DEAD"),
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



# Train the model with the selected methods with PRS + PSA + Gleason.
fit_modelcg <- list()
predcg      <- list()


for(i in 1:length(mlmethods)){
  tryCatch({
    fit_modelcg[[i]] <- train(as.factor(condition) ~ PRS+PSA+age,
                             data       = trainc.data,
                             method     = mlmethods[i],
                             trControl  = ctrl,
                             preProcess = c("center","scale"),
                             tuneLength = 20,
                             metric     = "Accuracy")
    predcg[[i]]      <- confusionMatrix(predict(fit_modelc[[i]],
                                               newdata = valc.data),
                                       as.factor(valc.data$condition)
    )
    print(paste0("Prediction done for model ", mlmethods[i]))
  }, error = function(e) {
    predcg[[i]] <- e$message
    print(paste0("Error for model ", mlmethods[i]))
  })
  gc()
}

# Predict aggressive cancer on the test data and calculate the AUC for each model.
predcg_roc <- list()
rocg_lda  <- list()

for(id in 1:length(mlmethods)){
  predcg_roc[[id]] <- predict(fit_modelcg[[id]], newdata = valc.data, type = "raw")
  predcg_roc[[id]] <- ifelse(predcg_roc[[id]] == "ALIVE", 0, 1)
  rocg_lda[[id]]  <- roc(as.factor(valc.data$condition),
                       predc_roc[[id]],
                       levels    = c("ALIVE", "DEAD"),
                       auc       = T,
                       ci        = T,
                       direction = "<")
}


# Get the sensitivity, specificity, accuracy, kappa, and AUC for the models.
comparisons_cg <- data.frame(Models   = mlmethods,
                          se       = rep(NA, length(mlmethods)),
                          sp       = rep(NA, length(mlmethods)),
                          accuracy = rep(NA, length(mlmethods)),
                          kappa    = rep(NA, length(mlmethods)),
                          auc      = rep(NA, length(mlmethods)))

for(id in 1:length(mlmethods)){
  comparisons_cg$se[id]       <- predc[[id]]$byClass[1]
  comparisons_cg$sp[id]       <- predc[[id]]$byClass[2]
  comparisons_cg$accuracy[id] <- predc[[id]]$overall[1]
  comparisons_cg$auc[id]      <- roc_lda[[1]]$auc
  comparisons_cg$kappa[id]    <- predc[[id]]$overall[2]
}



## Neural network ##
set.seed(123)
# Split the data into Training: 70%, and Validation: 30%.
trainc              <- createDataPartition(PRS_datac$pheno,
                                           times = 1,
                                           p     = 0.7,
                                           list  = F)
trainc.data         <- PRS_datac[trainc,]
trainc.data$pheno   <- as.factor(trainc.data$pheno)
valc.data           <- PRS_datac[-trainc,]
valc.data$pheno     <- as.factor(valc.data$pheno)

# Neural network method with 70% training and 30% validation split.
nc <- list()
for(i in 1:3){
  nc[[i]] <- tryCatch({neuralnet(pheno ~ PRS+PSA+Gleason,
                                     data          = trainc.data,
                                     hidden        = i,
                                     algorithm     = "rprop+",
                                     rep           = 5,
                                     err.fct       = "ce",
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

reps <- c(5,2,2)

# Compute the predicted values for ALIVE (0) and DEAD (1) on the validation data.
  outc  <- list()
  tabc  <- list()
  cm_nn <- list()
  auc   <- list()

  for(id in 1:3){
    outc[[id]]     <- predict(nc[[id]], valc.data, rep = reps[id])
    outc[[id]][,1] <- ifelse(outc[[id]][,1] > 0.5, 0, 1)
    tabc[[id]]     <- table(outc[[id]][,1], valc.data$pheno)
    cm_nn[[id]]    <- confusionMatrix(tabc[[id]])
    auc[[id]]      <- roc(as.factor(valc.data$condition),
                          outc[[id]][,1],
                          levels    = c("ALIVE", "DEAD"),
                          auc       = T,
                          ci        = T, direction = "<")
  }

# Get the sensitivity, specificity, accuracy, kappa, and AUC for the model.
  comparisons_nc <- data.frame(Models = paste0("ANN-", 1:3, " hidden layer"),
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




## Neural network (PRS+PSA+age)##

# Neural network method with 70% training and 30% validation split.
nc2 <- list()
for(i in 1:3){
  nc2[[i]] <- tryCatch({neuralnet(pheno ~ PRS+PSA+age,
                                 data          = trainc.data,
                                 hidden        = i,
                                 algorithm     = "rprop+",
                                 rep           = 5,
                                 err.fct       = "ce",
                                 linear.output = F,
                                 stepmax       = 1500000)},
                      error = function(e) {
                        message("length zero: ", e$message)
                        return(NULL)})
}

# Select the best rep value based on the lowest error rate.
for(i in 1:5) {
png(file = sprintf("hl%s.png", i))
plot(nc2[[1]], rep = i)
dev.off()
}

reps2 <- c(3,2,2)

# Compute the predicted values for ALIVE (0) and DEAD (1) on the validation data.
outc2  <- list()
tabc2  <- list()
cm_nn2 <- list()
auc2   <- list()

for(id in 1:3){
  outc2[[id]]     <- predict(nc2[[id]], valc.data, rep = reps2[id])
  outc2[[id]][,1] <- ifelse(outc2[[id]][,1] > 0.5, 0, 1)
  tabc2[[id]]     <- table(outc2[[id]][,1], valc.data$pheno)
  cm_nn2[[id]]    <- confusionMatrix(tabc2[[id]])
  auc2[[id]]      <- roc(as.factor(valc.data$condition),
                        outc2[[id]][,1],
                        levels    = c("ALIVE", "DEAD"),
                        auc       = T,
                        ci        = T, direction = "<")
}

# Get the sensitivity, specificity, accuracy, kappa, and AUC for the model.
comparisons_nc2 <- data.frame(Models = paste0("ANN-", 1:3, " hidden layer"),
                             se        = rep(NA, length(reps)),
                             sp        = rep(NA, length(reps)),
                             accuracy  = rep(NA, length(reps)),
                             kappa     = rep(NA, length(reps)),
                             auc       = rep(NA, length(reps)))

for(id in 1:3){
  comparisons_nc2$se[id]       <- cm_nn2[[id]]$byClass[1]
  comparisons_nc2$sp[id]       <- cm_nn2[[id]]$byClass[2]
  comparisons_nc2$accuracy[id] <- cm_nn2[[id]]$overall[1]
  comparisons_nc2$auc[id]      <- auc2[[id]]$auc
  comparisons_nc2$kappa[id]    <- cm_nn2[[id]]$overall[2]
}


## Neural network (PRS+PSA+Gleason+CAPRA)##

# Neural network method with 70% training and 30% validation split.
ncc <- list()
for(i in 1:3){
  ncc[[i]] <- tryCatch({neuralnet(pheno ~ PRS+PSA+Gleason+CAPRA,
                                     data          = trainc.data,
                                     hidden        = i,
                                     algorithm     = "rprop+",
                                     rep           = 5,
                                     err.fct       = "ce",
                                     linear.output = F,
                                     stepmax       = 1500000)},
                        error = function(e) {
                          message("length zero: ", e$message)
                          return(NULL)})
}

# Select the best rep value based on the lowest error rate.
for(i in 1:5) {
png(file = sprintf("PRS-PSA-gleason-capra%s.png", i))
plot(ncc[[1]], rep = i)
dev.off()
}

reps <- c(2,2,3)

# Compute the predicted values for ALIVE (0) and DEAD (1) on the validation data.
  outcc  <- list()
  tabcc  <- list()
  ccm_nn <- list()
  aucc   <- list()

  for(id in 1:3){
    outcc[[id]]     <- predict(ncc[[id]], valc.data, rep = reps[id])
    outcc[[id]][,1] <- ifelse(outcc[[id]][,1] > 0.5, 0, 1)
    tabcc[[id]]     <- table(outcc[[id]][,1], valc.data$pheno)
    ccm_nn[[id]]    <- confusionMatrix(tabcc[[id]])
    aucc[[id]]      <- roc(as.factor(valc.data$condition),
                          outc[[id]][,1],
                          levels    = c("ALIVE", "DEAD"),
                          auc       = T,
                          ci        = T, direction = "<")
  }

# Get the sensitivity, specificity, accuracy, kappa, and AUC for the model.
  comparisons_ncc <- data.frame(Models = paste0("ANN-", 1:3, " hidden layer"),
                               se        = rep(NA, length(reps)),
                               sp        = rep(NA, length(reps)),
                               accuracy  = rep(NA, length(reps)),
                               kappa     = rep(NA, length(reps)),
                               auc       = rep(NA, length(reps)))

  for(id in 1:3){
    comparisons_ncc$se[id]       <- ccm_nn[[id]]$byClass[1]
    comparisons_ncc$sp[id]       <- ccm_nn[[id]]$byClass[2]
    comparisons_ncc$accuracy[id] <- ccm_nn[[id]]$overall[1]
    comparisons_ncc$auc[id]      <- aucc[[id]]$auc
    comparisons_ncc$kappa[id]    <- ccm_nn[[id]]$overall[2]
  }


  ## Neural network (CAPRA)##

  # Split the data into Training: 70%, and Validation: 30%.

  # Neural network method with 70% training and 30% validation split.
  ncapra <- list()
  for(i in 1:3){
    ncapra[[i]] <- tryCatch({neuralnet(pheno ~ CAPRA,
                                       data          = trainc.data,
                                       hidden        = i,
                                       algorithm     = "rprop+",
                                       rep           = 5,
                                       err.fct       = "ce",
                                       linear.output = F,
                                       stepmax       = 1500000)},
                          error = function(e) {
                            message("length zero: ", e$message)
                            return(NULL)})
  }

  # Select the best rep value based on the lowest error rate.
  for(i in 1:5) {
  png(file = sprintf("capra%s.png", i))
  plot(ncapra[[1]], rep = i)
  dev.off()
  }

  reps <- c(3,3,2)

  # Compute the predicted values for ALIVE (0) and DEAD (1) on the validation data.
    outcapra  <- list()
    tabcapra  <- list()
    capram_nn <- list()
    aucapra   <- list()

    for(id in 1:3){
      outcapra[[id]]     <- predict(ncapra[[id]], valc.data, rep = reps[id])
      outcapra[[id]][,1] <- ifelse(outcapra[[id]][,1] > 0.5, 0, 1)
      tabcapra[[id]]     <- table(outcapra[[id]][,1], valc.data$pheno)
      capram_nn[[id]]    <- confusionMatrix(tabcapra[[id]])
      aucapra[[id]]      <- roc(as.factor(valc.data$condition),
                            outc[[id]][,1],
                            levels    = c("ALIVE", "DEAD"),
                            auc       = T,
                            ci        = T, direction = "<")
    }

  # Get the sensitivity, specificity, accuracy, kappa, and AUC for the model.
    comparisons_ncapra <- data.frame(Models = paste0("ANN-", 1:3, " hidden layer"),
                                 se        = rep(NA, length(reps)),
                                 sp        = rep(NA, length(reps)),
                                 accuracy  = rep(NA, length(reps)),
                                 kappa     = rep(NA, length(reps)),
                                 auc       = rep(NA, length(reps)))

    for(id in 1:3){
      comparisons_ncapra$se[id]       <- capram_nn[[id]]$byClass[1]
      comparisons_ncapra$sp[id]       <- capram_nn[[id]]$byClass[2]
      comparisons_ncapra$accuracy[id] <- capram_nn[[id]]$overall[1]
      comparisons_ncapra$auc[id]      <- aucapra[[id]]$auc
      comparisons_ncapra$kappa[id]    <- capram_nn[[id]]$overall[2]
    }
