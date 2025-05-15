#!/usr/bin/env Rscript
# R version 4.1.0
## Survival analysis using PSA and age to predict prostate cancer aggressiveness.

libraries <- c("arm", "survival", "survminer")

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
data_all$pheno       <- as.factor(ifelse(data_all$Status_2011 == "ALIVE", 0, 1))
data_all             <- subset(data_all, !is.na(CAPRA))


# Estimated survival probability from a Multi-variate Cox model (PSA, Gleason and CAPRA).
surv_data     <- data_all[, c("patientid", "Status_2011", "Age", 
                          "survival2011", "PSA", "LABNO", 
                          "Gleason", "pheno", "CAPRA")]

set.seed(150)
## Multivariate cox regression model.
res.cox_all2 <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PSA+CAPRA+Gleason, 
                      data = surv_data)
multi_cox2   <- summary(res.cox_all2)


# Test for the proportional hazard (PH) assumption.
PH1c <- cox.zph(res.cox_all2) 


# Uni-variate Cox model with PSA as the only feature, as it did not violate the PH assumption.
Feature       <- "PSA"
univ_formulas <- sapply(Feature,
                        function(x) as.formula(paste("Surv(survival2011, as.numeric(pheno))~", 
                                                     x)))

univ_models   <- lapply(univ_formulas, function(x){coxph(x, data = surv_data)})

## Extract the relevant info. 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         wald.test<-signif(x$wald["test"], digits=3)
                         beta<-signif(x$coef[1], digits=3);#coeficient beta
                         HR <-signif(x$coef[2], digits=3);#exp(beta) or hazard ration
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, " - ", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("Beta", "HR (95% CI)", "Wald.test ", "p.value")
                         return(res)
                       })

res.coxuni_all           <- t(as.data.frame(univ_results, check.names = F))
colnames(res.coxuni_all) <- c("Beta", "HR (95% CI)", "Wald test", "P-value")
res.coxuni_all           <- as.data.frame(res.coxuni_all)


# Estimated survival probability from Kaplan-Meier method with PSA level used as a predictor
surv_p   <- survfit(Surv(survival2011, as.numeric(pheno)) ~ PSA, 
                    data = surv_factor)

pvalue_p <- survdiff(Surv(survival2011, as.numeric(pheno)) ~ PSA, 
                     data = surv_factor)


q(save = "no")




