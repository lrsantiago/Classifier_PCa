#!/usr/bin/env Rscript
# R version 4.1.0
## Survival analysis using PRS to predict prostate cancer aggressiveness.

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


setwd("controlvscases/lassosum")

# Prepare the data set.
lasscc                     <- read.table("lass_selection.txt", header = T)
lasscc$IID                 <- gsub(pattern = "(.+\\/)|(_.+)", x = lasscc$IID, replacement = "")
surv_data_prs              <- lasscc
surv_data_prs$condition    <- master$Status_2011[match(surv_data_prs$IID, master$LABNO)]
surv_data_prs$CAPRA        <- master$CAPRA[match(surv_data_prs$IID, master$LABNO)]
surv_data_prs$pheno        <- ifelse(surv_data_prs$condition == "ALIVE", 0, 1)
surv_data_prs$pheno        <- as.factor(surv_data_prs$pheno)
surv_data_prs$survival2011 <- master$survival2011[match(surv_data_prs$IID, master$LABNO)] 
surv_data_prs$Age          <- master$Age[match(surv_data_prs$IID, master$LABNO)] 
surv_data_prs$Age          <- as.integer(surv_data_prs$Age)
surv_data_prs$survival2011 <- as.integer(surv_data_prs$survival2011)

set.seed(149)
## Multivariate cox regression model.
res.cox_all <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PRS+PSA+Gleason+CAPRA, data = surv_data_prs)
multi_cox   <- summary(res.cox_all)

# Testing the proportional hazard assumption
PH1b   <- cox.zph(res.cox_all) 

# Test if removing one of the variables would affect the model. 
test1    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PRS+PSA+Gleason+CAPRA, data = surv_data_prs)
test2    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PRS+PSA+Gleason, data = surv_data_prs)
test3    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PRS+PSA+CAPRA, data = surv_data_prs)
test4    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PRS+Gleason+CAPRA, data = surv_data_prs)
test5    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PSA+Gleason+CAPRA, data = surv_data_prs)
test6    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PRS+PSA, data = surv_data_prs)
test7    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PRS+Gleason, data = surv_data_prs)
test8    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PRS+CAPRA, data = surv_data_prs)
test9    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PSA+CAPRA, data = surv_data_prs)
test10    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PSA+Gleason, data = surv_data_prs)
test11    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ Gleason+CAPRA, data = surv_data_prs)
test12    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PRS, data = surv_data_prs)
test13    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PSA, data = surv_data_prs)
test14    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ Gleason, data = surv_data_prs)
test15    <- coxph(Surv(survival2011, as.numeric(pheno)) ~ CAPRA, data = surv_data_prs)

res.test <- anova(test1, test2, test3, test4, test5, test6, 
                  test7, test8, test9, test10, test11, 
                  test12, test13, test14, test15, test = "Chisq")

# Survival with cox proportional hazard model with PRS as predictor. 
res.coxuni       <- coxph(Surv(survival2011, as.numeric(pheno)) ~ PRS, data = surv_data_prs)
prs_res          <- summary(res.coxuni)
p.value          <- signif(prs_res$wald["pvalue"], digits=3)
wald.test        <- signif(prs_res$wald["test"], digits=3)
beta             <- signif(prs_res$coef[1], digits=3)
HR               <- signif(prs_res$coef[2], digits=3)
HR.confint.lower <- signif(prs_res$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(prs_res$conf.int[,"upper .95"],2)
HR               <- paste0(HR, " (", HR.confint.lower, " - ", HR.confint.upper, ")") 
res              <- c(beta, HR, wald.test, p.value)
res.cox          <- as.data.frame(t(as.data.frame(res, check.names = F)))
colnames(res.cox)<- c("Beta", "HR (95% CI)", "Wald_test", "P_value")

PH1   <- cox.zph(res.coxuni) 

# Survival with Kaplan-Meier method.
surv_factor_prs      <- surv_data_prs
surv_factor_prs1     <- subset(surv_factor_prs, pheno == 0)
surv_factor_prs2     <- subset(surv_factor_prs, pheno == 1)
surv_factor_prs1$PRS <- mean(surv_factor_prs1$PRS)
surv_factor_prs2$PRS <- mean(surv_factor_prs2$PRS)
surv_factor_prs      <- rbind(surv_factor_prs1, surv_factor_prs2)
surv_factor_prs$PRS  <- as.factor(surv_factor_prs$PRS)
surv_prs             <- survfit(Surv(survival2011, as.numeric(pheno)) ~ PRS, 
                                data = surv_factor_prs)
prs_diff             <- survdiff(Surv(survival2011, as.numeric(pheno)) ~ PRS, 
                                 data = surv_factor_prs)

manwhitney <- wilcox.test(subset(surv_data_prs, 
                                 condition == "ALIVE")$PRS, 
                          subset(surv_data_prs, 
                                 condition == "DEAD")$PRS)

q(save = "no")
