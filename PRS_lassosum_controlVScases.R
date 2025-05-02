#!/usr/bin/env Rscript
# R version 4.0.0

setwd("/wroking_dir/controlvscases/lassosum")


libraries <- c("lassosum", "data.table", "methods", 
               "magrittr", "parallel", "ggplot2", 
               "aod", "pROC", "ggsignif")

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

cl <- makeCluster(20)
system("mkdir -p controlvscases/lassosum")
setwd("controlvscases/lassosum")

# Load the GWAS base data converted to OR (exp(BETA)), which is the 1000G european prostate cancer data.
sum_stat <- "../../GWAS_data/base_data.or"
bfile    <- "../target_data_selection"

# Load the European population data (hg38).
ld.file <- "EUR.hg38"

# Load the summary statistics reference base data.
ss <- fread(sum_stat)

# Remove P-values = 0, which causes errors in the transformation.
ss <- ss[!'P' == 0]

# Transform the P-values into correlation
cor <- p2cor(p    = ss$P, 
             n    = ss$N, 
             sign = log(ss$OR))

# Run the lassosum pipeline
out <- lassosum.pipeline(cor        = cor,
                         chr        = ss$CHR,
                         pos        = ss$BP,
                         A1         = ss$A1,
                         A2         = ss$A2,
                         ref.bfile  = bfile,
                         test.bfile = bfile,
                         LDblocks   = ld.file,
                         cluster    = cl)

# Load the covariate file for control vs cases.
cov           <- fread("../target_data_sign.cov")
cov           <- as.data.frame(cov)
colnames(cov) <- c('FID', 'IID', 'age', 
                   'Gleason', 'PSA', 
                   'PC1', 'PC2')
covariate     <- cov[, c('FID', 'IID', 
                         'age', 'PC1', 'PC2')]

# Do not specify phenotypic data, as the function will use the .fam file automatically.
target_data  <- validate(out, 
                         test.bfile        = bfile, 
                         covar             = covariate, 
                         exclude.ambiguous = T, 
                         cluster           = cl)

# Create a data.frame with the best lambda, shrinkage parameter, best-fit PRS and the maximum R2, which says how much of phenotypic variation the best-fit PRS explains.
lass_res <- data.frame('R2'          = max(target_data$validation.table$value) ^ 2, 
                       'Best_lambda' = target_data$best.lambda, 
                       'Shrinkage'   = target_data$best.s)

# Get the best PRS results and merge them with the covariates.
lass              <- merge(target_data$results.table, cov, by = c('FID', 'IID'))
colnames(lass)[4] <- "PRS"

write.table(lass, 'lass_selection.txt', 
            quote     = F, 
            row.names = F, 
            col.names = T, 
            sep       = '\t')

write.table(lass_res, 'lass_lambda_selection.txt', 
            quote     = F, 
            row.names = F, 
            col.names = T, 
            sep       = '\t')

q(save = 'no')