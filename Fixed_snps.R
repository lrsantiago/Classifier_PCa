#!/usr/bin/env Rscript 
# R version 4.0.0

library('VariantAnnotation')
library('ggplot2')
library('rcompanion')
library('data.table')

setwd('variant_calling')

#Load vcf and the list of case and control samples.
vcf            <- readVcf(file = 'all_imputed.vcf.gz', genome = 'hg38')
master         <- read.csv('../TAPG_TURP_DATA.csv')
master         <- subset(master, pheno == 1)
samples_all    <- read.table('sample_names.txt')
samples_all$V2 <- gsub(pattern = '(.+/)|(_.+)', x = samples_all$V1, replacement = '')
samples_all$V3 <- ifelse(samples_all$V2 %in% master$LABNO, 0, 1)
case_sample    <- subset(samples_all, V3 == 1)
case_sample    <- data.frame('V1' = case_sample$V1)
control_sample <- subset(samples_all, V3 == 0)
control_sample <- data.frame('V1' = control_sample$V1)

# Extract the names for cases and controls only.
total_samples <- samples(header(vcf))
case          <- subset(total_samples, 
                        total_samples %in% case_sample$V1)
control       <- subset(total_samples, 
                        total_samples %in% control_sample$V1)

# Get the genotypes for all individuals within each group.
gt         <- geno(vcf)$GT
case_gt    <- gt[, case]
control_gt <- gt[, control]

# Extract SNPs which are fixed in either cases or controls.
fixed_case    <- apply(case_gt, 1, function(x) all(x == x[1]))
fixed_control <- apply(control_gt, 1, function(x) all(x == x[1]))

# Extract the fixed SNPs in cases that are found in controls.
case_diff_control <- sapply(1:nrow(case_gt), 
                            function(i) case_gt[i,1] != control_gt[i,1])

# Generate the new vcf files for SNPs fixed in cases, fixed in controls, and found in cases and missing in controls.
writeVcf(vcf[case_diff_control, ], 
         filename = 'cases_diffs.vcf')

writeVcf(vcf[fixed_case, ], 
         filename = 'fixed_case.vcf')

writeVcf(vcf[fixed_control, ], 
         filename = 'fixed_control.vcf')
