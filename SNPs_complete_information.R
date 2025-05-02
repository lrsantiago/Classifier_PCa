#!/usr/bin/env Rscript
# R version 4.0.0

library("rcompanion")
library("data.table")

setwd('~/NGS-store/leandro/TAPG_TURP/variant_call/merged_vcfs/noGleason')

# Load the sample information.
samples_all <- read.table('sample_names.txt')

# Load the results from Plink logistic regression combined with all annotated SNPs and saved in bed format.
all_filtered_annotated            <- read.delim("bedfiles/all_imputed_complete.bed", 
                                                header = F)

colnames(all_filtered_annotated)  <- c('chr', 'start', 'end', 
                                       'SNP', 'REF', 'ALT', 'V7', 
                                       'INFO', 'OR', 'L95', 'U95', 'P')

# Get the minimum allele frequency (MAF) and imputation score information from the loaded bed file.
all_filtered_annotated$MAF        <- gsub(".+MAF=|;.+", 
                                          replacement = "", 
                                          all_filtered_annotated$INFO)
all_filtered_annotated$IMPUTATION <- gsub(".+INFO=|;.+", 
                                          replacement = "", 
                                          all_filtered_annotated$INFO)

# Get only the INFO column from the loaded bed file and remove the first part.
all_info <- as.character(all_filtered_annotated$INFO)
all_info <- gsub(".+;ANN=[A-Z]+\\|", replacement = "", all_info)

# Split all content between '|' into different parts in separated lists in a list object.
all_info <- strsplit(x = all_info, split = "|", fixed = T)

# From the above list, get only position, effect and gene information.
snp_position   <- sapply(all_info, function(x) x[1])
snp_effect     <- sapply(all_info, function(x) x[2])
snp_gene       <- sapply(all_info, function(x) x[3])
snp_geneNoexon <- gsub(pattern = '(exon-)|(\\..+)', x = snp_gene, replacement = '')

# From the exon IDs get the gene IDs.
# Get the exon IDs.
exonid <- data.frame(exonID = snp_geneNoexon[grep(pattern = '[A-Z]+_.+', x = snp_geneNoexon)])

# Remove duplicates.
exonid_uniq <- subset(exonid, !duplicated(exonID))

# Export the unique exon IDs and submit it to a biojava-based script to retrieve gene IDs from exon IDs.
write.table(exonid_uniq, file = 'ids.txt', quote = F, col.names = F, row.names = F)
system(command = paste('sh', 'get_geneid.sh', sep = ' '))

# Load the results from the above command and replace exon IDs by gene IDs.
exonid_geneid <- read.delim('exon_geneid.txt', header = F)

for (id in 1:nrow(exonid_geneid)) {
  snp_geneNoexon <- gsub(pattern = exonid_geneid$V1[id], 
                         x = snp_geneNoexon, 
                         replacement = exonid_geneid$V2[id])
}


# Get all relevant SNP information from the loaded bed file.
chromosome  <- all_filtered_annotated$chr
position    <- all_filtered_annotated$start
REF         <- all_filtered_annotated$REF
ALT         <- all_filtered_annotated$ALT
OR          <- round(all_filtered_annotated$OR, digits = 3)
L95CI       <- round(all_filtered_annotated$L95, digits = 3)
U95CI       <- round(all_filtered_annotated$U95, digits = 3)
MAF         <- round(as.numeric(all_filtered_annotated$MAF), digits = 3)
N           <- rep(nrow(samples_all), nrow(all_filtered_annotated))

# Get the rs IDs.
rsID    <- read.delim("bedfiles/all_imputedrsID.bed", header = F)
rsID$V4 <- gsub(pattern = '.+;', x = rsID$V3, replacement = '')
rsID$V3 <- gsub(pattern = ';.+', x = rsID$V3, replacement = '')

# Generate a data frame with all information combined.
snp_data      <- data.frame("Chromosome"       = chromosome,
                            "Position_in_chr"  = position,
                            "Variant"          = all_filtered_annotated$SNP,
                            "variant_effect"   = snp_effect,
                            "REF"               = REF,
                            "ALT"               = ALT,
                            "MAF"              = MAF,
                            "Position_in_gene" = snp_position,
                            "Genes_affected"   = snp_geneNoexon,
                            "OR"               = OR,
                            "L95CI"            = L95CI,
                            "U95CI"            = U95CI,
                            "N"                = N,
                            "P"                = all_filtered_annotated$P,
                            "INFO"             = all_filtered_annotated$IMPUTATION)

snp_data$rsID <- rsID$V4[match(snp_data$Variant, rsID$V3)]

# Replace the genes in the Genes_affected column that did not match the genes found in the rs IDs list.
DbSnp              <- fread('bedfiles/DbSnp.bed', header = F, nThread = 30)
snp_data$genesrsID <- DbSnp$V4[match(snp_data$rsID, DbSnp$V3)]
snp_data$genesrsID <- gsub(pattern = '(.+GENEINFO=)|(:.+)', x = snp_data$genesrsID, replacement = '')


for (id in 1:nrow(snp_data)) {
  snp_data$match[id] <- ifelse(as.character(snp_data$Genes_affected[id]) == as.character(snp_data$genesrsID[id]), 'yes', 'no')
}

snpmatch             <- subset(snp_data, match == 'yes')
snpNOmatch           <- subset(snp_data, match == 'no')
snpNOmatch$genesrsID <- gsub(pattern = '(.+GENEINFO=)|(:.+)', x = snpNOmatch$genesrsID, replacement = '')
snp_data             <- snp_data[, 1:16]

# After a manual check of the difference between the listed genes in the DbSnp file and the snpEff annotation, load and merge the files.
write.csv(snpNOmatch, 'snpNOmatch.csv')
snpNOmatch     <- read.csv('snpNOmatch.csv')
snp_data_temp  <- rbind(snpmatch[,1:16], snpNOmatch)
snpmatch2      <- subset(snp_data, !(Variant %in% snp_data_temp$Variant))
snpmatch2$INFO <- NA
snpmatch2$rsID <- NA
snp_data_final <- rbind(snp_data_temp, snpmatch2)

# After running the vcf_fixed_snps.R script, load the results.
private_control               <- read.table('control_privateSNPs.txt')
snp_data_final$privateControl <- ifelse(snp_data_final$Variant %in% private_control$V1, 'yes', 'no')
fixed_control                 <- read.table('fixed_control.txt')
snp_data_final$fixedINcontrol <- ifelse(snp_data_final$Variant %in% fixed_control$V1, 'yes', 'no')

# Keep the variants that are not common between cases and controls.
cases_diff               <- read.table('cases_diffs.txt')
snp_data_sign            <- subset(snp_data_final, P < 0.05)
snp_data_sign$cases_diff <- ifelse(snp_data_sign$Variant %in% cases_diff$V1, 'yes', 'no')
snps                     <- subset(snp_data_sign, cases_diff == 'yes')

# Keep only the variants found in the clinvar database.
clinvar          <- read.delim('bedfiles/clinvar.bed', header = F)
clinvar$ids      <- paste(clinvar[,1], clinvar[,2], clinvar[,5], clinvar[,6], sep = ':')
snp_sign         <- subset(snp_data, P < 0.05)
snp_sign_clinvar <- snp_sign[which(snp_sign$Variant %in% clinvar$ids), ]
snp_clinvar      <- clinvar[which(clinvar$ids %in% snp_sign$Variant), ]

# Combine the that are not common between cases and controls and the ones found in the ClinVar database and use them for the polygenic score derivation.
selected_snps    <- c(snp_clinvar$ids, snps$Variant)
selected_snps    <- data.frame('V1' = unique(c(snp_clinvar$ids, snps$Variant)))

# Export everything as csv and txt files.
write.csv(snp_data_final, file = 'all_imputed_complete.csv')
write.table(snp_data_final, file = 'all_imputed_complete.txt', 
            quote = F, 
            row.names = F, 
            col.names = T, 
            sep = '\t')
write.table(selected_snps, 'snp_for_PRS.txt', 
            quote = F, 
            row.names = F, 
            col.names = F, 
            sep = '\t')
