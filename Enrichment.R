#!/usr/bin/env Rscript

library('topGO')

system('mkdir -p variant_calling/Enrichment')
setwd('variant_calling/Enrichment')

genelist  <- read.table('../../GWAS_data/mart_export.txt', 
                        header = T)

# Load the file downloaded from Ensembl, containing gene symbols and GO IDs (saved as mart_exportgo.txt).
goid_list <- read.delim('mart_exportgo.txt', header = T, sep = ',')
goid_list <- goid_list[, -1]
colnames(goid_list) <- c('HGNC_symbol', 'goid')

# Assign all GO terms to its respective gene.
mymerge <- function(x) {
  all_in_one <- paste(unlist(x), sep=",", collapse=",")
  split_term <- unlist(strsplit(all_in_one, split=","))
  return(paste(unique(split_term), sep=",", collapse=","))
}

goid_list     <- aggregate(goid_list[-1], by = list(goid_list$HGNC_symbol), mymerge)
genelist$goid <- goid_list$goid[match(genelist$hgnc_symbol, goid_list$Group.1)]
rm(goid_list)

# Load the file from the SNPs_complete_information.R script saved in csv format.
snp_data <- read.csv('../all_imputed_complete.csv')

# Assign gene symbols and GO IDs to the snp_data data frame.
for (row in 1:nrow(snp_data)) {
  for (row2 in 1:nrow(genelist)) {
    if(snp_data$Chromosome[row] == genelist$chromosome_name[row2] & 
       snp_data$Position_in_chr[row] >= genelist$start_position[row2] & 
       snp_data$Position_in_chr[row] <= genelist$end_position[row2]) {
       snp_data$goid[row]  <- as.character(genelist$goid[row2])
    }
  }
  if(row2 %% 1000 == 0) message(row)
}

# Remove items without GO terms assigned.
cases_goid <- subset(snp_data, !is.na(goid))

# ENRICHMENT ANALYSIS-----------------------------------------------------
# Create a list with GO terms and their respective gene symbols:
go_id   <- strsplit(x = as.character(cases_goid$goid), split = ",")
gene2go <- setNames(object = go_id, nm = cases_goid$Genes_affected)

# Create a vector with p-values and gene symbols as names.
topgo_pvalues        <- cases_goid$P
names(topgo_pvalues) <- cases_goid$Genes_affected

# Rank the p-values lowest to highest.
topgo_pvalues <- sort(x = topgo_pvalues, decreasing = F)

# Create a function to highlight only the significant variants.
topDiffGenes_fisher <- function(allScore) {
  return(allScore < 0.05)
}

# Create the TopGO object for each ontology category.
topgoBP <- new("topGOdata", 
               ontology = 'BP', 
               allGenes = topgo_pvalues, 
               geneSel = topDiffGenes_fisher, 
               nodeSize = 20, 
               annot = annFUN.gene2GO, 
               gene2GO = gene2go)
topgoMF <- new("topGOdata", 
               ontology = 'MF', 
               allGenes = topgo_pvalues, 
               geneSel = topDiffGenes_fisher, 
               nodeSize = 20, 
               annot = annFUN.gene2GO, 
               gene2GO = gene2go)
topgoCC <- new("topGOdata", 
               ontology = 'CC', 
               allGenes = topgo_pvalues, 
               geneSel = topDiffGenes_fisher, 
               nodeSize = 20, 
               annot = annFUN.gene2GO, 
               gene2GO = gene2go)


# Run the tests for each ontology category for Fisher and KS stats.
result_topgoBP_F  <- runTest(object = topgoBP, algorithm = "weight01", statistic = "fisher")
result_topgoMF_F  <- runTest(object = topgoMF, algorithm = "weight01", statistic = "fisher")
result_topgoCC_F  <- runTest(object = topgoCC, algorithm = "weight01", statistic = "fisher")
result_topgoBP_KS <- runTest(object = topgoBP, algorithm = "weight01", statistic = "ks")
result_topgoMF_KS <- runTest(object = topgoMF, algorithm = "weight01", statistic = "ks")
result_topgoCC_KS <- runTest(object = topgoCC, algorithm = "weight01", statistic = "ks")


# Get the GO terms listed in each TopGO object for each ontology category.
allGOBP <- usedGO(object = topgoBP)
allGOMF <- usedGO(object = topgoMF)
allGOCC <- usedGO(object = topgoCC)

# Get the results table for each ontology category. 
topGOResults_BP <- GenTable(topgoBP, 
                            Fisher = result_topgoBP_F, 
                            KS = result_topgoBP_KS, 
                            orderBy ="Fisher", 
                            topNodes = length(allGOBP), 
                            numChar  = 1000)
topGOResults_MF <- GenTable(topgoMF, 
                            Fisher = result_topgoMF_F, 
                            KS = result_topgoMF_KS, 
                            orderBy ="Fisher", 
                            topNodes = length(allGOMF), 
                            numChar  = 1000)
topGOResults_CC <- GenTable(topgoCC, 
                            Fisher = result_topgoCC_F, 
                            KS = result_topgoCC_KS, 
                            orderBy ="Fisher", 
                            topNodes = length(allGOCC), 
                            numChar  = 1000)

# Keep only the significant results according to KS test for each ontology category.
topGOResults_BP$Fisher <- as.numeric(topGOResults_BP$Fisher)
topGOResults_BP$KS     <- as.numeric(topGOResults_BP$KS)
topGOResults_BP_sig    <- subset(topGOResults_BP, KS < 0.05)
topGOResults_MF$Fisher <- as.numeric(topGOResults_MF$Fisher)
topGOResults_MF$KS     <- as.numeric(topGOResults_MF$KS)
topGOResults_MF_sig    <- subset(topGOResults_MF, KS < 0.05)
topGOResults_CC$Fisher <- as.numeric(topGOResults_CC$Fisher)
topGOResults_CC$KS     <- as.numeric(topGOResults_CC$KS)
topGOResults_CC_sig    <- subset(topGOResults_CC, KS < 0.05)

# Combined the results from each ontology category and remove the duplicates.
topGOResults_BP_MF_CC_sig <- rbind(topGOResults_BP_sig, topGOResults_MF_sig, topGOResults_CC_sig)
topGOResults_BP_MF_CC_sig <- subset(topGOResults_BP_MF_CC_sig, !duplicated(GO.ID))

# Export all results.
write.table(topGOResults_BP_sig, file = 'topGOResults_BP_sig.txt', 
            quote = F, 
            row.names = F, 
            col.names = T, 
            sep = '\t')
write.table(topGOResults_MF_sig, file = 'topGOResults_MF_sig.txt', 
            quote = F, 
            row.names = F, 
            col.names = T, 
            sep = '\t')
write.table(topGOResults_CC_sig, file = 'topGOResults_CC_sig.txt', 
            quote = F, 
            row.names = F, 
            col.names = T, 
            sep = '\t')
write.table(topGOResults_BP_MF_CC_sig, file = 'topGOResults_BP_MF_CC_sig.txt', 
            quote = F, 
            row.names = F,
            col.names = T, 
            sep = '\t')

# Export three inputs for pathway analysis:
# For all SNPs.
write.table(snp_data$Genes_affected, 'input_pathway_all.txt', 
            quote = F, 
            row.names = F, 
            col.names = F)

# For the significant SNPs.
write.table(subset(snp_data, P < 0.05)$Genes_affected, 'input_pathway_sign.txt', 
            quote = F, 
            row.names = F, 
            col.names = F)

# For the SNPs used to calculate the Polygenic score.
repanel     <- read.delim('../snp_for_PRS.txt', header = F)
snp_repanel <- subset(snp_data, Variant %in% repanel$V1 & P < 0.05)

write.table(snp_repanel$Genes_affected, 'input_pathway_selection.txt', 
            quote = F, 
            row.names = F, 
            col.names = F)

q(save = 'no')