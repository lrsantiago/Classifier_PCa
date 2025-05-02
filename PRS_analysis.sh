#!binsh

module load bcftools/1.8
module load plink/1.9-170906
module load anaconda3/2020.02
conda activate r_env_new # Environment with R v4.0.0 where all required packages were installed.
pip3 install CrossMap

NSLOTS=20

# Summary stats from ebi can be found here:
#ftp:ftp.ebi.ac.ukpubdatabasesgwassummary_statistics

## Analysis based on the tutorial accessed at https:choishingwan.github.ioPRS-Tutorial.
# Base data from https:doi.org10.1038s41588-018-0142-8.
# Description of the meta-analysis info is described in the file GWAS_MetaAnalysis_results_summary_public.docx
cd GWAS_data
wget http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip
wget http://practical.icr.ac.uk/blog/wp-content/uploads/docs/OncoArray/GWAS_MetaAnalysis_results_summary_public.docx

# Rearrange data with the following colnames:
#CHR	BP	SNP	A1	A2	N	SE	P	OR	INFO	MAF. A1 is the minor allele and A2 major allele, this is decided by their frequency, they are not ref and alt alleles.
# MAF has no data (.) as all are below 0.01 already. 142233 is the N column (number of individuals).
# A column called ID was added to match different data frames.
# All done in R, then liftOver the file from build 37 to the build 38 version of the human genome.

R -e "

library(data.table)
setDTthreads(20)

base                 <- fread('GWAS_datameta_v3_onco_euro_overall_ChrAll_1_release.txt',
                              header = T,
                              nThread = 20)
base$N               <- rep(142233, nrow(base))
base$MAF             <- rep('.', nrow(base))
base_final           <- base[, c('Chr', 'position', 'MarkerName',
                                 'Allele1', 'Allele2', 'N',
                                 'StdErr', 'Pvalue', 'Effect',
                                 'OncoArray_imputation_r2', 'MAF')]
base_final$Allele1   <- toupper(base_final$Allele1)
base_final$Allele2   <- toupper(base_final$Allele2)
base_final$MarkerName<- gsub(pattern     = '_',
                             x           = base_final$MarkerName,
                             replacement = ':')
colnames(base_final) <- c('CHR','BP','SNP',
                          'A1','A2','N','SE',
                          'P','BETA','INFO','MAF')
base_final$CHR       <- gsub('23',
                             x           =  base_final$CHR,
                             replacement = 'X')
base_final$ID        <- paste('ID', 1:nrow(base_final), sep = '-')

write.table(base_final, 'base_data_hg19.txt',
            quote = F,
            col.names = T,
            row.names = F,
            sep = '\t')

q(save = 'no')

"

awk \
-v OFS='\t' '{print $1, $2, $2, $3, $12}' \
base_data_hg19.txt | \
sed 1d > base_data_hg19.bed

CrossMap bed GRCh37_to_GRCh38.chain.gz base_data_hg19.bed base_data_hg38.bed

R -e "

  library(data.table)
  setDTthreads(20)

  base_final       <- fread('base_data_hg19.txt', header = T)
  hg38             <- fread('base_data_hg38.bed')
  hg38$V4          <- paste(hg38$V1, hg38$V2, gsub(pattern = '.+:[0-9]+:',
                                                   x = hg38$V4,
                                                   replacement = ''),
                            sep = ':')
  base_final$match <- hg38$V5[match(base_final$ID, hg38$V5)]
  base_hg38        <- subset(base_final, !is.na(match))
  base_hg38$SNP    <- hg38$V4[match(base_hg38$match, hg38$V5)]
  base_hg38$BP     <- hg38$V2[match(base_hg38$match, hg38$V5)]

  write.table(base_hg38[, c('CHR','BP','SNP',
                            'A1','A2','N','SE',
                            'P','BETA','INFO','MAF')],
              'base_data_hg38.txt',
              quote.    = F,
              col.names = T,
              row.names = F,
              sep       = '\t')


  q(save = 'no')

"

rm -f *.bed

# Remove duplicated SNPs.
awk '{print $3}' base_data_hg38.txt | sort | uniq -d > duplicated.snp
cat base_data_hg38.txt | grep -vf duplicated.snp > base_data_uniq.txt
rm -f duplicated.snp

# Remove ambiguous SNPs.
cat base_data_uniq.txt | \
awk '!( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")) {print}' \
> base_data.txt

# Convert BETA to OR and remove imputation scores below 0.6.
R -e "
  library(data.table)
  setDTthreads(20)

  dat     <- fread('base_data.txt', header = T)
  dat$OR  <- exp(dat$BETA)
  dat$SNP <- paste0('chr', dat$SNP)
  dat$MAF <- 0.02
  dat     <- subset(dat, INFO > 0.6)

  write.table(dat[, c('CHR', 'BP', 'SNP',
                    'A1', 'A2', 'N', 'SE',
                    'P', 'BETA', 'INFO', 'MAF')],
            'base_data.Transformed',
            quote     = F,
            row.names = F,
            col.names = T,
            sep       = '\t')

  write.table(dat[, c('CHR', 'BP', 'SNP',
                      'A1', 'A2', 'N', 'SE',
                      'P', 'INFO', 'MAF', 'OR')],
              'base_data.or',
              quote.    = F,
              row.names = F,
              col.names = T,
              sep       = '\t')


  q(save = 'no')

"

rm -f base_data_uniq.txt


# Process the target data (vcf with significant SNPs).
# Create a data set with cancer (indolent and agressive prostate cancer samples together) and healthy samples as controls (90 samples of males from the 1000G project used as reference for imputation).
# Remove female samples from the 1000G samples. Download all samples from the 1000G project.
wget http:ftp.1000genomes.ebi.ac.ukvol1ftpdata_collections1000G_2504_high_coverage20130606_g1k_3202_samples_ped_population.txt

R -e "
  # Load the list of samples from the 1000G.
  pop <- read.table('GWAS_data20130606_g1k_3202_samples_ped_population.txt',
                    header = T)

  # Keep only males from the European population.
  eur <- subset(pop, Sex == 1 & Superpopulation == 'EUR')

  write.table(eur$SampleID, 'GWAS_data1000G_male_eur_samples.txt',
            col.names = F,
            row.names = F,
            quote     = F)

  q(save = 'no')

"

# Combine all files from the 1000G project into one vcf.
ls 1000GP_filtered_chr*.gz > file

bcftools concat --threads ${NSLOTS} -f file -Ov | \
bgzip > 1000G_temp.vcf.gz

bcftools index --threads ${NSLOTS} -t 1000G_temp.vcf.gz

bcftools view --threads ${NSLOTS} \
-S 1000G_male_eur_samples.txt 1000G_temp.vcf.gz \
-Ov | \
bgzip > 1000G.vcf.gz

bcftools index --threads ${NSLOTS} -t 1000G.vcf.gz
rm -f file
rm -f 1000G_temp*


## Process healthy samples data.
mv 1000G.vcf* ../healthy_vs_all

cd ../healthy_vs_all

## QC the control data.
#Convert vcf to plink input format.
plink \
--vcf 1000G.vcf.gz \
--allow-extra-chr \
--double-id \
--allow-no-sex \
-no-parents \
--make-bed \
--recode \
--out 1000G

# Assign phenotype info (controls = 1) to the .fam file.
R -e "

fam         <- read.table('1000G.fam')
fam$V6      <- 1
fam$V5      <- 1

write.table(fam, file = '1000G.fam',
            quote     = F,
            row.names = F,
            col.names = F,
            sep       = '\t')

q(save = 'no')

"

# Update the .ped file adding phenotype and sex info.
awk '!($1=$2=$3=$4=$5=$6="")' 1000G.ped > temp2
rm -f 1000G.ped
paste 1000G.fam temp2 > 1000G.ped

rm -f tem*

# Get the MAF and heterozygosity files.
# Remove highly correlated SNPs and generate heterozygosity file.
plink \
--bfile 1000G
--fam 1000G.fam \
--allow-extra-chr \
--hwe 1e-6 \
--maf 0.01 \
--write-snplist \
--out 1000G

plink \
--bfile 1000G \
--allow-extra-chr \
--keep 1000G.fam \
--extract 1000G.snplist \
--indep-pairwise 200 50 0.25 \
--out 1000G

plink \
--bfile 1000G \
--allow-extra-chr \
--extract 1000G.prune.in \
--keep 1000G.fam \
--het \
--out 1000G

# Remove individuals with F coefficients that are more than 3 sd units from the mean.
R -e "

dat2   <- read.table('1000G.het', header = T)
m2     <- mean(dat2$F)
s2     <- sd(dat2$F)
valid2 <- subset(dat2, F <= m2+3*s2 & F >= m2-3*s2)

write.table(valid2[,c(1,2)], '1000G.valid.sample',
            quote.    = F,
            row.names = F)

q(save = 'no')

"

# QC the data.
plink \
--bfile 1000G \
--allow-extra-chr \
--keep 1000G.valid.sample \
--make-bed \
--extract 1000G.snplist \
--out 1000G_QC

# Generate the population stratification data to use as covariate.
plink \
--bfile 1000G_QC \
--indep-pairwise 200 50 0.25 \
--out 1000G_QC

# Calculate the first 5 PCs
plink \
--bfile 1000G_QC \
--extract 1000G_QC.prune.in \
--pca 5 \
--out 1000G_QC


R -e "
fam    <- read.table('1000G_QC.fam')
fam$V5 <- 1
fam$V6 <- 1

write.table(fam, file = '1000G_QC.fam',
            quote     = F,
            row.names = F,
            col.names = F,
            sep       = '\t')

q(save = 'no')

"
module unload bcftools/1.8
module load bcftools/1.4

awk '{print $1}' 1000G_QC.fam > samples.txt
awk '{print $2}' 1000G_QC.bim > snp.txt
bcftools view \
--threads ${NSLOTS} \
-S samples.txt 1000G.vcf.gz \
-Ov | \
bcftools view \
--threads ${NSLOTS} \
-i 'ID=@snp.txt' - | \
bgzip > 1000G_QC.vcf.gz

bcftools index --threads ${NSLOTS} -t 1000G_QC.vcf.gz

rm -f 1000G.*

module unload bcftools/1.4
module load bcftools/1.8

# From the 1000G genome control samples, keep only genes where SNPs from cancer samples are locaded.
# Download the list of genes from Biomart (mart_export.txt).

R-e "
library('data.table')

genelist           <- fread('mart_export.txt')
colnames(genelist) <- c('id', 'transcript', 'start',
                        'end', 'exon', 'chr', 'gene')

# Keep enes without duplicates.
snp_data      <- read.csv('../variant_calling/all_imputed_complete.csv')
snp_data_uniq <- subset(snp_data, !duplicated(Genes_affected))

# keep only chromosomes.
snp_genelist     <- genelist[-grep('(C.+)|(K.+)|(G.+)|(M.+)', genelist$chr),]
snp_genelist     <- subset(snp_genelist, !duplicated(gene))
snp_genelist$chr <- paste0('chr', snp_genelist$chr)

# Remove those found in the biomart file.
remaining              <- snp_data_uniq[-which(snp_data_uniq$Genes_affected %in%
                                                 snp_genelist$gene),]
remaining$start        <- remaining$Position_in_chr - 1000
remaining$end          <- remaining$Position_in_chr + 1000
colnames(remaining)[2] <- 'chr'
bed                    <- rbind(snp_genelist[, c('chr', 'start', 'end')],
                                remaining[, c('chr', 'start', 'end')])

write.table(bed, 'regions_keep.bed',
            quote     = F,
            row.names = F,
            col.names = F,
            sep       = '\t')

q(save = ' no')

"

# Sort the bed file.
bedtools sort -chrThenSizeA -i regions_keep.bed > regions_keep_sorted.bed
rm -f regions_keep.bed

bedtools intersect \
-wa \
-a 1000G_QC.vcf.gz \
-b regions_keep_sorted.bed \
-header | \
bgzip > 1000G_QC_sel.vcf.gz

bcftools index --threads ${NSLOTS} -t 1000G_QC_sel.vcf.gz

# Merge the 1000G genome control samples to the target imputed samples.
cat ../variant_calling/all_imputed.assoc.logistic | \
sed 1d | \
awk '$12<"0.05"' | \
awk '{print $2}' > logistic.snplist

module unload bcftools/1.8
module load bcftools/1.4

bcftools view \
--threads ${NSLOTS} \
-i 'ID=@logistic.snplist' ../variant_calling/all_imputed.vcf.gz | \
bgzip > ../variant_calling/all_imputed_sign.vcf.gz

bcftools index --threads ${NSLOTS} -t ../variant_calling/all_imputed_sign.vcf.gz

module unload bcftools/1.4
module load bcftools/1.8

bcftools merge \
--threads ${NSLOTS}  \
../variant_calling/all_imputed_sign.vcf.gz \
1000G_QC_sel.vcf.gz \
-Ov | \
bgzip > 1000G_all_filteredQC.vcf.gz

bcftools index --threads ${NSLOTS} -t 1000G_all_filteredQC.vcf.gz

rm -f logistic.snplist

# Do a PCA.
bcftools query \
-f '%CHROM\t%POS[\t%GT]\n' \
1000G_all_filteredQC.vcf.gz \
> snp_matrix_all.txt

bcftools query \
-l 1000G_all_filteredQC.vcf.gz \
> samples_all.txt

# calculate the PCs to use as covariates.
plink \
--vcf 1000G_all_filteredQC.vcf.gz \
--allow-extra-chr \
--double-id \
--allow-no-sex \
-no-parents \
--make-bed \
--recode \
--out target_data_sign

plink \
--bfile \
target_data_sign \
--pca 2 \
--out target_data_sign


# Create a file containing the population stratification data for the cancer and healthy samples.
R -e '

samples       <- read.table("samples.txt", header = F)
fam           <- read.table("target_data_sign.fam",
                            header = F)
colnames(fam) <- c("FID","IID", "V3",
                   "V4", "sex", "Phenotype")
fam$sex       <- 1
fam$Phenotype <- ifelse(fam$FID %in% samples$V1, 1, 2)

write.table(fam[, c("FID","IID", "V3", "V4",
                    "sex", "Phenotype")],
            "target_data_sign.fam",
            quote     = F,
            row.names = F,
            col.names = F)

write.table(fam[, c("FID","IID", "Phenotype")],
            "target_data_sign.pheno",
            quote     = F,
            row.names = F,
            col.names = T)

q(save = "no")

'

# Update the .ped file adding phenotype and sex info.
awk '!($1=$2=$3=$4=$5=$6="")' target_data_sign.ped > temp
rm -f target_data_sign.ped
paste target_data_sign.fam temp > target_data_sign.ped

rm -f temp

module unload bcftools/1.8
module load bcftools/1.4

# Subset for the strongest SNPs only.
# snps_selection.snplist is a file containing 12 SNPs found in the clinvar database, which has the strongest power for the PRS.
bcftools view \
--threads ${NSLOTS} \
-i 'ID=@../variant_calling/snp_for_PRS.txt' \
1000G_all_filteredQC.vcf.gz \
-Ov | \
bgzip > snps_selection.vcf.gz

bcftools index --threads ${NSLOTS} -t snps_selection.vcf.gz

plink \
--vcf snps_selection.vcf.gz \
--allow-extra-chr \
--double-id \
--allow-no-sex \
-no-parents \
--make-bed \
--recode \
--out target_data_selection

cp target_data_sign.fam target_data_selection.fam
cp target_data_sign.pheno target_data_selection.pheno

# Update the .ped file adding info on phenotype and sex.
awk '!($1=$2=$3=$4=$5=$6="")' target_data_selection.ped > temp
rm -f target_data_selection.ped
paste target_data_selection.fam temp > target_data_selection.ped
rm -f temp*

## Process Prostate cancer samples.
mkdir -p ../controlvscases
cd ../controlvscases

cp ../variant_calling/all_imputed_sign.vcf.* .

bcftools view \
--threads ${NSLOTS} \
-i 'ID=@../variant_calling/snp_for_PRS.txt' \
all_imputed_sign.vcf.gz \
-Ov | \
bgzip > snps_selection.vcf.gz

bcftools index --threads ${NSLOTS} -t snps_selection.vcf.gz

plink \
--vcf all_imputed_sign.vcf.gz \
--allow-extra-chr \
--double-id \
--allow-no-sex \
-no-parents \
--make-bed \
--recode \
--out target_data_sign

plink \
--vcf snps_selection.vcf.gz \
--allow-extra-chr \
--double-id \
--allow-no-sex \
-no-parents \
--make-bed \
--recode \
--out target_data_selection


# Get the population stratification data to use as covariate.
plink \
--bfile target_data_sign \
--pca 2 \
--out target_data_sign

rm -f target_data_sign.eigenval


R -e "

master   <- read.csv('../TAPG_TURP_DATA.csv')
master   <- subset(master, !is.na(Age) & LABNO
                   !='' & !duplicated(LABNO))
data     <- read.table('target_data_sign.fam',
                       header = F)
data$V5  <- 1
data$id  <- gsub(pattern = '(.+)|(_.+)',
                 x           = data$V1,
                 replacement = '')
master   <- subset(master, LABNO %in% data$id)
data$age     <- as.integer(master$Age[match(data$id, master$LABNO)])
data$Gleason <- master$Gleason[match(data$id, master$LABNO)]
data$PSA     <- master$PSA[match(data$id, master$LABNO)]
pc <- read.table('target_data_sign.eigenvec')
data <- merge(data, pc, by = c('V1','V2'))
colnames(data) <- c('V1', 'V2', 'V3',
                    'V4', 'V5', 'V6',
                    'id','age', 'Gleason',
                    'PSA', 'PC1', 'PC2')
data1    <- read.table('target_data_selection.fam',
                       header = F)
data1$V5 <- 1
control  <- read.table('../variant_calling/control_sample.txt')
data$V6  <- ifelse(data$V1 %in% control$V1, '1', '2')
data1$V6 <- ifelse(data1$V1 %in% control$V1, '1', '2')

write.table(data[, c('V1', 'V2', 'V3',
                     'V4', 'V5', 'V6')],
            'target_data_sign.fam',
            quote     = F,
            row.names = F,
            col.names = F)
write.table(data1,'target_data_selection.fam',
            quote     = F,
            row.names = F,
            col.names = F)
write.table(data[, c('V1', 'V2', 'age',
                     'Gleason', 'PSA',
                     'PC1', 'PC2')],
            file      = 'target_data_sign.cov',
            quote     = F,
            row.names = F,
            col.names = F,
            sep       = '\t')


q(save = 'no')

"

# Update the .ped file adding info on phenotype and sex.
awk '!($1=$2=$3=$4=$5=$6="")' target_data_sign.ped > temp
awk '!($1=$2=$3=$4=$5=$6="")' target_data_selection.ped > temp2
rm -f target_data*.ped
paste target_data_sign.fam temp > target_data_sign.ped
paste target_data_selection.fam temp2 > target_data_selection.ped

rm -f tem*
