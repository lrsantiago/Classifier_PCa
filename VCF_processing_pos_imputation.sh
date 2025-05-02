!/bin/sh

module load java/1.8.0_152-oracle
module load bedtools/1.8
mmodule load bcftools/1.8
module load plink/1.9-170906
module load anaconda3/2020.02
conda activate py36

NSLOTS=20

cd variant_calling

# Merge all imputed files into one and recalculate AF and MAF.
ls ../imputation/all_filtered_QC_imputedINFO*vcf.gz > temp
bcftools concat \
--threads ${NSLOTS} \
-f temp \
-Ov | \
bgzip > imputed1.vcf.gz

bcftools index --threads ${NSLOTS} -t imputed1.vcf.gz

# Re-calculate allele frequency and compute Impute2-like INFO score.
bcftools +fill-tags imputed1.vcf.gz -Ov -o imputed2.vcf -- -t AF,MAF

bgzip imputed2.vcf

bcftools index --threads ${NSLOTS} -t imputed2.vcf.gz

bcftools view \
--threads ${NSLOTS} \
-e 'INFO/INFO<0.7' \
imputed2.vcf.gz \
-Ov | \
bgzip > imputed3.vcf.gz

bcftools index --threads ${NSLOTS} -t imputed3.vcf.gz

# Retrieve the variants that were not imputed.
bcftools view \
--threads ${NSLOTS} \
-H imputed3.vcf.gz | \
awk '{print $1}' | \
sort | \
uniq > chromosome.temp

bcftools view \
--threads  ${NSLOTS} \
-H all_filtered_QC.vcf.gz | \
awk '{print $1}' | \
sort | \
uniq > chromosome2.temp

diff chromosome.temp chromosome2.temp | \
grep '>' | \
sed 's|> ||' \
> sel_chromosomes.temp

bcftools view \
--threads ${NSLOTS} \
-r chr11,chr14,chr15,chr16,chr18,chr20,chr21,chr22,chr4,chr8,chr9,chrX \
all_filtered_QC.vcf.gz \
-Ov | \
bgzip > other_snps.vcf.gz

bcftools index --threads ${NSLOTS} -t other_snps.vcf.gz

bcftools concat \
--threads ${NSLOTS} \
imputed3.vcf.gz \
other_snps.vcf.gz -Ov | \
bgzip > all_imputed.vcf.gz

bcftools index --threads ${NSLOTS} -t all_imputed.vcf.gz

rm -f *emp
rm -f other_snps.vcf*
rm -f imputed*

# Convert vcfs to plink.
plink \
--vcf all_imputed.vcf.gz \
--allow-extra-chr \
--double-id \
--allow-no-sex \
-no-parents \
--make-bed \
--recode \
--out all_imputed

# Assing phenotype info (cases = 2, controls = 1) to the .fam file and covariates info in the .cov file.

R -e "

master      <- read.csv('../TAPG_TURP_DATA.csv')
master      <- subset(master, !is.na(Age) &
                        LABNO !='' &
                        !duplicated(LABNO))
fam         <- read.table('all_imputed.fam')
fam$id      <- gsub(pattern = '(.+/)|(_.+)',
                    x           = fam$V1,
                    replacement = '')
master      <- subset(master, LABNO %in% fam$id)
fam$age     <- as.integer(master$Age[match(fam$id,
                                           master$LABNO)])
fam$Gleason <- master$Gleason[match(fam$id,
                                    master$LABNO)]
fam$PSA     <- master$PSA[match(fam$id,
                                master$LABNO)]
fam$V5      <- rep(1, nrow(fam))
fam$V6      <- master$pheno[match(fam$id,
                                  master$LABNO)]

write.table(fam[, c('V1', 'V2', 'V3',
                    'V4', 'V5', 'V6')],
            file      = 'all_imputed.fam',
            quote     = F,
            row.names = F,
            col.names = F,
            sep       = '\t')

write.table(fam[, c('V1', 'V2', 'age')],
            file      = 'all_imputed.cov',
            quote     = F,
            row.names = F,
            col.names = F,
            sep       = '\t')


q(save = 'no')

"

# Update the .ped file adding phenotype and sex info.
awk '!($1=$2=$3=$4=$5=$6="")' all_imputed.ped > temp
rm -f all_imputed.ped
paste all_imputed.fam temp > all_imputed.ped

rm -f temp

# Remove highly correlated SNPs and generate heterozygosity file.
plink \
--bfile all_imputed \
--fam all_imputed.fam \
--allow-extra-chr \
--hwe 1e-6 \
--maf 0.01 \
--write-snplist \
--out all_imputed

plink \
--bfile all_imputed \
--allow-extra-chr \
--keep all_imputed.fam \
--extract all_imputed.snplist \
--indep-pairwise 200 50 0.25 \
--out all_imputed

plink \
--bfile all_imputed \
--allow-extra-chr \
--extract all_imputed.prune.in \
--keep all_imputed.fam \
--het \
--out all_imputed

# Remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean.
R -e "

dat2   <- read.table('all_imputed.het', header = T)
m2     <- mean(dat2$F)
s2     <- sd(dat2$F)
valid2 <- subset(dat2, F <= m2+3*s2 & F >= m2-3*s2)

write.table(valid2[,c(1, 2)],
            'all_imputed.valid.sample',
            quote     = F,
            row.names = F)

q(save = 'no')

"

plink \
--bfile all_imputed \
--allow-extra-chr \
--keep all_imputed.valid.sample \
--make-bed \
--extract all_imputed.snplist \
--out all_imputedQC

R -e "

fam     <- read.table('all_imputedQC.fam')
cov     <- read.table('all_imputed.cov')
cov_sub <- subset(cov, V1 %in% fam$V1)

write.table(cov_sub,
            file      = 'all_imputedQC.cov',
            quote     = F,
            row.names = F,
            col.names = F,
            sep       = '\t')

q(save = 'no')

"

# Association test between cases and controls with plink (logistic).
plink \
--bfile all_imputedQC \
--fam all_imputedQC.fam \
--allow-extra-chr \
--logistic hide-covar \
--ci 0.95 \
--covar all_imputedQC.cov \
--out all_imputed

# Get the allele frequency and MAF.
plink --bfile all_imputedQC --allow-extra-chr --freq --out all_imputed

# keep only SNPs listed in the all_imputed.assoc.logistic
module unload bcftools/1.8
module load bcftools/1.4

awk '{print $2}' all_imputed.assoc.logistic | sed 1d > snps.txt
cat all_imputedQC.fam | awk '{print $1}' > samples.txt

bcftools view \
--threads ${NSLOTS} \
-S samples.txt \
all_imputed.vcf.gz \
-Ov | \
bcftools view \
--threads ${NSLOTS} \
-i 'ID=@snps.txt'  - | \
bgzip > all_imputed2.vcf.gz

bcftools index --threads ${NSLOTS} -t all_imputed2.vcf.gz

rm -f snps.txt
rm -f samples.txt

module unload bcftools/1.4
module load bcftools/1.8

# Rename chromosomes according to GRCh38 version.
bcftools annotate \
--threads ${NSLOTS} \
--rename-chrs chrNames_chrHg38.txt \
all_imputed2.vcf.gz \
-Ov | \
bgzip > all_imputed3.vcf.gz

bcftools index --threads ${NSLOTS} -t all_imputed3.vcf.gz

cd ../snpEff

# Run snpEff
java -Xmx40g -jar snpEff.jar ann \
-lof -s all_imputed.html \
-csvStats all_imputed.csv \
-v GRCh38.p13 \
../variant_calling/all_imputed3.vcf.gz > \
../variant_calling/all_imputed4.vcf

bgzip ../variant_calling/all_imputed4.vcf

bcftools index --threads ${NSLOTS} -t ../variant_calling/all_imputed4.vcf.gz

bcftools stats \
--threads ${NSLOTS} \
-F ../reference/GCF_000001405.39_GRCh38.p13_genomic.fa \
-s - \
../variant_calling/all_imputed4.vcf.gz \
> ../variant_calling/summary_stats/all_imputed.vcf.gz.stats

rtg-tools-3.11/rtg vcfstats \
../variant_calling/all_imputed4.vcf.gz \
> ../variant_calling/summary_stats/all_imputed.vcf.gz.sstats

mkdir -p ../variant_calling/summary_stats/plots/all_imputed

plot-vcfstats \
-p ../variant_calling/summary_stats/plots/all_imputed \
../variant_calling/summary_stats/all_imputed.vcf.gz.stats

cd ../variant_calling/summary_stats/plots/all_imputed

python3.6 plot.py

cd
cd variant_calling


# Rename chromosomes back to what it was.
bcftools annotate \
--threads ${NSLOTS}  \
--rename-chrs NC_chr_names.txt \
all_imputed4.vcf.gz \
-Ov | \
bgzip > all_imputed.vcf.gz

bcftools index --threads ${NSLOTS} -t all_imputed.vcf.gz

rm -f all_imputed*.vcf*

cd ../snpEff

# Create a vcf with the rsIDs.
java -Xmx40g -jar SnpSift.jar annotate \
-id cos_clin_dbsnp.vcf.gz \
../variant_calling/all_imputed.vcf.gz \
> ../variant_calling/all_imputedrsID.vcf

bgzip ../variant_calling/all_imputedrsID.vcf

bcftools index --threads ${NSLOTS} -t ../variant_calling/all_imputedrsID.vcf.gz

bcftools view \
--threads ${NSLOTS} \
-H ../variant_calling/all_imputedrsID.vcf.gz | \
awk \
-v OFS='\t' '{print $1, $2, $3}' \
> ../variant_calling/all_imputedrsID.bed

cd ../variant_calling

# Create a bed file from the all_imputed vcf and add columns with OR and P.
bcftools view \
--threads ${NSLOTS} \
-H all_imputed.vcf.gz | \
awk -v OFS='\t' '{print $1, $2, $2, $3, $4, $5, $6, $8}' \
> all_impute.bed

R -e "

bed     <- read.table('all_impute.bed')
log     <- read.table('all_imputed.assoc.logistic',
                      header = T)
log$CHR <- paste0('chr', log$CHR)
log$CHR <- gsub(pattern = 'chr23',
                x           = log$CHR,
                replacement = 'chrX')
bed     <- bed[order(bed$V1, bed$V2), ]
log     <- log[order(log$CHR, log$BP), ]

bed[, c('V9', 'V10',
        'V11', 'V12')] <- log[, c('OR', 'L95',
                                  'U95', 'P')]


write.table(bed,
            'all_imputed_complete.bed',
            quote     = F,
            row.names = F,
            col.names = F,
            sep       = '\t')


q(save = 'no')

"

rm -f all_impute.bed

# Get the sample names.
bcftools query -l all_imputed.vcf.gz > sample_names.txt

# Get the fixed differences in cases and the variants that are different between cases and controls.
Rscript scripts/vcf_fixed_snp.R

bgzip fixed_case.vcf

bgzip fixed_control.vcf

bgzip cases_diffs.vcf

bcftools index --threads ${NSLOTS} -t fixed_case.vcf.gz

bcftools index --threads ${NSLOTS} -t fixed_control.vcf.gz

bcftools index --threads ${NSLOTS} -t cases_diffs.vcf.gz

bcftools view \
--threads ${NSLOTS} \
-H cases_diffs.vcf.gz | \
awk \
-v OFS='\t' '{print $1, $2, $2, $3, $4, $5, $6, $7, $8}' \
> cases_diffs.bed

bedtools intersect \
-wo \
-a cases_diffs.bed \
-b all_imputed_complete.bed \
> diff_match.bed

cat diff_match.bed | \
awk \
-v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $15}' \
> variant_diff.bed

bcftools view \
--threads ${NSLOTS} \
-H  cases_diffs.vcf.gz | \
awk '{print $3}' > cases_diffs.txt

bcftools view \
--threads ${NSLOTS} \
-H fixed_case.vcf.gz | \
awk '{print $3}' > fixed_case.txt

bcftools view \
--threads ${NSLOTS} \
-H fixed_control.vcf.gz | \
awk '{print $3}' > fixed_control.txt

rm -f cases_diffs.bed
rm -f diff_match.bed


# Get only variant private to either control or cases.
bcftools view \
--threads ${NSLOTS} \
-X \
-S case_sample.txt \
all_imputed.vcf.gz \
-Ov | \
bgzip > case.privat.vcf.gz

bcftools view \
--threads ${NSLOTS} \
-X \
-S control_sample.txt \
all_imputed.vcf.gz \
-Ov | \
bgzip > control.privat.vcf.gz

bcftools index --threads ${NSLOTS} -t case.privat.vcf.gz

bcftools index --threads ${NSLOTS} -t control.privat.vcf.gz

# Get the variants present in A but not in B (0000.vcf) and present in B but not in A (0001.vcf)
bcftools isec \
control.privat.vcf.gz \
case.privat.vcf.gz \
  -p \
  -n \
  -1 \
  -c all

mv 0000.vcf control_privateSNPs.vcf

bgzip control_privateSNPs.vcf

bcftools index -t control_privateSNPs.vcf.gz

bcftools view \
-H control_privateSNPs.vcf.gz | \
awk '{print $3}' > control_privateSNPs.txt

mv 0001.vcf case_privateSNPs.vcf

bgzip case_privateSNPs.vcf

bcftools index --threads ${NSLOTS} -t case_privateSNPs.vcf.gz

bcftools view \
--threads ${NSLOTS} \
-H case_privateSNPs.vcf.gz | \
awk '{print $3}' > case_privateSNPs.txt

rm -f *privat.vcf*

#  Create a file with chr, position and genotypes to use afterwards.
bcftools query all_imputed.vcf.gz -f '%CHROM\t%POS[\t%GT]\n' > snp_matrix.txt

# create DbSnp, cosmic and clinvar bed files.
bcftools view \
--threads ${NSLOTS} \
-H ../snpEff/DbSnp-All.vcf.gz | \
awk \
-v OFS='\t' '{print $1, $2, $4, $8}' \
> DbSnp.temp

cat DbSnp.temp | \
awk '{print $4}'| \
sed 's|.*GENEINFO=||' | \
sed 's|:.*||' \
> DbSnp.temp2

cat DbSnp.temp | \
awk '{print $1"\t"$2"\t"$3}' \
> DbSnp.temp1

paste DbSnp.temp1 DbSnp.temp2 > DbSnp.bed

rm -f D*.tem*

# Get the SNPs that overlap between all samples and cosmic and clinvar.
bcftools view \
--threads ${NSLOTS} \
-H ../snpEff/clinvar_chr.vcf.gz | \
awk \
-v OFS='\t' '{print $1, $2, $2, $3, $4, $5, $6, $7, $8}' \
> clinvar.bed

bcftools view \
--threads ${NSLOTS} \
-H ../snpEff/rs_COSMIC.vcf.gz | \
awk \
-v OFS='\t' '{print $1, $2, $2, $3, $4, $5, $6, $7, $8}' \
> cosmic.bed

bcftools view \
--threads ${NSLOTS} \
-H ../snpEff/clinvar_chr.vcf.gz | \
awk \
-v OFS='\t' '{print $1, $2, $4, $5}' \
> clinvar.txt

# Get the SNPs that overlap between all samples and cosmic and clinvar.
bedtools intersect \
-wo \
-a all_imputed_complete.bed \
-b clinvar.bed \
> clinvar_imputed.bed

bedtools intersect \
-wo \
-a all_imputed_complete.bed \
-b cosmic.bed \
> cosmic_imputed.bed
