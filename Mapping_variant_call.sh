#!/bin/sh

module load java/1.8.0_152-oracle
module load bwa/0.7.17
module load samtools/1.8
module load vcftools/0.1.16
module load bedtools/1.8
mmodule load bcftools/1.8
module load anaconda3/2020.02
conda activate r_env_new # Environment with R v4.0.0 where all required packages were installed.

NSLOTS=20

mkdir -p variant_calling
mkdir -p bwa_indexes/
mkdir -p reference/
mkdir -p variant_calling/summary_stats
mkdir -p variant_calling/plots
mkdir -p variant_calling/raw_vcfs

# Index the genome as required for variant calling.
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -P reference/
gunzip reference/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
mv reference/GCF_000001405.39_GRCh38.p13_genomic.fna reference/GCF_000001405.39_GRCh38.p13_genomic.fa
samtools faidx reference/GCF_000001405.39_GRCh38.p13_genomic.fa

# Check if bam files are single or paired ends.
samtools view -@ ${NSLOTS} \
-c \
-f 1 \
merged/T_sorted_merged_rawlib.bam

# Generate bwa indexes.
bwa index \
-p bwa_indexes/GRCh38.p13 \
-a bwtsw \
reference/GCF_000001405.39_GRCh38.p13_genomic.fa

# Align all bam files to the reference, convert aligned sam to bam, sort, then select only the uniquely mapped and perform the variant call.
for samples in $(ls merged/*bam | sed 's|.*/||' | sed 's|_.*||')
do
  bwa aln \
  -t ${NSLOTS} \
  bwa_indexes/GRCh38.p13 \
  -b0 merged/${samples}*.bam \
  > merged/${samples}_mapped.sai

  bwa samse \
  -f merged/${samples}_mapped.sam \
  bwa_indexes/GRCh38.p13 \
  merged/${samples}_mapped.sai \
  merged/${samples}*.bam

  samtools view \
  -@ ${NSLOTS} \
  -bS merged/${samples}_mapped.sam \
  > merged/${samples}_mapped.bam

  samtools sort \
  -@ ${NSLOTS} \
  merged/${samples}_mapped.bam \
  -O BAM \
  -o merged/${samples}_mapped_sorted.bam

  samtools view \
  -@ ${NSLOTS} \
  -q 10 -b \
  merged/${samples}_mapped_sorted.bam \
  > merged/${samples}_mapped_sorted_uniq.bam

  samtools index \
  -@ ${NSLOTS} \
  merged/${samples}_mapped_sorted_uniq.bam

  samtools mpileup \
  -C50 \
  -ugp \
  -f reference/GCF_000001405.39_GRCh38.p13_genomic.fa \
  merged/${samples}_mapped_sorted_uniq.bam \
  -o variant_calling/raw_vcfs/${samples}_raw.bcf

  bcftools call \
  --threads ${NSLOTS} \
  -mv variant_calling/raw_vcfs/${samples}_raw.bcf \
  -O v \
  -o variant_calling/raw_vcfs/${samples}_raw.vcf

  bgzip variant_calling/raw_vcfs/${samples}_raw.vcf

  bcftools index \
  --threads ${NSLOTS} \
  -t variant_calling/raw_vcfs/${samples}_raw.vcf.gz

  bcftools view \
  --threads ${NSLOTS} \
  variant_calling/raw_vcfs/${samples}_raw.vcf.gz \
  -r NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10,NC_000020.11,NC_000021.9,NC_000022.11,NC_000023.11,NC_000024.10 \
  -Ov | \
  bgzip > variant_calling/raw_vcfs/${samples}_chr.vcf

  bcftools index \
  --threads ${NSLOTS} \
  -t variant_calling/raw_vcfs/${samples}_chr.vcf.gz

  bcftools stats \
  --threads ${NSLOTS} \
  -F reference/GCF_000001405.39_GRCh38.p13_genomic.fa \
  -s - variant_calling/raw_vcfs/${samples}_chr.vcf.gz \
  > variant_calling/summary_stats/${samples}_chr.vcf.gz.stats

  mkdir -p variant_calling/plots/${samples}_unfiltered
  plot-vcfstats \
  -p variant_calling/plots/${samples}_unfiltered \
  variant_calling/summary_stats/${samples}_chr.vcf.gz.stats

  cd variant_calling/plots/${samples}_unfiltered

  # Add double-slash in the line 73, then run the below command afterwards.
  python3.6 plot.py

  cd
  rtg-tools-3.11/rtg vcfstats \
  variant_calling/raw_vcfs/${samples}_chr.vcf.gz \
  > variant_calling/summary_stats/${samples}_chr_unfiltered.vcf.sstats

  echo $samples >> file_list.txt
done

rm -f merged/*.sa*
rm -f merged/*mapped.bam

# Create a file with the summary stats for every unfiltered vcf sample.
R -e "

files <- read.table('file_list.txt')

snps <- list()
for (i in 1:nrow(files)) {
  snps[[i]] <- read.delim(paste('variant_calling/summary_stats/',
                                files$V1[i], '_chr_unfiltered.vcf.sstats',
                                sep = ''),
                          header = F)
  temp        <- snps[[i]]
  files$V2[i] <- as.numeric(gsub(pattern     = '.+: ',
                                 x           = temp$V1[4],
                                 replacement =  ''))
  files$V3[i] <- as.numeric(gsub(pattern     = '.+: ',
                                 x           = temp$V1[6],
                                 replacement =  ''))
  files$V4[i] <- as.numeric(gsub(pattern     = '.+: ',
                                 x           = temp$V1[7],
                                 replacement =  ''))
}

colnames(files) <- c('Sample_ID', 'SNPs',
                     'Insertions', 'Deletions')

write.table(files,
            file      = 'variant_calling/unfiltered_variants.txt',
            quote     = F,
            row.names = F,
            col.names = T,
            sep       = '\t')


q(save = 'no')

"

rm -f file_list.txt

# Merge the raw vcfs.
ls variant_calling/raw_vcfs/*_raw.vcf.gz > raw_vcf_path.txt

bcftools merge \
--threads ${NSLOTS} \
-l raw_vcf_path.txt \
-Ov | \
bgzip > variant_calling/all_raw.vcf.gz

bcftools index --threads ${NSLOTS} -t variant_calling/all_raw.vcf.gz

# Merge all unfiltered vcfs and add chr:position:REF:ALT as the SNP ID.
ls variant_calling/raw_vcfs/*_chr.vcf.gz > unfiltered_vcf_path.txt

# Rename chr names by chr{number} and remove chrY.
bcftools merge \
--threads ${NSLOTS} \
-l unfiltered_vcf_path.txt \
-Ov | \
bgzip > variant_calling/all_unfiltered_raw.vcf.gz

bcftools index \
--threads ${NSLOTS} \
-t variant_calling/all_unfiltered_raw.vcf.gz

bcftools annotate \
--threads ${NSLOTS} \
--set-id '%CHROM\:%POS\:%REF\:%ALT' \
variant_calling/all_unfiltered_raw.vcf.gz \
-Ov | \
bgzip > variant_calling/all_unfiltered_1.vcf.gz

bcftools index --threads ${NSLOTS} -t variant_calling/all_unfiltered_1.vcf.gz

# Keep only the chromosomes.
bcftools view \
--threads ${NSLOTS} \
-r NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10,NC_000020.11,NC_000021.9,NC_000022.11,NC_000023.11 \
variant_calling/all_unfiltered_1.vcf.gz \
-Ov | \
bgzip > variant_calling/all_unfiltered.vcf.gz

bcftools index --threads ${NSLOTS} -t variant_calling/all_unfiltered.vcf.gz

rm -f raw_vcf_path.txt
rm -f unfiltered_vcf_path.txt
rm -f variant_calling/all_unfiltered_*


## Check the SNP effects of agressive prostate cancer (cases) on Human genes with snpEff.
# Download the snpEff database for the human genome GRCh38.p13 (hg38).
wget wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm -f snpEff_latest_core.zip

cd snpEff

# Building a database for the GRCh38.p13 version (NCBI Annotation Release:	109)
mkdir -p data
mkdir -p data/GRCh38.p13

cd data/GRCh38.p13

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz
mv GCF_000001405.39_GRCh38.p13_genomic.fna sequences.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_protein.faa.gz
gunzip GCF_000001405.39_GRCh38.p13_protein.faa.gz
mv GCF_000001405.39_GRCh38.p13_protein.faa protein.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
gunzip GCF_000001405.39_GRCh38.p13_genomic.gff.gz
mv GCF_000001405.39_GRCh38.p13_genomic.gff genes.gff

cd ../..

# Retrieve the DbSnp IDs.
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz

# Replace chromosome numbers by chr names.
bcftools annotate \
--threads ${NSLOTS} \
--rename-chrs number_chrnames.txt \
00-All.vcf.gz \
-O v | \
bgzip > DbSnp-All.vcf.gz

bcftools index --threads ${NSLOTS} -t DbSnp-All.vcf.gz

rm -f 00-All.vcf.gz

# Download the vcfs from clinvar, CGC and cosmic for comparison.
#cosmic
wget ftp://ftp.ncbi.nih.gov/snp/others/rs_COSMIC.vcf.gz
mv rs_COSMIC.vcf.gz rs_COSMIC0.vcf.gz
bcftools index --threads ${NSLOTS} -t rs_COSMIC0.vcf.gz
bcftools view \
--threads ${NSLOTS} \
-r NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10,NC_000020.11,NC_000021.9,NC_000022.11,NC_000023.11 \
clinvar.vcf.gz \
rs_COSMIC0.vcf.gz \
-Ov | \
bgzip > rs_COSMIC1.vcf.gz

bcftools index --threads ${NSLOTS} -t rs_COSMIC1.vcf.gz

bcftools annotate \
--threads ${NSLOTS} \
--rename-chrs NC_chr_names.txt \
rs_COSMIC1.vcf.gz \
-Ov | \
bgzip > rs_COSMIC.vcf.gz

bcftools index --threads ${NSLOTS} -t rs_COSMIC.vcf.gz

rm -f rs_COSMIC0.vcf*
rm -f rs_COSMIC1.vcf*

#CGC
wget ftp://ftp.ncbi.nih.gov/snp/others/snp_icgc.vcf.gz
bcftools index --threads ${NSLOTS} -t snp_icgc.vcf.gz

# clinVar. remove chr Y and MT and replace chr numbers for chr names.
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

chrs=$(echo {1..22} X | tr ' ' ',')

bcftools view \
--threads ${NSLOTS} \
-t $chrs clinvar.vcf.gz \
-Ov | \
bgzip > clinvarx.vcf.gz

bcftools index --threads ${NSLOTS} -t clinvarx.vcf.gz

bcftools annotate \
--threads ${NSLOTS} \
--rename-chrs number_chrnames.txt \
clinvarx.vcf.gz \
-Ov | \
bgzip > clinvar_chr.vcf.gz

bcftools index --threads ${NSLOTS} -t clinvar_chr.vcf.gz

rm -f clinvar.vcf.g*
rm -f clinvarx.vcf.g*

# Merge all data sets into one to get the highest number of SNP IDs
bcftools merge \
--threads ${NSLOTS} \
--force-samples DbSnp-All.vcf.gz \
clinvar_chr.vcf.gz \
rs_COSMIC.vcf.gz \
-Ov | \
bgzip > cos_clin_dbsnp.vcf.gz

bcftools index --threads ${NSLOTS} -t cos_clin_dbsnp.vcf.gz

bcftools view \
--threads ${NSLOTS} \
-H cos_clin_dbsnp.vcf.gz | \
awk \
-v OFS='\t' '{print $1, $2, $2, $3, $4, $5}' \
> cos_clin_dbsnp.bed

# Build the databases.
java -Xmx40g -jar snpEff.jar build -gff3 -v GRCh38.p13

cd

# Remove variants with quality below 20, depth above 100 and below 10, Read Position Bias (RPB, bigger is better), allele count below 3, minor allele count below 1 and minor allele frequency below 0.01 (i.e. removing rare variants at 1% cutoff).
bcftools filter \
--threads ${NSLOTS} \
-g 3 \
-G 10 \
-s PASS \
-e '%QUAL<20 || DP>100 || DP<10 || RPB<0.1 || AC<3' \
variant_calling/all_unfiltered.vcf.gz  \
-Ov | \
bgzip > variant_calling/all_filtered_temp.vcf.gz

bcftools index --threads ${NSLOTS} -t variant_calling/all_filtered_temp.vcf.gz

bcftools +fill-tags \
--threads ${NSLOTS} \
variant_calling/all_filtered_temp.vcf.gz \
-Ov \
-o variant_calling/all_filtered_1.vcf -- \
-t AF,MAF

bgzip variant_calling/all_filtered_1.vcf

bcftools index --threads ${NSLOTS} -t variant_calling/all_filtered_1.vcf.gz

# Remove duplicated variants.
bcftools norm \
--threads ${NSLOTS} \
-d none \
variant_calling/all_filtered_1.vcf.gz \
-Ov | \
bgzip > variant_calling/all_filtered_2.vcf.gz

bcftools index --threads ${NSLOTS} -t variant_calling/all_filtered_2.vcf.gz

cd snpEff

# Run snpEff.
java -Xmx160g -jar snpEff.jar ann \
-lof \
-s ../variant_calling/all_unfiltered_2.html \
-csvStats ../variant_calling/all_unfiltered_2.csv \
-v GRCh38.p13 \
../variant_calling/all_filtered_2.vcf.gz \
> ../variant_calling/all_filtered.vcf

bgzip ../variant_calling/all_filtered.vcf

bcftools index --threads ${NSLOTS} -t ../variant_calling/all_filtered.vcf.gz

bcftools stats \
--threads ${NSLOTS} \
-F ../reference/GCF_000001405.39_GRCh38.p13_genomic.fa \
-s - \
../variant_calling/all_filtered.vcf.gz \
> ../variant_calling/summary_stats/all_filtered.vcf.gz.stats

../rtg-tools-3.11/rtg vcfstats \
../variant_calling/all_filtered.vcf.gz \
> ../variant_calling/summary_stats/all_filtered.vcf.gz.sstats

rm -f ../variant_calling/all_filtered_*

cd

# Remove all sites with missing genotype above 50%, minor allele frequency below 1% and keep only biallelic sites.
vcftools \
--gzvcf variant_calling/all_filtered.vcf.gz \
--maf 0.01 \
--mac 1 \
--max-missing 0.5 \
--recode \
--recode-INFO-all \
--out variant_calling/all_filtered

mv variant_calling/all_filtered.recode.vcf variant_calling/all_filtered_QC1.vcf

bgzip variant_calling/all_filtered_QC1.vcf

bcftools index --threads ${NSLOTS} -t variant_calling/all_filtered_QC1.vcf.gz

bcftools view \
--threads ${NSLOTS} \
-m 2 \
-M 2 \
variant_calling/all_filtered_QC1.vcf.gz \
-Ov | \
bgzip > variant_calling/all_filtered_QC0.vcf.gz

bcftools index --threads ${NSLOTS} -t variant_calling/all_filtered_QC0.vcf.gz

bcftools stats \
--threads ${NSLOTS} \
-F reference/GCF_000001405.39_GRCh38.p13_genomic.fa \
-s - \
variant_calling/all_filtered_QC0.vcf.gz \
> variant_calling/summary_stats/all_filtered_QC.vcf.gz.stats

rtg-tools-3.11/rtg vcfstats \
variant_calling/all_filtered_QC0.vcf.gz \
> variant_calling/summary_stats/all_filtered_QC.vcf.gz.sstats

mkdir -p variant_calling/plots/all_filtered_QC

plot-vcfstats \
-p variant_calling/plots/all_filtered_QC \
variant_calling/summary_stats/all_filtered_QC.vcf.gz.stats

cd variant_calling/plots/all_filtered_QC

# add double-slash in the line 73.
python3.6 plot.py

cd

# Keep only the 163 samples and remane all chromosomes to chr{number,X}
bcftools view \
--threads ${NSLOTS} \
-S gleason_samples.txt \
variant_calling/all_filtered_QC0.vcf.gz \
-Ov | \
bgzip > variant_calling/all_filtered_QC00.vcf.gz

bcftools index --threads ${NSLOTS} -t variant_calling/all_filtered_QC00.vcf.gz

bcftools annotate \
--threads ${NSLOTS} \
--rename-chrs chrID_chrNames.txt \
variant_calling/all_filtered_QC00.vcf.gz \
-Ov | \
bcftools annotate \
--threads ${NSLOTS} \
--set-id '%CHROM\:%POS\:%REF\:%ALT' \
-Ov | \
bgzip > variant_calling/all_filtered_QC.vcf.gz

bcftools index --threads ${NSLOTS} -t variant_calling/all_filtered_QC.vcf.gz

rm -f variant_calling/all_filtered_QC0*
rm -f variant_calling/all_filtered_QC1*
