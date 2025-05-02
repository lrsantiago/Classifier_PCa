#!/bin/sh

module load bcftools/1.4
module load plink/1.9-170906
module load anaconda3/2020.02
module load R/4.0.0
conda activate py36 #python 3.6

NSLOTS=20

mkdir -p imputation
cd imputation

# Download the genetic map files and the software for phasing with Eagle.
wget https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz
tar -xfvz Eagle_v2.4.1.tar.gz

wget https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz \
-P Eagle_v2.4.1/tables

# Split files by chromosomes.
for CHR in {1..23}
do
    zcat Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    | grep ^${CHR} | \
    sed '1ichr position COMBINED_rate(cM/Mb) Genetic_Map(cM)' \
    > Eagle_v2.4.1/eagle_chr${CHR}_b38.map
done

# Replace 23 by X.
cat Eagle_v2.4.1/eagle_chr23_b38.map | \
sed 's|^23|X|' > Eagle_v2.4.1/eagle_chrX_hg38_temp.map

rm -f Eagle_v2.4.1/eagle_chr23_b38.map

mv Eagle_v2.4.1/eagle_chrX_hg38_temp.map \
Eagle_v2.4.1/eagle_chr23_b38.map


# Download genetic map files for imputation with Beagle.
mkdir -p Beagle
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip -P Beagle
unzip Beagle/plink.GRCh38.map.zip

# Change chromosome notation from numbers to chr(number), ex: 12 by chr12 and X, Y, by chr23 and chr24.
mv Beagle/plink.chrX.GRCh38.map Beagle/plink.chr23.GRCh38.map

for CHR in {1..23}
do
    cat Beagle/plink.chr${CHR}.GRCh38.map \
    | sed 's|^|chr|' \
    > Beagle/beagle_chr${CHR}_b38.map
done

rm -f Beagle/plink.GRCh38*
rm -f Beagle/plink.chr*
rm -f Beagle/*temp.map
rm -f Beagle/README.txt

# Download Imputation reference panel (hg38).
mkdir -p ../GWAS_data

cd ../GWAS_data

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{{1..22},X}.filtered.shapeit2-duohmm-phased.vcf.{gz,gz.tbi}

# Checksum the files.
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/phased-manifest_July2021.tsv

cat phased-manifest_July2021.tsv | \
sort -k 1 \
> phased-manifest_July2021.tsv_v2.txt

rm -f phased-manifest_July2021.tsv

md5sum CCDG* | \
sort -k 2 \
> md5.files

diff -ys <(awk '{print $3}' phased-manifest_July2021.tsv_v2.txt) \
<(awk '{print $3}' md5.files)

# Data processing.
for CHR in $(ls CCDG*.gz | sed 's|.*/CCDG.*_||' | sed 's|\..*||')
do
    # Rename chromosomes according to the reference genome (GRCh38).
    bcftools annotate \
    --threads ${NSLOTS} \
    --rename-chrs chr_NC_names.txt CCDG*${CHR}.filtered*.gz \
    -Ov | bgzip > CCDG_${CHR}_GRCh38.genotypes.vcf.gz
    bcftools index \
    --threads ${NSLOTS} \
    -t CCDG_${CHR}_GRCh38.genotypes.vcf.gz
    # Remove rare variants (singletons and doubletons) by setting AC threshold.
    # Split multiallelic sites to biallelic records. Align the variants to reference genome.
    # Remove multiallelic records. Remove sites containing missing data.
    bcftools view \
    --threads ${NSLOTS} \
    -e 'INFO/AC<3 | INFO/AN-INFO/AC<3' CCDG_${CHR}_GRCh38.genotypes.vcf.gz \
    -Ov | \
    bcftools norm \
    --threads ${NSLOTS} \
    -m \
    -any \
    -Ov | \
    bcftools norm \
    --threads ${NSLOTS} \
    -f ../reference/GCF_000001405.39_GRCh38.p13_genomic.fa \
    -d none \
    -Ov | \
    bcftools view \
    --threads ${NSLOTS} \
    -m 2 \
    -M 2 \
    -Ov | \
    bcftools view \
    --threads ${NSLOTS} \
    -g ^miss \
    -Ov | \
    bgzip > 1000GP_${CHR}_temp.vcf.gz

    bcftools index \
    --threads ${NSLOTS} \
    -t 1000GP_${CHR}_temp.vcf.gz
    # Rename chromosomes to chr(chromosome number), ex: NC_00001 to chr1.
    bcftools annotate \
    --threads ${NSLOTS} \
    --rename-chrs NC_chr_names.txt 1000GP_${CHR}_temp.vcf.gz \
    -Ov | \
    bcftools annotate \
    --threads ${NSLOTS} \
    --set-id '%CHROM\:%POS\:%REF\:%ALT' \
    -Ov | \
    bgzip > 1000GP_${CHR}.vcf.gz

    bcftools index --threads ${NSLOTS} -t 1000GP_${CHR}.vcf.gz
done

rm -f CCDG_chr*
rm -f 1000*temp.vcf*

# Change the chrX ploidy to phased diploid.
# Firstly create a ploidy.txt file containing space-separated CHROM,FROM,TO,SEX,PLOIDY.
echo "chrX 1 156040895 M 2" > ploidy.txt

bcftools +fixploidy 1000GP_chrX.vcf.gz \
-Ov -- \
-p ploidy.txt | \
sed 's#0/0#0\|0#g;s#1/1#1\|1#g' | \
bcftools view -Ov | bgzip > 1000GP_chr23.vcf.gz

bcftools index --threads ${NSLOTS} -t 1000GP_chr23.vcf.gz

rm -f 1000GP_chrX*
rm -f ploidy.txt

# Keep only somatic mutations.
for CHR in $(ls 1000*.gz | sed 's|.*/1000GP_||' | sed 's|.v.*||')
do
  java -Xmx40g -jar SnpSift.jar annotate \
  -id cos_clin_dbsnp.vcf.gz 1000GP_${CHR}.vcf.gz > 1000GP_${CHR}rsid.vcf
  bgzip 1000GP_${CHR}rsid.vcf
  bcftools index --threads ${NSLOTS} -t 1000GP_${CHR}rsid.vcf.gz
  bcftools view \
  -H 1000GP_${CHR}rsid.vcf.gz | \
  awk '{print $3}' | \
  grep 'rs' > ids${CHR}.txt
  bcftools view \
  -i ID=@ids${CHR}.txt 1000GP_${CHR}rsid.vcf.gz \
  -Ov | \
  bgzip > 1000GP_${CHR}rs.vcf.gz
  bcftools index --threads ${NSLOTS} -t 1000GP_${CHR}rs.vcf.gz
done

rm -f ids*.txt
rm -f 1000GP*

# Rename the SNP IDs column with chr:pos:ref:alt and remove duplicates.
for CHR in $(ls 1000*rs*.gz | sed 's|.*/1000GP_||' | sed 's|rs.*||')
do
  bcftools annotate \
  --threads ${NSLOTS} \
  --set-id '%CHROM\:%POS\:%REF\:%ALT' 1000GP_${CHR}rs.vcf.gz \
  -Ov | \
  bgzip > 1000GP_${CHR}new.vcf.gz
  bcftools index --threads ${NSLOTS} -t 1000GP_${CHR}new.vcf.gz
  bcftools query \
  -f '%ID\n' 1000GP_${CHR}new.vcf.gz | \
  sort | \
  uniq -d > 1000GP_${CHR}.dup_id
    if [[ -s 1000GP_${CHR}.dup_id ]]; then
    	bcftools view \
      -e ID=@1000GP_${CHR}.dup_id 1000GP_${CHR}new.vcf.gz \
      -Ov | \
      bgzip > 1000GP_filtered_${CHR}.vcf.gz
    else
    	mv 1000GP_${CHR}new.vcf.gz  1000GP_filtered_${CHR}.vcf.gz
    fi
  bcftools index --threads ${NSLOTS} -t 1000GP_filtered_${CHR}.vcf.gz
done

rm -f 1000GP_chr*
rm -f 1000GP*.dup_id

# Create a binary reference panel with Beagle.
wget https://faculty.washington.edu/browning/beagle/beagle.29May21.d6d.jar -P ../imputation/Beagle
wget https://faculty.washington.edu/browning/beagle/bref3.29May21.d6d.jar -P ../imputation/Beagle

# Create a freq file with CHR,SNP,REF,ALT,AF
echo -e 'CHR\tSNP\tREF\tALT\tAF' > 1000GP_imputation_all.frq

for CHR in $(ls 1000*.gz | sed 's|.*1000.*_||' | sed 's|.v.*||')
do
    bcftools query \
    -f '%CHROM\t%ID\t%REF\t%ALT\t%INFO/AF\n' 1000GP_filtered_${CHR}.vcf.gz \
    >> 1000GP_imputation_all.frq
    zcat 1000GP_filtered_${CHR}.vcf.gz | \
    java -Xmx50g -jar ../imputation/Beagle/bref3.29May21.d6d.jar > 1000GP_filtered_${CHR}.bref3
done

# Create a file with sample IDs.
bcftools query -l 1000GP_filtered_chr1.vcf.gz > 1000GP_sample_IDs.txt

# Generate a frequency file.
echo -e 'CHR\tSNP\tREF\tALT\tAF' > all_filtered_QC.frq

# Append info from the VCF file with the significant SNPs to the allele frequency file.
cp ../variant_calling/all_filtered_QC.vcf.gz .
bcftools query -f '%CHROM\t%ID\t%REF\t%ALT\t%INFO/AF\n' all_filtered_QC.vcf.gz \
>> all_filtered_QC.frq

# Plot the graphs to compare the frequency between the target data and the panel data (1000G data) for imputation.
DATASET=all_filtered_QC
PANEL_FRQ=1000GP_imputation_all.frq

# The R script will create the exclude.txt file with the discordant alleles.
Rscript --no-save ../scripts/plot_AF.R ${DATASET}.frq ${DATASET} ${PANEL_FRQ} 0.2 5

# Remove the discordant alleles and the ones not present in the panel.
sort -V all_filtered_QC_exclude.txt > all_filtered_QC_exclusion1.txt
sort -V all_filtered_QC_nonpanel_exclude.txt > all_filtered_QC_exclusion2.txt

# Remove the variants.
bcftools view \
--threads ${NSLOTS} \
-e 'ID=@all_filtered_QC_exclusion1.txt' all_filtered_QC.vcf.gz \
-Ov | \
bgzip > all_filtered_QC_exclusion1.vcf.gz

bcftools index --threads ${NSLOTS} -t all_filtered_QC_exclusion1.vcf.gz

bcftools view \
--threads ${NSLOTS} \
-e 'ID=@all_filtered_QC_exclusion2.txt' all_filtered_QC_exclusion1.vcf.gz \
-Ov | \
bgzip > all_filtered_QC_for_phasing.vcf.gz

bcftools index --threads ${NSLOTS} -t all_filtered_QC_for_phasing.vcf.gz

rm -f all_filtered_QC_exclusion*

# Run phasing for each chromosome.
for CHR in {1..23}
do

    if [ "$CHR" == "23" ]; then
          chrname=chrX
    else
          chrname=chr$CHR
    fi

    ../imputation/Eagle_v2.4.1/eagle \
    --vcf all_filtered_QC_for_phasing.vcf.gz \
    --chrom ${chrname} \
    --geneticMapFile Eagle_v2.4.1/eagle_chr${CHR}_b38.map \
    --numThreads=${NSLOTS} \
    --Kpbwt=20000 \
    --outPrefix all_filtered_QC_for_imputation_chr${CHR}
done

# Run imputation for each chromosome
cd

for CHR in $(ls *for_im*.gz | sed 's|.*_||' | sed 's|\..*||')
do
    java -jar ../imputation/Beagle/beagle.29May21.d6d.jar \
    gt=all_filtered_QC_for_imputation_${CHR}.vcf.gz \
    ref=1000GP_filtered_${CHR}.bref3 \
    map=../imputation/Beagle/beagle_${CHR}_b38.map \
    out=all_filtered_QC_imputed_${CHR} \
    nthreads=${NSLOTS} \
    impute=true \
    gp=true
    # Compress it again with bgzip.
    gunzip all_filtered_QC_imputed_${CHR}.vcf.gz
    bgzip all_filtered_QC_imputed_${CHR}.vcf
    bcftools index \
    --threads ${NSLOTS} \
    -t all_filtered_QC_imputed_${CHR}.vcf.gz
done

rm -f all_filtered*.frq
rm -f all_filtered_QC_f*
rm -f all_filtered_QC.v*
rm -f ../imputation/Beagle/*.bref3

module unload bcftools/1.4
module load bcftools/1.8

# Re-calculate and add INFO (imputation) values for each chromosome.
for CHR in $(ls *.gz | sed 's|.*_||' | sed 's|\..*||')
do
    # Re-calculate allele frequency and compute Impute2-like INFO score.
    bcftools +fill-tags \
    all_filtered_QC_imputed_${CHR}.vcf.gz \
    -Ov -- \
    -t AF,MAF | \
    bcftools +impute-info -Ov -o all_filtered_QC_imputedINFO_${CHR}.vcf
    bgzip all_filtered_QC_imputedINFO_${CHR}.vcf
    bcftools index --threads ${NSLOTS} -t all_filtered_QC_imputedINFO_${CHR}.vcf.gz
done

rm -f all_filtered_QC_imputed_c*

# Generate an allele frequency file for each chr.
for CHR in $(ls all_filtered_QC_imputedINFO*.gz | sed 's|.*_||' | sed 's|\..*||')
do
    echo -e 'CHR\tSNP\tREF\tALT\tAF\tINFO\tAF_GROUP' > all_filtered_QC_group_${CHR}.txt

    # Query the required fields and add frequency group (1, 2 or 3) as the last column.
    # $5 is AF values and $7 AF groups.
    bcftools query \
    -f '%CHROM\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/INFO\t-\n' \
    all_filtered_QC_imputedINFO_${CHR}.vcf.gz | \
    awk \
    -v \
    OFS="\t" '{if ($5>=0.01 && $5<=0.99) $7=1; else if(($5>=0.001 && $5<0.01) || ($5<=0.999 && $5>0.99)) $7=2; else $7=3} { print $1, $2, $3, $4, $5, $6, $7 }' \
    >> all_filtered_QC_group_${CHR}.txt
done


# Plot results to check the distribution.
mkdir -p ../imputation/plots

for CHR in $(ls *.txt | sed 's|.*_||' | sed 's|\..*||')
do
    Rscript \
    --no-save \
    ../scripts/plot_INFO_AF.R \
    all_filtered_QC_group_${CHR}.txt \
    ../imputation/plots/all_filtered_QC_plot_chr${CHR} \
    1000GP_imputation_all.frq
done

# Combine the plots per chromosome into a single pdf file.
convert $(ls ../imputation/plots/all_filtered_QC_plot*.png | sort -V) ../imputation/plots/all_filtered_QC_plot.pdf

rm -f  all_filtered_QC_group*.txt
rm -f  1000GP_imputation_all.frq

mv all_filtered_QC* ../imputation
