#!/bin/sh

module unload java
module load java/1.8.0_152-oracle
module load samtools/1.8
module load FastQC/0.11.9
module load R/4.1.1

NSLOTS=20
mkdir -p merged
cd merged

# Retrieve sample IDs from the IonXpress bam file.
ls /store/ArchiveIonTorrent/archivedReports/*TAPG_Prostate* | \
grep '/store/' | \
sed 's|:||' > TAPG_TURP_folders.txt

R -e "

TAPG_TURP_folders <- read.table('TAPG_TURP_folders.txt')
TAPG_TURP_files   <- c()

for (id in TAPG_TURP_folders$V1) {
  temp            <- paste0(id, '/', list.files(paste0(id, '/'), pattern = '.*.bam'))
  TAPG_TURP_files <- c(TAPG_TURP_files, temp)
}

TAPG_TURP_files <- grep(pattern = '.+.bai',
                        x       = TAPG_TURP_files,
                        invert  = T,
                        value   = T)
TAPG_TURP_files <- grep(pattern = '.+.sam',
                        x       = TAPG_TURP_files,
                        invert  = T,
                        value   = T)

write.table(TAPG_TURP_files, 'TAPG_TURP_files.txt',
            quote = F,
            col.names = F,
            row.names = F)

q(save='no')

"

for id in $(cat TAPG_TURP_files.txt)
do
 samtools view \
 -H $id | \
 grep @RG | \
 awk '{print $10}' | \
 sed 's|SM:||' | \
 uniq >> TAPG_TURP_ids.txt
done

R -e "

files  <- read.table('TAPG_TURP_files.txt')
ids    <- read.table('TAPG_TURP_ids.txt')
merged <- cbind(files, ids)

write.csv(merged, file = 'Samples_info.csv', sep = '\t')

q(save='no')

"

# Sort all files.
for name in $(cat Samples_info.csv | sed 1d | sed 's|.*store|/store|' | sed 's|.bam.*||')
do
  samtools sort \
  -@ ${NSLOTS} \
  -n ${name}.bam \
  -O BAM \
  -o ${name}_sorted.bam

done

R -e "

samples    <- read.csv('Samples_info.csv')
samples$V1 <- gsub(pattern     = '.bam',
                   x           = samples$V1,
                   replacement = '_sorted.bam')

# Get the files with repeats, so they can be merged.
test       <- unique(subset(samples, duplicated(V1.1))$V1.1)

# Take the files with no repeats and rename them as "_merged.bam".
uniq       <- subset(samples, !duplicated(V1.1))

for(id in test){
  temp <- subset(samples,
               V1.1 %in% id)$V1;
system(paste0('samtools merge -@ 20 -n ',
              '-O BAM ',
              'merged_',
              id,
              '_merged.bam ',
              ifelse(is.na(ifelse(length(temp[1]) < 0, '', as.character(temp[1]))), '',
                     ifelse(length(temp[1]) < 0, '', as.character(temp[1]))), ' ',
              ifelse(is.na(ifelse(length(temp[2]) < 0, '', as.character(temp[2]))), '',
                     ifelse(length(temp[2]) < 0, '', as.character(temp[2]))), ' ',
              ifelse(is.na(ifelse(length(temp[3]) < 0, '', as.character(temp[3]))), '',
                     ifelse(length(temp[3]) < 0, '', as.character(temp[3]))), ' ',
              ifelse(is.na(ifelse(length(temp[4]) < 0, '', as.character(temp[4]))), '',
                     ifelse(length(temp[4]) < 0, '', as.character(temp[4]))), ' ',
              ifelse(is.na(ifelse(length(temp[5]) < 0, '', as.character(temp[5]))), '',
                     ifelse(length(temp[5]) < 0, '', as.character(temp[5]))), ' ',
              ifelse(is.na(ifelse(length(temp[6]) < 0, '', as.character(temp[6]))), '',
                     ifelse(length(temp[6]) < 0, '', as.character(temp[6]))), ' ',
              ifelse(is.na(ifelse(length(temp[7]) < 0, '', as.character(temp[7]))), '',
                     ifelse(length(temp[7]) < 0, '', as.character(temp[7]))), ' ',
              ifelse(is.na(ifelse(length(temp[8]) < 0, '', as.character(temp[8]))), '',
                     ifelse(length(temp[8]) < 0, '', as.character(temp[8]))), ' ',
              ifelse(is.na(ifelse(length(temp[9]) < 0, '', as.character(temp[9]))), '',
                     ifelse(length(temp[9]) < 0, '', as.character(temp[9]))), ' ',
              ifelse(is.na(ifelse(length(temp[10]) < 0, '', as.character(temp[10]))), '',
                     ifelse(length(temp[10]) < 0, '', as.character(temp[10]))), ' ',
              ifelse(is.na(ifelse(length(temp[11]) < 0, '', as.character(temp[11]))), '',
                     ifelse(length(temp[11]) < 0, '', as.character(temp[11])))))
}

for(id in 1:nrow(uniq)){
  system(paste0('mv ',
              uniq$V1[id],
              ' ',
              'merged_',
              uniq$V1.1[id],
              '_merged.bam'))

}

q(save='no')

"

# Remove the sorted bam files.
for id in $(cat TAPG_TURP_folders.txt)
do
    rm -f ${id}/*_sorted.bam
done

rm -f TAPG_TURP*.txt
rm -f Samples_info.csv

# Check the quality of the samples.
for samples in $(ls *bam | sed 's|.*/||' | sed 's|_.*||')
do
  fastqc -t ${NSLOTS} -f bam ${samples}*.bam -o QC
done

# Keep the reads with quality score above 20.
samtools view -b -q ${NSLOTS} ${samples}_merged.bam -O BAM -o ${samples}_merged.bam
