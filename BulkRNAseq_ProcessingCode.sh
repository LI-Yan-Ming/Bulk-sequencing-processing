srun --tasks=1 --cpus-per-task=2 --mem=210G --time=120:00:00 --pty /bin/bash


##### fastqc
module load fastqc/0.12.1
module load java/jdk-18.0.1.1

cd /project/lemaire/14_bulkRNAseq/1_fastq/
for R in *.fq.gz
do
   fastqc -o ../2_fastqc/ $R
done



##### mapping by star
module load star/2.7.10b

# make index
cd /project/lemaire/0_ref/human/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v29.annotation.gtf.gz

STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir /project/lemaire/0_ref/human/hg38_StarIndex \
     --genomeFastaFiles /project/lemaire/0_ref/human/GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile /project/lemaire/0_ref/human/gencode.v29.annotation.gtf \
     --sjdbOverhang 100

# Run alignment
cd /project/lemaire/14_bulkRNAseq/1_fastq/
for R1 in *_1.fq.gz
do
   R2=${R1//_1.fq.gz/_2.fq.gz}
   Out=${R1//_1.fq.gz/}
   STAR --runThreadN 10 \
        --genomeDir /project/lemaire/0_ref/human/hg38_StarIndex/ \
        --sjdbGTFfile /project/lemaire/0_ref/human/gencode.v29.annotation.gtf \
        --readFilesIn $R1 $R2 \
        --readFilesCommand zcat \
        --outFilterType BySJout \
        --outFilterMultimapNmax 10 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000\
        --alignMatesGapMax 1000000 \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode TranscriptomeSAM \
        --outFileNamePrefix ../3_star/$Out
done

#--outSAMattributes NH HI AS NM MD \
#--outSAMunmapped Within \


##### Quantify gene/Transcript using RSEM
module load rsem/1.3.3
module load perl/5.32.1

# prepare reference
cd /project/lemaire/0_ref/human/
rsem-prepare-reference --gtf /project/lemaire/0_ref/human/gencode.v29.annotation.gtf \
		               --star \
					   --star-path /share/apps/star/2.7.10b/bin\
		               /project/lemaire/0_ref/human/GRCh38.primary_assembly.genome.fa \
		               /project/lemaire/0_ref/human/hg38_rsem_ref

# calculate expression
cd /project/lemaire/14_bulkRNAseq/3_star/
for R1 in *Aligned.toTranscriptome.out.bam
do
   Out=${R1//Aligned.sortedByCoord.out.bam/}
   rsem-calculate-expression --bam \
                             --estimate-rspd  \
                             --calc-ci \
                             --no-bam-output \
                             --seed 12345 \
                             -p 8 \
                             --ci-memory 1024 \
                             --paired-end \
                             --forward-prob 0.5 \
                             /project/lemaire/14_bulkRNAseq/3_star/$R1 \
                             /project/lemaire/0_ref/human/hg38_rsem_ref/hg38_rsem_ref \
                             /project/lemaire/14_bulkRNAseq/4_RSEM/$Out
done



