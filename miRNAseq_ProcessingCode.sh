srun --tasks=1 --cpus-per-task=2 --mem=120G --time=120:00:00 --pty /bin/bash


##### fastqc
module load fastqc/0.12.1
module load java/jdk-18.0.1.1

cd /project/lemaire/12_miRNAseq/1_fastq/
for R in *.fq.gz
do
   fastqc -o ../2_fastqc/ $R
done

for R in *.clean.fa.gz
do
   fastqc -o ../2_fastqc/ $R
done


######## trimmomatic
module load java/jdk-18.0.1.1
module load trimmomatic/0.39

cd /project/lemaire/12_miRNAseq/1_fastq/
for R in *.fq.gz
do
   Out=${R//.fq.gz/_trimed.fq.gz}

   java -jar /share/apps/trimmomatic/0.39/trimmomatic-0.39.jar SE -phred33 -threads 32 $R ../3_trimed/$Out ILLUMINACLIP:/share/apps/trimmomatic/0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:10
done





##### fastqc
cd /project/lemaire/12_miRNAseq/3_trimed/
for R in *_trimed.fq.gz
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
cd /project/lemaire/12_miRNAseq/3_trimed/
for R in *_trimed.fq.gz
do
   Out=${R//_trimed.fq.gz/}
   STAR --runThreadN 10 \
        --genomeDir /project/lemaire/0_ref/human/hg38_StarIndex/ \
        --sjdbGTFfile /project/lemaire/0_ref/human/gencode.v29.annotation.gtf \
        --readFilesIn $R \
        --readFilesCommand zcat \
        --alignEndsType EndToEnd \
        --outFilterMismatchNmax 1 \
        --outFilterMultimapScoreRange 0 \
        --quantMode TranscriptomeSAM GeneCounts \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 10 \
        --outSAMunmapped Within \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --outFilterMatchNmin 16 \
        --alignSJDBoverhangMin 1000 \
        --alignIntronMax 1 \
        --outWigType wiggle \
        --outWigStrand Stranded \
        --outWigNorm RPM \
        --outFileNamePrefix ../4_star/$Out
done


##### miRDeep2 processing
module load anaconda3/2019.10

conda create --prefix /project/lemaire/1_software/miRDeep2
source activate /project/lemaire/1_software/miRDeep2
conda install bioconda::mirdeep2

## bowtie index
bowtie-build /project/lemaire/0_ref/human/GRCh38.primary_assembly.genome.fa /project/lemaire/0_ref/human/GRCh38.primary_assembly_BowtieIndex.genome.fa

## collapse reads
cd /project/lemaire/12_miRNAseq/3_trimed/
for R in *_trimed.fq.gz
do
   R1=${R//_trimed.fq.gz/.fastq}
   Out1=${R//_trimed.fq.gz/_reads_collapsed.fa}
   Out2=${R//_trimed.fq.gz/_reads_vs_refdb.arf}
   zcat $R > $R1
   mapper.pl $R1 -e -h -i -j -l 18 -m -p /project/lemaire/0_ref/human/GRCh38.primary_assembly_BowtieIndex.genome.fa -s ../5_miRDeep2/$Out1 -t ../5_miRDeep2/$Out2 -v -o 4
done

## identify known and novel miRNA
cd /project/lemaire/12_miRNAseq/5_miRDeep2/
grep -A1 '>hsa-' hairpin.fa > hsa_hairpin.fa
grep -A1 '>hsa-' mature.fa > hsa_mature.fa
grep -A1 '>mmu-' mature.fa > mmu_mature.fa

for file in hsa_mature.fa mmu_mature.fa hsa_hairpin.fa
do
    awk '/^>/ {print $1; next} {print}' $file > ${file%.fa}_fixed.fa
done
#sed -n '170,180p' hsa_mature_fixed.fa
sed '/^--/d' hsa_mature_fixed.fa > hsa_mature_fixed2.fa
sed '/^--/d' mmu_mature_fixed.fa > mmu_mature_fixed2.fa
sed '/^--/d' hsa_hairpin_fixed.fa > hsa_hairpin_fixed2.fa

for R1 in *_reads_collapsed.fa
do
   R2=${R1//_reads_collapsed.fa/_reads_vs_refdb.arf}
   Out=${R1//_reads_collapsed.fa/_report.log}
   miRDeep2.pl $R1 /project/lemaire/0_ref/human/GRCh38.primary_assembly.genome.fixed.fa $R2 hsa_mature_fixed2.fa mmu_mature_fixed2.fa hsa_hairpin_fixed2.fa -t hsa 2>$Out
done

conda deactivate







##### prepare anno file for downstream analysis
cd /project/lemaire/0_ref/human/
awk '{ if($3 == "gene") print $1,$3,$4,$5,$10,$12,$14}' gencode.v29.annotation.gtf > gencode.v29.annotation.gene.txt
sed 's/\"//g' gencode.v29.annotation.gene.txt > gencode.v29.annotation.gene2.txt
sed 's/\;//g' gencode.v29.annotation.gene2.txt > gencode.v29.annotation.gene3.txt

