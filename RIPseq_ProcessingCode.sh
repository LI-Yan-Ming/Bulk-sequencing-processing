srun --tasks=1 --cpus-per-task=2 --mem=120G --time=120:00:00 --pty /bin/bash


##### fastqc
module load fastqc/0.12.1
module load java/jdk-18.0.1.1

cd /project/lemaire/7_RIPseq/1_fastq/
for R in *.fq.gz
do
   fastqc -o ../2_fastqc/ $R
done


######## trimmomatic remove Illumina Universal Adapter
module load java/jdk-18.0.1.1
module load trimmomatic/0.39

cd /project/lemaire/7_RIPseq/1_fastq/
for R1 in *_1.fq.gz
do
   R2=${R1//_1.fq.gz/_2.fq.gz}
   R1paired=${R1//.fq.gz/_paired.fq.gz}
   R1unpaired=${R1//.fq.gz/_unpaired.fq.gz}	
   R2paired=${R2//.fq.gz/_paired.fq.gz}
   R2unpaired=${R2//.fq.gz/_unpaired.fq.gz}	
   java -jar /share/apps/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 -threads 32 $R1 $R2 ../3_trimmomatic/$R1paired ../3_trimmomatic/$R1unpaired ../3_trimmomatic/$R2paired ../3_trimmomatic/$R2unpaired ILLUMINACLIP:/share/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done


##### fastqc of trimed data
module load fastqc/0.12.1
module load java/jdk-18.0.1.1

cd /project/lemaire/7_RIPseq/3_trimmomatic/
for R in *_paired.fq.gz
do
   fastqc -o ../2_fastqc/ $R
done


##### mapping by star
module load star/2.7.10b

# make index
cd /project/lemaire/0_ref/human/
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
# gunzip GRCh38.primary_assembly.genome.fa.gz
# gunzip gencode.v29.annotation.gtf.gz


cd /project/lemaire/0_ref/human/
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir hg38_StarIndex_forCLIPseq --sjdbGTFfile ENCFF159KBI.gtf --genomeFastaFiles GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta --sjdbOverhang 99

# Run alignment
cd /project/lemaire/7_RIPseq/3_trimmomatic/

for R1 in *_1_paired.fq.gz
do
   R2=${R1//_1_paired.fq.gz/_2_paired.fq.gz}
   Out=${R1//_1_paired.fq.gz/}
   STAR --alignEndsType EndToEnd \
        --genomeDir /project/lemaire/0_ref/human/hg38_StarIndex_forCLIPseq \
        --genomeLoad NoSharedMemory \
        --outBAMcompression 10 \
        --outFileNamePrefix ../4_star/$Out \
        --outFilterMultimapNmax 1 \
        --outFilterMultimapScoreRange 1 \
        --outFilterScoreMin 10 \
        --outFilterType BySJout \
        --outReadsUnmapped Fastx \
        --outSAMattrRGline ID:foo \
        --outSAMattributes All \
        --outSAMmode Full \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outStd Log \
        --readFilesIn $R1 $R2 \
        --readFilesCommand zcat \
        --runMode alignReads \
        --runThreadN 8
done


##### samtools and bedtools
##### sort, remove deplicate (not avaliable for single end), remove repeat
module load samtools/1.15.1
module load bedtools/2.31.1

cd /project/lemaire/7_RIPseq/4_star/
for R1 in *.bam
do
   Out1=${R1//.sam/.bam}
   Out2=${R1//.bam/_fixmate.bam}
   Out3=${R1//.bam/_sorted.bam}
   Out4=${R1//.bam/_sorted_markdup.bam}
   Out5=${R1//.bam/_sorted_rmBlackList.bam}
   OutTxt=${R1//.bam/.samtools.out.txt}
   
   samtools view -S -b $R1 > ../5_samtools/$Out1
   samtools fixmate -m ../5_samtools/$Out1 ../5_samtools/$Out2 # 
   samtools sort ../5_samtools/$Out2 -o ../5_samtools/$Out3
   samtools markdup -r -f ../5_samtools/$OutTxt ../5_samtools/$Out3 ../5_samtools/$Out4
   bedtools intersect -v -abam ../5_samtools/$Out4 -b /project/lemaire/0_ref/human/hg38_rmsk.bed > ../5_samtools/$Out5 ## remove repeat elements
   samtools index ../5_samtools/$Out5
done



####### macs2 call peaks of pair-end reads
module load  anaconda3/2019.10

cd /project/lemaire/7_RIPseq/5_samtools/

macs2 callpeak -t CCCP1_*_rmBlackList.bam -c Input_CCCP_*_rmBlackList.bam -f BAMPE -g 3e8 -n WT_CCCP_1 -B --nomodel --shift -100 --extsize 200 --keep-dup all --outdir /project/lemaire/7_RIPseq/7_macs2/ -q 0.01

macs2 callpeak -t CCCP2_*_rmBlackList.bam -c Input_CCCP_*_rmBlackList.bam -f BAMPE -g 3e8 -n WT_CCCP_2 -B --nomodel --shift -100 --extsize 200 --keep-dup all --outdir /project/lemaire/7_RIPseq/7_macs2/ -q 0.01

macs2 callpeak -t Control1_*_rmBlackList.bam -c Input_control_*_rmBlackList.bam -f BAMPE -g 3e8 -n WT_Ctrl_1 -B --nomodel --shift -100 --extsize 200 --keep-dup all --outdir /project/lemaire/7_RIPseq/7_macs2/ -q 0.01

macs2 callpeak -t OE_TFAM1_*_rmBlackList.bam -c Input_OE_TFAM_CKDL240035418-1A_*_rmBlackList.bam -f BAMPE -g 3e8 -n OETFAM_Ctrl_1 -B --nomodel --shift -100 --extsize 200 --keep-dup all --outdir /project/lemaire/7_RIPseq/7_macs2/ -q 0.01

macs2 callpeak -t OE_TFAM2_*_rmBlackList.bam -c Input_OE_TFAM_CKDL240035418-1A_*_rmBlackList.bam -f BAMPE -g 3e8 -n OETFAM_Ctrl_2 -B --nomodel --shift -100 --extsize 200 --keep-dup all --outdir /project/lemaire/7_RIPseq/7_macs2/ -q 0.01

macs2 callpeak -t OE_TFAM_CCCP1_*_rmBlackList.bam -c Input_OE_TFAM_CCP_*_rmBlackList.bam -f BAMPE -g 3e8 -n OETFAM_CCCP_1 -B --nomodel --shift -100 --extsize 200 --keep-dup all --outdir /project/lemaire/7_RIPseq/7_macs2/ -q 0.01

macs2 callpeak -t OE_TFAM_CCCP2_*_rmBlackList.bam -c Input_OE_TFAM_CCP_*_rmBlackList.bam -f BAMPE -g 3e8 -n OETFAM_CCCP_2 -B --nomodel --shift -100 --extsize 200 --keep-dup all --outdir /project/lemaire/7_RIPseq/7_macs2/ -q 0.01


######## deeptools covert to bigwig
module load  anaconda3/2019.10

cd /project/lemaire/7_RIPseq/5_samtools/
for D in *_rmBlackList.bam
do
    Out1=${D//.bam/.SeqDepthNorm.bw}
	bamCoverage --bam $D -o ../6_deepTools/$Out1 --binSize 10 --effectiveGenomeSize 2913022398 --normalizeUsing RPGC --ignoreDuplicates --extendReads -p max
done

plotFingerprint -b *_rmBlackList.bam --labels WT_CCCP_1 WT_CCCP_2 WT_Ctrl_1 Input_WT_CCCP Input_WT_Ctrl Input_OETFAM_CCCP Input_OETFAM_Ctrl OETFAM_Ctrl_1 OETFAM_Ctrl_2 OETFAM_CCCP_1 OETFAM_CCCP_2 --minMappingQuality 30 --skipZeros --region 1 --numberOfSamples 50000 -T "Fingerprints of different samples" --plotFile fingerprints_chr1_SMC.png --outRawCounts fingerprints_chr1.tab


######### quantify reads pf peaks use featureCounts

featureCounts -T 8 -a peaks.bed -F SAF -o counts.txt sample.bam


##### ploting heatmap of all H3K27ac pleaks (or H3K27Me3 peaks)
module load bedtools/2.31.1
cd /project/lemaire/7_RIPseq/7_macs2/

cat *.narrowPeak > TFAM_RIP_peaks.bed
bedtools sort -i TFAM_RIP_peaks.bed > TFAM_RIP_peaks_sorted.bed
bedtools merge -i TFAM_RIP_peaks_sorted.bed > TFAM_RIP_peaks_sorted_merge.bed

cd /project/lemaire/7_RIPseq/6_deepTools/
computeMatrix reference-point --referencePoint center -S *.SeqDepthNorm.bw -R ../7_macs2/TFAM_RIP_peaks_sorted_merge.bed -a 3000 -b 3000 -o matrix_TFAM_RIP_peaks_Summit.mat.gz --skipZeros --sortRegions descend --sortUsing max --outFileSortedRegions regions_TFAM_RIP_peaks_Summit.bed --outFileNameMatrix matrix_TFAM_RIP_peaks_Summit.txt
plotHeatmap -m matrix_TFAM_RIP_peaks_Summit.mat.gz -out Heatmap_on_TFAM_RIP_peaks_Summit.pdf





