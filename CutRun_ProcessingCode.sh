srun --tasks=1 --cpus-per-task=2 --mem=120G --time=120:00:00 --pty /bin/bash

##### fastqc
module load fastqc/0.12.1
module load java/jdk-18.0.1.1

cd /project/lemaire/11_CutAndRun/1_fastq/
for R in *.fq.gz
do
   fastqc -o ../3_fastqc/ $R
done

####### trimmomatic
module load java/jdk-18.0.1.1
module load trimmomatic/0.39

cd /project/lemaire/11_CutAndRun/1_fastq/
for R1 in *_1.fq.gz
do
   R2=${R1//_1.fq.gz/_2.fq.gz}
   R1paired=${R1//.fq.gz/_paired.fq.gz}
   R1unpaired=${R1//.fq.gz/_unpaired.fq.gz}	
   R2paired=${R2//.fq.gz/_paired.fq.gz}
   R2unpaired=${R2//.fq.gz/_unpaired.fq.gz}	
   java -jar /share/apps/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 -threads 32 $R1 $R2 ../2_trimed/$R1paired ../2_trimed/$R1unpaired ../2_trimed/$R2paired ../2_trimed/$R2unpaired ILLUMINACLIP:/share/apps/trimmomatic/0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

##### fastqc of trimed data
module load fastqc/0.12.1
module load java/jdk-18.0.1.1

cd /project/lemaire/11_CutAndRun/2_trimed/
for R in *_paired.fq.gz
do
   fastqc -o ../3_fastqc/ $R
done

##### bowtie2 alignment
module load bowtie2/2.5.1
#module load python/3.4.1

cd /project/lemaire/11_CutAndRun/2_trimed/

for R1 in *_1_paired.fq.gz
do
   R2=${R1//_1_paired.fq.gz/_2_paired.fq.gz}
   Out=${R1//_1_paired.fq.gz/.sam}
   bowtie2 --local -x /project/lemaire/0_ref/human/indexBowtie2_hg38/indexBowtie2_hg38 --fr -1 $R1 -2 $R2 -S ../4_alignment/$Out
done


##### samtools and bedtools
module load samtools/1.15.1
module load bedtools/2.31.1

cd /project/lemaire/11_CutAndRun/4_alignment/
for R1 in *.sam
do
   Out1=${R1//.sam/.bam}
   Out2=${R1//.sam/_fixmate.bam}
   Out3=${R1//.sam/_fixmate_sorted.bam}
   Out4=${R1//.sam/_fixmate_sorted_markdup.bam}
   Out5=${R1//.sam/_fixmate_sorted_markdup_rmBlackList.bam}
   OutTxt=${R1//.sam/.samtools.out.txt}
   
   samtools view -S -b $R1 > ../5_samtools/$Out1
   samtools fixmate -m ../5_samtools/$Out1 ../5_samtools/$Out2 # 
   samtools sort ../5_samtools/$Out2 -o ../5_samtools/$Out3
   samtools markdup -r -f ../5_samtools/$OutTxt ../5_samtools/$Out3 ../5_samtools/$Out4
   bedtools intersect -v -abam ../5_samtools/$Out4 -b ../hg38-blacklist.v2.bed > ../5_samtools/$Out5
   samtools index ../5_samtools/$Out5
done

### need to combine the sam/bam files for sample sample
cd /project/lemaire/11_CutAndRun/5_samtools/
samtools merge MAC_H3K4M3.bam MAC_H3K4M3_*.bam
samtools merge MAC_IgG.bam MAC_IgG_*.bam
samtools merge MAC_IRF3C_H3K27AC.bam MAC_IRF3C_H3K27AC_*.bam
samtools merge MAC_IRF3C_H3K27M3.bam MAC_IRF3C_H3K27M3_CKDL240017754-1A_H3VJMDSXC_L2.bam MAC_IRF3C_H3K27m3_CKDL240020738-1A_HHCG2DSXC_L3.bam
samtools merge MAC_IRF3G_H3K27M3.bam MAC_IRF3G_H3K27M3_CKDL240017755-1A_H3VJMDSXC_L2.bam MAC_IRF3G_H3K27m3_CKDL240020739-1A_HHCG2DSXC_L3.bam
samtools merge MAC_STINGC_H3K27AC.bam MAC_STINGC_H3K27AC_*.bam
samtools merge MAC_STINGC_H3K27M3.bam MAC_STINGC_H3K27M3_CKDL240017756-1A_H3VJMDSXC_L2.bam MAC_STINGC_H3K27m3_CKDL240020742-1A_HHCG2DSXC_L3.bam
samtools merge MAC_STINGG_H3K27AC.bam MAC_STINGG_H3K27AC_*.bam
samtools merge MAC_STINGG_H3K27M3.bam MAC_STINGG_H3K27M3_CKDL240017757-1A_H3VJMDSXC_L2.bam MAC_STINGG_H3K27m3_CKDL240020743-1A_HHCG2DSXC_L3.bam
samtools merge MAC_WTC_H3K27AC.bam MAC_WTC_H3K27AC_*.bam
samtools merge MAC_WTC_H3K27M3.bam MAC_WTC_H3K27M3_CKDL240017752-1A_H3VJMDSXC_L2.bam MAC_WTC_H3K27m3_CKDL240020736-1A_HHCG2DSXC_L3.bam
samtools merge MAC_WTC_IRF3.bam MAC_WTC_IRF3_*.bam
samtools merge MAC_WTG_H3K27AC.bam MAC_WTG_H3K27AC_*.bam
samtools merge MAC_WTG_H3K27M3.bam MAC_WTG_H3K27M3_CKDL240017753-1A_H3VJMDSXC_L2.bam MAC_WTG_H3K27m3_CKDL240020737-1A_HHCG2DSXC_L3.bam
samtools merge MAC_WTG_IRF3.bam MAC_WTG_IRF3_*.bam

for R1 in *.bam
do
   Out2=${R1//.bam/_fixmate.bam}
   Out3=${R1//.bam/_fixmate_sorted.bam}
   Out4=${R1//.bam/_fixmate_sorted_markdup.bam}
   Out5=${R1//.bam/_fixmate_sorted_markdup_rmBlackList.bam}
   OutTxt=${R1//.bam/.samtools.out.txt}
   
   samtools fixmate -m $R1 $Out2 # 
   samtools sort $Out2 -o $Out3
   samtools markdup -r -f $OutTxt $Out3 $Out4
   bedtools intersect -v -abam $Out4 -b ../hg38-blacklist.v2.bed > $Out5
   samtools index $Out5
done



####### macs2 call peaks
module load  anaconda3/2019.10

cd /project/lemaire/11_CutAndRun/5_samtools/
for D in *_fixmate_sorted_markdup_rmBlackList.bam
do
   Out1=${D//_fixmate_sorted_markdup_rmBlackList.bam/_rmBL}
   
   macs2 callpeak -t $D -c TFs/MAC_IgG_CKDL240017750-1A_H3VJMDSXC_L2_fixmate_sorted_markdup_rmBlackList.bam -f BAMPE -g hs -n $Out1 -B --slocal 1500 --outdir ../7_macs2/ -q 0.05 --broad
done

cd /project/lemaire/11_CutAndRun/5_samtools/TFs
for D in *_fixmate_sorted_markdup_rmBlackList.bam
do
   Out1=${D//_fixmate_sorted_markdup_rmBlackList.bam/_rmBL}
   
   macs2 callpeak -t $D -c MAC_IgG_CKDL240017750-1A_H3VJMDSXC_L2_fixmate_sorted_markdup_rmBlackList.bam -f BAMPE -g hs -n $Out1 -B --slocal 1500 --outdir /project/lemaire/11_CutAndRun/7_macs2/ -q 0.05 #--broad
done


macs2 callpeak -t SMC_Ct_H3K4M3_CKDL240008316-1A_H5MHCDSXC_L2_fixmate_sorted_markdup_rmBlackList.bam -c SMC_Ct_IgG_CKDL240008317-1A_H5MHCDSXC_L2_fixmate_sorted_markdup_rmBlackList.bam -f BAMPE -g hs -n M_c_H3K4Me3_rmBL -B --slocal 1500 --outdir ../7_macs2/ -q 0.05 --broad


##### deeptools covert to bigwig
module load  anaconda3/2019.10

cd /project/lemaire/11_CutAndRun/5_samtools/
for D in *_fixmate_sorted_markdup_rmBlackList.bam
do
    Out1=${D//.bam/.SeqDepthNorm.bw}
	bamCoverage --bam $D -o ../6_deepTools/$Out1 --binSize 10 --effectiveGenomeSize 2913022398 --normalizeUsing RPGC --ignoreDuplicates --ignoreForNormalization chrM --extendReads -p max
done

plotFingerprint -b *_fixmate_sorted_markdup_rmBlackList.bam --labels SMC_cp_TFAM_ba SMC_cp_TFAM_cst SMC_Ct_H3K4M3 SMC_Ct_IgG SMC_Ct_IRF3_a SMC_Ct_IRF3 SMC_Ct_MEF2C SMC_Ct_TFAM_ba SMC_Ct_TFAM_cst SMC_H2_TFAM_ba SMC_H2_TFAM_cst SMC_HT_IRF3 --minMappingQuality 30 --skipZeros --region 1 --numberOfSamples 50000 -T "Fingerprints of different samples" --plotFile fingerprints_chr1_SMC.png --outRawCounts fingerprints_chr1.tab


bamCoverage --bam M_c_IRF3_b_fixmate_sorted_markdup_rmBlackList.bam -o ../6_deepTools/M_c_IRF3_b_fixmate_sorted_markdup_rmBlackList.SeqDepthNorm.bw --binSize 10 --normalizeTo1x 2913022398 --normalizeUsingRPKM --ignoreDuplicates --ignoreForNormalization chrM --extendReads -p max

#ploting
cd /project/lemaire/11_CutAndRun/6_deepTools/
computeMatrix reference-point --referencePoint center -S SMC_Ct_IRF3_a_CKDL240008319-1A_H5MHCDSXC_L2_fixmate_sorted_markdup_rmBlackList.SeqDepthNorm.bw SMC_Ct_IRF3_CKDL240008326-1A_H5MHCDSXC_L2_fixmate_sorted_markdup_rmBlackList.SeqDepthNorm.bw SMC_HT_IRF3_CKDL240008327-1A_H5MHCDSXC_L2_fixmate_sorted_markdup_rmBlackList.SeqDepthNorm.bw -R SMC_HT_IRF3_CKDL240008327-1A_H5MHCDSXC_L2_rmBL_peaks.narrowPeak -a 3000 -b 3000 -o matrix_SMC_HtDNA_IRF3_Summit.mat.gz --skipZeros --sortRegions descend --sortUsing max --outFileSortedRegions regions_SMC_HtDNA_IRF3_Summit.bed --outFileNameMatrix matrix_SMC_HtDNA_IRF3_Summit.txt
plotHeatmap -m matrix_SMC_HtDNA_IRF3_Summit.mat.gz -out Heatmap_on_SMC_HtDNA_IRF3_SummitSummit.pdf


######################## creat my own environment for bowtie2
module load anaconda3/2019.10
conda create --prefix /project/lemaire/1_software/bowtie2
source activate /project/lemaire/1_software/bowtie2
conda install -c bioconda bowtie2

srun --tasks=1 --cpus-per-task=2 --mem=210G --time=72:00:00 --pty /bin/bash
module load anaconda3/2019.10
source activate /project/lemaire/1_software/bowtie2


conda create --prefix /project/lemaire/1_software/macs2
source activate /project/lemaire/1_software/macs2
conda install -c bioconda macs2

conda create --prefix /project/lemaire/1_software/deepTools
source activate /project/lemaire/1_software/deepTools

conda create --prefix /project/lemaire/1_software/cellDancer python==3.7.6



