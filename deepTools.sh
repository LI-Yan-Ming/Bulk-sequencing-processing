#!/bin/bash
#PBS -m e 
#PBS -M yanming.li@bcm.edu
#PBS -k oe
#PBS -l nodes=1:ppn=2
#PBS -l vmem=64GB
#PBS -l walltime=60:00:00

module load anaconda2/4.2.0 deepTools/2.4.2 
cd /project/lemaire/13_Rongmo/5_samtools/

for D in *_fixmate_sorted_markdup_rmBlackList.bam
do
    Out1=${D//.bam/.SeqDepthNorm.bw}
	bamCoverage --bam $D -o ../7_deepTools/$Out1 --binSize 10 --normalizeTo1x 2913022398 --normalizeUsingRPKM --ignoreDuplicates --ignoreForNormalization chrM --extendReads -p max
done



cd /project/lemaire/13_Rongmo/7_deepTools/

plotFingerprint -b ../5_samtools/ZV-61037_S45_L001_fixmate_sorted_markdup.bam ../5_samtools/Aligned_HUS-58562_S6_L002.sorted.bam --labels SMC_H2O2_MEF2C SMC_Ctrl_input --minMappingQuality 30 --skipZeros --region 1 --numberOfSamples 50000 -T "Fingerprints of different samples" --plotFile fingerprints_chr1.png --outRawCounts fingerprints_chr1.tab

plotCoverage -b ../5_samtools/ZV-61037_S45_L001_fixmate_sorted_markdup.bam ../5_samtools/Aligned_HUS-58562_S6_L002.sorted.bam --plotFile chr19_coverage -n 1000000 --plotTitle "chr19_coverage_MEF2C" --outRawCounts coverage.tab --ignoreDuplicates --minMappingQuality 10 --region 19

#computeGCBias -b ../5_samtools/*sorted.bam --effectiveGenomeSize 2913022398 --genome genome.2bit -l 200 -freq freq_test.txt --region 19 --biasPlot chr19_gc.png

bamCoverage --bam ../5_samtools/ZV-61037_S45_L001_fixmate_sorted_markdup.bam -o SMC_H2O2_MEF2C.SeqDepthNorm.bw --binSize 10 --normalizeTo1x 2913022398 --normalizeUsingRPKM --ignoreDuplicates --ignoreForNormalization chrM --extendReads -p max


## plot heatmap
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -R ../9_bedtools/TGF_down_peak_skinny_uniq.bed -S ../7_deepTools/Round1/*_rmBlackList.SeqDepthNorm.bw --skipZeros -o matrix1_TGFpeakDown_center.gz --outFileSortedRegions TGFpeakDown_allSample.bed

plotHeatmap -m matrix1_TGFpeakDown_center.gz -out TGFpeakDown_Heatmap1.png






	
