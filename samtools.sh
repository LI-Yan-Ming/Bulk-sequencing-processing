#!/bin/bash
#PBS -m e 
#PBS -M yanming.li@bcm.edu
#PBS -k oe
#PBS -l nodes=1:ppn=2
#PBS -l vmem=64GB
#PBS -l walltime=72:00:00

module load samtools/1.15.1
cd /project/lemaire/13_Rongmo/4_alignment/

for R1 in *.sam
do
   Out1=${R1//.sam/.bam}
   Out2=${R1//.sam/_fixmate.bam}
   Out3=${R1//.sam/_fixmate_sorted.bam}
   Out4=${R1//.sam/_fixmate_sorted_markdup.bam}
   OutTxt=${R1//.sam/.samtools.out.txt}
   
   samtools view -S -b $R1 > ../5_samtools/$Out1
   samtools fixmate -m ../5_samtools/$Out1 ../5_samtools/$Out2
   samtools sort ../5_samtools/$Out2 -o ../5_samtools/$Out3
   samtools markdup -r -f ../5_samtools/$OutTxt ../5_samtools/$Out3 ../5_samtools/$Out4
   samtools index ../5_samtools/$Out4
done


module load bedtools/2.29.2
module load samtools/1.15.1
cd /project/lemaire/13_Rongmo/5_samtools/

for R in *_fixmate_sorted_markdup.bam
do
   Out=${R//.bam/_rmBlackList.bam}
   bedtools intersect -v -abam $R -b /project/lemaire/11_CutAndRun/hg38-blacklist.v2.bed > $Out
   samtools index $Out
done



