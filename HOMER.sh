#!/bin/bash
#PBS -m e 
#PBS -M yanming.li@bcm.edu
#PBS -k oe
#PBS -l nodes=1:ppn=2
#PBS -l vmem=64GB
#PBS -l walltime=60:00:00

module load homer/4.11
module load R/4.0.5q()

cd /project/lemaire/13_Rongmo/5_samtools/

for R in *_fixmate_sorted_markdup_rmBlackList.bam
do
   Dic=${R//_fixmate_sorted_markdup_rmBlackList.bam/_rmBlackList}
   makeTagDirectory ../8_Homer/$Dic/ $R     
done

cd /project/lemaire/13_Rongmo/8_Homer/

getDifferentialPeaksReplicates.pl -t TGF1_S4_L006_rmBlackList/ TGF2_S5_L006_rmBlackList/ TGF3_S6_L006_rmBlackList/ -b Scr1_S1_L006_rmBlackList/ Scr2_S2_L006_rmBlackList/ Scr3_S3_L006_rmBlackList/ > DiffPeak_TGFvsSCR.txt

getDifferentialPeaksReplicates.pl -t siPDK4-1_S7_L006_rmBlackList/ siPDK4-2_S8_L006_rmBlackList/ siPDK4-3_S9_L006_rmBlackList/ -b Scr1_S1_L006_rmBlackList/ Scr2_S2_L006_rmBlackList/ Scr3_S3_L006_rmBlackList/ > DiffPeak_siPDK4vsSCR.txt




