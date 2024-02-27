#!/bin/bash
#PBS -m e 
#PBS -M yanming.li@bcm.edu
#PBS -k oe
#PBS -l nodes=1:ppn=2
#PBS -l vmem=64GB
#PBS -l walltime=60:00:00

cd /project/lemaire/13_Rongmo/5_samtools/

module load macs/2.1.0

for R1 in *_fixmate_sorted_markdup.bam
do
   R2=${R1//_fixmate_sorted_markdup.bam/_MACS2}
   macs2 callpeak -t $R1 -f BAMPE -g hs -n $R2 -B --slocal 1500 --outdir ../6_MACS2/ -q 0.05 
done

for R1 in *_fixmate_sorted_markdup_rmBlackList.bam
do
   R2=${R1//_fixmate_sorted_markdup_rmBlackList.bam/_rmBlackList_MACS2}
   macs2 callpeak -t $R1 -f BAMPE -g hs -n $R2 -B --slocal 1500 --outdir ../6_MACS2/ -q 0.05 
done

