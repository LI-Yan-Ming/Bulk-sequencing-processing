#!/bin/bash
#PBS -m e 
#PBS -M yanming.li@bcm.edu
#PBS -k oe
#PBS -l nodes=1:ppn=2
#PBS -l vmem=64GB
#PBS -l walltime=72:00:00

module load bowtie2/2.4.1
module load python/3.4.1

cd /project/lemaire/13_Rongmo/2_trimed/

for R1 in *_R1_001_paired.fastq.gz
do
   R2=${R1//R1_001_paired.fastq.gz/R2_001_paired.fastq.gz}
   Out=${R1//_R1_001_paired.fastq.gz/.sam}
   bowtie2 --local -x /project/lemaire/0_ref/human/indexBowtie2_hg38/indexBowtie2_hg38 --fr -1 $R1 -2 $R2 -S ../4_alignment/$Out
done


