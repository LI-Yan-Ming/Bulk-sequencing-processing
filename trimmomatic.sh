#!/bin/bash
#PBS -m e
#PBS -M yanming.li@bcm.edu
#PBS -k oe
#PBS -l nodes=1:ppn=2
#PBS -l vmem=60GB
#PBS -l walltime=12:00:00

cd /project/lemaire/13_Rongmo/1_fastq/
module load java/jdk-18.0.1.1
module load trimmomatic/0.39


for R1 in *_R1_001.fastq.gz
do
   R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz}
   R1paired=${R1//.fastq.gz/_paired.fastq.gz}
   R1unpaired=${R1//.fastq.gz/_unpaired.fastq.gz}	
   R2paired=${R2//.fastq.gz/_paired.fastq.gz}
   R2unpaired=${R2//.fastq.gz/_unpaired.fastq.gz}	
   java -jar /share/apps/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 -threads 32 $R1 $R2 ../2_trimed/$R1paired ../2_trimed/$R1unpaired ../2_trimed/$R2paired ../2_trimed/$R2unpaired ILLUMINACLIP:/share/apps/trimmomatic/0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

## specify minimal length to be 29
for R1 in *_R1_001.fastq.gz
do
   R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz}
   R1paired=${R1//.fastq.gz/_min29_paired.fastq.gz}
   R1unpaired=${R1//.fastq.gz/_min29_unpaired.fastq.gz}	
   R2paired=${R2//.fastq.gz/_min29_paired.fastq.gz}
   R2unpaired=${R2//.fastq.gz/_min29_unpaired.fastq.gz}	
   java -jar /share/apps/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 -threads 32 $R1 $R2 ../2_trimed/$R1paired ../2_trimed/$R1unpaired ../2_trimed/$R2paired ../2_trimed/$R2unpaired ILLUMINACLIP:/share/apps/trimmomatic/0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:29
done


