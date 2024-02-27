#!/bin/bash
#PBS -m e 
#PBS -M yanming.li@bcm.edu
#PBS -k oe
#PBS -l nodes=1:ppn=2
#PBS -l vmem=60GB
#PBS -l walltime=24:00:00

#!/bin/bash
#
#SBATCH -J MITgcm_exp1
#SBATCH –n 48
#SBATCH –p edr
#SBATCH –t 2-12:00:00  # format is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --mem_per_cpu=60G


cd /project/lemaire/13_Rongmo/2_trimed/
module load fastqc/0.12.1
java/1.8.0u71

for R in *_paired.fastq.gz
do
   fastqc -o ../3_fastqc/ $R
done

cd /project/lemaire/13_Rongmo/1_fastq/
for R in *.fastq.gz
do
   fastqc -o ../3_fastqc/ $R
done
