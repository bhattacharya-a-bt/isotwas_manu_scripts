#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=20G,h_rt=40:00:00,highp
#$ -o /u/scratch/a/abtbhatt/PEC_Fetal_Models/outFiles
#$ -t 1-32
#$ -pe shared 1
#$ -tc 16

. /u/local/Modules/default/init/modules.sh
module load R
cd /u/scratch/a/abtbhatt
Rscript train_PEC_Fetal.R --index ${SGE_TASK_ID} --number 1000