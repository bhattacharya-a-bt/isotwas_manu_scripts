#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=30G,h_rt=48:00:00,highp
#$ -o /u/scratch/a/abtbhatt/PEC_Fetal_Models/outFiles
#$ -t 1-72
#$ -pe shared 1

. /u/local/Modules/default/init/modules.sh
module load intel/2022.1.1; module load R/4.2.1
cd /u/project/gandalm/abtbhatt
Rscript trainLCWeights_isotwas.R --index ${SGE_TASK_ID} --number 200