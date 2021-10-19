#!/bin/bash
#SBATCH --account="gbru_fy21_tomato_ralstonia"
#SBATCH --job-name=keio
#SBATCH --out=keio_%A_%a.log
#SBATCH --time=14-00:00:00
#SBATCH --array=1-58
#SBATCH -p atlas
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem 250G


date;hostname;pwd

# load the required programs
# hpg1-compute
#module load conda
#source activate keio


RUN=${SLURM_ARRAY_TASK_ID}
echo "My Slurm RUN_ID: '${RUN}'"
echo "My TMPDIR IS: " $TMPDIR


infile=$(ls keio_output/*.csv | sed -n ${RUN}p)
echo "$infile"

time Rscript keio.R $infile
