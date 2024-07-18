#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --account=putnamlab
#SBATCH --error=GOMWU.%J.stderr
#SBATCH --output=GOMWU.%J.stdout
#SBATCH --exclusive
#SBATCH --mem=250GB

cd /data/putnamlab/zdellaert/Pdam-TagSeq/go_MWU
module load R/4.2.2-foss-2022b
Rscript GO_MWU_all_mod.R