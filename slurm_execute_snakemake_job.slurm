#!/bin/bash
#SBATCH -J snakemake
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -o %N.%J.out # STDOUT
#SBATCH -e %N.%j.err # STDERR
#SBATCH --time=48:00:00

ml load Apptainer
ml load anaconda3
conda activate snakemakeconda activate snakemake

snakemake -j 1 --unlock
snakemake --profile config/snakemake_profile/slurm "$@"
