#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=6-20:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks=20
#SBATCH --output=outputs/primates/01_filter_genetrees.log

module load apps/R

Rscript scripts/immune_allometry/01_filter_genetrees.r
