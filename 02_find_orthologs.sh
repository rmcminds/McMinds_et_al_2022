#!/bin/bash
#SBATCH --job-name=02_find_orthologs
#SBATCH --mem=170G
#SBATCH --time=6-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks=24
#SBATCH --output=outputs/primates_20230304/02_find_orthologs.log

## extract longest isoforms into a fasta
## run orthofinder with all our inferred transcripts but also all the ensembl cds files - only keep genes that have at least one human ensembl identifier in their family tree, so that annotations can be applied to entire tree. maybe collapse all within-species clusters of 'genes' into a single gene with multiple transcripts, before identifying putative '1:1' orthologs. shouldn't matter for our case if theres a duplication; we want to know if there are more transcripts. duplications that aren't 'monophyletic' will still be a problem.

## only run salmon on transcripts that have human ensembl family members. don't worry about potential pooling until at the tx2gene step for tximport
