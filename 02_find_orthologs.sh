#!/bin/bash
#SBATCH --job-name=02_find_orthologs
#SBATCH --mem=175G
#SBATCH --time=7-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks-per-node=24
#SBATCH --output=outputs/primates_20230304/02_find_orthologs.log

out_dir=outputs/primates_20230304/02_find_orthologs/

mkdir -p ${out_dir}/transcripts

cp outputs/primates_20230304/01_find_transcripts/longest/* ${out_dir}/transcripts

##download ensembl human cds and add to transcripts folder
wget -e robots=off -r -N -l1 -nd  https://ftp.ensembl.org/pub/release-109/fasta/Homo_sapiens/cds/    Homo_sapiens.GRCh38.cds.all.fa.gz -O outputs/primates_20230304/02_find_orthologs/transcripts/Homo_sapiens.GRCh38.cds.all.fa.gz

zcat outputs/primates_20230304/02_find_orthologs/transcripts/Homo_sapiens.GRCh38.cds.all.fa.gz > outputs/primates_20230304/02_find_orthologs/transcripts/Homo_sapiens_ensembl.fa
rm outputs/primates_20230304/02_find_orthologs/transcripts/Homo_sapiens.GRCh38.cds.all.fa.gz

module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda activate orthofinder

orthofinder -t 24 -d -M msa -A mafft -T iqtree \
  -s raw_data/20221215_primate_allometry/primates_ensemblDup.newick \
  -f outputs/primates_20230304/02_find_orthologs/transcripts \
  -o outputs/primates_20230304/02_find_orthologs/of_out

## run orthofinder with all our inferred transcripts but also all the ensembl cds files - only keep genes that have at least one human ensembl identifier in their family tree, so that annotations can be applied to entire tree. maybe collapse all within-species clusters of 'genes' into a single gene with multiple transcripts, before identifying putative '1:1' orthologs. shouldn't matter for our case if theres a duplication; we want to know if there are more transcripts. duplications that aren't 'monophyletic' will still be a problem.

