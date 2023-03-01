#!/bin/bash
#SBATCH --job-name=00_references
#SBATCH --mem=170G
#SBATCH --time=6-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks=20
#SBATCH --output=outputs/primates_20230224/00_references.log

mkdir -p outputs/primates_20230224/00_references
cd outputs/primates_20230224/00_references

## download ensembl gene trees as flatfile
wget -N https://ftp.ensembl.org/pub/release-109/emf/ensembl-compara/homologies/Compara.109.protein_default.nh.emf.gz

## extract newick strings from flatfile for later manipulation
zcat Compara.109.protein_default.nh.emf.gz | grep '^(' > Compara.109.protein_default.newick

module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda activate salmon
for spec in callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii; do

  ## download transcriptome reference
  wget -e robots=off -r -N -l1 -nd -A '*.cdna.all.fa.gz' https://ftp.ensembl.org/pub/release-109/fasta/${spec}/cdna/

  ## build transcriptome index
  salmon index --index ${spec}_index --transcripts ${spec^}.*.cdna.all.fa.gz

done
