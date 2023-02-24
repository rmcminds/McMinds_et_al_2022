#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=20:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks=20
#SBATCH --output=outputs/primates_20230224/00_references.log

mkdir -p outputs/primates_20230224/00_references
cd outputs/primates_20230224/00_references

wget https://ftp.ensembl.org/pub/release-109/emf/ensembl-compara/homologies/Compara.109.protein_default.nh.emf.gz

zcat Compara.109.protein_default.nh.emf.gz | grep '^(' > Compara.109.protein_default.newick

module purge
module load apps/hisat2/2.1.0
for spec in callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii; do

  ## download genome
  wget -e robots=off -r -l1 -nd -A https://ftp.ensembl.org/pub/release-109/fasta/${spec}/dna/*.dna.toplevel.fa.gz ./

  ## build genome index
  hisat2-build -p 20 <(zcat ${spec^}.*.fa.gz) ${spec}_index

done
