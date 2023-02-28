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
module load apps/hisat2/2.1.0
for spec in callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii; do

  ## download genome.
  ## genomes seem to need to be unzipped for hisat2 indexing but can be stored zipped after that
  ## if a 'primary assembly' file exists, it means the 'toplevel' files contain alternate chromosomes which will artificially influence multimapping rates, so we want to prioritize primary assemblies if they exist
  wget -e robots=off -r -N -l1 -nd -A '*.dna.primary_assembly.fa.gz,*.dna.toplevel.fa.gz' https://ftp.ensembl.org/pub/release-109/fasta/${spec}/dna/
  if [ -f  ${spec^}.*.primary_assembly.fa.gz ]; then
    zcat ${spec^}.*.primary_assembly.fa.gz > tmp.fa
  else
    zcat ${spec^}.*.toplevel.fa.gz > tmp.fa
  fi
  
  ## download genome annotations for later stringtie calculations. looks like stringtie needs this unzipped
  wget -e robots=off -r -N -l1 -nd -A '*.109.gff3.gz' https://ftp.ensembl.org/pub/release-109/gff3/${spec}/
  gunzip ${spec^}.*.gff3.gz

  ## build genome index
  hisat2-build -p 20 tmp.fa ${spec}_index
  
  rm tmp.fa

done
