#!/bin/bash
#SBATCH --job-name=02_find_orthologs
#SBATCH --mem=175G
#SBATCH --time=7-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks-per-node=24
#SBATCH --output=outputs/primates_20230304/02_find_orthologs.log

wd=$(pwd)
out_dir=outputs/primates_20230304/02_find_orthologs/

mkdir -p ${out_dir}/peptides

## copy seqs into directory with clean filenames
for i in outputs/primates_20230304/01_find_transcripts/*_longest_peptide_per_gene.pep; do
  newname=$(basename ${i})
  newname=${newname/_longest_peptide_per_gene.pep/.fa}
  newname=${newname^}
  cp ${i} ${out_dir}/peptides/${newname}
done

## download ensembl human cds and add to peptides folder
mkdir ${out_dir}/ensembl
cd ${out_dir}/ensembl

wget -N https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
wget -N https://ftp.ensembl.org/pub/release-109/emf/ensembl-compara/homologies/Compara.109.protein_default.aa.fasta.gz

cd ${wd}

## keep only the reference peptides that ensembl uses for their gene trees. Could also just not use the Homo_sapiens file and extract from alignment, but then would want to remove alignment dashes
awk 'NR == FNR {
       a[$1]++
       next
     }
     NR != FNR {
       split($1, b, ".")
       if(">"b[1] in a) {
         sub(/\n$/, "")
         print ">"$0
       }
     }' <(zcat ${out_dir}/ensembl/Compara.109.protein_default.aa.fasta.gz | grep 'ENSP[[:digit:]]') RS='(^|\n)>' <(zcat ${out_dir}/ensembl/Homo_sapiens.GRCh38.pep.all.fa.gz) > ${out_dir}/peptides/Homo_ensembl.fa

module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda activate orthofinder

## run orthofinder with all transdecoder longest isoform peptides, plus ensembl human peptides for annotation

orthofinder -t 24 -M msa -A mafft -T iqtree \
  -s raw_data/20221215_primate_allometry/primates_ensemblDup.newick \
  -f ${out_dir}/peptides \
  -o ${out_dir}/of_out
