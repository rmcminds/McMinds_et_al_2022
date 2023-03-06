#!/bin/bash
#SBATCH --job-name=02_find_orthologs
#SBATCH --mem=175G
#SBATCH --time=7-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks-per-node=24
#SBATCH --output=outputs/primates_20230304/02_find_orthologs.log

out_dir=outputs/primates_20230304/02_find_orthologs/

mkdir -p ${out_dir}/peptides

## make sure this is 1 per gene, not 1 per transcript
for i in outputs/primates_20230304/01_find_transcripts/*_longest_peptide_per_gene.pep; do
newname=$(basename ${i})
newname=${i/_longest_peptide_per_gene.pep/.fa}
newname=${i^}
cp ${i} ${out_dir}/peptides/${newname}
done

##download ensembl human cds and add to transcripts folder
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O ${out_dir}/peptides/Homo_sapiens.GRCh38.pep.all.fa.gz

zcat ${out_dir}/peptides/Homo_sapiens.GRCh38.pep.all.fa.gz > ${out_dir}/peptides/Homo_sapiens_ensembl.pep
rm ${out_dir}/peptides/Homo_sapiens.GRCh38.pep.all.fa.gz

module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda activate orthofinder

## run orthofinder with all transdecoder longest isoform peptides, plus ensembl human peptides for annotation

orthofinder -t 24 -M msa -A mafft -T 'iqtree -nt' \
  -s raw_data/20221215_primate_allometry/primates_ensemblDup.newick \
  -f ${out_dir}/peptides \
  -o ${out_dir}/of_out
