#!/bin/bash
#SBATCH --job-name=02_find_orthologs
#SBATCH --mem=175G
#SBATCH --time=7-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks-per-node=24
#SBATCH --output=outputs/primates_20230314_mixed/02_find_orthologs.log

wd=$(pwd)
out_dir=outputs/primates_20230314_mixed/02_find_orthologs/

mkdir -p ${out_dir}/peptides

## copy seqs into directory with clean filenames
for i in outputs/primates_20230314_mixed/01_find_transcripts/*_longest_peptide_per_gene.pep; do
  newname=$(basename ${i})
  newname=${newname/_longest_peptide_per_gene.pep/.fa}
  newname=${newname^}
  cp ${i} ${out_dir}/peptides/${newname}
done

## download ensembl peptides and add to peptides folder
mkdir ${out_dir}/ensembl
cd ${out_dir}/ensembl

wget -N https://ftp.ensembl.org/pub/release-109/emf/ensembl-compara/homologies/Compara.109.protein_default.aa.fasta.gz

species=(callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii)
spestrings=('ENSCJAP[[:digit:]]' 'ENSP[[:digit:]]' 'ENSMMUP[[:digit:]]' 'ENSMICP[[:digit:]]' 'ENSPANP[[:digit:]]' 'ENSPPYP[[:digit:]]')

for i in 0..5; do

  spec=${species[$i]}
  
  ## download
  wget -e robots=off -r -N -l1 -nd -A '*.pep.all.fa.gz' https://ftp.ensembl.org/pub/release-109/fasta/${spec}/pep/
  wget -e robots=off -r -N -l1 -nd -A '*.cdna.all.fa.gz' https://ftp.ensembl.org/pub/release-109/fasta/${spec}/cdna/
  zcat ${spec^}.*.cdna.all.fa.gz > ${spec}_transcripts.fa

  ## keep only the reference peptides that ensembl uses for their gene trees.
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
       }' <(zcat Compara.109.protein_default.aa.fasta.gz | grep "${spestrings[$i]}") RS='(^|\n)>' <(zcat ${spec^}.*.pep.all.fa.gz) > ${wd}/${out_dir}/peptides/${spec^}.fa

done

cd ${wd}

module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda activate orthofinder

## run orthofinder with all transdecoder longest isoform peptides, plus ensembl peptides

orthofinder -t 24 \
  -s raw_data/20221215_primate_allometry/primates.newick \
  -f ${out_dir}/peptides \
  -o ${out_dir}/of_out


