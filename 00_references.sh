#!/bin/bash
#SBATCH --job-name=00_references
#SBATCH --mem=170G
#SBATCH --time=6-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks=24
#SBATCH --output=outputs/primates_20230304/00_references.log

mkdir -p outputs/primates_20230304/00_references
cd outputs/primates_20230304/00_references

(

    ncbispec=(daubentonia_madagascariensis lemur_catta sapajus_appella)
    ncbiacc=(GCA_023783475.1 GCF_020740605.2 GCF_009761245.1)

    module purge
    module load hub.apps/anaconda3
    source activate ncbi_datasets

    datasets download genome accession ${ncbiacc[@]} --dehydrated --no-progressbar --include genome

    unzip ncbi_dataset.zip
    rm ncbi_dataset.zip

    datasets rehydrate --gzip --directory ./

    conda deactivate
    module purge
    module load apps/hisat2/2.1.0
    for i in 0 1 2; do

      (
          filename=$(basename ncbi_dataset/data/${ncbiacc[$i]}/*.fna.gz)

          cp ncbi_dataset/data/${ncbiacc[$i]}/*.fna.gz ${ncbispec[$i]^}_${filename/fna.gz/fa.gz}

          ## genomes seem to need to be unzipped for hisat2 indexing but can be stored zipped after that
          zcat ${ncbispec[$i]^}_${filename/fna.gz/fa.gz} > ${ncbispec[$i]}_tmp.fa

          ## build genome index
          hisat2-build -p 5 ${ncbispec[$i]}_tmp.fa ${ncbispec[$i]}_index
          
          rm ${ncbispec[$i]}_tmp.fa
      ) &
      
    done
    
    wait

    rm -r ncbi_dataset
    
) &

module purge
module load apps/hisat2/2.1.0
for spec in callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii; do

  (
      ## download genome.
      ## if a 'primary assembly' file exists, it means the 'toplevel' files contain alternate chromosomes which will artificially influence multimapping rates, so we want to prioritize primary assemblies if they exist
      wget -e robots=off -r -N -l1 -nd -A '*.dna.primary_assembly.fa.gz,*.dna.toplevel.fa.gz' https://ftp.ensembl.org/pub/release-109/fasta/${spec}/dna/
      if [ -f  ${spec^}.*.primary_assembly.fa.gz ]; then
        rm ${spec^}.*.toplevel.fa.gz
        zcat ${spec^}.*.primary_assembly.fa.gz > ${spec}_tmp.fa
      else
        rm ${spec^}.*.primary_assembly.fa.gz
        zcat ${spec^}.*.toplevel.fa.gz > ${spec}_tmp.fa
      fi
      
      ## build genome index
      hisat2-build -p 5 ${spec}_tmp.fa ${spec}_index
      
      rm ${spec}_tmp.fa
  ) &

done
