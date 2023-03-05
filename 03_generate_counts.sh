#!/bin/bash
#SBATCH --job-name=03_generate_counts
#SBATCH --mem=170G
#SBATCH --time=6-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks=20
#SBATCH --output=outputs/primates_20230304/01_generate_counts/logs/03_generate_counts_%a.log
#SBATCH --array=0-8

species=(callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii daubentonia_madagascariensis lemur_catta sapajus_appella)

spec=${species[$SLURM_ARRAY_TASK_ID]}

in_dir=outputs/primates_20230304/01_find_transcripts
out_dir=outputs/primates_20230304/03_generate_counts

fwds=(raw_data/20221215_primate_allometry/fastqs/${spec}*_R1_001.fastq.gz)

## run orthofinder with all our inferred transcripts but also all the ensembl cds files - only keep genes that have at least one human ensembl identifier in their family tree, so that annotations can be applied to entire tree. maybe collapse all within-species clusters of 'genes' into a single gene with multiple transcripts, before identifying putative '1:1' orthologs. shouldn't matter for our case if theres a duplication; we want to know if there are more transcripts. duplications that aren't 'monophyletic' will still be a problem.

## only run salmon on transcripts that have human ensembl family members. don't worry about potential pooling until at the tx2gene step for tximport. use reference genomes in indexing for 'decoys' (might be particularly imiportant since the reference transcripts will be incomplete, so we wouldn't want things looking like they align when they actually came from parts of the genome that aren't in the ref transcripts)

## build transcriptome index
module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda activate salmon

salmon index --index ${out_dir}/${spec}_salmon --transcripts ${out_dir}/${spec}_transcripts.fasta

for fwd in ${fwds[@]}; do

    # extract sample name from file name
    rev=${fwd%_R1_001.fastq.gz}_R2_001.fastq.gz
    sample=$(basename ${fwd})
    sample=${sample#*_}
    sample=${sample#*_}
    for i in {1..4}; do
      sample=${sample%_*}
    done

    salmon quant --index ${out_dir}/${spec}_salmon \
      --threads 20 \
      --libType A \
      -1 ${fwd} \
      -2 ${rev} \
      --output ${out_dir}/${spec}_${sample}_salmon

done
