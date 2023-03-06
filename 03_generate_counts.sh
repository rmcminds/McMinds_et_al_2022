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

## run salmon on full stringtie transcripts to ensure all genome locations are available for mapping. then just filter the quantified genes to 1:1 orthologs including the reference human ensembl, identified via orthofinder (and thus implicilty filtering to only protein coding genes). perhaps collapse all putative gene duplicates within a single species because they could be different isoforms of the same gene, and our rough pipeline didn't pick that up?

## build transcriptome index
module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda activate salmon

salmon index --index ${out_dir}/${spec}_salmon_index --transcripts ${in_dir}/${spec}_transcripts.fasta

for fwd in ${fwds[@]}; do

    # extract sample name from file name
    rev=${fwd%_R1_001.fastq.gz}_R2_001.fastq.gz
    sample=$(basename ${fwd})
    sample=${sample#*_}
    sample=${sample#*_}
    for i in {1..4}; do
      sample=${sample%_*}
    done

    salmon quant --index ${out_dir}/${spec}_salmon_index \
      --threads 20 \
      --libType A \
      -1 ${fwd} \
      -2 ${rev} \
      --output ${out_dir}/${spec}_${sample}_salmon

done
