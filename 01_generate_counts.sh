#!/bin/bash
#SBATCH --job-name=01_generate_counts_ensembl
#SBATCH --mem=170G
#SBATCH --time=6-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks=20
#SBATCH --output=outputs/primates_20230308_ensembl/01_generate_counts/logs/01_generate_counts_%a.log
#SBATCH --array=0-5

species=(callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii)
spec=${species[$SLURM_ARRAY_TASK_ID]}

ref_dir=outputs/primates_20230308_ensembl/00_references
out_dir=outputs/primates_20230308_ensembl/01_generate_counts

module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda activate salmon
for fwd in raw_data/20221215_primate_allometry/fastqs/${spec}*_R1_001.fastq.gz; do

    # extract sample name from file name
    rev=${fwd%_R1_001.fastq.gz}_R2_001.fastq.gz
    sample=$(basename ${fwd})
    sample=${sample#*_}
    sample=${sample#*_}
    for i in {1..4}; do
      sample=${sample%_*}
    done

    salmon quant --index ${ref_dir}/${spec}_index \
      --threads 20 \
      --libType A \
      -1 ${fwd} \
      -2 ${rev} \
      --output ${out_dir}/${spec}_${sample}_salmon

done
