#!/bin/bash
#SBATCH --job-name=01_generate_counts
#SBATCH --mem=170G
#SBATCH --time=6-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks=20
#SBATCH --output=outputs/primates_20230304/01_generate_counts/logs/01_generate_counts_%a.log
#SBATCH --array=0-8

species=(callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii daubentonia_madagascariensis lemur_catta sapajus_appella)

spec=${species[$SLURM_ARRAY_TASK_ID]}

ref_dir=outputs/primates_20230304/00_references
out_dir=outputs/primates_20230304/01_generate_counts

fwds=(raw_data/20221215_primate_allometry/fastqs/${spec}*_R1_001.fastq.gz)
revs=(raw_data/20221215_primate_allometry/fastqs/${spec}*_R2_001.fastq.gz)

## map reads to genome
module purge
module load apps/hisat2/2.1.0
module load apps/samtools/1.3.1

mkdir ${out_dir}/${spec}
hisat2 -p 20 -x ${ref_dir}/${spec}_index \
  --dta-cufflinks \
  -1 $(IFS=,; echo "${fwds[*]}") \
  -2 $(IFS=,; echo "${revs[*]}") |
  samtools view -@ 20 -bS - |
  samtools sort -T ${out_dir}/${spec}/tmp -m 5G -@ 20 - > ${out_dir}/${spec}/${spec}.bam

## create transcriptome from reads and genome
module purge
module load apps/stringtie/1.3.4b

stringtie -p 20 \
  -l ${spec} \
  -o ${out_dir}/${spec}/${spec}_transcripts.gtf \
  ${out_dir}/${spec}/${spec}.bam

module purge
module load hub.apps/bedtools/2.30.0

zcat ${ref_dir}/${spec^}_*.fa.gz > ${out_dir}/${spec}/${spec^}_tmp.fa
head -2 ${out_dir}/${spec}/${spec}_transcripts.gtf > ${out_dir}/${spec}/${spec}_transcripts_named.gtf
paste $(tail -n +3 ${out_dir}/${spec}/${spec}_transcripts.gtf |
        awk -F ';' '{print $2}' |
        awk -F '"' '{print $2}') $(tail -n +3 ${out_dir}/${spec}/${spec}_transcripts.gtf |
                                   cut -f 2-) |
awk '$3 == "transcript"' ${out_dir}/${spec}/${spec}_transcripts.gtf >> ${out_dir}/${spec}/${spec}_transcripts_named.gtf

bedtools getfasta \
  -fi ${out_dir}/${spec}/${spec^}_tmp.fa \
  -bed ${out_dir}/${spec}/${spec}_transcripts_named.gtf \
  -fo ${out_dir}/${spec}_transcripts.fasta

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

    salmon quant --index ${ref_dir}/${spec}_salmon \
      --threads 20 \
      --libType A \
      -1 ${fwd} \
      -2 ${rev} \
      --output ${out_dir}/${spec}_${sample}_salmon

done
