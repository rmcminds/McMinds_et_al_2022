#!/bin/bash
#SBATCH --job-name=01_find_transcripts
#SBATCH --mem=175G
#SBATCH --time=6-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks-per-node=24
#SBATCH --output=outputs/primates_20230304/01_find_transcripts/logs/01_find_transcripts_%a.log
#SBATCH --array=0-8

species=(callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii daubentonia_madagascariensis lemur_catta sapajus_appella)

spec=${species[$SLURM_ARRAY_TASK_ID]}

ref_dir=outputs/primates_20230304/00_references
out_dir=outputs/primates_20230304/01_find_transcripts

fwds=(raw_data/20221215_primate_allometry/fastqs/${spec}*_R1_001.fastq.gz)
revs=(raw_data/20221215_primate_allometry/fastqs/${spec}*_R2_001.fastq.gz)

## map reads to genome
module purge
module load apps/hisat2/2.1.0
module load apps/samtools/1.3.1

mkdir -p work/tmp/
hisat2 -p 24 -x ${ref_dir}/${spec}_index \
  --no-discordant \
  --no-mixed \
  --dta-cufflinks \
  -1 $(IFS=,; echo "${fwds[*]}") \
  -2 $(IFS=,; echo "${revs[*]}") |
  samtools view -@ 4 -bS - |
  samtools sort -T work/tmp/${spec} -m 7G -@ 4 - > ${out_dir}/${spec}.bam

## create transcriptome from reads and genome
module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda activate stringtie_2.2.1

stringtie -p 24 \
  --conservative \
  -l ${spec} \
  -o ${out_dir}/${spec}_transcripts.gtf \
  ${out_dir}/${spec}.bam

## predict coding and protein sequences from transcripts
conda deactivate
conda activate transdecoder

gtf_genome_to_cdna_fasta.pl ${out_dir}/${spec}_transcripts.gtf <(zcat ${ref_dir}/${spec^}.*.fa.gz) > ${out_dir}/${spec}_transcripts.fa

gtf_to_alignment_gff3.pl ${out_dir}/${spec}_transcripts.gtf > ${out_dir}/${spec}_transcripts.gff3

TransDecoder.LongOrfs -t ${out_dir}/${spec}_transcripts.fa --output_dir ${out_dir}/${spec}_transdecoder
