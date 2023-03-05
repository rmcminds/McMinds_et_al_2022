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
## convert gtf to bed while extracting transcript names which can be linked to gene names. bedtools directly from gtf had a single nucleotide difference in start vs bedtools from bed with same numbers, so this is adjusted ($4-1)
head -2 ${out_dir}/${spec}/${spec}_transcripts.gtf > ${out_dir}/${spec}/${spec}_transcripts.bed
tail -n +2 ${out_dir}/${spec}/${spec}_transcripts.gtf | awk '$3 == "transcript" {split($0,a,";"); split(a[2],b,"\""); print $1"\t"$4-1"\t"$5"\t"b[2]}' >> ${out_dir}/${spec}/${spec}_transcripts.bed

## find longest trancsript of each gene
awk '' ${out_dir}/${spec}/${spec}_transcripts.bed > ${out_dir}/${spec}/${spec}_longest.txt

sort -k 4 ${out_dir}/${spec}/${spec}_transcripts.bed | awk 'BEGIN {curgene='initiate'; longest='initiate'; curlen=0} NR>2 {split($4,a,"."); gene=a[1]"."a[2]; len=$3-$2; if(gene != curgene) {curgene=gene; curlen=len; print longest; longest=$4} else if(len > curlen) {longest=$4;curlen=len}} END {print longest}' > ${out_dir}/${spec}/${spec}_longest.txt

bedtools getfasta \
  -name \
  -fi ${out_dir}/${spec}/${spec^}_tmp.fa \
  -bed ${out_dir}/${spec}/${spec}_transcripts.bed \
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
