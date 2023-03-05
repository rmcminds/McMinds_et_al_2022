#!/bin/bash
#SBATCH --job-name=01_find_transcripts
#SBATCH --mem=170G
#SBATCH --time=6-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks=20
#SBATCH --output=outputs/primates_20230304/01_generate_counts/logs/01_find_transcripts_%a.log
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
  --conservative \
  -l ${spec} \
  -o ${out_dir}/${spec}/${spec}_transcripts.gtf \
  ${out_dir}/${spec}/${spec}.bam

module purge
module load hub.apps/bedtools/2.30.0

zcat ${ref_dir}/${spec^}_*.fa.gz > ${out_dir}/${spec}/${spec^}_tmp.fa
## convert gtf to bed while extracting transcript names which can be linked to gene names. bedtools directly from gtf had a single nucleotide difference in start vs bedtools from bed with same numbers, so this is adjusted ($4-1)
head -2 ${out_dir}/${spec}/${spec}_transcripts.gtf > ${out_dir}/${spec}/${spec}_transcripts.bed
tail -n +2 ${out_dir}/${spec}/${spec}_transcripts.gtf | awk '$3 == "transcript" {split($0,a,";"); split(a[2],b,"\""); print $1"\t"$4-1"\t"$5"\t"b[2]}' >> ${out_dir}/${spec}/${spec}_transcripts.bed

bedtools getfasta \
  -name \
  -fi ${out_dir}/${spec}/${spec^}_tmp.fa \
  -bed ${out_dir}/${spec}/${spec}_transcripts.bed \
  -fo ${out_dir}/${spec}_transcripts.fasta

## find longest trancsript of each gene
tail -n +3 ${out_dir}/${spec}/${spec}_transcripts.bed | sort -V -k 4 | awk 'NR==1 {split($4,a,"."); curgene=a[1]"."a[2]; curlen=$3-$2; longest=$4} NR>1 {split($4,a,"."); gene=a[1]"."a[2]; len=$3-$2; if(gene != curgene) {curgene=gene; curlen=len; print longest; longest=$4} else if(len > curlen) {longest=$4;curlen=len}} END {print longest}' > ${out_dir}/${spec}/${spec}_longest.txt

awk 'NR==FNR {a[$1]++; next} $4 in a' ${out_dir}/${spec}/${spec}_longest.txt ${out_dir}/${spec}/${spec}_transcripts.bed > ${out_dir}/${spec}/${spec}_transcripts_longest.bed

bedtools getfasta \
  -name \
  -fi ${out_dir}/${spec}/${spec^}_tmp.fa \
  -bed ${out_dir}/${spec}/${spec}_transcripts_longest.bed \
  -fo ${out_dir}/${spec}_transcripts_longest.fasta

## extract longest isoforms into a fasta
## run orthofinder with all our inferred transcripts but also all the ensembl cds files - only keep genes that have at least one human ensembl identifier in their family tree, so that annotations can be applied to entire tree. maybe collapse all within-species clusters of 'genes' into a single gene with multiple transcripts, before identifying putative '1:1' orthologs. shouldn't matter for our case if theres a duplication; we want to know if there are more transcripts. duplications that aren't 'monophyletic' will still be a problem.

## only run salmon on transcripts that have human ensembl family members. don't worry about potential pooling until at the tx2gene step for tximport
