#!/bin/bash
#SBATCH --job-name=01_find_transcripts
#SBATCH --mem=175G
#SBATCH --time=6-00:00:00
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --ntasks-per-node=24
#SBATCH --output=outputs/primates_20230309_all/01_find_transcripts/logs/01_find_transcripts_%a.log
#SBATCH --array=0-8

species=(callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii daubentonia_madagascariensis lemur_catta sapajus_appella)

spec=${species[$SLURM_ARRAY_TASK_ID]}

ref_dir=outputs/primates_20230309_all/00_references
out_dir=outputs/primates_20230309_all/01_find_transcripts

mkdir -p work/tmp/

module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda activate stringtie_2.2.1

superreads.pl <(zcat raw_data/20221215_primate_allometry/fastqs/${spec}*_R1_001.fastq.gz) <(zcat raw_data/20221215_primate_allometry/fastqs/${spec}*_R2_001.fastq.gz) /shares/omicshub/apps/anaconda3/envs/masurca -l work/tmp/${spec}_superreads.fastq -u work/tmp/${spec}_unassembled_

## map reads to genome
module purge
module load apps/hisat2/2.1.0
module load apps/samtools/1.3.1

hisat2 -p 24 -x ${ref_dir}/${spec}_index \
  --no-discordant \
  --no-mixed \
  --dta-cufflinks \
  -1 work/tmp/${spec}_unassembled_R1.fq.gz \
  -2 work/tmp/${spec}_unassembled_R2.fq.gz \
  -U work/tmp/${spec}_superreads.fastq |
  samtools view -@ 4 -bS - > work/tmp/${spec}_unsorted.bam

samtools sort -T work/tmp/${spec} -m 7G -@ 4 -o ${out_dir}/${spec}.bam work/tmp/${spec}_unsorted.bam

rm work/tmp/${spec}_superreads.fastq work/tmp/${spec}_unassembled_* work/tmp/${spec}_unsorted.bam

## create transcriptome from reads and genome
module purge
module load hub.apps/anaconda3/2020.11
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda activate stringtie_2.2.1

stringtie -p 24 \
  -l ${spec} \
  -o ${out_dir}/${spec}_transcripts.gtf \
  ${out_dir}/${spec}.bam

## predict coding and protein sequences from transcripts
conda deactivate
conda activate transdecoder

gtf_genome_to_cdna_fasta.pl ${out_dir}/${spec}_transcripts.gtf <(zcat ${ref_dir}/${spec^}.*.fa.gz) > ${out_dir}/${spec}_transcripts.fa

gtf_to_alignment_gff3.pl ${out_dir}/${spec}_transcripts.gtf > ${out_dir}/${spec}_transcripts.gff3

TransDecoder.LongOrfs -S -t ${out_dir}/${spec}_transcripts.fa --output_dir ${out_dir}/${spec}_transdecoder

cd ${out_dir}
TransDecoder.Predict -t ${spec}_transcripts.fa --single_best_only --output_dir ${spec}_transdecoder

## find longest isoform per gene
grep '^>' ${spec}_transcripts.fa.transdecoder.pep |
  sort -V -k 1,1 |
    awk 'NR == 1 {
          split($5, a, ":")
          curlen = a[2]
          split($1, b, ".")
          curgene = b[1]"."b[2]
          sub(">", "", $1)
          longest = $1
        }
        NR > 1 {
          split($5, a, ":")
          curlen = a[2]
          split($1, b, ".")
          gene = b[1]"."b[2]
          if(gene != curgene) {
            curgene = gene
            curlen = len
            print longest
            sub(">", "", $1)
            longest = $1
          }
          else if(len > curlen) {
            sub(">", "", $1)
            longest = $1
            curlen = len
          }
        }
        END {
          print longest
        }' > ${spec}_longest_transcript_per_gene.txt

## filter fasta to only keep longest isoforms
awk 'NR == FNR {
       a[$1]++
       next
     }
     $1 in a {
       sub(/\n$/, "")
       print ">"$0
     }' ${spec}_longest_transcript_per_gene.txt RS='(^|\n)>' ${spec}_transcripts.fa.transdecoder.pep > ${spec}_longest_peptide_per_gene.pep
