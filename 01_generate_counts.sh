#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=100G
#SBATCH --job-name=01_generate_counts
#SBATCH --output=outputs/primates_20230224/01_generate_counts/logs/01_generate_counts_%a.log
#SBATCH --array=0-5

species=(callithrix_jacchus homo_sapiens macaca_mulatta microcebus_murinus papio_anubis pongo_abelii)
spec=${species[$SLURM_ARRAY_TASK_ID]}

ref_dir=outputs/primates_20230224/00_references
out_dir=outputs/primates_20230224/01_generate_counts

module purge
module load apps/samtools/1.3.1
module load apps/hisat2/2.1.0
module load apps/stringtie/1.3.4b

for fwd in raw_data/20221215_primate_allometry/fastqs/*_${spec}*_R1_001.fastq.gz; do

    rev=${fwd%_R1_001.fastq.gz}_R2_001.fastq.gz
    base=$(basename ${fwd})
    sample=${base%%_*}

    # align reads to reference genome with HISAT2, and sort and convert output to bam
    hisat2 -p 20 -x ${ref_dir}/${spec}_index --dta-cufflinks -1 ${fwd} -2 ${rev} --summary-file ${out_dir}/${sample}_summary.txt | samtools view -@ 20 -bS - | samtools sort -@ 20 - > ${out_dir}/${sample}.bam

    # calculate transcript abundances with StringTie
    stringtie -p 20 -B -G ${ref_dir}/${spec^}.*.gff3 -o ${out_dir}/${sample}_transcripts.gtf -A ${out_dir}/${sample}_abundances.tsv -l ${sample} ${out_dir}/${sample}.bam

done
