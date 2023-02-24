#!/bin/bash
#
#SBATCH --workdir=${workdir}/${species}/
#SBATCH --job-name=${sample}
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=120000
#SBATCH -t 0-72:00:00
#SBATCH -o ${sample}.out
#SBATCH -e ${sample}.err

## this is a script template based on a set of individual scripts written for each sample independently. To replicate original code, assign the appropriate directory and name variables for each sample and species
#species={species_name}
#genome_dir={reference_location}
#genome={reference_gtf_name}
#genome_fa={reference_fa_name}
#species_index={location_of_indexed_genome_output}
#in1={forward_reads_fastq}
#in2={reverse_reads_Fastq}

module purge
module load apps/samtools/1.3.1
module load apps/hisat2/2.1.0

# build index (only needed once per species)
hisat2-build -p $SLURM_NTASKS ${genome_dir} ${species_index}

mkdir ${sample}
cd ${sample}

# align reads to reference genome with HISAT2, and sort and convert output to bam
hisat2 -p 20 -x ${species_index} --dta-cufflinks -1 ${in1} -2 ${in2} --summary-file ${sample}_summary.txt | samtools view -@ 20 -bS - | samtools sort -@ 20 - > ${sample}.bam

# calculate transcript abundances with StringTie
module purge
module load apps/stringtie/1.3.4b

stringtie -p 20 -B -G ${genome_dir}/${genome} -o ${sample}_transcripts.gtf -A ${sample}_abundances.tsv -l ${sample} ${sample}.bam

