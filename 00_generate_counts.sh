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
module load apps/bowtie/2.3.4.1
module load apps/samtools/1.3.1
module load apps/hisat2/2.1.0
module load apps/stringtie/1.3.4b

# build index (only needed once per species)
hisat2-build -p $SLURM_NTASKS ${genome_dir} ${species_index}

# Setup Directory
mkdir ${sample}
cd ${sample}

# Start HISAT2
hisat2 -p 20 -x ${species_index} --dta-cufflinks -1 ${in1} -2 ${in2} -S ./${sample}.sam --summary-file ${sample}_summary.txt

# Start samtools
samtools view -@ 20 -bS ${sample}.sam > ${sample}_unsorted.bam
samtools sort -@ 20 ${sample}_unsorted.bam -o ${sample}.bam

rm *_unsorted.bam
rm *.sam

# Start Stringtie
stringtie -p 20 -G ${genome_dir}/${genome} -o ${sample}_transcripts.gtf -A ${sample}_abundances.tsv -l ${sample} ${sample}.bam

module purge
module load apps/subread/1.6.3

# Start featurecounts
featureCounts -T 20 -p -t gene -g gene_id -a ${genome_dir}/${genome} -G ${genome_dir}/${genome_fa} -J -o ${sample}_counts.txt ${sample}.bam

echo "Finish - `date`"

## run once after all samples have been processed independently:

mkdir inputs/reference_counts_cleaned
cp ../Seq/*/*_counts.txt inputs/reference_counts_cleaned

