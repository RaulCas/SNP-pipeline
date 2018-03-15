#!/bin/bash -l

# Requires sorted .bam file as imput (ie, RNAseq reads mapped with STAR)

#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=16G
#SBATCH --time=0-01:00:00     
#SBATCH --output=snps.stdout
#SBATCH --mail-type=ALL
#SBATCH --job-name="ASE"

module load picard
module load gatk
module load bowtie
module load samtools

cd $SLURM_SUBMIT_DIR

######## First add RG groups with picard ######

picard AddOrReplaceReadGroups \
I=/Aligned.sortedByCoord.out.bam \
O=/Aligned.sortedByCoord_RG.out.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20

 
# Index the sorted bam

samtools index Aligned.sortedByCoord_RG.out.bam

# Run ASEReadCounter

java -jar /opt/linux/centos/7.x/x86_64/pkgs/GATK/3.6/GenomeAnalysisTK.jar -R PleosPC15_2_Assembly_scaffolds.fasta -T ASEReadCounter -o $SLURM_SUBMIT_DIR/out.csv -I /Aligned.sortedByCoord_RG.out.bam -sites resultsVCFXMLsortedcleaned.vcf -U ALLOW_N_CIGAR_READS








