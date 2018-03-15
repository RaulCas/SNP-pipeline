#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=0-00:20:00     # 0 day and 20 minutes
#SBATCH --output=vcf.stdout
#SBATCH --mail-user=inaki.etxeberria@unavarra.es
#SBATCH --mail-type=ALL
#SBATCH --job-name="script1"

module load gatk
module load ncbi-blast
module load samtools
module load picard
module load gatk
module load vcftools

cd $PBS_O_WORKDIR

#This script runs several other scripts (over Python). It's feeded with .fasta and .gff/.gtf files and as a result gives a file with an SNP per line, indicating reference and alternative alelles and their number of reads.

#For more info, read each script's manual

#######  PARAMETERS  #######

path="shared/inaki/proyecto/"

# -1- Generating MultiFasta files
p1GFForGTF_IN=$path"PleosPC9_1_GeneModels_Filteredmodels2.gff"
p1FASTA_IN=$path"PleosPC9_1_Assembly_scaffolds.fasta"
p1MultiFASTA_OUT=$path"PC9_multiFASTA.fasta"

p2GFForGTF_IN=$path"PleosPC15_2_GeneModels_FilteredModels1.gtf"
p2FASTA_IN=$path"PleosPC15_2_Assembly_scaffolds.fasta"
p2MultiFASTA_OUT=$path"PC15_multiFASTA.fasta"

# -2- Generating BLAST file


pvalue=80
evalue=0.00000000000000000001
#evalue=80


# -3- Finding reciprocal genes' SNP

p1BLAST=$path"PC9xPC15.xml"
p2BLAST=$path"PC15xPC9.xml"
snpOUT=$path"resultsVCFXML.vcf"

# -4- Ordering results
xmlSorted=$path"resultsVCFXMLsorted.vcf"
xmlSortedCleaned=$path"resultsVCFXMLsortedcleaned.vcf"

#######  PROGRAM #######
#Do NOT change anything after this line)
echo "Program running"

python $path"multiFasta.py" $p1GFForGTF_IN $p1FASTA_IN $p1MultiFASTA_OUT

python $path"multiFasta.py" $p2GFForGTF_IN $p2FASTA_IN $p2MultiFASTA_OUT

makeblastdb -in $p1MultiFASTA_OUT -dbtype nucl
makeblastdb -in $p2MultiFASTA_OUT -dbtype nucl

blastn -outfmt 5 -perc_identity $pvalue -num_alignments 1 -evalue $evalue -query $p1MultiFASTA_OUT -db $p2MultiFASTA_OUT -out $p1BLAST
blastn -outfmt 5 -perc_identity $pvalue -num_alignments 1 -evalue $evalue -query $p2MultiFASTA_OUT -db $p1MultiFASTA_OUT -out $p2BLAST

#blastn -outfmt 5 -max_hsps 1 -perc_identity $pvalue -num_alignments 1 -evalue $evalue -query $p1MultiFASTA_OUT -subject $p2MultiFASTA_OUT -out $p1BLAST
#blastn -outfmt 5 -max_hsps 1 -perc_identity $pvalue -num_alignments 1 -evalue $evalue -query $p2MultiFASTA_OUT -subject $p1MultiFASTA_OUT -out $p2BLAST

python $path"snpsVCFXML.py" $p2BLAST $p1BLAST > $snpOUT

vcf-sort $snpOUT > $xmlSorted

python $path"deleteDuplicates.py" $xmlSorted > $xmlSortedCleaned

rm $path"resultsVCFsortedcleaned.vcf.idx"
rm $path"*.nhr"
rm $path"*.nin"
rm $path"*.nsq"

echo "Program ended"
echo "Final file is " $xmlSortedCleaned
