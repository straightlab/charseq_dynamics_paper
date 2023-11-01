#!/bin/bash

function usage {
    echo "usage: bams2chibam rna_in dna_in out rna_filter dna_filter picard_path rgid rglb rgsm"
    echo ""
    echo "Combines RNA and DNA alignment bam files into single "
    echo ""
    echo "  bridgeF      forward bridge sequence"
    echo "  bridgeF      reverse bridge sequence"
    echo "  fastq_filaname      a fasqt.gz file of reads"
    echo "  out_ub_filename     name of output file where all the reads with a single bridge are written"
    echo "  out_ub_filename     name of output file where the position and orientation of the unique bridges are written"
    echo "  out_mb_filename     name of output file where all the reads with multiple bridges are written"
    echo "EXAMPLES"
    echo "  bash bridge_count.sh ACCGGCGTCCAAG CTTGGACGCCGGT myreads.fastq ub.fastq.gz pos.csv mb.fastq.gz"
    exit 1
}

if [ $# -lt 2 ]; then
	usage
else
RNA_FILTER="$4"
DNA_FILTER="$5"
PICARD="$6"
# RGID="$7"
# RGLB="$8"
# RGSM="$9"
# first apply desired filtering, parallel branches for r and d
echo "Filtering RNA reads"
samtools view -b $RNA_FILTER "$1" > "$1".filtered.bam
echo "Filtering DNA reads"
samtools view -b $DNA_FILTER "$2" > "$2".filtered.bam
echo "Finding shared reads"
bedtools bamtobed -i "$1".filtered.bam | cut -f4 | sort | uniq > "$1".filtered.uniq.txt
bedtools bamtobed -i "$2".filtered.bam | cut -f4 | sort | uniq > "$2".filtered.uniq.txt
join "$1".filtered.uniq.txt "$2".filtered.uniq.txt > "$1".filtered.dual.txt
echo "Cleaning up intermediate files"
rm "$1".filtered.uniq.txt
rm "$2".filtered.uniq.txt

echo "Keeping only common reads"
java -Xmx4g -jar "$PICARD" FilterSamReads I="$1".filtered.bam O="$3".rna.bam READ_LIST_FILE="$1".filtered.dual.txt FILTER=includeReadList SORT_ORDER=queryname
java -Xmx4g -jar "$PICARD" FilterSamReads I="$2".filtered.bam O="$3".dna.bam READ_LIST_FILE="$1".filtered.dual.txt FILTER=includeReadList SORT_ORDER=queryname

echo "Creating pairs file"
python bams2pairs.py "$3".rna.bam "$3".dna.bam "$3".pairs

# echo "Adding flags to DNA reads"
# python flagreads <(java -Xmx4g -jar "$PICARD" FilterSamReads I="$2".filtered.bam O=/dev/stdout READ_LIST_FILE="$1".filtered.dual.txt FILTER=includeReadList SORT_ORDER=queryname
# ) "$3".dna.bam $6 $7 $8
#
# echo "Creating chibams files"
# java -Xmx4g -jar "$PICARD" MergeSamFiles I="$3".rna.bam I="$3".dna.bam O="$3".chi.bam SORT_ORDER=queryname ASSUME_SORTED=true
#
# echo "And finally... cleaning up temporary files"
# rm "$3".rna.bam
# rm "$3".dna.bam

fi
