#!/bin/bash
# cd "$1"
touch "$1"/"$2".stats
rm "$1"/"$2".stats

echo "Filtering pairs"
cat "$1"/final.pairs | grep "$2" | cut -d " " -f 1 > "$1"/"$2".pairs

echo "Making bams"
java -Xmx4g -jar /share/PI/astraigh/software/picard.jar FilterSamReads I="$1"/final.dna.bam O="$1"/"$2".dna.bam READ_LIST_FILE="$1"/"$2".pairs FILTER=includeReadList SORT_ORDER=coordinate

echo "Indexing"
samtools index "$1"/"$2".dna.bam

echo "Counting"
samtools view -c -q 10 "$1"/"$2".dna.bam >>"$1"/"$2".stats
# echo -e "" >>"$2".stats
samtools view -c -q 10 "$1"/"$2".dna.bam $3 >>"$1"/"$2".stats
# echo -e "" >>"$2".stats
samtools view -c -q 10 "$1"/"$2".dna.bam $3:$4 >>"$1"/"$2".stats


# 73820651 73852753

# 73733238 73933374

# XIST chrX 73720651-73952753
#73,820,651-73,851,091
