#!/bin/bash
export LC_ALL=C
genome_bam="$1"
tr_bam="$2"
tr_bam_out="$3"

echo "Finding shared reads"
COMMON_READS="${genome_bam}_EXCLUSIVE.txt"

join --check-order -v 1 <(samtools view "$genome_bam" | cut -f1 | uniq) <(samtools view "$tr_bam" | cut -f1 | uniq) >"$COMMON_READS"
#
echo "Filtering reads unique to genome"
picard FilterSamReads I="$genome_bam" O="$genome_bam"_exclusive.bam READ_LIST_FILE="$COMMON_READS" FILTER=includeReadList SORT_ORDER=queryname

echo "Merging reads unique to genome and reads with tr match"
picard MergeSamFiles I="$genome_bam"_exclusive.bam I="$tr_bam" O="$tr_bam_out" MERGE_SEQUENCE_DICTIONARIES=true SORT_ORDER=queryname

