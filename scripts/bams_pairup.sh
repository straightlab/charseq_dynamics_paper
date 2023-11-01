#!/bin/bash
export LC_ALL=C
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

  PICARD=""
  RNA_FILTER=""
  DNA_FILTER=""
  VAL=""
  while [[ $# -gt 0 ]]; do
      case "$1" in
  			-picard ) shift; PICARD="$1"; shift ;;
        -rna-filter ) shift; RNA_FILTER="$1"; shift ;;
  			-dna-filter ) shift; DNA_FILTER="$1"; shift ;;
  			* ) break ;;
  	esac
  done
  VAL=$(echo -n "$RNA_FILTER""$DNA_FILTER" | md5sum | cut -d " " -f 1)

echo "Finding shared reads"
COMMON_READS=$1"$VAL""_COMMON".txt
echo "join <(samtools view ""$RNA_FILTER"" ""$1"" | cut -f1 | sort | uniq) <(samtools view ""$DNA_FILTER"" ""$2"" | cut -f1 | sort | uniq)"
# bedtools bamtobed -i "$2" | cut -f4 | sort | uniq > "$2".uniq.txt
join --check-order <(samtools view ${RNA_FILTER} "$1" | cut -f1 | sort | uniq) <(samtools view ${DNA_FILTER} "$2" | cut -f1 | sort | uniq) >"$COMMON_READS"
#
echo "Keeping only common reads"
java -Xmx16g -jar "$PICARD" FilterSamReads I="$1" O="$3" READ_LIST_FILE="$COMMON_READS" FILTER=includeReadList SORT_ORDER=queryname
java -Xmx16g -jar "$PICARD" FilterSamReads I="$2" O="$4" READ_LIST_FILE="$COMMON_READS" FILTER=includeReadList SORT_ORDER=queryname

echo "Cleaning up temporary files"
rm "$COMMON_READS"

fi
