
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

#1 fq1, #2 fq2 --> use to pick up names, #3 bam1, #4 bam2, --> outputs bam1_dedup, bam2_dedup

if [ $# -lt 4 ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
  exit
fi

  PICARD="picard.jar"
  OUTSUFFIX=".dedup"
  OUT1=""
  OUT2=""
  SUBS=0
  XMX="3G"
  KEEPTEMP=0
  TRIMOPTIONS=""
  while [[ $# -gt 0 ]]; do
      case "$1" in
  			-picard ) shift; PICARD="$1"; shift ;;
        -out1 ) shift; OUT1="$1"; shift ;;
        -out2 ) shift; OUT2="$1"; shift ;;
        -keeptemp ) shift; KEEPTEMP=$1; shift ;;
        -subs ) shift; SUBS=$1; shift ;;
        -xmx ) shift; XMX=$1; shift ;;
        -dedup ) shift; DEDUP=$1; shift ;;
        -trimoptions ) shift; TRIMOPTIONS="$1"; shift ;;
  			* ) break ;;
  	esac
  done

  if [ $# -lt 4 ]; then
  	usage
    exit
  fi


IN1_PREFIX=""
OUT1_FQ=""
DEDUP_READS_ID=""
if [[ ${1: -3} == ".gz" ]]; then
  IN1_PREFIX2="${1%.*}"
  # IN1_EXT2=${1##*.}
  IN1_PREFIX="${IN1_PREFIX2%.*}"
  # IN1_EXT=${IN1_PREFIX2##*.}
else
  IN1_PREFIX="${1%.*}"
  # IN1_EXT=${1##*.}
fi
IN3_PREFIX="${3%.*}"
if [[ -z $OUT1 ]]; then
  OUT1_FQ="${IN1_PREFIX}${OUTSUFFIX}.fastq"
  DEDUP_READS_ID="${IN1_PREFIX}${OUTSUFFIX}.fastq.reads2keep.txt"
  OUT1="${IN3_PREFIX}${OUTSUFFIX}.bam"
else
  OUT1_PREFIX="${OUT1%.*}"
  OUT1_FQ="${OUT1_PREFIX}.fastq"
  DEDUP_READS_ID="${OUT1_PREFIX}.reads2keep.txt"
fi


IN2_PREFIX=""
OUT2_FQ=""
if [[ ${2: -3} == ".gz" ]]; then
  IN2_PREFIX2="${2%.*}"
  # IN1_EXT2=${1##*.}
  IN2_PREFIX="${IN2_PREFIX2%.*}"
  # IN1_EXT=${IN1_PREFIX2##*.}
else
  IN2_PREFIX="${2%.*}"
  # IN1_EXT=${1##*.}
fi
IN4_PREFIX="${4%.*}"
if [[ -z $OUT2 ]]; then
  OUT2_FQ="${IN2_PREFIX}${OUTSUFFIX}.fastq"
  OUT2="${IN4_PREFIX}${OUTSUFFIX}.bam"
else
  OUT2_PREFIX="${OUT2%.*}"
  OUT2_FQ="${OUT2_PREFIX}.fastq"
fi

echo "$OUT1_FQ"
echo "$OUT2_FQ"
echo "$DEDUP_READS_ID"
echo "$OUT1"
echo "$OUT2"

echo "Deduping paired fastq files";
echo "Running command : clumpify.sh in=${1} in2=${2} dedupe=t out=${OUT1_FQ} out2=${OUT2_FQ} subs=${SUBS} -Xmx${XMX}"
clumpify.sh in="$1" in2="$2" dedupe=t out="$OUT1_FQ" out2="$OUT2_FQ" subs=$SUBS -Xmx$XMX

if [[ -z "$TRIMOPTIONS" ]]; then
  echo "No trimming applied"
else
  echo "Trimming deduplicated reads"
  OUT1_FQ_2="${OUT1_PREFIX}.TMP.fastq"
  STATS="${OUT1_PREFIX}.trimstats.txt"
  bbduk.sh in="$OUT1_FQ" out="$OUT1_FQ_2" stats="$STATS" -Xmx${XMX} "$TRIMOPTIONS"
  mv "$OUT1_FQ_2" "$OUT1_FQ"
fi

echo "$TRIMOPTIONS"

echo "Finding read ID of duplicates to remove"
awk 'NR%4==1 && NF==1' "$OUT1_FQ" | cut -c2- >"$DEDUP_READS_ID"

# k=15 mink=8 hdist=1 hdist2=0  cardinality=t cardinalityout=t loglogk=25 minlength=15 ktrim=r rcomp=f literal=AAAAAAAAAAAAAAA stats=stats_.txt out=pass_.fq outm=fail_.fq overwrite=t
echo "Deduplicating in1"
java -Xmx$XMX -jar "$PICARD" FilterSamReads I="$3" O="$OUT1" READ_LIST_FILE="$DEDUP_READS_ID" FILTER=includeReadList SORT_ORDER=queryname
echo "Deduplicating in2"
java -Xmx$XMX -jar "$PICARD" FilterSamReads I="$4" O="$OUT2"  READ_LIST_FILE="$DEDUP_READS_ID" FILTER=includeReadList SORT_ORDER=queryname

if [[ $KEEPTEMP -lt 1 ]]; then
echo "Cleaning up temporary files"
rm "$DEDUP_READS_ID"
rm "$OUT1_FQ"
rm "$OUT2_FQ"
fi
