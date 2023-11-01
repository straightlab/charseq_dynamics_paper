#!/bin/bash
export LC_ALL=C
FQ="$1"
CONTAMINANTS_BAM="$2"
FQOUT="$3"
OUTDIR="$4"
THREADS=$5

OUTDIR_NOCON="${OUTDIR%/}.decon/"
OUTDIR_CON="${OUTDIR%/}.contaminants/"
OUTDIR_TMP="${OUTDIR%/}.decon/tmp"

mkdir -p "$OUTDIR_NOCON"
mkdir -p "$OUTDIR_CON"
mkdir -p "$OUTDIR_TMP"

FQ_IN=$(basename "$FQ")
FQ_CON="gzip >${OUTDIR_CON}${FQOUT}.fastq.gz"
FQ_NOCON="gzip >${OUTDIR_NOCON}${FQOUT}.fastq.gz"



cat /dev/null | gzip -c>"${OUTDIR_CON}${FQOUT}.fastq.gz"
cat /dev/null | gzip -c>"${OUTDIR_NOCON}${FQOUT}.fastq.gz"


#F=4 not set bit = mapped to con, means CONTAMINANT
echo "FQ_IN=${FQ} | FQ_CON=${FQ_CON} | FQ_NOCON=${FQ_NOCON} | CONTAMINANTS_BAM=${CONTAMINANTS_BAM}" 
join --check-order -t $'\t' -a 1 <(zcat "$FQ" | awk '{if (NR%4==1){print $1;} else{print $0;}}' | paste - - - - | sort -T "${OUTFOLDER}/tmp" -S1G --parallel=$THREADS -t $'\t' -k1,1) <(samtools view -F4 "$CONTAMINANTS_BAM" | cut -f 1-2 | sed -e 's/^/@/') | awk -F $'\t' -v out_con="$FQ_CON" -v out_nocon="$FQ_NOCON" 'BEGIN{OFS="\n"}{if (NF==5){print $1, $2, $3, $4 | out_con;} else {print $1, $2, $3, $4 | out_nocon;}}'

