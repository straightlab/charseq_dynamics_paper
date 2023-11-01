#!/bin/bash
export LC_ALL=C
REF="$1"
FQOUT="$2"
OUTFOLDER="$(dirname ${2})"
THREADS=$3
mkdir -p "$OUTFOLDER"/tmp
shift 3


join --check-order -t $'\t' <(cut -d $'\t' -f1-1 "$REF") <(gunzip -c "$@" | awk '{if (NR%4==1){print substr($1,2);} else{print $0;}}' | paste - - - - | sort -T "${OUTFOLDER}/tmp" -S1G --parallel=$THREADS -t $'\t' -k1,1) | awk -F $'\t' 'BEGIN{OFS="\n"}{print "@"$1, $2, $3, $4}' | gzip -c > "$FQOUT"

rm -rf "$OUTFOLDER"/tmp