#!/bin/bash
FASTA=0
NALIGN=3
OUTFMT="6 qseqid stitle qlen length pident qcovs"
while [[ $# -gt 0 ]]; do
    case "$1" in
			-f ) FASTA=1; shift ;;
			-n ) shift; NALIGN=$1; shift ;;
			-outformat ) shift; OUFMT="$1"; shift ;;
			* ) break ;;
	esac
done

x="$1"
x2=$(basename "$x")
y="${x2%.*}"
if [ $NALIGN -gt 0 ]; then
	join -a 1 -t $'\t' <(unzip -c "$1" "$y"/fastqc_data.txt | awk '/^>>Overrepresented/{flag=1;getline; next;}/>>END_MODULE/{flag=0}flag' | awk 'BEGIN{OFS=","}{printf("%03d,",NR); print $2, $3"\t"$1}') <(blastn -query <(unzip -c "$1" "$y"/fastqc_data.txt | awk '/^>>Overrepresented/{flag=1;getline; next;}/>>END_MODULE/{flag=0}flag' | awk 'BEGIN{OFS=","}{printf(">%03d,",NR); print $2, $3"\n"$1}') -db "$2" -num_alignments $NALIGN -outfmt "$OUTFMT")
else
	unzip -c "$1" "$y"/fastqc_data.txt | awk '/^>>Overrepresented/{flag=1;getline; next;}/>>END_MODULE/{flag=0}flag' | awk 'BEGIN{OFS=","}{printf("%03d,",NR); print $2, $3"\t"$1}'
fi


# join -a 1 t2.txt t1.txt
