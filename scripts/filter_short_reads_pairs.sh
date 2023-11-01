#!/bin/bash

LENGTH_THRESHOLD=$1
OUTDIR="$2"
RNA="$3"
DNA="$4"

OUT_LONG="${OUTDIR%/}/long/"
OUT_SHORTR="${OUTDIR%/}/shortR/"
OUT_SHORTD="${OUTDIR%/}/shortD/"
OUT_SHORTRD="${OUTDIR%/}/shortRD/"
OUT_STATS="${OUTDIR%/}/filtering.stats.txt"


mkdir -p "$OUT_LONG"
mkdir -p "$OUT_SHORTR"
mkdir -p "$OUT_SHORTD"
mkdir -p "$OUT_SHORTRD"

echo -e "total,long,shortR,shortD,shortRD" >"$OUT_STATS"

cat /dev/null | gzip -c>${OUT_LONG}dna.fastq.gz
cat /dev/null | gzip -c>${OUT_LONG}rna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTR}dna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTR}rna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTD}dna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTD}rna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTRD}dna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTRD}rna.fastq.gz

OUT_LONG_R="gzip >${OUT_LONG}rna.fastq.gz"
OUT_LONG_D="gzip >${OUT_LONG}dna.fastq.gz"
OUT_SHORTR_R="gzip >${OUT_SHORTR}rna.fastq.gz"
OUT_SHORTR_D="gzip >${OUT_SHORTR}dna.fastq.gz"
OUT_SHORTD_R="gzip >${OUT_SHORTD}rna.fastq.gz"
OUT_SHORTD_D="gzip >${OUT_SHORTD}dna.fastq.gz"
OUT_SHORTRD_R="gzip >${OUT_SHORTRD}rna.fastq.gz"
OUT_SHORTRD_D="gzip >${OUT_SHORTRD}dna.fastq.gz"


echo "RNA:${RNA} | DNA:${DNA} | l=${LENGTH_THRESHOLD}"
echo  "l=${LENGTH_THRESHOLD} -v long_r=${OUT_LONG_R} -v long_d=${OUT_LONG_D}  -v shortR_r=${OUT_SHORTR_R} -v shortR_d=${OUT_SHORTR_D}  -v shortD_r=${OUT_SHORTD_R} -v shortD_d=${OUT_SHORTD_D}  -v shortRD_r=${OUT_SHORTRD_R} -v shortRD_d=${OUT_SHORTRD_D} -v stats=${OUT_STATS}"

paste <(gunzip -c "$RNA") <(gunzip -c "$DNA") | awk -v l=$LENGTH_THRESHOLD -v long_r="$OUT_LONG_R" -v long_d="$OUT_LONG_D"  -v shortR_r="$OUT_SHORTR_R" -v shortR_d="$OUT_SHORTR_D"  -v shortD_r="$OUT_SHORTD_R" -v shortD_d="$OUT_SHORTD_D"  -v shortRD_r="$OUT_SHORTRD_R" -v shortRD_d="$OUT_SHORTRD_D" -v stats="$OUT_STATS"  -F $'\t' 'BEGIN{OFS="\n"; n=0; n_long=0; n_shortr=0; n_shortd=0; n_shortrd=0; print l}{
    x1=$1; x2=$2; getline; r1=$1; r2=$2; getline; getline; n+=1;
            if (length(r1)>=l){
                    if(length(r2)>=l){
                        n_long+=1;
                        print x1, r1, "+", $1 | long_r;
                        print x2, r2, "+", $2 | long_d;
                    }
                    else{
                        n_shortd+=1;
                        print x1, r1, "+", $1 | shortD_r;
                        print x2, r2, "+", $2 | shortD_d;
                    }
            }
            else{
                     if(length(r2)>=l){
                        n_shortr+=1;
                        print x1, r1, "+", $1 | shortR_r;
                        print x2, r2, "+", $2 | shortR_d;
                    }
                    else{
                        n_shortrd+=1;
                        print x1, r1, "+", $1 | shortRD_r;
                        print x2, r2, "+", $2 | shortRD_d;
                    }

            }
}END{print n","n_long","n_shortr","n_shortd","n_shortrd"\n" >>stats}'


