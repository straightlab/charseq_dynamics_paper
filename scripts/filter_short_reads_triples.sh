#!/bin/bash

LENGTH_THRESHOLD=$1
OUTDIR="$2"
RNA="$3"
DNA1="$4"
DNA2="$5"

OUT_LONG="${OUTDIR%/}/"
OUT_SHORTR="${OUTDIR%/}/shortR/"
OUT_SHORTD="${OUTDIR%/}/shortD/"
OUT_SHORTRD="${OUTDIR%/}/shortRD/"

mkdir -p "$OUT_LONG"
mkdir -p "$OUT_SHORTR"
mkdir -p "$OUT_SHORTD"
mkdir -p "$OUT_SHORTRD"


OUT_LONG_R="gzip >${OUT_LONG}rna.fastq.gz"
OUT_LONG_D="gzip >${OUT_LONG}dna.fastq.gz"
OUT_SHORTR_R="gzip >${OUT_SHORTR}rna.fastq.gz"
OUT_SHORTR_D="gzip >${OUT_SHORTR}dna.fastq.gz"
OUT_SHORTD_R="gzip >${OUT_SHORTD}rna.fastq.gz"
OUT_SHORTD_D="gzip >${OUT_SHORTD}dna.fastq.gz"
OUT_SHORTRD_R="gzip >${OUT_SHORTRD}rna.fastq.gz"
OUT_SHORTRD_D="gzip >${OUT_SHORTRD}dna.fastq.gz"


echo "RNA:${RNA} | DNA:${DNA} | l=${LENGTH_THRESHOLD}"
echo  "l=${LENGTH_THRESHOLD} -v long_r=${OUT_LONG_R} -v long_d=${OUT_LONG_D}  -v shortR_r=${OUT_SHORTR_R} -v shortR_d=${OUT_SHORTR_D}  -v shortD_r=${OUT_SHORTD_R} -v shortD_d=${OUT_SHORTD_D}  -v shortRD_r=${OUT_SHORTRD_R} -v shortRD_d=${OUT_SHORTRD_D}"

paste <(gunzip -c "$RNA_SHORT") <(gunzip -c "$DNA_LONG1") <(gunzip -c "$DNA_LONG2") | awk -v l=$LENGTH_THRESHOLD -v long_r="$OUT_LONG_R" -v long_d1="$OUT_LONG_D1" -v long_d2="$OUT_LONG_D2"  -v shortR_r="$OUT_SHORTR_R" -v shortR_d1="$OUT_SHORTR_D1" -v shortR_d2="$OUT_SHORTR_D2"  -v shortD_r="$OUT_SHORTD_R" -v shortD_d1="$OUT_SHORTD_D1" -v shortD_d2="$OUT_SHORTD_D2"  -v shortRD_r="$OUT_SHORTRD_R" -v shortRD_d1="$OUT_SHORTRD_D1" -v shortRD_d2="$OUT_SHORTRD_D2" -F $'\t' 'BEGIN{OFS="\n"; print l}{
    x1=$1; x2=$2; x3=$3; getline; r1=$1; r2=$2; r3=$3; getline; getline; 
            if (length(r1)>=l){
                    if(length(r2)>=l){
                        print x1, r1, "+", $1 | long_r;
                        print x2, r2, "+", $2 | long_d1;
                        print x3, r3, "+", $3 | long_d2;

                    }
                    else{
                        print x1, r1, "+", $1 | shortD_r;
                        print x2, r2, "+", $2 | shortD_d1;
                        print x3, r3, "+", $3 | shortD_d2;

                    }
            }
            else{

                    if(length(r2)>=l){
                        print x1, r1, "+", $1 | shortR_r;
                        print x2, r2, "+", $2 | shortR_d1;
                        print x3, r3, "+", $3 | shortR_d2;
                    }
                    else{
                        print x1, r1, "+", $1 | shortRD_r;
                        print x2, r2, "+", $2 | shortRD_d1;
                        print x3, r3, "+", $3 | shortRD_d2;
                    }

            }
}'


