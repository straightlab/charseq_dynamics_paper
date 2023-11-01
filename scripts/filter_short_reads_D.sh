#!/bin/bash

LENGTH_THRESHOLD=$1
OUTDIR="$2"
RNA="$3"
DNA1="$4"
DNA2="$5"

OUT_LONG="${OUTDIR%/}/long/D/"
OUT_SHORTR="${OUTDIR%/}/shortR/D/"
OUT_SHORTD="${OUTDIR%/}/shortD/D/"
OUT_SHORTRD="${OUTDIR%/}/shortRD/D/"
OUT_STATS="${OUTDIR%/}/D.filtering.stats.txt"

mkdir -p "$OUT_LONG"
mkdir -p "$OUT_SHORTR"
mkdir -p "$OUT_SHORTD"
mkdir -p "$OUT_SHORTRD"

echo -e "total,long,shortR,shortD,shortRD" >"$OUT_STATS"

cat /dev/null | gzip -c>${OUT_LONG}dna.1.fastq.gz
cat /dev/null | gzip -c>${OUT_LONG}dna.2.fastq.gz
cat /dev/null | gzip -c>${OUT_LONG}rna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTR}dna.1.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTR}dna.2.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTR}rna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTD}dna.1.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTD}dna.2.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTD}rna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTRD}dna.1.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTRD}dna.2.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTRD}rna.fastq.gz

OUT_LONG_R="gzip >${OUT_LONG}rna.fastq.gz"
OUT_LONG_D1="gzip >${OUT_LONG}dna.1.fastq.gz"
OUT_LONG_D2="gzip >${OUT_LONG}dna.2.fastq.gz"
OUT_SHORTR_R="gzip >${OUT_SHORTR}rna.fastq.gz"
OUT_SHORTR_D1="gzip >${OUT_SHORTR}dna.1.fastq.gz"
OUT_SHORTR_D2="gzip >${OUT_SHORTR}dna.2.fastq.gz"
OUT_SHORTD_R="gzip >${OUT_SHORTD}rna.fastq.gz"
OUT_SHORTD_D1="gzip >${OUT_SHORTD}dna.1.fastq.gz"
OUT_SHORTD_D2="gzip >${OUT_SHORTD}dna.2.fastq.gz"
OUT_SHORTRD_R="gzip >${OUT_SHORTRD}rna.fastq.gz"
OUT_SHORTRD_D1="gzip >${OUT_SHORTRD}dna.1.fastq.gz"
OUT_SHORTRD_D2="gzip >${OUT_SHORTRD}dna.2.fastq.gz"

echo "RNA:${RNA} | DNA1:${DNA1} | DNA2:${DNA2}"
echo  "-v l=${LENGTH_THRESHOLD} -v long_r=${OUT_LONG_R} -v long_d1=${OUT_LONG_D1} -v long_d2=${OUT_LONG_D2}  -v shortR_r=${OUT_SHORTR_R} -v shortR_d1=${OUT_SHORTR_D1} -v shortR_d2=${OUT_SHORTR_D2}  -v shortD_r=${OUT_SHORTD_R} -v shortD_d1=${OUT_SHORTD_D1} -v shortD_d2=${OUT_SHORTD_D2}  -v shortRD_r=${OUT_SHORTRD_R} -v shortRD_d1=${OUT_SHORTRD_D1} -v shortRD_d2=${OUT_SHORTRD_D2}"

paste <(gunzip -c "$RNA") <(gunzip -c "$DNA1") <(gunzip -c "$DNA2") | awk -v l=$LENGTH_THRESHOLD -v long_r="$OUT_LONG_R" -v long_d1="$OUT_LONG_D1" -v long_d2="$OUT_LONG_D2"  -v shortR_r="$OUT_SHORTR_R" -v shortR_d1="$OUT_SHORTR_D1" -v shortR_d2="$OUT_SHORTR_D2"  -v shortD_r="$OUT_SHORTD_R" -v shortD_d1="$OUT_SHORTD_D1" -v shortD_d2="$OUT_SHORTD_D2"  -v shortRD_r="$OUT_SHORTRD_R" -v shortRD_d1="$OUT_SHORTRD_D1" -v shortRD_d2="$OUT_SHORTRD_D2" -v stats="$OUT_STATS" -F $'\t' 'BEGIN{OFS="\n"; n=0; n_long=0; n_shortr=0; n_shortd=0; n_shortrd=0; print l}{
    x1=$1; x2=$2; x3=$3; getline; r1=$1; r2=$2; r3=$3; getline; getline; n+=1;
            if (length(r1)>=l){
                    if(length(r2)>=l){
                        n_long+=1;
                        print x1, r1, "+", $1 | long_r;
                        print x2, r2, "+", $2 | long_d1;
                        print x3, r3, "+", $3 | long_d2;
                    }
                    else{
                        n_shortd+=1;
                        print x1, r1, "+", $1 | shortD_r;
                        print x2, r2, "+", $2 | shortD_d1;
                        print x3, r3, "+", $3 | shortD_d2;
                    }
            }
            else{
                    if(length(r2)>=l){
                        n_shortr+=1;
                        print x1, r1, "+", $1 | shortR_r;
                        print x2, r2, "+", $2 | shortR_d1;
                        print x3, r3, "+", $3 | shortR_d2;
                    }
                    else{
                        n_shortrd+=1;
                        print x1, r1, "+", $1 | shortRD_r;
                        print x2, r2, "+", $2 | shortRD_d1;
                        print x3, r3, "+", $3 | shortRD_d2;
                    }

            }
}END{print n","n_long","n_shortr","n_shortd","n_shortrd"\n" >>stats}'

