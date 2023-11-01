#!/bin/bash

LENGTH_THRESHOLD=$1
OUTDIR="$2"
RNA1="$3"
RNA2="$4"
DNA="$5"

OUT_LONG="${OUTDIR%/}/long/R/"
OUT_SHORTR="${OUTDIR%/}/shortR/R/"
OUT_SHORTD="${OUTDIR%/}/shortD/R/"
OUT_SHORTRD="${OUTDIR%/}/shortRD/R/"
OUT_STATS="${OUTDIR%/}/R.filtering.stats.txt"


mkdir -p "$OUT_LONG"
mkdir -p "$OUT_SHORTR"
mkdir -p "$OUT_SHORTD"
mkdir -p "$OUT_SHORTRD"

echo -e "total,long,shortR,shortD,shortRD" >"$OUT_STATS"

cat /dev/null | gzip -c>${OUT_LONG}rna.1.fastq.gz
cat /dev/null | gzip -c>${OUT_LONG}rna.2.fastq.gz
cat /dev/null | gzip -c>${OUT_LONG}dna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTR}rna.1.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTR}rna.2.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTR}dna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTD}rna.1.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTD}rna.2.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTD}dna.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTRD}rna.1.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTRD}rna.2.fastq.gz
cat /dev/null | gzip -c>${OUT_SHORTRD}dna.fastq.gz



OUT_LONG_R1="gzip >${OUT_LONG}rna.1.fastq.gz"
OUT_LONG_R2="gzip >${OUT_LONG}rna.2.fastq.gz"
OUT_LONG_D="gzip >${OUT_LONG}dna.fastq.gz"
OUT_SHORTR_R1="gzip >${OUT_SHORTR}rna.1.fastq.gz"
OUT_SHORTR_R2="gzip >${OUT_SHORTR}rna.2.fastq.gz"
OUT_SHORTR_D="gzip >${OUT_SHORTR}dna.fastq.gz"
OUT_SHORTD_R1="gzip >${OUT_SHORTD}rna.1.fastq.gz"
OUT_SHORTD_R2="gzip >${OUT_SHORTD}rna.2.fastq.gz"
OUT_SHORTD_D="gzip >${OUT_SHORTD}dna.fastq.gz"
OUT_SHORTRD_R1="gzip >${OUT_SHORTRD}rna.1.fastq.gz"
OUT_SHORTRD_R2="gzip >${OUT_SHORTRD}rna.2.fastq.gz"
OUT_SHORTRD_D="gzip >${OUT_SHORTRD}dna.fastq.gz"


echo "RNA1:${RNA1} | RNA2:${RNA2} | DNA:${DNA}"
echo "-v l=${LENGTH_THRESHOLD} -v long_r1=${OUT_LONG_R1} -v long_r2=${OUT_LONG_R2} -v long_d=${OUT_LONG_D}  -v shortR_r1=${OUT_SHORTR_R1} -v shortR_r2=${OUT_SHORTR_R2} -v shortR_d=${OUT_SHORTR_D}  -v shortD_r1=${OUT_SHORTD_R1} -v shortD_r2=${OUT_SHORTD_R2} -v shortD_d=${OUT_SHORTD_D}  -v shortRD_r1=${OUT_SHORTRD_R1} -v shortRD_r2=${OUT_SHORTRD_R2} -v shortRD_d=${OUT_SHORTRD_D}"

paste <(gunzip -c "$RNA1") <(gunzip -c "$RNA2") <(gunzip -c "$DNA") | awk -v l=$LENGTH_THRESHOLD -v long_r1="$OUT_LONG_R1" -v long_r2="$OUT_LONG_R2" -v long_d="$OUT_LONG_D"  -v shortR_r1="$OUT_SHORTR_R1" -v shortR_r2="$OUT_SHORTR_R2" -v shortR_d="$OUT_SHORTR_D"  -v shortD_r1="$OUT_SHORTD_R1" -v shortD_r2="$OUT_SHORTD_R2" -v shortD_d="$OUT_SHORTD_D"  -v shortRD_r1="$OUT_SHORTRD_R1" -v shortRD_r2="$OUT_SHORTRD_R2" -v shortRD_d="$OUT_SHORTRD_D" -v stats="$OUT_STATS" -F $'\t' 'BEGIN{OFS="\n"; n=0; n_long=0; n_shortr=0; n_shortd=0; n_shortrd=0;print l}{
    x1=$1; x2=$2; x3=$3; getline; r1=$1; r2=$2; r3=$3; getline; getline; n+=1;
            if (length(r2)>=l){
                    if(length(r3)>=l){
                        n_long+=1;
                        print x1, r1, "+", $1 | long_r1;
                        print x2, r2, "+", $2 | long_r2;
                        print x3, r3, "+", $3 | long_d;

                    }
                    else{
                        n_shortd+=1;
                        print x1, r1, "+", $1 | shortD_r1;
                        print x2, r2, "+", $2 | shortD_r2;
                        print x3, r3, "+", $3 | shortD_d;

                    }
            }
            else{

                    if(length(r3)>=l){
                        n_shortr+=1;
                        print x1, r1, "+", $1 | shortR_r1;
                        print x2, r2, "+", $2 | shortR_r2;
                        print x3, r3, "+", $3 | shortR_d;
                    }
                    else{
                        n_shortrd+=1;
                        print x1, r1, "+", $1 | shortRD_r1;
                        print x2, r2, "+", $2 | shortRD_r2;
                        print x3, r3, "+", $3 | shortRD_d;
                    }

            }
}END{print n","n_long","n_shortr","n_shortd","n_shortrd"\n" >>stats}'


