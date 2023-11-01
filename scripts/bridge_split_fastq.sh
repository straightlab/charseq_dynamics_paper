#!/bin/bash

function usage {
    echo "usage: bridge_count bridgeF bridgeR fastqGZ_filename out_ub_filename out_pos_filename [out_mb_filename]"
    echo ""
    echo "Separates reads from a fastq.gz file based on the presence of a bridge sequence"
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
BRIDGEF=$1
BRIDGER=$2
FQ_IN=$3 #FastQ Input
FQ_UB=$4 #FastQ Unique Bridge -> output file
POS=$5 #Position of bridge in unique bridge -> output
echo "Simple bridge count with"
echo "Bridge forward: "$BRIDGEF
echo "Bridge reverse: "$BRIDGER

echo "Sending reads with single bridge to separate fastq"
zcat $FQ_IN | eval "awk 'BEGIN{OFS=\",\"}/^@/ {{i++; a=\$0; next;}} /$BRIDGER|$BRIDGEF/ {{noccF=gsub(/$BRIDGEF/,\"$BRIDGEF\"); noccR=gsub(/$BRIDGER/,\"$BRIDGER\"); if (noccF==0 && noccR==1) {{print i, \"0\", index(\$0,\"$BRIDGER\") > \"$POS\"; print a; print; getline; print; getline; print;}}; if (noccF==1 && noccR==0) {{print i, \"1\", index(\$0,\"$BRIDGEF\") > \"$POS\"; print a; print; getline; print; getline; print;}}}};'" | gzip -c > $FQ_UB

if [ ! -z $6 ]; then
    FQ_MB=$6 #FastQ Multiple Bridges -> output file
    echo "Sending reads with multiple bridges to separate fastq"
    zcat $FQ_IN | eval "awk '/^@/ {{a=\$0; next;}} /$BRIDGER|$BRIDGEF/ {{noccF=gsub(/$BRIDGEF/,\"$BRIDGEF\"); noccR=gsub(/$BRIDGER/,\"$BRIDGER\"); if (noccF+noccR>1) {{print a; print; getline; print; getline; print;}}}}'" | gzip -c > $FQ_MB
fi

fi
