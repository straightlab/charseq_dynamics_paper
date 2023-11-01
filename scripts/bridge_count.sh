#!/bin/bash

function usage {
    echo "usage: bridge_count bridgeF bridgeR [fastq_filename]"
    echo ""
    echo "Does a simple count of the total number of reads (N) and the number of reads where the bridge occurs exactly once in the forward/reverse direction (NbrF/NbrR), multiple times in the forward/reverse/both direction (NbrDupF/NbrDupR/NbrDupFR)"
    echo ""
    echo "  bridgeF      forward bridge sequence"
    echo "  bridgeR      reverse bridge sequence"
    echo "  fastq_filaname      a fasqt file of reads. If omitted, takes input from STDIN"
    echo ""
    echo "EXAMPLES"
    echo "  bash bridge_count.sh ACCGGCGTCCAAG CTTGGACGCCGGT myreads.fastq"
    echo "  zcat myreads.fastq.gz | bash bridge_count.sh ACCGGCGTCCAAG CTTGGACGCCGGT"
    exit 1
}

if [ $# -lt 2 ]; then
	usage
else
BRIDGEF=$1
BRIDGER=$2
echo "Simple bridge count with"
echo "Bridge forward: "$BRIDGEF
echo "Bridge reverse: "$BRIDGER
FILE=${3:--}


eval "awk '/^@/ {{i++; next;}} /$BRIDGEF|$BRIDGER/ {{noccF=gsub(/$BRIDGEF/,\"+\"); noccR=gsub(/$BRIDGER/,\"-\"); if (noccF+noccR==1) {{j+=noccF; k+=noccR;}} else {{if (noccF>0 && noccR>0) alldblFR++; else {{if (noccR>1) alldblR++; if (noccF>1) alldblF++}}; getline; getline;}}}} END {{print \"N=\",i+0, \"\tNbrF=\", j+0, \"\tNbrR=\", k+0, \"\tNbrDupFR=\", alldblFR+0, \"\tNbrDupF=\", alldblF+0, \"\tNbrDupR=\", alldblR+0}}'" "$FILE"
fi


eval "awk '/^@/ {{i++; next;}} /$BRIDGEF|$BRIDGER/ {{noccF=gsub(/$BRIDGEF/,\"+\"); noccR=gsub(/$BRIDGER/,\"-\"); if (noccF+noccR==1) {{j+=noccF; k+=noccR;}} else {{if (noccF>0 && noccR>0) alldblFR++; else {{if (noccR>1) alldblR++; if (noccF>1) alldblF++}}; getline; getline;}}}} END {{print \"N=\",i+0, \"\tNbrF=\", j+0, \"\tNbrR=\", k+0, \"\tNbrDupFR=\", alldblFR+0, \"\tNbrDupF=\", alldblF+0, \"\tNbrDupR=\", alldblR+0}}'" "$FILE"
fi
