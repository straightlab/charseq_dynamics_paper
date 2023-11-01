#!/bin/bash
outprefix="$1"
shift
nthreads=$1
shift

INFILESTR=${@}
INFILES=(${INFILESTR// / })
NFILES=${#INFILES[@]}

if [ $NFILES -eq 1 ]
then

    cp $1 $outprefix.pairs.gz
    pairix -f $outprefix.pairs.gz

else

    # header
    echo "Preparing header"
    gunzip -c ${INFILES[0]} | awk '(/^#/){print $0; next;}(/s+/){next;}{exit;}' | grep -v '^#command:'  > "${outprefix}_TMP.pairs"
    
    # unzipping to named pipes
    echo "Unzipping to named pipes"
    mydir=$(dirname "$outprefix")
    
    arg=''
    k=1
    for f in $INFILESTR
    do
    rm -f "${mydir}/pp.${k}"
    mkfifo "${mydir}/pp.${k}"
    arg="$arg ${mydir}/pp.${k}"
    gunzip -c "$f" | grep -v '^#' > "${mydir}/pp.${k}" &
    let "k++"
    done
    
    # merging
    echo "Merging"
    echo "sort -T ${mydir} -S3G --parallel=$nthreads -m -k2,2 -k4,4 -k3,3g -k5,5g ${arg} >> ${outprefix}_TMP.pairs" 
    sort -T "$mydir" -S3G --parallel=$nthreads -m -k2,2 -k4,4 -k3,3g -k5,5g $arg >> "${outprefix}_TMP.pairs"
    
    # compressing
    echo "Compressing merged file"
    bgzip --threads $nthreads -f "${outprefix}_TMP.pairs" && mv "${outprefix}_TMP.pairs.gz" "${outprefix}.pairs.gz"
    
    # indexing
    echo "Indexing"
    pairix -f "${outprefix}.pairs.gz"
    
    # clean up
    echo "Cleaning up"
    k=1
    for f in $INFILESTR
    do
    rm "${mydir}/pp.${k}"
    let "k++"
    done
fi

