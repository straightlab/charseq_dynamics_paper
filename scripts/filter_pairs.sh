
#!/bin/bash

function usage {
    echo "usage: bams2chibam rna_in dna_in out rna_filter dna_filter picard_path rgid rglb rgsm"
    echo ""
    echo "Combines RNA and DNA alignment bam files into single "
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

#1 fq1, #2 fq2 --> use to pick up names, #3 bam1, #4 bam2, --> outputs bam1_dedup, bam2_dedup
gunzip -c "$1" | awk 'BEGIN{Ntot=0; NmultiR=0; NmultiD=0; NbadR=0; NbadD=0; Nmulti=0; Nbad=0;}(/^#/){print; next;}{Ntot+=1; mv=0; if($8<2){NmultiR+=1; mv=1}; if ($9<2){NmultiD+=1; mv=1;} if($8<2 && $9<2){Nmulti+=1};  if($8>1 && $8<10){NbadR+=1; mv=1;}; if($9>1 && $9<10){NbadD+=1; mv=1;}; if($9>1 && $9<10 && $8>1 && $8<10){Nbad+=1}; if(mv==0){print $0};}END{printf("Ntot NmultiR NmultiD Nmulti NbadR NbadD Nbad\n%g %g %g %g %g %g %g\n", Ntot, NmultiR, NmultiD, Nmulti, NbadR, NbadD, Nbad) > "/dev/stderr"}'
