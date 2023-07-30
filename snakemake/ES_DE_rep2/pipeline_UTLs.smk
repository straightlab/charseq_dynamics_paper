# import relevant module
import os
# shell.prefix("PATH=$PATH:"+config["scripts_path"])
import yaml

# load the sample definition file if exists
if 'samples_def_file' in config:
    with open(config['samples_def_file']) as file:
        samples_config = yaml.load(file)
        # print(samples_config)
        print("Loading samples definition from "+config['samples_def_file'])
else:
    samples_config=config
    print('No sample definition file specified, loading samples definition from Snakemake configuration file')

# specify the working directory
workdir: config["workdir"]

# define samples to run
SAMPLES = samples_config["samples"]
for key, value in SAMPLES.items():
    SAMPLES[key]["PE1"]=srcdir(os.path.join(samples_config["data_root"],value["PE1"]))
    SAMPLES[key]["PE2"]=srcdir(os.path.join(samples_config["data_root"],value["PE2"]))

if 'run_samples' in config:
    if (len(config["run_samples"])==1) and (config["run_samples"][0]=='all'):
        SAMPLESID=list(SAMPLES.keys())
    else:
        SAMPLESID=config.get("run_samples",list(SAMPLES.keys()))
else:
    SAMPLESID=list(SAMPLES.keys())
print("Running samples :")
print(SAMPLESID)

# define pairing mode and types to run
PAIRING_TYPES=config.get("run_pairing_types",['exons'])
PAIRING_MODES=config.get("run_pairing_modes",['strngtieCharOnly_hg38'])

# define outputs
OUTPUTS_DEFAULT=[
    expand('{sampleID}/pairs/{pairing_mode}/{pair_type}/filtered/DNAq15-blacklisted_RNAunq/{dnarna}.bed.gz', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, pair_type=['exons'],dnarna=['dna','rna'])
    ]

OUTPUTS=config.get("outputs",OUTPUTS_DEFAULT)


wildcard_constraints:
    sampleID="[^/]+",
    rna_alignment_mode="[^/]+",
    alignment_type="[^/]+",
    rna_alignment_mode_parent="[^/]+",
    pairing_mode="[^/]+",
    bedORbam="bam|bed.gz",
    sense='F|R|U',
    pair_type='exons',
    splitOR5='split|5|3',
    strand='\+|-|=|-rev',
    strandStrict='\+|-|=',
    cistrancs='cis|trans|dnaambig'

### =========== RULES START ===========================
rule all:
    input: OUTPUTS

# rule prepare_input:
#     input:
#         PE1=lambda wildcards: SAMPLES[wildcards.sampleID]["PE1"],
#         PE2=lambda wildcards: SAMPLES[wildcards.sampleID]["PE2"]
#     output:
#         PE1='{sampleID}/raw/chimera.1.fastq.gz',
#         PE2='{sampleID}/raw/chimera.2.fastq.gz'
#     run:
#         nlines=4*config["nreads"]
#         if nlines<0:
#             print("Running pipeline on ALL reads")
#             shell("mkdir -p {wildcards.sampleID}")
#             shell("ln -s $(readlink -e {input.PE1}) {output.PE1}")
#             shell("ln -s $(readlink -e {input.PE2}) {output.PE2}")
#         else:
#             print("Running pipeline on subset of " + str(nlines) + " reads")
#             shell("mkdir -p {wildcards.sampleID}")
#             shell("set +o pipefail; pwd; zcat $(readlink -f {input.PE1}) | head -n " + str(nlines) + " | gzip > {output.PE1}")
#             shell("set +o pipefail; pwd; zcat $(readlink -f {input.PE2}) | head -n " + str(nlines) + " | gzip > {output.PE2}")



# # # === ALIGN SINGLE END DNA AND RNA READS TO GENOME ===
rule align_RNA_star_SE: #2g = align to genome
    input: lambda wildcards: '{sampleID}/{fq}'.format(sampleID=wildcards.sampleID, fq=config['rna_alignment_modes'][wildcards.rna_alignment_mode]['fq'])
    output:
        bam_transcriptome=temp('{sampleID}/alignments/rna/{rna_alignment_mode}/rna.transcriptome.all.bam'),
        bam_genome=temp('{sampleID}/alignments/rna/{rna_alignment_mode}/rna.genome.all.bam')
    params:
        star_config=lambda wildcards: config['star_configurations'].get(config["rna_alignment_modes"][wildcards.rna_alignment_mode]['config'], []).get('S',[]),
        xmx=config["max_mem"]
    threads: config["max_threads"]
    shell: """
STAR --runThreadN {threads} --readFilesIn {input} \
--outSAMattrRGline ID:{wildcards.sampleID} SM:{wildcards.sampleID} \
--quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted --readFilesCommand zcat --outFileNamePrefix {wildcards.sampleID}/alignments/rna/{wildcards.rna_alignment_mode}/rna. {params.star_config} && picard -Xmx{params.xmx} SortSam I={wildcards.sampleID}/alignments/rna/{wildcards.rna_alignment_mode}/rna.Aligned.out.bam O={output.bam_genome} SORT_ORDER=queryname && picard -Xmx{params.xmx} SortSam I={wildcards.sampleID}/alignments/rna/{wildcards.rna_alignment_mode}/rna.Aligned.toTranscriptome.out.bam O={output.bam_transcriptome} SORT_ORDER=queryname && rm {wildcards.sampleID}/alignments/rna/{wildcards.rna_alignment_mode}/rna.Aligned.toTranscriptome.out.bam {wildcards.sampleID}/alignments/rna/{wildcards.rna_alignment_mode}/rna.Aligned.out.bam
"""

rule filter_intergenic:
    input:
        bam_transcriptome='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.transcriptome.all.bam',
        bam_genome='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.genome.all.bam',
        bam_intergenic= lambda wildcards: '{sampleID}/alignments/rna/{rna_alignment_mode_parent}/bytype/rna.intergenic.bam'.format(sampleID=wildcards.sampleID, rna_alignment_mode_parent=config['rna_alignment_modes'][wildcards.rna_alignment_mode]['alignment_mode_parent'])
    output:
        bam_transcriptome='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.transcriptome.bam',
        bam_genome='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.genome.bam',
        reads2keep=temp('{sampleID}/alignments/rna/{rna_alignment_mode}/rna.intergenic.readidsToKeep.txt')
    threads: 4
    params:
        xmx=str(4*config["mem_per_cpu"])+"g"
    shell: """
samtools view {input.bam_intergenic} | cut -f1 >{output.reads2keep}; picard -Xmx{params.xmx} FilterSamReads I={input.bam_transcriptome} READ_LIST_FILE={output.reads2keep} FILTER=includeReadList SORT_ORDER=queryname O={output.bam_transcriptome};
picard -Xmx{params.xmx} FilterSamReads I={input.bam_genome} READ_LIST_FILE={output.reads2keep} FILTER=includeReadList SORT_ORDER=queryname O={output.bam_genome}
 """

rule tag_bam_rna:
    input: 
        bam_transcriptome='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.transcriptome.bam',
        bam_genome='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.genome.bam'
    output:
        tag_bam='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.tag.bam',
        notag_bam='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.notag.bam',
        stats='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.tag.stats',
        tag='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.tag.json'
    threads: 1
    params:
        exons_file=lambda wildcards: config['star_configurations'].get(config["rna_alignment_modes"][wildcards.rna_alignment_mode]['config'], [])['exons'],
        annot_file=lambda wildcards: config['star_configurations'].get(config["rna_alignment_modes"][wildcards.rna_alignment_mode]['config'], [])['annot']
    shell: """
tagtools tag --outnoannot {output.notag_bam} --annotfile {params.annot_file} --outstats {output.stats} --ambivout {output.tag} --strand 1 {input.bam_genome} {input.bam_transcriptome} {params.exons_file} {output.tag_bam}
"""


rule salmon_rna_bam:
    input:
        bam_transcriptome='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.transcriptome.bam'
    output:
        quantfile='{sampleID}/alignments/rna/{rna_alignment_mode}/salmon_quant/quant.sf'
    params:
        transcriptome_fa=lambda wildcards: config['star_configurations'].get(config["rna_alignment_modes"][wildcards.rna_alignment_mode]['config'], [])['transcriptome_fa']
    threads: 4
    shell:"""
salmon quant -t {params.transcriptome_fa} -l A -a {input.bam_transcriptome} -o $(dirname {output.quantfile}) --dumpEq
"""

rule resolve_tag_bam_rna:
    input: 
        tag_bam='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.tag.bam',
        quantfile='{sampleID}/alignments/rna/{rna_alignment_mode}/salmon_quant/quant.sf',
        tag='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.tag.json'
    output:
        resolved_tag_bam='{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.exons.bam'
    threads: 1
    params:
        annot_file=lambda wildcards: config['star_configurations'].get(config["rna_alignment_modes"][wildcards.rna_alignment_mode]['config'], [])['annot']
    shell: """
tagtools resolve --annotfile {params.annot_file} {input.tag_bam} {input.tag} {input.quantfile} {output.resolved_tag_bam}
"""



rule merge_rna_bams:   
    input: lambda wildcards: ['{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/{rna}.bam'.format(sampleID=wildcards.sampleID, rna_alignment_mode=wildcards.rna_alignment_mode, rna=r) for r in ['rna.exons','rna.introns', 'rna.intergenic', 'rna.novel']]
    output: '{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.all.bam'
    params:
        xmx=config["max_mem"]
    run:
        for ins in input:
            shell("samtools sort -n -o "+ins+"_samsort.bam "+ins)
        shell("samtools merge -f -n {output}_samsort.bam "+" ".join([ins2+"_samsort.bam" for ins2 in input]))
        shell("picard -Xmx{params.xmx} SortSam I={output}_samsort.bam O={output} SORT_ORDER=queryname && rm {output}_samsort.bam")
        for ins in input:
            shell("rm "+ins+"_samsort.bam ")

rule generate_pairs: 
    input: 
        rna= lambda wildcards: '{sampleID}/alignments/rna/{rna_alignment_mode}/{rna_bam}'.format(sampleID=wildcards.sampleID, rna_alignment_mode=config["pairing_modes"][wildcards.pairing_mode]['rna'], rna_bam=config["pairing_modes"][wildcards.pairing_mode]["alignment_types"][wildcards.alignment_type]),
        dna= lambda wildcards: '{sampleID}/alignments/dna/{dna_alignment_mode}/dna.bam'.format(sampleID=wildcards.sampleID, dna_alignment_mode=config["pairing_modes"][wildcards.pairing_mode]['dna'])
    output:
        rna='{sampleID}/pairs/{pairing_mode}/{alignment_type}/paired.rna.bam',
        dna='{sampleID}/pairs/{pairing_mode}/{alignment_type}/paired.dna.bam',
        pairs='{sampleID}/pairs/{pairing_mode}/{alignment_type}/rd.pairs', #readid sorted
        stats='{sampleID}/pairs/{pairing_mode}/{alignment_type}/rd.stats.txt',
        rna_link='{sampleID}/pairs/{pairing_mode}/{alignment_type}/rna.bam'
    params:
        pairing_options=lambda wildcards: config["pairing_modes"][wildcards.pairing_mode]["pairing_options"]
    threads: 1
    shell: """
chartools pairup -m prefix -s {output.stats} -o {output.pairs} --outpairedRNA {output.rna} --outpairedDNA {output.dna} {params.pairing_options} {input.rna} {input.dna}; ln -s $(readlink -e {input.rna}) {output.rna_link}
"""

rule index_pairs:
    input:
        pairs='{sampleID}/pairs/{pairing_mode}/{pair_type}/rd.pairs'
    output:
        compressed='{sampleID}/pairs/{pairing_mode}/{pair_type}/rd.indexed.pairs.gz'
    params:
        pairix_header=lambda wildcards: config["pairing_modes"][wildcards.pairing_mode]["pairix_header"]
    threads: 6
    shell: """
export LC_ALL=C; cat {params.pairix_header} <(sort -T $(dirname {output.compressed}) -S3G --parallel={threads} -k2,2 -k4,4 -k3,3n -k5,5n {input.pairs}) | bgzip --thread {threads} -c > {output.compressed}; pairix {output.compressed}
"""

#1=CHR_interaction_locus, 2=START_interaction_locu, 3=STOP_interaction_locus, 4=READID, 5=MAQ, 6=STRAND, 
#7=CHR_RNA, 8=START_RNA_3'end, 9=STOP_RNA_3'end, 10=unused, 11=unsed, 12=RNA ensembl ID (ENST# for exon, ENSG# for intron), 13=RNA_3'end, 14=RNA name, 15=RNA type, 16=ambivalence group, 17=RNA gene ensemble ID (ENSG#)

rule make_bed_dna:
    input:
        pairs='{sampleID}/pairs/{pairing_mode}/{pair_type}/rd.pairs'
    output:
        dna='{sampleID}/pairs/{pairing_mode}/{pair_type}/dna.bed.gz'
    threads: 2
    shell: """
export LC_ALL=C; cat {input.pairs} | awk -F $'\t' 'BEGIN{{OFS=FS}}{{same=0; $2=substr($2,3); flight=$5-$3; if ($7=="-"){{$5=$5-3;}}; if($2==$4){{same=1;}}; print $4, $5, $5+1, $1, $9, $7, $2, $3, $3+1, $8, $6, $18, $19, $20, $21, $25, $26, same, flight;}}' | sort -k1,1 -k2,2n | bgzip --thread {threads} -c > {output.dna}; tabix -0 -f -p bed {output.dna}
"""

rule filter_blacklist_bed_dna: #Q15, blacklisted, $25==1, no unassigned chr scaffolds
    input:'{sampleID}/pairs/{pairing_mode}/{pair_type}/dna.bed.gz'
    output: '{sampleID}/pairs/{pairing_mode}/{pair_type}/filtered/DNAq15-blacklisted_RNAunq/dna.bed.gz'
    params: 
        blacklist_bed=config['blacklist_bed']
    shell: """
export LC_ALL=C; bedtools intersect -sorted -wa -v -a <(zcat {input} | awk -F $'\t' 'BEGIN{{OFS=FS}}($16==1 && $5>15)') -b {params.blacklist_bed} | bgzip -c > {output}; tabix -0 -f -p bed {output} 
    """

#Dna print $4, $5, $5+1, $1, $9, $7, $2, $3, $3+1, $8, $6, $18, $19, $20, $21, $25, $26, same, flight;
#rna       $2, $3, $3+1, $1, $8, $6, $4, $5, $5+1, $9, $7, $18, $19, $20, $21, $25, $26, same, flight
# convert: $7, $8, $9, $4, $10, $11, $1, $2, $3, $5, $6, $12, $13, $14, $15, $16, $17, $18, $19
rule dnabed_to_rna:
    input:'{sampleID}/{bed_folder}/dna.bed.gz'
    output:'{sampleID}/{bed_folder}/rna.bed.gz'
    shell: """
export LC_ALL=C; gunzip -c {input} | awk -F $'\t' 'BEGIN{{OFS=FS}}{{print $7, $8, $9, $4, $10, $11, $1, $2, $3, $5, $6, $12, $13, $14, $15, $16, $17, $18, $19}}' | sort -k1,1 -k2,2n | bgzip --thread {threads} -c > {output}; tabix -0 -f -p bed {output}
    """

#1=CHR_interaction_locus, 2=START_interaction_locu, 3=STOP_interaction_locus, 4=READID, 5=MAQ, 6=STRAND, 
#7=CHR_RNA, 8=START_RNA_3'end, 9=STOP_RNA_3'end, 10=unused, 11=unsed, 12=RNA ensembl ID (ENST# for exon, ENSG# for intron), 13=RNA_3'end, 14=RNA name, 15=RNA type, 16=ambivalence group, 17=RNA gene ensemble ID (ENSG#)


rule gene_counts_table:
    input: '{sampleID}/pairs/{pairing_mode}/{pair_type}/rd.pairs'
    output: 
        complete='{sampleID}/pairs/{pairing_mode}/{pair_type}/stats/rnaTable.stats.txt',
        simple='{sampleID}/pairs/{pairing_mode}/{pair_type}/stats/rnaTableSimple.stats.txt'
    shell: """
export LC_ALL=C; awk -F $'\t' 'BEGIN{{OFS=" "}}{{sure=-1;istrans=-1; if($25==1){{sure=1;}}else{{if ($22==1){{sure=0;}};}}; if(($9>15) && (substr($2,3)==$4)){{istrans=0;}}else{{if($9>15){{istrans=1;}};}}; {{print $18, sure, istrans;}}}}' {input} | tee >(sort | uniq -c | sort -k1,1nr | sed -e 's/ *//' >{output.complete}) | cut -d " " -f1,1 |sort | uniq -c | sort -k1,1nr | sed -e 's/ *//' >{output.simple}
"""

rule gene_counts_summary:
    input: '{sampleID}/pairs/{pairing_mode}/{pair_type}/stats/rnaTable.stats.txt'
    params: 
        txdb=config['annotations']['txdb'],
        genedb=config['annotations']['genedb'],
        chrdb=config['annotations']['chrdb']
    output: 
        all_parquet='{sampleID}/pairs/{pairing_mode}/{pair_type}/stats/gene_expression_multimapRNA.parquet',
        all_xlsx='{sampleID}/pairs/{pairing_mode}/{pair_type}/stats/gene_expression_multimapRNA.xlsx',
        unq_parquet='{sampleID}/pairs/{pairing_mode}/{pair_type}/stats/gene_expression_NOmultimapRNA.parquet',
        unq_xlsx='{sampleID}/pairs/{pairing_mode}/{pair_type}/stats/gene_expression_NOmultimapRNA.xlsx'
    threads: 1
    shell: """
chartools genestats3 {input} {params.chrdb} {params.txdb} {params.genedb} "$(dirname {output.all_parquet})"/gene_expression_
"""

rule make_bed_rna:
    input:
        pairs='{sampleID}/pairs/{pairing_mode}/{pair_type}/rd.pairs'
    output:
        rna='{sampleID}/pairs/{pairing_mode}/{pair_type}/rna.bed.gz'
    threads: 2
    shell: """
export LC_ALL=C; cat {input.pairs} | awk -F $'\t' 'BEGIN{{OFS=FS}}{{same=0; $2=substr($2,3); flight=$5-$3; if ($7=="-"){{$5=$5-3;}}; if($2==$4){{same=1;}}; print $2, $3, $3+1, $1, $8, $6, $4, $5, $5+1, $9, $7, $18, $19, $20, $21, $25, $26, same, flight;}}' | sort -k1,1 -k2,2n | bgzip --thread {threads} -c > {output.rna}; tabix -0 -f -p bed {output.rna}
"""

rule make_bw:
    input: '{sampleID}/{bed}.{bedOrBam}'
    output: 
        bg=temp('{sampleID}/{bed}.{bedORbam}.{sense}.bedgraph'),
        bw='{sampleID}/{bed}.{bedORbam}.{sense}.bw'
    params:
        chrdb=config['annotations']['chrdb']
    threads: 2
    run:
        if bedOrBam=="bam":
            if (sense=="F"):
                shell("export LC_ALL=C; bedtools coverage -ibam {input} -strand + -bg | sort -k1,1 -k2,2n >{output.bg}")
            elif (sense=="R"):
                shell("export LC_ALL=C; bedtools coverage -ibam {input} -strand - -bg | sort -k1,1 -k2,2n >{output.bg}")
            else:
                shell("export LC_ALL=C; bedtools coverage -ibam {input} -bg | sort -k1,1 -k2,2n >{output.bg}")
        else:
            if (sense=="F"):
                shell("export LC_ALL=C; bedtools coverage -i <(gunzip -c {input}) -g {params.chrdb} -strand + -bg | sort -k1,1 -k2,2n >{output.bg}")
            elif (sense=="R"):
                shell("export LC_ALL=C; bedtools coverage -i <(gunzip -c {input}) -g {params.chrdb} -strand - -bg | sort -k1,1 -k2,2n >{output.bg}")
            else:
                shell("export LC_ALL=C; bedtools coverage -i <(gunzip -c {input}) -g {params.chrdb} -bg | sort -k1,1 -k2,2n >{output.bg}")
        print("Running bedGraphToBigwig")
        shell("bedGraphToBigwig {output.bg} {params.chrdb} {output.bw}")