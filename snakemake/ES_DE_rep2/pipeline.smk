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
PAIRING_TYPES=config.get("run_pairing_types",['all','novel','intergenic','exons','introns'])
PAIRING_MODES=config.get("run_pairing_modes",['gencondeV29_hg38'])

# define outputs
OUTPUTS_DEFAULT=[
    expand('{sampleID}/pairs/{pairing_mode}/{pair_type}/stats/gene_expression_NOmultimapRNA.parquet', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, pair_type=['all']), 
    expand('{sampleID}/pairs/{pairing_mode}/{pair_type}/rd.pairs', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, pair_type=PAIRING_TYPES),
    expand('{sampleID}/pairs/{pairing_mode}/{pair_type}/rd.indexed.pairs.gz', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, pair_type=PAIRING_TYPES),
    expand('{sampleID}/pairs/{pairing_mode}/{pair_type}/{dnarna}.bed.gz', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, pair_type=['all'],dnarna=['dna','rna']),
    expand('{sampleID}/pairs/{pairing_mode}/{pair_type}/filtered/DNAq15-blacklisted_RNAunq/{dnarna}.bed.gz', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, pair_type=['all'],dnarna=['dna','rna']),
    # expand('{sampleID}/pairs/gencondeV29_hg38/{type}/{rd}.bed.gz', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, type=PAIRING_TYPES, rd=['rna','dna']),
    expand('{sampleID}/alignments/{salmon_alignment_mode}/stats/gene_expression_salmon.parquet', sampleID=SAMPLESID,salmon_alignment_mode=['salmon_gencodeV29']),
    expand('{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/{dnarna}.{cistrans}.{splitOR5}.{strand}.bw', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, dnarna=['dna','rna'], alignment_type=['all','introns','exons','intergenic'], strand=['-rev','+'], splitOR5=['5','split'], cistrans=['cis','trans']),
    expand('{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/{dnarna}.{cistrans}.{splitOR5}.{strand}.bw', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, dnarna=['rna'], alignment_type=['all','intergenic'], strand=['=','+','-rev'], splitOR5=['3'], cistrans=['cis','trans']),
    expand('{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/{dnarna}.{cistrans}.{splitOR5}.{strand}.bw', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, dnarna=['rna','dna'], alignment_type=['all','intergenic'], strand=['+','-rev'], splitOR5=['3','split'], cistrans=['dnaambig']),
    expand('{sampleID}/pairs/{pairing_mode}/{pair_type}/rd.simple.pairs',sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, pair_type=['exons','introns']),
    expand('{sampleID}/alignments/rna/star_gencodeV29/noalign.stats.txt', sampleID=SAMPLESID)
    # expand('{sampleID}/pairs/{pairing_mode}/{pair_type}/{rd}.bed.gz.{strand}.bw', sampleID=SAMPLESID, pairing_mode=PAIRING_MODES, pair_type=['all','novel','intergenic','exons','introns'], rd=['rna','dna'], strand=["F","R","U"])
    ]

OUTPUTS=config.get("outputs",OUTPUTS_DEFAULT)

# INTERNAL DEFINITIONS, don't mess with it
readtype_SE=['F.dna','F.rna','R.dna','R.rna','0.xxx']
readtype_PE=['00.xxx.1', '00.xxx.2', 'F0.dna.1', 'F0.dna.2', 'F0.rna.1', 'FF.dna.1', 'FF.dna.2', 'FF.rna.1', 'FF.rna.2', 'FM.dna.1', 'FM.rna.1', 'FM.xxx.2', 'FR.dna.1', 'FR.dna.2', 'FR.rna.1', 'FR.rna.2', 'M0.xxx.1', 'M0.xxx.2', 'MM.xxx.1', 'MM.xxx.2', '0R.dna', '0R.rna.1', '0R.rna.2', 'RM.dna.1', 'RM.rna.1', 'RM.xxx.2', 'RR.dna.1',  'RR.dna.2', 'RR.rna.1', 'RR.rna.2']

wildcard_constraints:
    sampleID="[^/]+",
    chimera_typeP="allPairs|unmergedPairs_pear",
    chimera_typeS="[^/]+",
    filter="long|shortR|shortD|shortRD|unfiltered",
    mates_pairing_mode="[^/]+",
    rna_alignment_mode="[^/]+",
    dna_alignment_mode="[^/]+",
    pairing_mode="[^/]+",
    bedORbam="bam|bed.gz",
    sense='F|R|U',
    pair_type='all|novel|intergenic|exons|introns',
    splitOR5='split|5|3',
    strand='\+|-|=|-rev',
    strandStrict='\+|-|=',
    cistrancs='cis|trans|dnaambig'
    

### =========== RULES START ===========================
rule all:
    input: OUTPUTS

# === LINK INPUT FILE INTO THE PIPELINE RUNNING FOLDER ===
# a symlink to the raw data is created in the pipeline running folders {sampleID}. If the "nreads" option is set (not -1), then instead of a symlink, the first "nreads" of the raw data are COPIED in the {sampleID} folder.
rule prepare_input:
    input:
        PE1=lambda wildcards: SAMPLES[wildcards.sampleID]["PE1"],
        PE2=lambda wildcards: SAMPLES[wildcards.sampleID]["PE2"]
    output:
        PE1='{sampleID}/raw/chimera.1.fastq.gz',
        PE2='{sampleID}/raw/chimera.2.fastq.gz'
    run:
        nlines=4*config["nreads"]
        if nlines<0:
            print("Running pipeline on ALL reads")
            shell("mkdir -p {wildcards.sampleID}")
            shell("ln -s $(readlink -e {input.PE1}) {output.PE1}")
            shell("ln -s $(readlink -e {input.PE2}) {output.PE2}")
        else:
            print("Running pipeline on subset of " + str(nlines) + " reads")
            shell("mkdir -p {wildcards.sampleID}")
            shell("set +o pipefail; pwd; zcat $(readlink -f {input.PE1}) | head -n " + str(nlines) + " | gzip > {output.PE1}")
            shell("set +o pipefail; pwd; zcat $(readlink -f {input.PE2}) | head -n " + str(nlines) + " | gzip > {output.PE2}")


# one pass through the file to get bridge occurence statistics : count #reads, #reads with unique bridge in forward/reverse orientation, #reads with multiple bridges.
rule simple_bridge_count_PE:
    input:
        PE1='{sampleID}/{upstream}/chimera.1.fastq.gz',
        PE2='{sampleID}/{upstream}/chimera.2.fastq.gz'
    params:
        bridgef=lambda wildcards: samples_config["samples"][wildcards.sampleID].get("bridgeF",samples_config["bridgeF"]),
        bridger=lambda wildcards: samples_config["samples"][wildcards.sampleID].get("bridgeR",samples_config["bridgeR"]),
        L=config['read_length']
    output: '{sampleID}/{upstream}/bridge.stats.PE.txt'
    shell: """
debridge.jl {params.bridgef} {params.bridger} {wildcards.sampleID}/{wildcards.upstream}/ {input.PE1} {input.PE2} -l {params.L} -d 0 -p 0 -e 0 -v && mv {wildcards.sampleID}/{wildcards.upstream}/summary.PE.txt {output}
"""

rule simple_bridge_count_SE:
    input: '{sampleID}/{upstream}/chimera.fastq.gz'
    params:
        bridgef=lambda wildcards: samples_config["samples"][wildcards.sampleID].get("bridgeF",samples_config["bridgeF"]),
        bridger=lambda wildcards: samples_config["samples"][wildcards.sampleID].get("bridgeR",samples_config["bridgeR"]),
        L=config['read_length']
    output: '{sampleID}/{upstream}/bridge.stats.txt'
    shell: """
debridge.jl {params.bridgef} {params.bridger} {wildcards.sampleID}/{wildcards.upstream}/ {input} -l {params.L} -d 0 -p 0 -e 0 -v -s && mv {wildcards.sampleID}/{wildcards.upstream}/summary.PE.txt {output}
"""

rule deduplicate:
    input:
        PE1='{sampleID}/raw/chimera.1.fastq.gz',
        PE2='{sampleID}/raw/chimera.2.fastq.gz'
    output:
        PE1=temp('{sampleID}/chimeras/deduped/chimera.1.fastq.gz'),
        PE2=temp('{sampleID}/chimeras/deduped/chimera.2.fastq.gz')
    #benchmark: 'benchmarks/{sampleID}_preprocessing_2_deduplicate.tsv'
    log: '{sampleID}/chimeras/deduped/deduping.log'
    params:
        xmx=config["max_mem"]
    threads: config["max_threads"]
    shell: """
clumpify.sh in1={input.PE1} in2={input.PE2} out1={output.PE1} out2={output.PE2} dedupe=t subs=0 reorder=f overwrite=t -Xmx{params.xmx} -eoom 2>{log}
"""

# trims adapters
rule trim:
    input: 
        PE1='{sampleID}/chimeras/deduped/chimera.1.fastq.gz',
        PE2='{sampleID}/chimeras/deduped/chimera.2.fastq.gz'
    output: 
        PE1='{sampleID}/chimeras/deduped.trimmed/allPairs/chimera.1.fastq.gz',
        PE2='{sampleID}/chimeras/deduped.trimmed/allPairs/chimera.2.fastq.gz',
        readthrough='{sampleID}/chimeras/deduped.trimmed/readthrough/chimera.fastq.gz',
   
    log: '{sampleID}/chimeras/deduped.trimmed/trimming.log'
    params:
        trimmomatic_options=config["trimmomatic_options"],
        trimmomatic_filepath=srcdir(config["trimmomatic_filepath"]),
        trimfasta=lambda wildcards: SAMPLES[wildcards.sampleID]["clip"],
        xmx=config["max_mem"]
    threads: config["max_threads"]
    shell: """
java -XX:ParallelGCThreads={threads} -Xmx{params.xmx} -jar {params.trimmomatic_filepath} PE -threads {threads} -phred33 {input.PE1} {input.PE2} {output.PE1} {output.readthrough}_U1.fastq.gz {output.PE2} {output.readthrough}_U2.fastq.gz ILLUMINACLIP:{params.trimfasta}:{params.trimmomatic_options} 2>{log} && gzip -c <(gunzip -c {output.readthrough}_U1.fastq.gz) <(gunzip -c {output.readthrough}_U2.fastq.gz) > {output.readthrough} && rm {output.readthrough}_U1.fastq.gz && rm {output.readthrough}_U2.fastq.gz
"""

rule merge_mates_pear:
    input:
        PE1='{sampleID}/chimeras/deduped.trimmed/allPairs/chimera.1.fastq.gz',
        PE2='{sampleID}/chimeras/deduped.trimmed/allPairs/chimera.2.fastq.gz'
    output:
        M='{sampleID}/chimeras/deduped.trimmed/mergedPairs_pear/chimera.fastq.gz',
        UM1='{sampleID}/chimeras/deduped.trimmed/unmergedPairs_pear/chimera.1.fastq.gz',
        UM2='{sampleID}/chimeras/deduped.trimmed/unmergedPairs_pear/chimera.2.fastq.gz',
        dropped='{sampleID}/chimeras/deduped.trimmed/unmergedPairs_pear/dropped.fastq.gz',
        stats='{sampleID}/chimeras/deduped.trimmed/mergedPairs_pear/merging_pear.stats.txt'
    log: '{sampleID}/chimeras/deduped.trimmed/pear.log'
    params:
        pear_options=config["pear_options"]
    threads: config["max_threads"]
    shell: """
pear -f {input.PE1} -r {input.PE2} -o "$(dirname {output.M})"/pear {params.pear_options} -j {threads} 1>{output.stats} 2>{log} && (pigz -f -p {threads} "$(dirname {output.M})"/pear.assembled.fastq ; pigz -f -p {threads} "$(dirname {output.M})"/pear.unassembled.forward.fastq ; pigz -f -p {threads} "$(dirname {output.M})"/pear.unassembled.reverse.fastq; pigz -f -p {threads} "$(dirname {output.M})"/pear.discarded.fastq) && (mv "$(dirname {output.M})"/pear.assembled.fastq.gz {output.M} ; mv "$(dirname {output.M})"/pear.unassembled.forward.fastq.gz {output.UM1} ; mv "$(dirname {output.M})"/pear.unassembled.reverse.fastq.gz {output.UM2} ; mv "$(dirname {output.M})"/pear.discarded.fastq.gz {output.dropped}) 
"""

# DEBRIDGE
rule debridge_SE:
# Split the cDNA/DNA chimeric single-end reads into appropriate files
# TODO : describe this better!
    input: # chimera type: allPairs, mergedPairs, unmergebalePairs, mergeablePairs, readthrough
        chimera='{sampleID}/chimeras/deduped.trimmed/{chimera_typeS}/chimera.fastq.gz',
    params:
        bridgef=lambda wildcards: samples_config["samples"][wildcards.sampleID].get("bridgeF",samples_config["bridgeF"]),
        bridger=lambda wildcards: samples_config["samples"][wildcards.sampleID].get("bridgeR",samples_config["bridgeR"]),
        L=config['read_length']
    output: ['{sampleID}/split_chimeras/_debridging/{chimera_typeS}/' + x + '.fastq.gz' for x in readtype_SE]
    shell: """
debridge.jl {params.bridgef} {params.bridger} {wildcards.sampleID}/split_chimeras/_debridging/{wildcards.chimera_typeS}/ {input} -l {params.L} -s -d 2 -p 2 -e 0 -r -v >/dev/null && mv {wildcards.sampleID}/split_chimeras/_debridging/{wildcards.chimera_typeS}/summary.SE.txt {wildcards.sampleID}/split_chimeras/_debridging/{wildcards.chimera_typeS}/bridge.stats.txt

"""

# DEBRIDGE Paire-End reads
rule debridge_PE:
# Split the cDNA/DNA chimeric single-end reads into appropriate files
# TODO : describe this better!
    input: 
        PE1='{sampleID}/chimeras/deduped.trimmed/{chimera_typeP}/chimera.1.fastq.gz',
        PE2='{sampleID}/chimeras/deduped.trimmed/{chimera_typeP}/chimera.2.fastq.gz'
    params:
        bridgef=lambda wildcards: samples_config["samples"][wildcards.sampleID].get("bridgeF",samples_config["bridgeF"]),
        bridger=lambda wildcards: samples_config["samples"][wildcards.sampleID].get("bridgeR",samples_config["bridgeR"]),
        L=config['read_length']
    output:
        fq=['{sampleID}/split_chimeras/_debridging/{chimera_typeP}/' + x + '.fastq.gz' for x in readtype_PE]
    shell: """
debridge.jl {params.bridgef} {params.bridger} {wildcards.sampleID}/split_chimeras/_debridging/{wildcards.chimera_typeP}/ {input.PE1} {input.PE2} -l {params.L} -d 2 -p 2 -e 0 -r -v >/dev/null && mv {wildcards.sampleID}/split_chimeras/_debridging/{wildcards.chimera_typeP}/summary.PE.txt {wildcards.sampleID}/split_chimeras/_debridging/{wildcards.chimera_typeP}/bridge.stats.txt
    """
# chimera group: P(all PE), S (all SE), M (merged), U(unmergeable)
#
# === Makes a fastq .rna and fastq.dna master file ===

def get_mates_pairing_config(wildcards):
    p=config["mates_pairing_modes"][wildcards.mates_pairing_mode]
    rna=['{sampleID}/{base}/{fq}.fastq.gz'.format(sampleID=wildcards.sampleID, base=p['base'], fq=fq) for fq in p['rna']]
    dna=['{sampleID}/{base}/{fq}.fastq.gz'.format(sampleID=wildcards.sampleID, base=p['base'], fq=fq) for fq in p['dna']]

    return {'RNA':rna, 'DNA':dna}

rule consolidate_debridged_reads:
    input: unpack(get_mates_pairing_config)
    output: 
        RNA='{sampleID}/split_chimeras/{mates_pairing_mode}/unfiltered/rna.fastq.gz',
        DNA='{sampleID}/split_chimeras/{mates_pairing_mode}/unfiltered/dna.fastq.gz'
    run:
        allinputs=input
        allinputsRNA=" ".join(allinputs['RNA'])
        allinputsDNA=" ".join(allinputs['DNA'])
        shellCMD_RNA="zcat "+ allinputsRNA +" | pigz -c -p {threads} > {output.RNA}"
        shellCMD_DNA="zcat "+ allinputsDNA +" | pigz -c -p {threads} > {output.DNA}"
        print(shellCMD_RNA)
        shell(shellCMD_RNA)
        print(shellCMD_DNA)
        shell(shellCMD_DNA)
       
rule remove_short_reads_SE: #sends the good reads as paired somewhere, the short ones somewhere else
    input:
        RNA='{sampleID}/split_chimeras/{mates_pairing_mode}/unfiltered/rna.fastq.gz',
        DNA='{sampleID}/split_chimeras/{mates_pairing_mode}/unfiltered/dna.fastq.gz'
    output:
        fqs=expand('{{sampleID}}/split_chimeras/{{mates_pairing_mode}}/{filter}/{fq}.fastq.gz',filter=['long','shortR','shortD','shortRD'],fq=['rna','dna']),
        stats='{sampleID}/split_chimeras/{mates_pairing_mode}/filtering.stats.txt'
    params:
        cutoff_len_RNA=config["cutoff_len_RNA"]
    shell: """
filter_short_reads_pairs.sh {params.cutoff_len_RNA} {wildcards.sampleID}/split_chimeras/{wildcards.mates_pairing_mode} {input.RNA} {input.DNA}
"""

rule decon_RNA_SE:
    input: 
        rna='{sampleID}/split_chimeras/{mates_pairing_mode}/{filter}/rna.fastq.gz',
    output: 
        decon='{sampleID}/split_chimeras/{mates_pairing_mode}/{filter}.decon/rna.fastq.gz',
        contaminants='{sampleID}/split_chimeras/{mates_pairing_mode}/{filter}.contaminants/rna.fastq.gz',
        contaminants_bam='{sampleID}/split_chimeras/{mates_pairing_mode}/{filter}.contaminants/rna.bam',
        stats='{sampleID}/split_chimeras/{mates_pairing_mode}/{filter}.contaminants/rna.stats.txt'
    threads: config["max_threads"]
    params:
        decon_config=config["decon_config"]
    shell: """
bowtie2 -p{threads} {params.decon_config} --rg-id {wildcards.sampleID}:{wildcards.mates_pairing_mode}/{wildcards.filter}/rna --rg SM:{wildcards.sampleID} -U {input} 2>{output.stats} | picard SortSam I=/dev/stdin O={output.contaminants_bam} SORT_ORDER=queryname && decon_reads_pairs.sh {input.rna} {output.contaminants_bam} rna {wildcards.sampleID}/split_chimeras/{wildcards.mates_pairing_mode}/{wildcards.filter} {threads}
"""

rule decon_DNA_SE:
    input: 
        contaminants_bam='{sampleID}/split_chimeras/{mates_pairing_mode}/{filter}.contaminants/rna.bam',
        dna='{sampleID}/split_chimeras/{mates_pairing_mode}/{filter}/dna.fastq.gz'
    output: 
        decon='{sampleID}/split_chimeras/{mates_pairing_mode}/{filter}.decon/dna.fastq.gz',
        contaminants='{sampleID}/split_chimeras/{mates_pairing_mode}/{filter}.contaminants/dna.fastq.gz'
    threads: config["max_threads"]
    shell: """
decon_reads_pairs.sh {input.dna} {input.contaminants_bam} dna {wildcards.sampleID}/split_chimeras/{wildcards.mates_pairing_mode}/{wildcards.filter} {threads} 
"""

# finalize prealignment
rule revcomp_fastq:
    input: ['{sampleID}/{fq2revcomp}.fastq.gz']
    output: ['{sampleID}/{fq2revcomp}.revcomp.fastq.gz']
    threads: 1
    shell: """
revcompFQ.jl {input} {output}
    """

rule align_DNA_bowtie_SE:
    input: lambda wildcards: '{sampleID}/{fq}'.format(sampleID=wildcards.sampleID, fq=config['dna_alignment_modes'][wildcards.dna_alignment_mode]['fq'])
    params:
        bowtie_config=lambda wildcards: config['bowtie_configurations'].get(config["dna_alignment_modes"][wildcards.dna_alignment_mode]['config'], []).get('S',[])
    output:
        bam= '{sampleID}/alignments/dna/{dna_alignment_mode}/dna.bam',
        stats='{sampleID}/alignments/dna/{dna_alignment_mode}/dna.stats.txt'
    threads: config["max_threads"]
    shell: """
bowtie2 -p{threads} {params.bowtie_config} --rg-id {wildcards.sampleID} --rg SM:{wildcards.sampleID} -U {input} 2>{output.stats} | picard SortSam I=/dev/stdin O={output.bam} SORT_ORDER=queryname
"""

# # # === ALIGN SINGLE END DNA AND RNA READS TO GENOME ===
rule align_RNA_star_SE: #2g = align to genome
    input: lambda wildcards: '{sampleID}/{fq}'.format(sampleID=wildcards.sampleID, fq=config['rna_alignment_modes'][wildcards.rna_alignment_mode]['fq'])
    output:
        bam_transcriptome='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.transcriptome.bam',
        bam_genome='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.genome.bam'
    params:
        star_config=lambda wildcards: config['star_configurations'].get(config["rna_alignment_modes"][wildcards.rna_alignment_mode]['config'], []).get('S',[]),
        xmx=config["max_mem"]
    threads: config["max_threads"]
    shell: """
STAR --runThreadN {threads} --readFilesIn {input} \
--outSAMattrRGline ID:{wildcards.sampleID} SM:{wildcards.sampleID} \
--quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted --readFilesCommand zcat --outFileNamePrefix {wildcards.sampleID}/alignments/rna/{wildcards.rna_alignment_mode}/rna. {params.star_config} && picard -Xmx{params.xmx} SortSam I={wildcards.sampleID}/alignments/rna/{wildcards.rna_alignment_mode}/rna.Aligned.out.bam O={output.bam_genome} SORT_ORDER=queryname && picard -Xmx{params.xmx} SortSam I={wildcards.sampleID}/alignments/rna/{wildcards.rna_alignment_mode}/rna.Aligned.toTranscriptome.out.bam O={output.bam_transcriptome} SORT_ORDER=queryname && rm {wildcards.sampleID}/alignments/rna/{wildcards.rna_alignment_mode}/rna.Aligned.toTranscriptome.out.bam {wildcards.sampleID}/alignments/rna/{wildcards.rna_alignment_mode}/rna.Aligned.out.bam
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

rule geneannotate_noannot_bam:
    input: 
        bam_genes= lambda wildcards: '{sampleID}/alignments/rna/{rna_alignment_mode_genes}/rna.tag.bam'.format(sampleID=wildcards.sampleID, rna_alignment_mode_genes=config['rna_alignment_modes'][wildcards.rna_alignment_mode]['alignment_mode_genebodies']),
        bam_genes_tag= lambda wildcards: '{sampleID}/alignments/rna/{rna_alignment_mode_genes}/rna.tag.json'.format(sampleID=wildcards.sampleID, rna_alignment_mode_genes=config['rna_alignment_modes'][wildcards.rna_alignment_mode]['alignment_mode_genebodies']),
        bam_withannot= '{sampleID}/alignments/rna/{rna_alignment_mode}/rna.tag.bam'
        
    output:
        bam_genes='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.genes.tag.bam',
        bam_genes_tag= '{sampleID}/alignments/rna/{rna_alignment_mode}/rna.genes.tag.json',
        exons_id=temp('{sampleID}/alignments/rna/{rna_alignment_mode}/rna.tag.readids.txt')
    threads: 4
    params:
        xmx=str(4*config["mem_per_cpu"])+"g"
    shell: """
cp {input.bam_genes_tag} {output.bam_genes_tag}; samtools view {input.bam_withannot} | cut -f1 >{output.exons_id}; picard -Xmx{params.xmx} FilterSamReads I={input.bam_genes} READ_LIST_FILE={output.exons_id} FILTER=excludeReadList SORT_ORDER=queryname O={output.bam_genes}
 """

rule make_PREintergenic_bam:
    input: 
        bam_intergenic= lambda wildcards: '{sampleID}/alignments/rna/{rna_alignment_mode_genes}/rna.notag.bam'.format(sampleID=wildcards.sampleID, rna_alignment_mode_genes=config['rna_alignment_modes'][wildcards.rna_alignment_mode]['alignment_mode_genebodies']),
        exons_id= '{sampleID}/alignments/rna/{rna_alignment_mode}/rna.tag.readids.txt'
    output: temp('{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.PREintergenic.bam')
    threads: 4
    params:
        xmx=str(4*config["mem_per_cpu"])+"g"
    shell: """
picard -Xmx{params.xmx} FilterSamReads I={input.bam_intergenic} READ_LIST_FILE={input.exons_id} FILTER=excludeReadList SORT_ORDER=queryname O={output}
"""

rule make_intergenic_bam:
    input: '{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.PREintergenic.bam'
    output:
        intergenic='{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.intergenic.bam',
        ambiguous=temp('{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.ambiguous.bam')
    shell: """
samtools view -F4 -hb {input} >{output.intergenic}; samtools view -f4 -hb {input} >{output.ambiguous}
"""

rule make_PREnovel_bam:
    input:
        ambiguous='{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.ambiguous.bam',
        bam_genome='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.genome.bam'
    output:
        PREnovel=temp('{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.PREnovel.bam')
    params:
        xmx=str(4*config["mem_per_cpu"])+"g"
    shell: """
samtools view {input.ambiguous} | cut -f1 >$(dirname {input.ambiguous})/rna.ambiguous.readids.txt; picard -Xmx{params.xmx} FilterSamReads I={input.bam_genome} READ_LIST_FILE=$(dirname {input.ambiguous})/rna.ambiguous.readids.txt FILTER=includeReadList SORT_ORDER=queryname O={output.PREnovel}; rm $(dirname {input.ambiguous})/rna.ambiguous.readids.txt
"""

rule make_novel_bam:
    input: '{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.PREnovel.bam'
    output:
        novel='{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.novel.bam',
        unmapped='{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.unmapped.bam',
        toomany='{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.toomany.bam'
    shell: """
samtools view -F4 -F256 -bh {input} >{output.novel}; samtools view -f4 -h {input} |
 grep -v "uT:A:3" | samtools view -bh >{output.unmapped}; cat <(samtools view -f4 -H {input}) <(samtools view -f4 {input} | grep "uT:A:3") | samtools view -bh >{output.toomany}
"""

rule count_noalign:
    input:
        unmapped='{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.unmapped.bam',
        toomany='{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.toomany.bam'
    output: '{sampleID}/alignments/rna/{rna_alignment_mode}/noalign.stats.txt'
    shell: """
x2=$(/home/groups/astraigh/software/samtools-1.10/bin/samtools view -c {input.unmapped});  x1=$(/home/groups/astraigh/software/samtools-1.10/bin/samtools view -c {input.toomany}); echo -e "toomany,unmapped\n${{x1}},${{x2}}\n" >{output}
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

rule resolve_tag_bam_rna_genebodies:
    input: 
        tag_bam='{sampleID}/alignments/rna/{rna_alignment_mode}/rna.genes.tag.bam',
        quantfile=lambda wildcards: '{sampleID}/alignments/rna/{rna_alignment_mode_genes}/salmon_quant/quant.sf'.format(sampleID=wildcards.sampleID, rna_alignment_mode_genes=config['rna_alignment_modes'][wildcards.rna_alignment_mode]['alignment_mode_genebodies']),
        tag=lambda wildcards: '{sampleID}/alignments/rna/{rna_alignment_mode_genes}/rna.tag.json'.format(sampleID=wildcards.sampleID, rna_alignment_mode_genes=config['rna_alignment_modes'][wildcards.rna_alignment_mode]['alignment_mode_genebodies'])
    output:
        resolved_tag_bam='{sampleID}/alignments/rna/{rna_alignment_mode}/bytype/rna.introns.bam'
    threads: 1
    params:
        annot_file=lambda wildcards: config['star_configurations'].get(config["rna_alignment_modes"].get(config["rna_alignment_modes"][wildcards.rna_alignment_mode]['alignment_mode_genebodies'],[])['config'], [])['annot']
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

rule split_cistrans:
    input: 
        rna='{sampleID}/pairs/{pairing_mode}/{alignment_type}/paired.rna.bam',
        dna='{sampleID}/pairs/{pairing_mode}/{alignment_type}/paired.dna.bam'
    output:
        rna_cis='{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/rna.cis.bam',
        dna_cis='{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/dna.cis.bam',
        rna_trans='{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/rna.trans.bam',
        dna_trans='{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/dna.trans.bam',
        rna_ambig='{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/rna.ambiguous.bam',
        dna_ambig='{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/dna.ambiguous.bam'
    shell: """
chartools cistrans --rcis {output.rna_cis} --rtrans {output.rna_trans} --dcis {output.dna_cis} --dtrans {output.dna_trans} --rambig {output.rna_ambig} --dambig {output.dna_ambig} --qrna 255 --qdna 15 {input.rna} {input.dna}
"""


rule split_ambiguous:
    input: 
        rna='{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/rna.ambiguous.bam',
        dna='{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/dna.ambiguous.bam'
    output:
        rna_cis='{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/rna.dnaambig.bam',
        dna_cis='{sampleID}/pairs/{pairing_mode}/{alignment_type}/filtered/cistrans_Q255-Q15/dna.dnaambig.bam'
    shell: """
chartools deambig --rout {output.rna_cis} --dout {output.dna_cis} --qrna 255 {input.rna} {input.dna}
"""


rule make_bg_coverage:
    input: '{sampleID}/{bamfile}.sorted.bam'
    output: temp('{sampleID}/{bamfile}.{splitOR5}.{strandStrict}.bg')
    params: 
        strand_arg=lambda wildcards: "" if (wildcards['strandStrict']=="=") else ("-strand "+wildcards['strandStrict'])
    threads: 1
    shell: """
export LC_ALL="C"; bedtools genomecov -bg -{wildcards.splitOR5} {params.strand_arg} -ibam {input} | sort -k1,1 -k2,2n >{output}
"""

rule reverse_bg:
    input: '{sampleID}/{bamfile}.{splitOR5}.-.bg'
    output: temp('{sampleID}/{bamfile}.{splitOR5}.-rev.bg')
    threads: 1
    shell: """
awk -F $'\t' 'BEGIN{{OFS=FS}}{{print $1, $2, $3, -$4}}' {input} >{output}
"""

rule make_bw_coverage:
    input: '{sampleID}/{bedgraph}.bg'
    output: '{sampleID}/{bedgraph}.bw'
    params:
        chrNameLength=config['chromosomesFile']
    threads: 1
    shell: """
bedGraphToBigWig {input} {params.chrNameLength} {output}
"""

rule make_simple_pairs:
    input:
        rna='{sampleID}/pairs/{pairing_mode}/{alignment_type}/paired.rna.bam',
        dna='{sampleID}/pairs/{pairing_mode}/{alignment_type}/paired.dna.bam'
    output: '{sampleID}/pairs/{pairing_mode}/{alignment_type}/rd.simple.pairs'
    params:
        chrNameLength=config['chromosomesFile']
    threads: 1
    shell: """
chartools pairup_simple {input.rna} {input.dna} {output} {params.chrNameLength} --qrna 255 --qdna 15
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

rule gene_counts_table_fromDNA:
    input: '{sampleID}/pairs/{pairing_mode}/{pair_type}/{bed_folder}/dna.bed.gz'
    output: 
        complete='{sampleID}/pairs/{pairing_mode}/{pair_type}/{bed_folder}/stats/rnaTable.stats.txt',
        simple='{sampleID}/pairs/{pairing_mode}/{pair_type}/{bed_folder}/stats/rnaTableSimple.stats.txt'
    shell: """
export LC_ALL=C; awk -F $'\t' 'BEGIN{{OFS=" "}}{{sure=-1;istrans=-1; if($25==1){{sure=1;}}else{{if ($22==1){{sure=0;}};}}; if(($9>15) && (substr($2,3)==$4)){{istrans=0;}}else{{if($9>15){{istrans=1;}};}}; {{print $18, sure, istrans;}}}}' {input} | tee >(sort | uniq -c | sort -k1,1nr | sed -e 's/ *//' >{output.complete}) | cut -d " " -f1,1 |sort | uniq -c | sort -k1,1nr | sed -e 's/ *//' >{output.simple}
"""

rule gene_counts_summary:
    input: '{sampleID}/pairs/{pairing_mode}/{pair_type2}/stats/rnaTable.stats.txt'
    params: 
        txdb=config['annotations']['txdb'],
        genedb=config['annotations']['genedb'],
        chrdb=config['annotations']['chrdb']
    output: 
        all_parquet='{sampleID}/pairs/{pairing_mode}/{pair_type2}/stats/gene_expression_multimapRNA.parquet',
        all_xlsx='{sampleID}/pairs/{pairing_mode}/{pair_type2}/stats/gene_expression_multimapRNA.xlsx',
        unq_parquet='{sampleID}/pairs/{pairing_mode}/{pair_type2}/stats/gene_expression_NOmultimapRNA.parquet',
        unq_xlsx='{sampleID}/pairs/{pairing_mode}/{pair_type2}/stats/gene_expression_NOmultimapRNA.xlsx'
    threads: 1
    shell: """
chartools genestats3 {input} {params.chrdb} {params.txdb} {params.genedb} "$(dirname {output.all_parquet})"/gene_expression_
"""

def get_salmon_input(wildcards):
    fq1=config['salmon_alignment_modes'][wildcards.salmon_alignment_mode].get("fq1","")
    fq2=config['salmon_alignment_modes'][wildcards.salmon_alignment_mode].get("fq2","")
    if (len(fq1)>0) and (len(fq2)>0):
        return {'fq1': wildcards.sampleID+"/"+fq1, 'fq2': wildcards.sampleID+"/"+fq2}
    else:
        fq=config['salmon_alignment_modes'][wildcards.salmon_alignment_mode].get("fq","")
        return {'fq': wildcards.sampleID+"/"+fq}

rule pseudoalign_salmon:
    input: unpack(get_salmon_input)
    output:
        salmonquant='{sampleID}/alignments/{salmon_alignment_mode}/quant.sf'
    threads: config["max_threads"]
    params:
        index=lambda wildcards: config['salmon_alignment_modes'][wildcards.salmon_alignment_mode]['salmon_index'],
        salmon_params=lambda wildcards: config['salmon_alignment_modes'][wildcards.salmon_alignment_mode]['params']
    run:
        if len(input)>1:
            shell("salmon quant -p {threads} -i {params.index} {params.salmon_params} -1 {input.fq1} -2 {input.fq2} -o $(dirname {output.salmonquant})")
        else:
            shell("salmon quant -p {threads} -i {params.index} {params.salmon_params} -r {input.fq} -o $(dirname {output.salmonquant})")

rule salmon_summary:
    input: '{sampleID}/{upstream}/quant.sf'
    params:
        txdb=config['annotations']['txdb'],
        genedb=config['annotations']['genedb'],
        chrdb=config['annotations']['chrdb']
    output: 
        all_csv='{sampleID}/{upstream}/stats/gene_expression_salmon.csv',
        all_xlsx='{sampleID}/{upstream}/stats/gene_expression_salmon.xlsx',
        all_parquet='{sampleID}/{upstream}/stats/gene_expression_salmon.parquet'
    shell: """
chartools salmonstats {input} {params.chrdb} {params.txdb} {params.genedb} "$(dirname {output.all_csv})"/gene_expression_salmon
"""

# rule make_bed_rna:
#     input:
#         pairs='{sampleID}/pairs/{pairing_mode}/{pair_type}/rd.pairs'
#     output:
#         rna='{sampleID}/pairs/{pairing_mode}/{pair_type}/rna.bed.gz'
#     threads: 2
#     shell: """
# export LC_ALL=C; cat {input.pairs} | awk -F $'\t' 'BEGIN{{OFS=FS}}{{same=0; $2=substr($2,3); flight=$5-$3; if ($7=="-"){{$5=$5-3;}}; if($2==$4){{same=1;}}; print $2, $3, $3+1, $1, $8, $6, $4, $5, $5+1, $9, $7, $18, $19, $20, $21, $25, $26, same, flight;}}' | sort -k1,1 -k2,2n | bgzip --thread {threads} -c > {output.rna}; tabix -0 -f -p bed {output.rna}
# """

rule sort_bam:
    input: '{sampleID}/{bam}.bam'
    output: 
        bam='{sampleID}/{bam}.sorted.bam',
        index='{sampleID}/{bam}.sorted.bam.bai'
    threads: 2
    shell: """
samtools sort --threads {threads} -o {output.bam} {input}; samtools index {output.bam}
"""

rule make_bw:
    input: '{sampleID}/{bed}.{bedORbam}'
    output: 
        bg=temp('{sampleID}/{bed}.{bedORbam}.{sense}.bedgraph'),
        bw='{sampleID}/{bed}.{bedORbam}.{sense}.bw'
    params:
        chrdb=config['annotations']['chrdb']
    threads: 2
    run:
        if wildcards.bedORbam=="bam":
            if (wildcards.sense=="F"):
                shell("export LC_ALL=C; bedtools genomecov -ibam {input} -strand + -bg | sort -k1,1 -k2,2n >{output.bg}")
            elif (wildcards.sense=="R"):
                shell("export LC_ALL=C; bedtools genomecov -ibam {input} -strand - -bg | sort -k1,1 -k2,2n >{output.bg}")
            else:
                shell("export LC_ALL=C; bedtools genomecov -ibam {input} -bg | sort -k1,1 -k2,2n >{output.bg}")
        else:
            if (wildcards.sense=="F"):
                shell("export LC_ALL=C; bedtools genomecov -i <(gunzip -c {input}) -g {params.chrdb} -strand + -bg | sort -k1,1 -k2,2n >{output.bg}")
            elif (wildcards.sense=="R"):
                shell("export LC_ALL=C; bedtools genomecov -i <(gunzip -c {input}) -g {params.chrdb} -strand - -bg | sort -k1,1 -k2,2n >{output.bg}")
            else:
                shell("export LC_ALL=C; bedtools genomecov -i <(gunzip -c {input}) -g {params.chrdb} -bg | sort -k1,1 -k2,2n >{output.bg}")
        print("Running bedGraphToBigWig")
        shell("bedGraphToBigWig {output.bg} {params.chrdb} {output.bw}")