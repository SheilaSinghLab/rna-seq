from os.path import join
import os
from snakemake.io import glob_wildcards, expand
import re
import datetime

#Defining all relevant paths for analysis
FASTQ_DIR = 'data/raw_fastq/'
CLEAN_DIR = 'data/clean_fastq/'
ALIGNED = 'data/aligned/'
QC_DIR = 'data/reports/'
PATH_TO_GENOME = 'data/ref/Ensembl86.sjdbOverhang49'
REF_FLAT = 'data/ref/annotations/Homo_sapiens.GRCh38.Ensembl86.ref_flat.tsv'

#Regular expression to find patterns of file sample names 
pattern = '([A-Z][a-z]+)_(SS[0-9]{2})_[A-Z]*?_(L[0-9]{3})_(R[0-9])_[0-9]{3}'
FILES = glob_wildcards(join(FASTQ_DIR,'{sample}.fastq.gz')).sample
pattern = ('([A-Z][a-z]+)_(SS[0-9]{2})_[A-Z]*?_(L[0-9]{3})_(R[0-9])_[0-9]{3}')
today = datetime.date.today()
todaystr = today.isoformat()
#Generating SAMPLES python dictionary with keys=PatientNum, vals=Another dictionary where (Keys = Lane number, Val=file name)
SAMPLES = {}
for f in FILES:
    m = re.search(pattern, f)
    if m:
        try:
            SAMPLES[m.group(2)][m.group(3)] = f
        except KeyError:
            SAMPLES[m.group(2)] = {}
            SAMPLES[m.group(2)][m.group(3)]= f

#Rules that don't need to be submitted to a slurm job (Job manager for computer clusters)
localrules: merge_counts, all, qc_all, metrics_all

try:
    normMethod = config['nm'].lower()
    if normMethod != "tmm" and normMethod != "rpkm":
        normMethod = "tmm"
except KeyError:
    normMethod = "tmm"

#Current final output is a count_matrix of all the samples in a given raw folder
rule all:
    input: 
        counts = expand('data/counts/count_matrix_norm_{date}.txt', date = todaystr)

#Optional path to generate just quality control metrics
rule qc_all:
    input: expand(join(QC_DIR,'{files}_fastqc.html'), files=FILES)

rule quality_control:
    input: join(FASTQ_DIR,'{files}.fastq.gz')
    output: 
        html=join(QC_DIR, '{files}_fastqc.html'),
        z=join(QC_DIR, '{files}_fastqc.zip')
    threads: 4
    shell:"fastqc {input} --outdir ./data/reports"

rule trim:
    input:  join(FASTQ_DIR,'{sample}.fastq.gz')
    output: join(CLEAN_DIR,'{sample}_clean.fq')
    resources: cpus=4, mem_mb=4000
    shell: 
        "module load bbmap\n"
        "bbduk.sh -Xmx1g in={input} out={output} ref=adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=10 minlen=32 tpe tbo\n"
        "module unload bbmap"
rule star_align:
    input:
        fq_l1 =lambda wildcards: join(CLEAN_DIR, SAMPLES[wildcards.sample]['L001'])+ '_clean.fq',
        fq_l2 =lambda wildcards: join(CLEAN_DIR, SAMPLES[wildcards.sample]['L002']) + '_clean.fq',
        fa=PATH_TO_GENOME
    output: 
        bam = join(ALIGNED,'{sample}','{sample}_Aligned.sortedByCoord.out.bam'),
        counts =join(ALIGNED,'{sample}',"{sample}_ReadsPerGene.out.tab")
    resources: cpus=20, time_min=60, mem_mb=35000
    #log:join(ALIGNED, 'SS01','star.map.log')
    params:fq = lambda wildcards: ",".join(expand(join(CLEAN_DIR,'{sample}_clean.fq'), sample=SAMPLES))
    shell: 
        'module load star\n'
        'STAR --runThreadN 20'
        ' --genomeDir {input.fa}' 
        ' --readFilesIn {input.fq_l1},{input.fq_l2}'
        ' --outSAMstrandField intronMotif'
        ' --outSAMtype BAM SortedByCoordinate'
        ' --outFileNamePrefix ./data/aligned/{wildcards.sample}/{wildcards.sample}_'
        ' --quantMode GeneCounts'
        ' --sjdbGTFfile ./data/ref/annotations/Homo_sapiens.GRCh38.Ensembl86.gtf\n'
        'module unload star'

rule merge_counts:
    input:
        counts=expand(join(ALIGNED,'{sample}','{sample}_ReadsPerGene.out.tab'), sample=sorted(list(SAMPLES.keys())))
    output:expand('data/counts/count_matrix_{date}.txt', date=todaystr)
    shell: "python scripts/merge_counts.py {output} {input}"

rule rna_metrics:
    input:
        bam=join(ALIGNED,'{sample}','{sample}_Aligned.sortedByCoord.out.bam'),
        ref=REF_FLAT
    output:join(ALIGNED,'{sample}','{sample}_RNA_Metrics')
    threads:4
    shell:
        '''java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics \
        I={input.bam} \
        O={output} \
        REF_FLAT={input.ref} \
        STRAND=FIRST_READ_TRANSCRIPTION_STRAND'''

rule metrics_all:
    input: expand(join(ALIGNED,'{sample}','{sample}_RNA_Metrics'), sample=sorted(list(SAMPLES.keys())))

# rule that takes all counts and produces normalized counts by TPM
rule normalize:
    input:
         counts = expand('data/counts/count_matrix_{date}.txt',date=todaystr),
         rFile = 'scripts/normalize.R'
    output: expand('data/counts/count_matrix_norm_{date}.txt',date=todaystr)
    shell:
        '''
        module load gcc/7.3.0 r/3.6.0 r-bundle-bioconductor/3.9
        echo normalizing {input.counts} to create {output}
        Rscript {input.rFile} {normMethod} {input.counts} {output}
        module purge
        '''
