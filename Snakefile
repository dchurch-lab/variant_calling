"""
Author: Katarzyna Kedzierska
Affiliation: WCHG, Unieversity of Oxford
Aim: A simple Snakemake workflow for processing WES/WGS data.
Date: February 2020
Run: snakemake --use-conda -j NUM_OF_THREADS
"""

#TODO:
# automatize config.yaml creation
# automatize pattern recognition
# BUG: 'set +u'   in variants_annotated to avoid:
#      bash: THEME_SHOW_CLOCK: unbound variable
#      bash: RBFU_RUBY_VERSION: unbound variable

from glob import glob
import re
import os
from snakemake.logging import logger

configfile: "config.yaml"
SAMPLES_DICT = dict()
with open(config['SAMPLE_SHEET'], "r+") as fil:
    next(fil) # skip the header
    for lin in fil.readlines():
        row = lin.strip("\n").split("\t")
        sample_id = row[0]
        sample_name = row[1]
        if sample_name in SAMPLES_DICT.keys():
            SAMPLES_DICT[sample_name].append(sample_id)
        else:
            SAMPLES_DICT[sample_name] = [sample_id]


SAMPLES = list(SAMPLES_DICT.keys())
SAMPLE_IDS = [sample_id for sample in SAMPLES_DICT.values()
                for sample_id in sample]

rule all:
    input:
        # Filtered alignment
        #expand(os.path.join(config['FILTERED'], '{sample}.bam{ext}'),
        #       sample = SAMPLES,
        #       ext = ['', '.bai']),
        # Annotated variants
        expand(os.path.join(config['ANNOTATED_VARIANTS'], '{sample}.vcf.gz'),
              sample = SAMPLES),

       # multiqc report
       # "multiqc_report.html"

    message:
        '\n#################### WES/WGS pipeline #####################\n'
        'Running all necessary rules to produce complete output.\n'
        '############################################################'

rule annotate_variants:
    """Annotate variants with VEP
    """
    input:
        os.path.join(config['FILTERED_VARIANTS'], '{sample}.vcf.gz')
    output:
        os.path.join(config['ANNOTATED_VARIANTS'], '{sample}.vcf.gz')
    params:
        cache = config['VEP_CACHE']
    resources:
        mem_mb = 4096
    log:
        os.path.join(config['LOGS'], 'vep', '{sample}.log')
    conda:
        config['CONDA_VEP']
    message:
        '\n####################### Variant annotation ######################\n'
        'Annotate variants and saving to file: {output}\n'
        '############################################################'
    shell:
        'set +u; '
        'vep -i {input}'
        ' -o stdout'
        ' --cache --dir_cache {params.cache}'
        ' --force_overwrite --vcf --offline'
        ' --canonical --symbol --af 2> {log} |'
        ' gzip > {output}'

rule variant_filtering:
    """Filter variants with default filtering mode (mutect2)
    """
    input:
        variants = os.path.join(config['VARIANTS'], '{sample}.vcf.gz'),
        stats = os.path.join(config['VARIANTS'], '{sample}.vcf.gz.stats'),
        index = os.path.join(config['VARIANTS'], '{sample}.vcf.gz.tbi')
    output:
        variants = os.path.join(config['FILTERED_VARIANTS'], '{sample}.vcf.gz'),
        stats = os.path.join(config['FILTERED_VARIANTS'],
                             '{sample}.vcf.gz.filteringStats.tsv')
    params:
        ref = config['REFERENCE']
    resources:
        mem_mb = 4096
    log:
        os.path.join(config['LOGS'], 'mutect2', 'filtering_{sample}.log')
    conda:
        config['CONDA_MUTECT2']
    message:
        '\n################### Variant filtering ######################\n'
        'Filtering variants and saving to file: {output}\n'
        '############################################################'
    shell:
        'gatk --java-options \'-Xmx{resources.mem_mb}m\' FilterMutectCalls'
        ' -R {params.ref}'
        ' -V {input.variants} '
        ' -O {output.variants} 2> {log}'

#print('config["FILTERED"]'); #print(config['FILTERED'])

rule variant_calling:
    """Call variants with Mutect2 single sample mode
    """
    input:
        os.path.join(config['FILTERED'], '{sample}.bam')
    output:
        variants = os.path.join(config['VARIANTS'], '{sample}.vcf.gz'),
        stats = os.path.join(config['VARIANTS'], '{sample}.vcf.gz.stats'),
        index = os.path.join(config['VARIANTS'], '{sample}.vcf.gz.tbi')
    params:
        ref = config['REFERENCE']
    resources:
        mem_mb = 4096
    log:
        os.path.join(config['LOGS'], 'mutect2', '{sample}.log')
    conda:
        config['CONDA_MUTECT2']
    message:
        '\n##################### Variant calling ######################\n'
        'Variant calling and creating: {output.variants}\n'
        '############################################################'
    shell:
        'gatk --java-options \'-Xmx{resources.mem_mb}m\' Mutect2'
        ' -R {params.ref}'
        ' -I {input}'
        ' -O {output.variants} 2> {log}'

rule filter:
    """Clean up alignments
    Flags to filter out (-F):
      read unmapped (0x4 = 4)
      mate unmapped (0x8 = 8)
      not primary alignment (0x100 = 256)
      read fails platform/vendor quality checks (0x200 = 512)
      read is PCR or optical duplicate (0x400 = 1024)
      supplementary alignment (0x800 = 2048)

    Flags to keep (-f):
      read paired (0x1 = 1)
      read mapped in proper pair (0x2 = 2)

    Also, filtering out blackilsted regions.
    """
    input:
        bam = os.path.join(config['ALIGNMENT'], '{sample}_sorted_md.bam'),
        index = os.path.join(config['ALIGNMENT'], '{sample}_sorted_md.bam.bai')
    params:
        blacklist = config['BLACKLIST']
    output:
        bam = os.path.join(config['FILTERED'], '{sample}.bam'),
        index = os.path.join(config['FILTERED'], '{sample}.bam.bai')
    threads:
        8
    conda:
        config['CONDA_ALIGNMENT']
    message:
        '\n######################### Filtering ########################\n'
        'Filtering alignment followed by indexing\n'
        '{output.bam}\n'
        '{output.index}\n'
        '############################################################'
    shell:
        'samtools view -b -h -f 3 -F 3852 -@ {threads} {input.bam} '
        '$(for i in $(echo $(seq 22) X); do echo chr$i; done | xargs) | '
        'bedtools intersect -v -abam stdin -b {params.blacklist} > {output.bam}; '
        'samtools index {output.bam}'

rule mark_duplicates:
    input:
        lambda wildcards: os.path.join(config['ALIGNMENT'],
                                       SAMPLES_DICT[wildcards.sample][0] +
                                        '_sorted.bam')
    output:
        bam = os.path.join(config['ALIGNMENT'], '{sample}_sorted_md.bam'),
        index = os.path.join(config['ALIGNMENT'], '{sample}_sorted_md.bam.bai'),
        flagstat_metrics = os.path.join(config['ALIGNMENT_QUAL'],
                                        'flagstat_{sample}.txt'),
        md_metrics = os.path.join(config['ALIGNMENT_QUAL'],
                                 'DuplicationMetrics_{sample}.txt')
    params:
        tmp_dir = config['ALIGNMENT']
    resources:
        mem_mb = 4096
    message:
        '\n##################### Marking Duplicates ###################\n'
        'Running Picard MarkDuplicates followed by indexing\n'
        '{output.bam}\n'
        '{output.index}\n'
        '############################################################'
    log:
        os.path.join(config['LOGS'], 'MarkDuplicates', '{sample}.log')
    conda:
        config['CONDA_ALIGNMENT']
    shell:
        'picard -Xmx{resources.mem_mb}m -Xms{resources.mem_mb}m MarkDuplicates'
        ' I={input} O={output.bam} ASSUME_SORTED=true'
        ' METRICS_FILE={output.md_metrics}'
        ' TMP_DIR={params.tmp_dir} 2> {log};'
        ' samtools index {output.bam};'
        ' samtools flagstat {output.bam} > {output.flagstat_metrics}'

rule alignment:
    input:
        ref = config['REFERENCE'],
        forw = os.path.join(config['TRIMMED'], '{sample_id}_1_val_1.fq.gz'),
        rev = os.path.join(config['TRIMMED'], '{sample_id}_2_val_2.fq.gz'),
        amb = config['REFERENCE'] + '.amb',
        ann = config['REFERENCE'] + '.ann',
        bwt = config['REFERENCE'] + '.bwt',
        pac = config['REFERENCE'] + '.pac',
        sa = config['REFERENCE'] + '.sa'
    output:
        os.path.join(config['ALIGNMENT'], '{sample_id}_sorted.bam')
    params:
        sort_tmp = os.path.join(config['ALIGNMENT'], '{sample_id}_sort_tmp'),
        platform = config['SEQUENCING_PLATFORM'],
        sample_sheet = config['SAMPLE_SHEET']
    resources:
        mem_mb = 10240
    threads:
        4
    log:
        bwa = os.path.join(config['LOGS'], 'bwa', '{sample_id}.log'),
        samtools = os.path.join(config['LOGS'],
                                'samtools', 'sort_{sample_id}.log')
    conda:
        config['CONDA_ALIGNMENT']
    message:
        '\n######################### Mapping ##########################\n'
        'Running bwa_mem follow by sorting to produce:\n'
        '{output}\n'
        '############################################################'
    shell:
        'orgnl_id="$(basename {input.forw} | rev | cut -f4- -d\'_\' | rev)";'
        'sample_name=$(grep "$orgnl_id" {params.sample_sheet} | cut -f2); '
        'run_barcode=$(echo "$orgnl_id" | cut -f2 -d\'_\'); '
        'sequencing_centre=$(echo "$orgnl_id" | cut -f1 -d\'_\'); '
        'bwa mem -t {threads} {input.ref} {input.forw} {input.rev} '
        '-R "@RG\\tID:$orgnl_id\\t'
        'SM:$sample_name\\tPU:$run_barcode\\t'
        'CN:$sequencing_centre\\tPL:{params.platform}" '
        '2> {log.bwa} | '
        'samtools sort -m 1g -@ {threads} -O bam -T {params.sort_tmp} '
        '-o {output} - 2> {log.samtools}'

rule build_index:
    input:
        config['REFERENCE']
    output:
        config['REFERENCE'] + '.amb',
        config['REFERENCE'] + '.ann',
        config['REFERENCE'] + '.bwt',
        config['REFERENCE'] + '.pac',
        config['REFERENCE'] + '.sa'
    message:
        '\n######################### Indexing ##########################\n'
        'BWA index not found, running bwa index command:\n'
        'bwa index -a bwstw {input}\n'
        '############################################################'
    log:
        os.path.join(config["LOGS"], 'indexing.log')
    conda:
        config['CONDA_ALIGNMENT']
    shell:
        'bwa index -a bwtsw {input} 2> {log}'

rule trimming:
    input:
        os.path.join(config['READS'], '{sample_id}_1.fastq.gz'),
        os.path.join(config['READS'], '{sample_id}_2.fastq.gz')
    output:
        os.path.join(config['TRIMMED'], '{sample_id}_1_val_1.fq.gz'),
        os.path.join(config['TRIMMED'], '{sample_id}_2_val_2.fq.gz'),
        os.path.join(config['TRIMMED'],
                     '{sample_id}_1.fastq.gz_trimming_report.txt'),
        os.path.join(config['TRIMMED'],
                     '{sample_id}_2.fastq.gz_trimming_report.txt')
    params:
        qc_outdir = config["FASTQC"],
        outdir = config['TRIMMED']
    log:
        os.path.join(config["LOGS"], 'trim_galore', '{sample_id}.log')
    conda:
        config['CONDA_QUALITY']
    message:
        '\n######################### Trimming #########################\n'
        'Running trim_galore on:\n'
        '{input}\n'
        '############################################################'
    shell:
        'trim_galore --fastqc_args "--outdir {params.qc_outdir}" --gzip '
        '-o {params.outdir} --paired {input} &> {log}'

rule fastqc:
    input:
        os.path.join(config['READS'], '{sample_id}_1.fastq.gz'),
        os.path.join(config['READS'], '{sample_id}_2.fastq.gz')
    output:
        os.path.join(config['FASTQC'], '{sample_id}_1_fastqc.zip'),
        os.path.join(config['FASTQC'], '{sample_id}_2_fastqc.zip')
    params:
        outdir = config["FASTQC"]
    log:
        os.path.join(config["LOGS"], 'fastqc', '{sample_id}.log')
    conda:
        config['CONDA_QUALITY']
    shell:
        'fastqc {input} -o {params.outdir} &> {log}'

rule multiqc:
    input:
        # FASTQC output for RAW reads
        expand(os.path.join(config['FASTQC'],
                            '{sample_id}_{read}_fastqc.zip'),
               sample_id = SAMPLE_IDS,
               read = ['1', '2']),

        # Trimming reports
        expand(os.path.join(config['TRIMMED'],
                            '{sample_id}_{read}.fastq.gz_trimming_report.txt'),
               sample_id = SAMPLE_IDS,
               read = ['1', '2']),

        # Mark Duplicates metrcis
        expand(os.path.join(config['ALIGNMENT_QUAL'],
                            'DuplicationMetrics_{sample}.txt'),
               sample = SAMPLES),

        # Alignment Stats
        expand(os.path.join(config['ALIGNMENT_QUAL'],
                            'flagstat_{sample}.txt'),
               sample = SAMPLES)

    output:
        "multiqc_report.html"
    conda:
        config['CONDA_QUALITY']
    log:
        os.path.join(config["LOGS"], 'multiqc.log')
    message:
        '\n######################### Multiqc ##########################\n'
        'Running multiqc on all intermediate\n'
        'quality checks.\n'
        '############################################################'
    shell:
        'multiqc --force . --ignore .snakemake/ &> {log}'
