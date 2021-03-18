"""
Author: Kasia Kedzierska
Affiliation: WHG, Unieversity of Oxford
Aim: A simple Snakemake workflow for variant calling from single sample.
Date: February 2020 -
Run: snakemake --use-conda -j NUM_OF_THREADS
"""

#TODO:
# automatize config.yaml creation
# automatize pattern recognition
# BUG: 'set +u'   in variants_annotated to avoid:
#      bash: THEME_SHOW_CLOCK: unbound variable
#      bash: RBFU_RUBY_VERSION: unbound variable

# had to add '../' to conda env files in smk files
# asthe path to conda envs is interpretted relative to file

from glob import glob
import re
import os
from snakemake.logging import logger

configfile: "config.yaml"
SAMPLES = []
with open(config['SAMPLE_SHEET'], "r+") as sample_sheet:
    for line in sample_sheet:
      if not line.startswith('#'): # skip the header
        row = line.strip("\n").split("\t")
        sample = row[0]
        SAMPLES.append(sample)

rule all:
    input:
        # Filtered alignment
        expand(os.path.join(config['FILTERED'], '{sample}.bam{ext}'),
               sample = SAMPLES,
               ext = ['', '.bai']),
        # Annotated variants
        expand(os.path.join(config['ANNOTATED_VARIANTS'], '{sample}.vcf.gz'),
              sample = SAMPLES),
        # multiqc report
        "multiqc_report.html"

    message:
        '\n################# Variant calling pipeline #################\n'
        'Running all necessary rules to produce complete output.\n'
        '############################################################'

rule multiqc:
    input:
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


include: "utils/align.smk"
include: "utils/process_reads.smk"
include: "utils/variant_calling.smk"
