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

