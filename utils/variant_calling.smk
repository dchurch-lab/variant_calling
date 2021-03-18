localrules: install_vep_cache

rule call_variants:
    """Call variants with Mutect2 single sample mode
    """
    input:
        os.path.join(config['FILTERED'], '{sample}.bam')
    output:
        variants = os.path.join(config['VARIANTS'], '{sample}.vcf.gz'),
        stats = os.path.join(config['VARIANTS'], '{sample}.vcf.gz.stats'),
        index = os.path.join(config['VARIANTS'], '{sample}.vcf.gz.tbi')
    params:
        ref = config['REFERENCE'],
        germline = config['MUTECT2_GERMLINE'],
        pon = config['MUTECT2_PON']
    resources:
        mem_mb = 4096
    log:
        os.path.join(config['LOGS'], 'mutect2', '{sample}.log')
    conda:
        '../' + config['CONDA_MUTECT2']
    message:
        '\n##################### Variant calling ######################\n'
        'Variant calling and creating: {output.variants}\n'
        '############################################################'
    shell:
        'gatk --java-options \'-Xmx{resources.mem_mb}m\' Mutect2'
        ' -R {params.ref}'
        ' -I {input}'
        ' --germline-resource {params.germline}'
        ' --panel-of-normals {params.pon}'
        ' -O {output.variants} 2> {log}'

rule filter_variants:
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
        '../' + config['CONDA_MUTECT2']
    message:
        '\n################### Variant filtering ######################\n'
        'Filtering variants and saving to file: {output}\n'
        '############################################################'
    shell:
        'gatk --java-options \'-Xmx{resources.mem_mb}m\' FilterMutectCalls'
        ' -R {params.ref}'
        ' -V {input.variants} '
        ' -O {output.variants} 2> {log}'

rule install_vep_cache:
    """
    If run for the first time and the VEP cache doesn't exist, download
    and convert the cache.
    """
    output:
        # dummy file to make sure we can check if the cache has been downloaded
        os.path.join(config['VEP_CACHE'], 
                     "_".join([config['VEP_SPECIES'], config['VEP_VERSION'], 
                              config['VEP_ASSEMBLY'], "info.txt"]))
    params:
        cache = config['VEP_CACHE'],
        species = config['VEP_SPECIES'],
        assembly = config['VEP_ASSEMBLY'],
        version = config['VEP_VERSION']
    log:
        os.path.join(config['LOGS'], 'install_vep_cache.log')
    conda:
        '../' + config['CONDA_VEP']
    message:
        '\n################### Installing VEP Cache ###################\n'
        'Installing VEP cache for VEP {params.version}:\n'
        'organism: {params.species}; assembly: {params.assembly}.\n'
        '############################################################'
    shell:
        #'rm -frv {params.cache}/{params.species}/{params.version}_{params.assembly};'
        'vep_install --AUTO cf'
        ' --SPECIES {params.species}'
        ' --ASSEMBLY {params.assembly}'
        ' --VERSION {params.version}'
        ' --CACHEDIR {params.cache}'
        ' --CONVERT;'
        ' cp {params.cache}/{params.species}/{params.version}_{params.assembly}/info.txt {output}'

rule annotate_variants:
    """Annotate variants with VEP
    """
    input:
        vcf = os.path.join(config['FILTERED_VARIANTS'], '{sample}.vcf.gz'),
        # this is the path for a dummy file to check if we have the proper cache
        cache_info = os.path.join(config['VEP_CACHE'], 
                                  "_".join([config['VEP_SPECIES'], config['VEP_VERSION'], 
                                            config['VEP_ASSEMBLY'], "info.txt"]))
    output:
        os.path.join(config['ANNOTATED_VARIANTS'], '{sample}.vcf.gz')
    resources:
        mem_mb = 4096
    params:
        cache = config['VEP_CACHE']
    log:
        os.path.join(config['LOGS'], 'vep', '{sample}.log')
    conda:
        '../' + config['CONDA_VEP']
    message:
        '\n##################### Variant annotation ####################\n'
        'Annotate variants and saving to file: {output}\n'
        '############################################################'
    shell:
        'set +u; '
        'vep -i {input.vcf}'
        ' -o stdout'
        ' --cache --dir_cache {params.cache}'
        ' --force_overwrite --vcf --offline'
        ' --canonical --symbol --af 2> {log} |'
        ' gzip > {output}'

