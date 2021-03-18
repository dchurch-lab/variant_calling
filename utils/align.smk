rule align:
    input:
        ref = config['REFERENCE'],
        forw = os.path.join(config['TRIMMED'], '{sample}_1_val_1.fq.gz'),
        rev = os.path.join(config['TRIMMED'], '{sample}_2_val_2.fq.gz'),
        amb = config['REFERENCE'] + '.amb',
        ann = config['REFERENCE'] + '.ann',
        bwt = config['REFERENCE'] + '.bwt',
        pac = config['REFERENCE'] + '.pac',
        sa = config['REFERENCE'] + '.sa'
    output:
        os.path.join(config['ALIGNMENT'], '{sample}_sorted.bam')
    params:
        sort_tmp = os.path.join(config['ALIGNMENT'], '{sample}_sort_tmp'),
        platform = config['SEQUENCING_PLATFORM'],
        sample_sheet = config['SAMPLE_SHEET']
    resources:
        mem_mb = 10240
    threads:
        4
    log:
        bwa = os.path.join(config['LOGS'], 'bwa', '{sample}.log'),
        samtools = os.path.join(config['LOGS'],
                                'samtools', 'sort_{sample}.log')
    conda:
        '../' + config['CONDA_ALIGNMENT']
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
        '../' + config['CONDA_ALIGNMENT']
    shell:
        'bwa index -a bwtsw {input} 2> {log}'

rule mark_duplicates:
    input:
        os.path.join(config['ALIGNMENT'], '{sample}_sorted.bam')
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
        '../' + config['CONDA_ALIGNMENT']
    shell:
        'picard -Xmx{resources.mem_mb}m -Xms{resources.mem_mb}m MarkDuplicates'
        ' I={input} O={output.bam} ASSUME_SORTED=true'
        ' METRICS_FILE={output.md_metrics}'
        ' TMP_DIR={params.tmp_dir} 2> {log};'
        ' samtools index {output.bam};'
        ' samtools flagstat {output.bam} > {output.flagstat_metrics}'

rule filter_alignment:
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
        '../' + config['CONDA_ALIGNMENT']
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
