rule fastqc:
    input:
        os.path.join(config['READS'], '{sample}_1.fastq.gz'),
        os.path.join(config['READS'], '{sample}_2.fastq.gz')
    output:
        os.path.join(config['FASTQC'], '{sample}_1_fastqc.zip'),
        os.path.join(config['FASTQC'], '{sample}_2_fastqc.zip')
    params:
        outdir = config["FASTQC"]
    log:
        os.path.join(config["LOGS"], 'fastqc', '{sample}.log')
    conda:
        config['CONDA_QUALITY']
    shell:
        'fastqc {input} -o {params.outdir} &> {log}'

rule trimming:
    input:
        os.path.join(config['READS'], '{sample}_1.fastq.gz'),
        os.path.join(config['READS'], '{sample}_2.fastq.gz')
    output:
        os.path.join(config['TRIMMED'], '{sample}_1_val_1.fq.gz'),
        os.path.join(config['TRIMMED'], '{sample}_2_val_2.fq.gz'),
        os.path.join(config['TRIMMED'],
                     '{sample}_1.fastq.gz_trimming_report.txt'),
        os.path.join(config['TRIMMED'],
                     '{sample}_2.fastq.gz_trimming_report.txt')
    params:
        qc_outdir = config["FASTQC"],
        outdir = config['TRIMMED']
    log:
        os.path.join(config["LOGS"], 'trim_galore', '{sample}.log')
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
