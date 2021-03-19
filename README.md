# WES pipeline

This directory hosts the pipeline for WES.

[![DOI](https://zenodo.org/badge/327455301.svg)](https://zenodo.org/badge/latestdoi/327455301)


## Steps of the pipeline

![Outline for a pipeline](dag.png)

1. Quality control of the reads.
2. Reads trimming, followed by QC.
3. Aligning the reads with bwa.
4. Marking duplicate alignments.
5. Quality based filtering.
6. Calling variants.
7. Filtering variants.
8. Annotating variants.

## Prerequisites

This pipleine relies on **conda** and obviosuly **snakemake**. All the necessary software is stored in the environments, therefore there is no need of installing them on their own; they will be downloaded and installed first time the pipeline is run.

### Sample sheet
The pipeline requires the sample_sheet.tsv to have the sample_id in first column (the ID of the fastq file) and second column to store the name of the sample that one whishes to call their samples. At marking duplicates step the sample naming switches from sample_ids to sample_names.

### Raw reads
The pipeline starts with raw reads which need to be stored in read/raw. The naming assumes sample_id_[12]{1}.fastq.gz.

## Usage

Best way of setting up this pipeline is to copy the content of this directory into the directory where one wants to store the results of the pipeline.

I also recommend first the dry run:

```bash
snakemake -np --use-conda
```

then, running it on one sample by running this:

```bash
head -n2 sample_sheet.tsv > sample_sheet_one.tsv
sed 's/sample_sheet.tsv/sample_sheet_one.tsv/g' config.yaml > tmp && mv tmp config.yaml
snakemake -np --use-conda
snakemake -p --use-conda
```

When successful, run the whole pipeline (*-j 16* is an option to run the pipeline with 16 cores, remember to change accordingly):

```bash
rm -f sample_sheet_one.tsv
sed 's/sample_sheet_one.tsv/sample_sheet.tsv/g' config.yaml > tmp && mv tmp config.yaml
snakemake -p --use-conda -j 16
```
