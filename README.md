# TeloMosa


![workflow](figure/TeloMosa_Schematic_diagram.svg)


TeloMosa is a bioinformatic tool for measuring telomere length and identifying TelC–TelG mosaic telomeric repeats (mTRs) from Oxford Nanopore long-read sequencing data.

## Requirements

1. Python >=3.8.0
2. minimap2 (v2.24)
3. samtools (v1.13)

## Workflows

1. Extract telomeric reads

   Use basecalled nanopore fastq.gz file as input, reads with >=2 TelC or TelG were extracted and saved as output.fq.gz

```bash
python3 extract_telomeric_reads_ONT_DNA.py \
    input.fastq.gz \
    -o output \
    --motif TTAGGG \
    --repeat 2 \
    --approx \
    --max-mismatch 2 \
    --step 6
```

Output:

```text
output.fq.gz
```
   Aligen extracted reads on T2T-CHM13 reference genome, and retention those mapped on end 500 kb.
   
```bash
#!/bin/bash

BED=./CHM13_chrEnd_500kb.bed # This file can be found in this repository (TeloMosa/material)

## CHM13v2.mmi were required
minimap2 \
    -ax map-ont \
    -t 28 \
    CHM13_v2.mmi \
   output.fq.gz \
    > output.sam

samtools view \
    -S -b -h \
    output.sam \
    > output.bam

samtools sort \
    -@ 16 \
    output.bam \
    -o output.sorted.bam

samtools index output.sorted.bam

samtools view \
    -b \
    -L ${BED} \
    output.sorted.bam \
    > output.2telo.500kb.bam

samtools index output.2telo.500kb.bam

samtools fastq \
    output.2telo.500kb.bam | \
    gzip -c > output_2telo_500kb.fastq.gz
```
2.2 Measure Telomere Length

Telomeric repeat length was estimated using a sliding-window approach. A 12-bp window, corresponding to two consecutive telomeric repeats, was scanned across each read with a step size of 6 bp. Both G-rich (TTAGGG) and C-rich (CCCTAA) telomeric repeats were considered. Up to two mismatches were allowed within each 12-bp window. Adjacent telomeric windows were merged into continuous telomeric blocks, and only blocks containing at least four consecutive telomeric repeats were retained. The lengths of all retained telomeric blocks within a read were summed and reported as the estimated telomeric repeat length of that read.

```bash
python3 TeloMoas_length_calculation.py \
    output_2telo_500kb.fastq.gz \
    -o output
```

Output:

```text
output.telomeric_repeat_length.tsv
```








