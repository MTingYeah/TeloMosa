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
   
