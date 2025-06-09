# mtb-mixed-infection-pipeline

mtb-mixed-infection-pipeline is a reproducible bioinformatics workflow for detecting mixed *Mycobacterium tuberculosis* infections from whole-genome sequencing (WGS) data. It integrates **Snippy**, **FreeBayes** (pooled-discrete mode), and a patched version of **MixInfect2.R** to process raw sequencing reads through alignment, joint variant calling, & statistical analysis. This pipeline outputs multi-sample VCF files with key genotype information suitable for identifying mixed-strain infections.

---

[![Conda](https://img.shields.io/conda/vn/conda-forge/r-mclust.svg)](https://anaconda.org/conda-forge/r-mclust)
[![R Version](https://img.shields.io/badge/R-4.0+-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)


---

## Reference genome
- **Organism**: *Mycobacterium tuberculosis* H37Rv  
- **GenBank accession**: `AL123456.3`

---

##  Running Snippy

### 1. Create conda environment for Snippy

```bash
conda env create -f snippy_env.yaml
conda activate snippy-env
```

### 2. Run Snippy on multiple samples

#### `sample.list.txt`

**Format:** `SampleID<TAB>Forward_Read<TAB>Reverse_Read`

```
TN106985	~/TN106985_1.fastq.gz	~/TN106985_2.fastq.gz
TN106727	~/TN106727_1.fastq.gz	~/TN106727_2.fastq.gz
TN106925	~/TN106925_1.fastq.gz	~/TN106925_2.fastq.gz
TN106439	~/TN106439_1.fastq.gz	~/TN106439_2.fastq.gz
```

```bash

./run_snippy_mtb.sh -i ~/sample.list.txt -r ~/AL123456_MTB_H37Rv.fasta -t 8

```

---

## FreeBayes joint-genotyping (pooled-discrete)

### 1. Prepare BAM files

```bash

samtools indexed sample bams '**snps.bam**' are available from Snippy run on multiple samples

```

### 2. Run FreeBayes

```bash
freebayes -f reference.fasta \
  -b sample1.bam -b sample2.bam -b sample3.bam \
  --pooled-discrete \
  --use-best-n-alleles 4 \
  --genotype-qualities \
  --min-alternate-fraction 0.01 \
  --min-alternate-count 2 \
  > multisample.vcf
```

Tip: You can generate an index with:

 ```bash
samtools faidx AL123456.fasta

 ```

---

## Detect mixed tuberculosis infections with MixInfect2

### 1. Create R environment
```bash
conda env create -f r-mixinfect.yaml
conda activate r-mixinfect

# If needed

conda install -c conda-forge r-data.table r-optparse r-ggplot2 r-mclust "r-base>=4.0" "icu=73.2"

```

### 2. Run MixInfect2

This script analyzes multi-strain Mycobacterium tuberculosis infections from VCF files using variant clustering and SNP-based heuristics.

```bash
(r-mixinfect):$ Rscript ~/MixInfect2/MixInfect2.R \
  --VCFfile ~/multisample.vcf \
  --prefix output \
  --maskFile ~/MixInfect2/MaskedRegions.csv \
  --minQual 10 \
  --useFilter FALSE
```

### Outputs

- `output_MixSampleSummary.csv`: Summary of sample classifications.
- `output_BICvalues.csv`: BIC scores & inferred strain counts.

---

## Loaded R packages

```
Loading required package: mclust
Package 'mclust' version 6.1
Type 'citation("mclust")' for citing this R package in publications.
Loading required package: stringr
Loading required package: optparse
Loading required package: foreach
Loading required package: doParallel
Loading required package: iterators
Loading required package: parallel
```

---

## Processing summary

```
Total variants read from VCF:          590,051
After removing indels:                 576,927
After applying QUAL > 10:                1,121
After removing spanning variants (*):    1,121
After masking regions:                    884
```

---

## Output summary table

| SampleName | Mix.Non.mix | hSNPs | Total.SNPs | Proportion.hSNPs_totalSNPs | No.strains | Major.strain.proportion |
|------------|-------------|-------|------------|-----------------------------|------------|--------------------------|
| TN106985   | Mix         |   459 |        517 | 88.78%                      | 2          | 0.78                     |
| TN106727   | Mix         |   477 |        505 | 94.46%                      | 2          | 0.71                     |
| TN106925   | Non-mix     |     6 |        149 | 4.03%                       | 1          | NA                       |
| TN106439   | Non-mix     |     9 |        252 | 3.57%                       | 1          | NA                       |

---
### License

This project is licensed under the MIT License.


### References
- [MixInfect2 GitHub](https://github.com/bensobkowiak/MixInfect2)
- [Snippy Documentation](https://github.com/tseemann/snippy)
- [FreeBayes GitHub](https://github.com/freebayes/freebayes)


