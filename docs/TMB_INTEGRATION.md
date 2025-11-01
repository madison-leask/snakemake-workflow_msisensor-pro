# TMB Integration Guide

## Overview

This document provides guidance on integrating Tumor Mutational Burden (TMB) calculation into the MSIsensor-pro workflow.

## Background

**Tumor Mutational Burden (TMB)** is the total number of somatic mutations per coding area of a tumor genome, typically measured in mutations per megabase (mut/Mb). It is an emerging biomarker for predicting response to immune checkpoint inhibitor therapy.

**Microsatellite Instability (MSI)** and TMB are complementary biomarkers:
- MSI-High tumors often have high TMB
- Both predict immunotherapy response
- They can be analyzed independently or together

## Current State

The MSIsensor-pro workflow currently provides:
- ✅ MSI detection and classification
- ✅ Analysis notebooks for MSI results
- ✅ Framework for combined MSI/TMB analysis
- ⚠️ TMB calculation: **NOT YET IMPLEMENTED** (simulated in notebooks for demonstration)

## TMB Calculation Requirements

To add real TMB calculation, you need:

### 1. Variant Calling
Generate somatic variant calls in VCF format using a variant caller such as:
- **Mutect2** (GATK) - recommended for tumor-only and tumor-normal
- **VarScan2** - good for tumor-normal pairs
- **Strelka2** - fast, accurate for tumor-normal
- **LoFreq** - good for low-frequency variants

### 2. Variant Filtering
Filter variants to include only:
- High-quality somatic mutations (not germline)
- Non-synonymous mutations in coding regions
- Remove common polymorphisms (filter against dbSNP, gnomAD)
- Apply quality filters (depth, VAF, strand bias, etc.)

### 3. Callable Genome Size
Determine the callable/sequenceable genome size:
- **Whole Exome Sequencing (WES)**: typically 30-50 Mb
- **Whole Genome Sequencing (WGS)**: typically ~3000 Mb
- **Targeted panels**: varies by panel (e.g., 1-2 Mb)

### 4. TMB Calculation
```
TMB = (Number of somatic mutations) / (Callable genome size in Mb)
```

## Implementation Plan

### Step 1: Add Variant Calling Rule

Create a new rule file `workflow/rules/variant_calling.smk`:

```python
rule mutect2_calling:
    input:
        tumor="results/recal/{sample}.bam",
        normal="results/recal/{normal}.bam",  # if tumor-normal mode
        ref="resources/{genome_version}.fasta",
        ref_dict="resources/{genome_version}.dict",
    output:
        vcf="results/variants/{sample}.mutect2.vcf.gz",
        stats="results/variants/{sample}.mutect2.vcf.gz.stats",
    log:
        "logs/variants/{sample}.mutect2.log"
    params:
        extra="--germline-resource resources/af-only-gnomad.vcf.gz",
    threads: 4
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk Mutect2 \\
            -R {input.ref} \\
            -I {input.tumor} \\
            -tumor {wildcards.sample} \\
            -O {output.vcf} \\
            {params.extra} \\
            > {log} 2>&1
        """
```

### Step 2: Add Variant Filtering Rule

```python
rule filter_variants:
    input:
        vcf="results/variants/{sample}.mutect2.vcf.gz",
        ref="resources/{genome_version}.fasta",
    output:
        filtered_vcf="results/variants/{sample}.filtered.vcf.gz",
    log:
        "logs/variants/{sample}.filter.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk FilterMutectCalls \\
            -R {input.ref} \\
            -V {input.vcf} \\
            -O {output.filtered_vcf} \\
            > {log} 2>&1
        """
```

### Step 3: Add TMB Calculation Script

Create `workflow/scripts/calculate_tmb.py`:

```python
#!/usr/bin/env python3
"""Calculate Tumor Mutational Burden from VCF file."""

import sys
import gzip
from pathlib import Path

def count_somatic_mutations(vcf_file):
    """Count PASS somatic mutations in VCF."""
    count = 0
    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            # Check FILTER field (7th column, 0-indexed as 6)
            if fields[6] == 'PASS':
                count += 1
    return count

def calculate_tmb(vcf_file, callable_mb, output_file):
    """Calculate TMB and write to file."""
    mutation_count = count_somatic_mutations(vcf_file)
    tmb = mutation_count / callable_mb
    
    with open(output_file, 'w') as f:
        f.write("sample\ttotal_mutations\tcallable_mb\ttmb\n")
        sample_name = Path(vcf_file).stem
        f.write(f"{sample_name}\t{mutation_count}\t{callable_mb}\t{tmb:.2f}\n")
    
    return tmb

if __name__ == '__main__':
    vcf = snakemake.input.vcf
    callable_mb = float(snakemake.params.callable_mb)
    output = snakemake.output.tmb
    
    tmb = calculate_tmb(vcf, callable_mb, output)
    print(f"TMB calculated: {tmb:.2f} mutations/Mb")
```

### Step 4: Add TMB Calculation Rule

```python
rule calculate_tmb:
    input:
        vcf="results/variants/{sample}.filtered.vcf.gz",
    output:
        tmb="results/tmb/{sample}.tmb.tsv",
    params:
        callable_mb=30.0,  # Adjust based on sequencing type (WES/WGS/panel)
    log:
        "logs/tmb/{sample}.tmb.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/calculate_tmb.py"
```

### Step 5: Add TMB Merging Rule

```python
rule merge_tmb_results:
    input:
        tmb_files=expand(
            "results/tmb/{sample}.tmb.tsv",
            sample=get_tumor_samples()
        ),
    output:
        merged="results/tmb.all_samples.tsv",
    log:
        "logs/tmb.merge.log"
    run:
        import pandas as pd
        dfs = [pd.read_csv(f, sep='\\t') for f in input.tmb_files]
        merged = pd.concat(dfs, ignore_index=True)
        merged.to_csv(output.merged, sep='\\t', index=False)
```

### Step 6: Create GATK Environment

Create `workflow/envs/gatk.yaml`:

```yaml
name: gatk
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - gatk4>=4.3
  - samtools>=1.15
  - bcftools>=1.15
```

### Step 7: Update Workflow to Include TMB

In `workflow/Snakefile`, add:

```python
include: "rules/variant_calling.smk"
```

Update the `get_final_output()` function in `workflow/rules/common.smk`:

```python
def get_final_output(wildcards):
    final_output = []
    
    # MSI results
    final_output.extend([
        f"results/{workflow_mode}.{genome_version}.all_samples.tsv"
    ])
    
    # TMB results (if variant calling enabled)
    if config.get("enable_tmb", False):
        final_output.extend([
            "results/tmb.all_samples.tsv"
        ])
    
    return final_output
```

### Step 8: Update Configuration

Add to `config/config.yaml`:

```yaml
# TMB calculation settings (optional)
enable_tmb: true
tmb:
  # Callable genome size in megabases
  # - WES (whole exome): ~30-50 Mb
  # - WGS (whole genome): ~3000 Mb
  # - Custom panel: calculate from BED file
  callable_mb: 35.0
  
  # Variant caller to use
  caller: mutect2  # options: mutect2, varscan, strelka
  
  # Germline resource for filtering (optional)
  germline_resource: "resources/af-only-gnomad.vcf.gz"
```

## Testing TMB Integration

1. **Dry run**: Test the workflow without running
   ```bash
   snakemake --dry-run --cores 1
   ```

2. **Small test**: Run on a subset of data
   ```bash
   snakemake --cores 4 --sdm conda results/tmb/test_sample.tmb.tsv
   ```

3. **Full workflow**: Run complete analysis
   ```bash
   snakemake --cores 8 --sdm conda
   ```

## Updating the Notebooks

Once TMB is calculated:

1. Update `workflow/notebooks/tmb_msi_combined_analysis.ipynb`
2. Replace simulated TMB loading with:
   ```python
   # Load real TMB results
   tmb_file = paths['results'] / 'tmb.all_samples.tsv'
   tmb_results = pd.read_csv(tmb_file, sep='\\t')
   
   # Merge with MSI results
   combined = msi_results.merge(tmb_results, on='sample')
   ```

3. Remove simulation code and warnings about simulated data

## References

- [GATK Best Practices for Somatic Variant Calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731)
- [TMB as a Biomarker (Yarchoan et al., 2017)](https://doi.org/10.1056/NEJMoa1609279)
- [FDA Approval of Pembrolizumab for TMB-H](https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approves-pembrolizumab-adults-and-children-tmb-h-solid-tumors)

## Support

For questions or issues with TMB integration, please open an issue in the repository.
