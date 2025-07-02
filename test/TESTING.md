# SmORFinder Testing Guide

This guide explains how to test the custom SmORFinder command to ensure it produces identical results to the standard command (except for headers).

## Overview

The testing setup compares two SmORFinder runs:
1. **Standard run**: `smorf single genome.fna` (uses Prodigal internally)
2. **Custom run**: `smorf custom genome.fna proteins.faa genes.gff` (uses pre-predicted genes)

Both should produce identical protein sequences and results, differing only in the locus tag headers.

## Setup

### 1. Install Dependencies

```bash
# Install SmORFinder
pip install -e .

# Install system dependencies
sudo apt-get install prodigal hmmer

### 2. Prepare Test Data

Create a directory with test genome FASTA files:

```bash
mkdir test_genomes
# Add your .fna files here
```

Run the setup script to generate .faa and .gff files:

```bash
python test/setup_test_data.py --genome-dir test_genomes --output-dir genome_files
```

This will:
- Copy genome files to `genome_files/`
- Run Prodigal on each genome
- Generate corresponding .faa and .gff files

## Running Tests

### Manual Testing

```bash
# Run standard SmORFinder
smorf single genome_files/genome1.fna --outdir test_output_standard

# Run custom SmORFinder
smorf custom genome_files/genome1.fna genome_files/genome1.faa genome_files/genome1.gff --outdir test_output_custom

# Compare results
python test/compare_results.py test_output_standard test_output_custom
```

### Automated Testing (GitHub Actions)

The `.github/workflows/test.yml` file automatically runs tests on:
- Push to main/develop branches
- Pull requests

The workflow:
1. Sets up the environment
2. Downloads test data
3. Runs both standard and custom commands
4. Compares results
5. Uploads test artifacts

## Expected Results

### What Should Be Identical
- **Protein sequences**: The actual amino acid sequences should be identical
- **Number of proteins**: Both runs should identify the same number of small proteins
- **Model predictions**: DSN1 and DSN2 scores should be identical
- **HMM results**: SmORFam annotations should be identical

### What Will Differ
- **Locus tag headers**: The custom command preserves your original headers
- **File organization**: Internal file names may differ

## Test Output

The comparison script will output:

```
SmORFinder Results Comparison
========================================

Comparing Protein sequences (.faa)...
✓ Protein sequences (.faa): Sequences are identical (42 proteins)

Comparing Detailed results (.tsv)...
✓ Detailed results (.tsv): Data is identical (42 entries)

========================================
✓ ALL TESTS PASSED: Results are identical except for headers
```

## Troubleshooting

### Common Issues

1. **Different number of proteins**
   - Check that your .faa file contains all proteins (not just ≤50aa)
   - SmORFinder will filter to ≤51aa internally

2. **Missing files**
   - Ensure Prodigal generated both .faa and .gff files
   - Check that GFF IDs match FASTA headers

3. **Sequence differences**
   - Verify that your .faa file contains the same proteins as Prodigal would predict
   - Check for any preprocessing or filtering you may have applied

### Debug Mode

To see detailed output, run with verbose logging:

```bash
smorf single test_data/genome1.fna --outdir test_output_standard --force
smorf custom test_data/genome1.fna test_data/genome1.faa test_data/genome1.gff --outdir test_output_custom --force
```

## Adding More Test Genomes

1. Add genome FASTA files to `test_genomes/`
2. Run the setup script again
3. Update the GitHub Actions workflow if needed

## File Format Requirements

### Genome FASTA (.fna)
- Standard FASTA format
- Contig/chromosome sequences

### Protein FASTA (.faa)
- Standard FASTA format
- Headers must match GFF IDs
- Can contain all proteins (SmORFinder will filter)

### GFF File (.gff)
- Standard GFF3 format
- Must have `ID=` field in attributes
- IDs must match FASTA headers
- Should contain CDS features

## Example Test Data Structure

```
test_data/
├── genome1.fna
├── genome1.faa
├── genome1.gff
├── genome2.fna
├── genome2.faa
└── genome2.gff
``` 