# SmORFinder Testing Guide

This guide explains how to test SmORFinder to ensure both workflows produce identical results.

## Overview

SmORFinder has two main workflows:

1. **Standard run**: `smorf single genome.fna` (runs Prodigal internally)
2. **Pre-called run**: `smorf pre-called genome.fna proteins.faa nucleotides.ffn genes.gff` (uses pre-predicted genes)

Both should produce identical protein sequences and results, differing only in the locus tag headers.

## Test Setup

### 1. Prepare Test Data

Run the setup script to generate test files:

```bash
python test/setup_test_data.py --genome-dir test_genomes --output-dir test_data
```

This will:
- Clean genome FASTA headers
- Run Prodigal on each genome
- Rename proteins to match SmORFinder's format
- Generate `.faa`, `.ffn`, and `.gff` files

### 2. Run Both Workflows

```bash
# Standard workflow (runs Prodigal internally)
smorf single test_data/genome1.fna --outdir test_output_standard --force

# Pre-called workflow (uses your Prodigal output)
smorf pre-called test_data/genome1.fna test_data/genome1.faa test_data/genome1.ffn test_data/genome1.gff --outdir test_output_pre_called --force
```

### 3. Compare Results

```bash
python test/compare_results.py test_output_standard test_output_pre_called
```

## Expected Results

The comparison should show:

- **Protein sequences**: The actual amino acid sequences should be identical
- **Number of proteins**: Both runs should identify the same number of small proteins

Note: Detailed results (TSV files) are not compared because probability scores can vary slightly between runs due to numerical precision differences, even with identical sequences.

## Example Output

```
Comparing results:
  Standard workflow: test_output_standard
  Pre-called workflow: test_output_pre_called

Files in standard directory (test_output_standard):
  - 160.faa
  - 160.ffn
  - 160.gff
  - 160.tsv

Files in pre-called directory (test_output_pre_called):
  - 160.faa
  - 160.ffn
  - 160.gff
  - 160.tsv

Comparing Protein sequences (.faa)...
✓ Protein sequences (.faa): Sequences are identical (42 proteins)

==================================================
✓ ALL TESTS PASSED: Protein sequences and counts are identical
```

## Troubleshooting

### Common Issues

1. **Different number of proteins**
   - Check that your .faa file contains all proteins (not just ≤50aa)
   - SmORFinder filters to ≤50aa internally
   - Verify that your .faa file contains the same proteins as Prodigal would predict

2. **Different protein sequences**
   - Ensure you're using the same Prodigal version and parameters
   - Check that your .faa and .gff files are from the same Prodigal run
   - Verify that protein IDs match between .faa and .gff files

3. **Missing files**
   - Ensure all required files exist: `.fna`, `.faa`, `.ffn`, `.gff`
   - Check file permissions and paths

### File Requirements

#### Genome FASTA (.fna)
- Must be a valid FASTA file
- Headers should be clean (no spaces or special characters)

#### Protein FASTA (.faa)
- Must be from Prodigal output
- Headers must match GFF IDs exactly
- Can contain all proteins (SmORFinder will filter)

#### Nucleotide FASTA (.ffn)
- Must be from Prodigal output
- Headers must match GFF IDs exactly
- Can contain all genes (SmORFinder will filter)

#### GFF File (.gff)
- Must be from Prodigal output
- Must have ID field in attributes column
- IDs must match protein and nucleotide headers

## Automated Testing

The GitHub Actions workflow automatically runs these tests on multiple genomes:

1. Downloads test genomes
2. Runs both workflows
3. Compares results
4. Reports any differences

Check the Actions tab in your repository to see test results. 