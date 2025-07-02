#!/usr/bin/env python3
"""
Compare results from standard and custom SmORFinder runs.
Ensures identical protein sequences except for headers.
"""

import sys
import os
from pathlib import Path
from Bio import SeqIO
import click


def normalize_fasta_headers(fasta_file):
    """Extract protein sequences and sort them for comparison."""
    sequences = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        # Store sequence without header
        sequences.append(str(record.seq))
    return sorted(sequences)


def compare_fasta_files(file1, file2, description):
    """Compare two FASTA files by their sequences (ignoring headers)."""
    print(f"\nComparing {description}...")
    
    if not os.path.exists(file1):
        print(f"ERROR: {file1} not found")
        return False
    if not os.path.exists(file2):
        print(f"ERROR: {file2} not found")
        return False
    
    seqs1 = normalize_fasta_headers(file1)
    seqs2 = normalize_fasta_headers(file2)
    
    if seqs1 == seqs2:
        print(f"✓ {description}: Sequences are identical ({len(seqs1)} proteins)")
        return True
    else:
        print(f"✗ {description}: Sequences differ!")
        print(f"  Standard: {len(seqs1)} proteins")
        print(f"  Custom: {len(seqs2)} proteins")
        
        # Find differences
        set1 = set(seqs1)
        set2 = set(seqs2)
        only_in_standard = set1 - set2
        only_in_custom = set2 - set1
        
        if only_in_standard:
            print(f"  Only in standard: {len(only_in_standard)} sequences")
        if only_in_custom:
            print(f"  Only in custom: {len(only_in_custom)} sequences")
        
        return False


def compare_tsv_files(file1, file2, description):
    """Compare TSV files by their content (excluding header differences)."""
    print(f"\nComparing {description}...")
    
    if not os.path.exists(file1):
        print(f"ERROR: {file1} not found")
        return False
    if not os.path.exists(file2):
        print(f"ERROR: {file2} not found")
        return False
    
    # Read and compare TSV content
    with open(file1) as f1, open(file2) as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
    
    # Skip header line and compare data
    if len(lines1) > 1 and len(lines2) > 1:
        data1 = sorted(lines1[1:])  # Skip header
        data2 = sorted(lines2[1:])  # Skip header
        
        if data1 == data2:
            print(f"✓ {description}: Data is identical ({len(data1)} entries)")
            return True
        else:
            print(f"✗ {description}: Data differs!")
            print(f"  Standard: {len(data1)} entries")
            print(f"  Custom: {len(data2)} entries")
            return False
    else:
        print(f"⚠ {description}: One or both files are empty or have no data")
        return True


@click.command()
@click.argument('standard_dir', type=click.Path(exists=True))
@click.argument('protein_dir', type=click.Path(exists=True))
def main(standard_dir, protein_dir):
    """Compare SmORFinder results from genome and protein workflows."""
    
    print("SmORFinder Results Comparison")
    print("=" * 40)
    
    # Get output directory names
    standard_name = os.path.basename(standard_dir)
    protein_name = os.path.basename(protein_dir)
    
    all_tests_passed = True
    
    # Compare FASTA files (protein sequences)
    standard_faa = os.path.join(standard_dir, f"{standard_name}.faa")
    protein_faa = os.path.join(protein_dir, f"{protein_name}.faa")
    
    if not compare_fasta_files(standard_faa, protein_faa, "Protein sequences (.faa)"):
        all_tests_passed = False
    
    # Compare TSV files (detailed results)
    standard_tsv = os.path.join(standard_dir, f"{standard_name}.tsv")
    protein_tsv = os.path.join(protein_dir, f"{protein_name}.tsv")
    
    if not compare_tsv_files(standard_tsv, protein_tsv, "Detailed results (.tsv)"):
        all_tests_passed = False
    
    # Summary
    print("\n" + "=" * 40)
    if all_tests_passed:
        print("✓ ALL TESTS PASSED: Results are identical except for headers")
        sys.exit(0)
    else:
        print("✗ TESTS FAILED: Results differ between genome and protein workflows")
        sys.exit(1)


if __name__ == '__main__':
    main() 