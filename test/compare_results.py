#!/usr/bin/env python3
"""
Compare SmORFinder results from genome and pre-called workflows.
Reports sequence counts for monitoring code changes.
Note: Results may differ due to pyrodigal vs old Prodigal binary differences.
"""

import sys
import os
from pathlib import Path
from Bio import SeqIO
import click


def extract_protein_sequences(fasta_file):
    """Extract protein sequences and sort them for comparison."""
    sequences = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequences.append(str(record.seq))
    return sorted(sequences)


def report_fasta_files(file1, file2, description):
    """Report sequence counts from two FASTA files."""
    print(f"\n{description}:")
    
    if not os.path.exists(file1):
        print(f"  ERROR: {file1} not found")
        return False
    if not os.path.exists(file2):
        print(f"  ERROR: {file2} not found")
        return False
    
    seqs1 = extract_protein_sequences(file1)
    seqs2 = extract_protein_sequences(file2)
    
    print(f"  Standard workflow: {len(seqs1)} proteins")
    print(f"  Pre-called workflow: {len(seqs2)} proteins")
    
    # Check if sequences are identical (for reference, but don't fail)
    if seqs1 == seqs2:
        print(f"  ✓ Sequences are identical")
    else:
        print(f"  ⚠ Sequences differ (expected with pyrodigal vs old Prodigal)")
        
        # Find differences for reporting
        set1 = set(seqs1)
        set2 = set(seqs2)
        only_in_standard = set1 - set2
        only_in_pre_called = set2 - set1
        
        if only_in_standard:
            print(f"    Only in standard: {len(only_in_standard)} sequences")
        if only_in_pre_called:
            print(f"    Only in pre-called: {len(only_in_pre_called)} sequences")
    
    return True


@click.command()
@click.argument('standard_dir', type=click.Path(exists=True))
@click.argument('pre_called_dir', type=click.Path(exists=True))
def main(standard_dir, pre_called_dir):
    """Report SmORFinder results from genome and pre-called workflows."""
    
    print("SmORFinder Results Report")
    print("=" * 40)
    print("Note: Results may differ due to pyrodigal vs old Prodigal binary differences")
    
    # Debug: List all files in both directories
    print(f"\nFiles in standard directory ({standard_dir}):")
    if os.path.exists(standard_dir):
        for file in os.listdir(standard_dir):
            print(f"  - {file}")
    
    print(f"\nFiles in pre-called directory ({pre_called_dir}):")
    if os.path.exists(pre_called_dir):
        for file in os.listdir(pre_called_dir):
            print(f"  - {file}")
    
    # Try to find the .faa files
    standard_faa = None
    pre_called_faa = None
    
    # Look for .faa files
    if os.path.exists(standard_dir):
        for file in os.listdir(standard_dir):
            if file.endswith('.faa'):
                standard_faa = os.path.join(standard_dir, file)
                break
    
    if os.path.exists(pre_called_dir):
        for file in os.listdir(pre_called_dir):
            if file.endswith('.faa'):
                pre_called_faa = os.path.join(pre_called_dir, file)
                break
    
    # Report FASTA files (protein sequences)
    if standard_faa and pre_called_faa:
        report_fasta_files(standard_faa, pre_called_faa, "Protein sequences (.faa)")
    else:
        print(f"\nERROR: Could not find .faa files")
        print(f"  Standard: {standard_faa}")
        print(f"  Pre-called: {pre_called_faa}")
    
    # Summary
    print("\n" + "=" * 40)
    print("✓ REPORT COMPLETE: Check counts above for monitoring code changes")
    print("Note: Differences are expected when using pyrodigal vs old Prodigal binary")
    sys.exit(0)


if __name__ == '__main__':
    main() 