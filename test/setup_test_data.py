#!/usr/bin/env python3
"""
Setup test data for SmORFinder testing.
Downloads test genomes and runs Prodigal to generate .faa and .gff files.
"""

import os
import subprocess
import click
from pathlib import Path


def run_prodigal(genome_file, output_prefix):
    """Run Prodigal on a genome file."""
    cmd = [
        'prodigal',
        '-i', genome_file,
        '-a', f'{output_prefix}.faa',
        '-o', f'{output_prefix}.gff',
        '-f', 'gff'
    ]
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running Prodigal: {result.stderr}")
        return False
    
    print(f"✓ Generated {output_prefix}.faa and {output_prefix}.gff")
    return True


@click.command()
@click.option('--genome-dir', default='test_genomes', help='Directory containing genome FASTA files')
@click.option('--output-dir', default='test_data', help='Output directory for processed files')
def main(genome_dir, output_dir):
    """Setup test data for SmORFinder testing."""
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Find all .fna files in genome directory
    genome_files = list(Path(genome_dir).glob('*.fna'))
    
    if not genome_files:
        print(f"No .fna files found in {genome_dir}")
        print("Please place your test genome FASTA files in the genome directory.")
        return
    
    print(f"Found {len(genome_files)} genome files")
    
    for genome_file in genome_files:
        genome_name = genome_file.stem
        output_prefix = os.path.join(output_dir, genome_name)
        
        print(f"\nProcessing {genome_name}...")
        
        # Copy genome file to output directory
        genome_output = f"{output_prefix}.fna"
        if not os.path.exists(genome_output):
            subprocess.run(['cp', str(genome_file), genome_output])
        
        # Run Prodigal
        if run_prodigal(genome_output, output_prefix):
            print(f"✓ Successfully processed {genome_name}")
        else:
            print(f"✗ Failed to process {genome_name}")
    
    print(f"\nTest data setup complete!")
    print(f"Files are in: {output_dir}")
    print(f"\nYou can now run tests with:")
    print(f"smorf single {output_dir}/genome1.fna --outdir test_output_standard")
    print(f"smorf custom {output_dir}/genome1.fna {output_dir}/genome1.faa {output_dir}/genome1.gff --outdir test_output_custom")


if __name__ == '__main__':
    main() 