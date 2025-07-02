#!/usr/bin/env python3
"""
Setup test genome files for SmORFinder testing.
Runs Prodigal on genome FASTA files and organizes them in test/genome_files/.
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


def verify_files(genome_name, test_dir):
    """Verify that all required files exist and are valid."""
    faa_file = test_dir / f"{genome_name}.faa"
    gff_file = test_dir / f"{genome_name}.gff"
    
    if not faa_file.exists():
        print(f"✗ Missing {faa_file}")
        return False
    
    if not gff_file.exists():
        print(f"✗ Missing {gff_file}")
        return False
    
    # Check that GFF has ID fields
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            if 'ID=' in line:
                print("✓ GFF file has ID fields")
                break
        else:
            print("✗ GFF file missing ID fields")
            return False
    
    print(f"✓ All files verified for {genome_name}")
    return True


@click.command()
@click.option('--input-dir', default='.', help='Directory containing genome FASTA files')
@click.option('--output-dir', default='test/genome_files', help='Output directory for test files')
@click.option('--genome-pattern', default='*.fna', help='Pattern to match genome files')
def main(input_dir, output_dir, genome_pattern):
    """Setup test genome files for SmORFinder testing."""
    
    # Create output directory
    test_dir = Path(output_dir)
    test_dir.mkdir(parents=True, exist_ok=True)
    
    # Find genome files
    input_path = Path(input_dir)
    genome_files = list(input_path.glob(genome_pattern))
    
    if not genome_files:
        print(f"No genome files found matching pattern '{genome_pattern}' in {input_dir}")
        print("Please place your test genome FASTA files in the input directory.")
        return
    
    print(f"Found {len(genome_files)} genome files")
    
    for genome_file in genome_files:
        genome_name = genome_file.stem
        output_prefix = test_dir / genome_name
        
        print(f"\nProcessing {genome_name}...")
        
        # Copy genome file to test directory
        genome_output = test_dir / f"{genome_name}.fna"
        if not genome_output.exists():
            import shutil
            shutil.copy(genome_file, genome_output)
            print(f"✓ Copied {genome_file} to {genome_output}")
        
        # Run Prodigal
        if run_prodigal(str(genome_output), str(output_prefix)):
            if verify_files(genome_name, test_dir):
                print(f"✓ Successfully processed {genome_name}")
            else:
                print(f"✗ Verification failed for {genome_name}")
        else:
            print(f"✗ Failed to process {genome_name}")
    
    print(f"\nTest genome setup complete!")
    print(f"Files are in: {test_dir}")
    print(f"\nYou can now run tests with:")
    print(f"smorf single {test_dir}/genome1.fna --outdir test_output_standard")
    print(f"smorf custom {test_dir}/genome1.fna {test_dir}/genome1.faa {test_dir}/genome1.gff --outdir test_output_custom")


if __name__ == '__main__':
    main() 