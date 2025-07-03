#!/usr/bin/env python3
"""
Setup test data for SmORFinder testing.
Downloads test genomes and runs Pyrodigal to generate .faa and .gff files.
"""

import os
import subprocess
import click
from pathlib import Path
from Bio import SeqIO
from random import choices
import string


def clean_fasta_headers(input_file, output_file):
    """Clean FASTA headers by removing everything after the first whitespace."""
    cleaned_records = []
    
    for record in SeqIO.parse(input_file, 'fasta'):
        # Keep only the part before the first whitespace
        clean_id = record.id.split()[0]
        record.id = clean_id
        record.name = clean_id
        record.description = clean_id
        cleaned_records.append(record)
    
    SeqIO.write(cleaned_records, output_file, 'fasta')
    print(f"✓ Cleaned headers in {output_file}")


def rename_proteins(faa_file, ffn_file, gff_file):
    """Rename proteins exactly like SmORFinder does - with random prefix and sequential numbering."""
    out_faa, out_ffn, out_gff = [], [], []
    prefix = ''.join(choices(string.ascii_uppercase, k=6))
    num = 0
    
    # Read all files
    faa_records = list(SeqIO.parse(faa_file, 'fasta'))
    ffn_records = list(SeqIO.parse(ffn_file, 'fasta'))
    
    with open(gff_file, 'r') as f:
        gff_lines = [line.strip() for line in f if not line.startswith('#') and len(line.strip()) > 10]
    
    # Rename each protein
    for faa_rec, ffn_rec, gff_line in zip(faa_records, ffn_records, gff_lines):
        num += 1
        new_id = f"{prefix}_{num}"

        # Update FASTA record
        faa_rec.id = new_id
        faa_rec.name = new_id
        desc = faa_rec.description.split()
        if len(desc) > 0:
            desc[-1] = f'ID={new_id};' + ';'.join(desc[-1].split(';')[1:])
            faa_rec.description = ' '.join(desc)
        out_faa.append(faa_rec)

        # Update FFN record
        ffn_rec.id = new_id
        ffn_rec.name = new_id
        desc = ffn_rec.description.split()
        if len(desc) > 0:
            desc[-1] = f'ID={new_id};' + ';'.join(desc[-1].split(';')[1:])
            ffn_rec.description = ' '.join(desc)
        out_ffn.append(ffn_rec)

        # Update GFF line
        gff_parts = gff_line.split('\t')
        if len(gff_parts) >= 9:
            gff_parts[-1] = f'ID={new_id};' + ';'.join(gff_parts[-1].split(';')[1:])
            out_gff.append('\t'.join(gff_parts))

    # Write updated files
    SeqIO.write(out_faa, faa_file, 'fasta')
    SeqIO.write(out_ffn, ffn_file, 'fasta')
    
    with open(gff_file, 'w') as outfile:
        for line in out_gff:
            outfile.write(line + '\n')
    
    print(f"✓ Renamed proteins with prefix {prefix} ({num} proteins)")


def run_pyrodigal(genome_file, output_prefix, meta=False):
    """Run Pyrodigal on a genome file using the exact same command as SmORFinder."""
    cmd = [
        'pyrodigal',
        '--min-gene', '15',
        '--max-overlap', '15',
        '-c',  # Closed ends
        '-f', 'gff',
        '-o', f'{output_prefix}.gff', 
        '-a', f'{output_prefix}.faa', 
        '-d', f'{output_prefix}.ffn',
        '-i', genome_file
    ]
    
    # Add metagenomic mode if specified
    if meta:
        cmd.append('-p')  # Metagenomic mode
        print(f"Running in metagenomic mode...")
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running Pyrodigal: {result.stderr}")
        return False
    
    print(f"✓ Generated {output_prefix}.faa, {output_prefix}.ffn, and {output_prefix}.gff")
    return True


def clean_gff_ids(gff_file):
    """Clean GFF file to match cleaned FASTA headers."""
    cleaned_lines = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                cleaned_lines.append(line)
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                # Clean the attributes column (last column)
                attributes = parts[8]
                if 'ID=' in attributes:
                    # Extract the ID and clean it
                    for attr in attributes.split(';'):
                        if attr.startswith('ID='):
                            old_id = attr.split('=')[1]
                            clean_id = old_id.split()[0]  # Remove everything after whitespace
                            new_attr = f"ID={clean_id}"
                            # Replace the old ID attribute
                            attributes = attributes.replace(attr, new_attr)
                            break
                
                parts[8] = attributes
                cleaned_lines.append('\t'.join(parts) + '\n')
            else:
                cleaned_lines.append(line)
    
    # Write back to the same file
    with open(gff_file, 'w') as f:
        f.writelines(cleaned_lines)
    
    print(f"✓ Cleaned GFF IDs in {gff_file}")


@click.command()
@click.option('--genome-dir', default='test_genomes', help='Directory containing genome FASTA files')
@click.option('--output-dir', default='test_data', help='Output directory for processed files')

def main(genome_dir, output_dir):
    """Setup test data for SmORFinder testing."""
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Find all .fna and .fa files in genome directory
    genome_files = []
    genome_path = Path(genome_dir)
    
    # Search for both .fna and .fa files
    genome_files.extend(list(genome_path.glob('*.fna')))
    genome_files.extend(list(genome_path.glob('*.fa')))
    
    # Remove duplicates (in case of multiple extensions)
    genome_files = list(set(genome_files))
    
    if not genome_files:
        print(f"No .fna or .fa files found in {genome_dir}")
        print("Please place your test genome FASTA files (.fna or .fa) in the genome directory.")
        return
    
    print(f"Found {len(genome_files)} genome files")
    
    for genome_file in genome_files:
        genome_name = genome_file.stem
        output_prefix = os.path.join(output_dir, genome_name)
        
        print(f"\nProcessing {genome_name}...")
        
        # Copy and clean genome file to output directory with .fna extension
        genome_output = f"{output_prefix}.fna"
        if not os.path.exists(genome_output):
            # Clean the genome FASTA headers
            clean_fasta_headers(str(genome_file), genome_output)
        else:
            print(f"✓ Genome file already exists: {genome_output}")
        
        # Run Pyrodigal (try single genome mode first, fallback to metagenomic)
        if run_pyrodigal(genome_output, output_prefix, meta=False):
            # Rename proteins like SmORFinder does
            rename_proteins(
                f"{output_prefix}.faa",
                f"{output_prefix}.ffn", 
                f"{output_prefix}.gff"
            )
            
            print(f"✓ Successfully processed {genome_name}")
        else:
            print(f"✗ Failed to process {genome_name} in single genome mode")
            print(f"Trying metagenomic mode...")
            
            # Try metagenomic mode as fallback
            if run_pyrodigal(genome_output, output_prefix, meta=True):
                # Rename proteins like SmORFinder does
                rename_proteins(
                    f"{output_prefix}.faa",
                    f"{output_prefix}.ffn", 
                    f"{output_prefix}.gff"
                )
                
                print(f"✓ Successfully processed {genome_name} in metagenomic mode")
            else:
                print(f"✗ Failed to process {genome_name} in both modes")
    
    print(f"\nTest data setup complete!")
    print(f"Files are in: {output_dir}")
    print(f"\nYou can now run tests with:")
    print(f"smorf single {output_dir}/genome1.fna --outdir test_output_standard")
    print(f"smorf pre-called {output_dir}/genome1.fna {output_dir}/genome1.faa {output_dir}/genome1.ffn {output_dir}/genome1.gff --outdir test_output_pre_called")


if __name__ == '__main__':
    main() 