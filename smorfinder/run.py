import sys
import click
from os import makedirs
from os.path import join, isdir
import shutil
from Bio import SeqIO
from smorfinder.prodigal import *
from smorfinder.hmmsearch import run_hmmsearch
from smorfinder.model import *
from smorfinder.finalize import _finalize, _finalize_pre_called


def _run(fasta, outdir, threads, prodigal_path, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff, mode):
    """Run SmORFinder on genome sequences using Prodigal for gene prediction."""
    tmp_dir = join(outdir, 'tmp')
    if force and isdir(outdir):
        shutil.rmtree(outdir)
    try:
        makedirs(tmp_dir)
    except FileExistsError:
        click.echo("Output directory exists, please delete or overwrite with --force")
        sys.exit(1)

    click.echo("Running Prodigal...")
    if mode == 'meta':
        run_prodigal_multithread(prodigal_path, fasta, tmp_dir, threads)
    else:
        run_prodigal(prodigal_path, fasta, tmp_dir, meta=False)

    rename_proteins(join(tmp_dir, 'prodigal.faa'), join(tmp_dir, 'prodigal.ffn'), join(tmp_dir, 'prodigal.gff'))
    click.echo("Filtering to only predicted genes less than or equal to 50 aa in length...")
    filter_prodigal_small_genes(tmp_dir)
    click.echo("Running HMMSEARCH...")
    run_hmmsearch(hmmsearch_path, smorf_hmm_path, join(tmp_dir, 'prodigal.small.faa'), tmp_dir)
    click.echo("Extracting nucleotide sequences...")
    names, fiveprime, orf, threeprime = extract_sequences(fasta, join(tmp_dir, 'prodigal.small.gff'))
    click.echo("Running deep learning models on predicted smORFs...")
    dsn1_predictions = run_model(fiveprime, orf, threeprime, dsn1_model_path)
    dsn2_predictions = run_model(fiveprime, orf, threeprime, dsn2_model_path)

    write_results_to_file(dsn1_predictions, dsn2_predictions, names, fiveprime, orf, threeprime, join(tmp_dir, 'model_predictions.tsv'))
    click.echo("Finalizing results...")
    _finalize(outdir, tmp_dir, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff)


def _run_pre_called(fasta, protein_faa, nucleotide_ffn, gff, outdir, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff):
    """Run SmORFinder on pre-called genes from Prodigal output (≤50aa, custom locus tags), requiring all Prodigal output files for context extraction."""
    tmp_dir = join(outdir, 'tmp')
    if force and isdir(outdir):
        shutil.rmtree(outdir)
    makedirs(tmp_dir, exist_ok=True)

    # Copy all Prodigal output files to tmp directory
    shutil.copy(protein_faa, join(tmp_dir, 'input.faa'))
    shutil.copy(nucleotide_ffn, join(tmp_dir, 'input.ffn'))
    shutil.copy(gff, join(tmp_dir, 'input.gff'))

    click.echo("Filtering to only predicted genes less than or equal to 50 aa in length...")
    filter_pre_called_genes(tmp_dir)

    click.echo("Running HMMSEARCH...")
    run_hmmsearch(hmmsearch_path, smorf_hmm_path, join(tmp_dir, 'input.small.faa'), tmp_dir)

    click.echo("Extracting nucleotide sequences for deep learning...")
    names, fiveprime, orf, threeprime = extract_sequences(fasta, join(tmp_dir, 'input.small.gff'))
    click.echo("Running deep learning models on predicted smORFs...")
    dsn1_predictions = run_model(fiveprime, orf, threeprime, dsn1_model_path)
    dsn2_predictions = run_model(fiveprime, orf, threeprime, dsn2_model_path)
    write_results_to_file(dsn1_predictions, dsn2_predictions, names, fiveprime, orf, threeprime, join(tmp_dir, 'model_predictions.tsv'))
    click.echo("Finalizing results...")
    _finalize_pre_called(outdir, tmp_dir, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff)


def filter_pre_called_genes(tmp_dir):
    """Filter pre-called gene predictions to only small genes (≤51 aa, ≤153 nt)."""
    # Filter protein sequences
    faa_recs = []
    for rec in SeqIO.parse(join(tmp_dir, 'input.faa'), 'fasta'):
        if len(rec.seq) <= 51:
            faa_recs.append(rec)
    SeqIO.write(faa_recs, join(tmp_dir, 'input.small.faa'), 'fasta')

    # Filter nucleotide sequences
    ffn_recs = []
    for rec in SeqIO.parse(join(tmp_dir, 'input.ffn'), 'fasta'):
        if len(rec.seq) <= 153:
            ffn_recs.append(rec)
    SeqIO.write(ffn_recs, join(tmp_dir, 'input.small.ffn'), 'fasta')

    # Filter GFF file
    outgff = open(join(tmp_dir, 'input.small.gff'), 'w')
    with open(join(tmp_dir, 'input.gff')) as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')

            if line[2] != 'CDS':
                continue

            if int(line[4]) - int(line[3]) + 1 <= 153:
                print('\t'.join(line), file=outgff)
    outgff.close()

