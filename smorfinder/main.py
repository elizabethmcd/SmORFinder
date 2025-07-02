import click
from smorfinder import *
from smorfinder.help import CustomHelp
from smorfinder.run import _run, _run_protein

@click.group(cls=CustomHelp)
def cli():
    """SmORFinder: Predict and annotate small protein sequences in genomic data
    
    SmORFinder uses deep learning models and HMM searches to identify small ORFs (≤50aa)
    from genomic sequences. Two main workflows are available:
    
    GENOME WORKFLOWS (use Prodigal for gene prediction):
    • single: Complete/draft genome assembly of a single species
    • meta: Metagenomic assembly (multi-threaded)
    
    PROTEIN WORKFLOW (use your own gene predictions):
    • protein: Pre-predicted proteins with custom locus tags
    """
    pass


@cli.command(short_help='Run SmORFinder on a complete or draft genome assembly (single species)', help_priority=1)
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--outdir', '-o', default='smorf_output', help='Output directory for results')
@click.option('--prodigal-path', '-pp', default=PRODIGAL_PATH, type=click.Path(exists=True), help='Path to Prodigal executable')
@click.option('--dsn1-model-path', '-d1p', default=DSN1_MODEL_PATH, type=click.Path(exists=True), help='Path to DSN1 model file')
@click.option('--dsn2-model-path', '-d2p', default=DSN2_MODEL_PATH, type=click.Path(exists=True), help='Path to DSN2 model file')
@click.option('--smorf-hmm-path', '-shp', default=SMORFHMM_PATH, type=click.Path(exists=True), help='Path to SmORF HMM database')
@click.option('--hmmsearch-path', '-hp', default=HMMSEARCH_PATH, type=click.Path(exists=True), help='Path to HMMER hmmsearch executable')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory")
@click.option('--dsn1-indiv-cutoff', '-idsn1', default=0.9999, help='DSN1 individual significance cutoff (0-1, default=0.9999)')
@click.option('--dsn2-indiv-cutoff', '-idsn2', default=0.9999, help='DSN2 individual significance cutoff (0-1, default=0.9999)')
@click.option('--phmm-indiv-cutoff', '-iphmm', default=1e-6, help='pHMM individual significance cutoff (default=1e-6)')
@click.option('--dsn1-overlap-cutoff', '-odsn1', default=0.5, help='DSN1 overlap significance cutoff (0-1, default=0.5)')
@click.option('--dsn2-overlap-cutoff', '-odsn2', default=0.5, help='DSN2 overlap significance cutoff (0-1, default=0.5)')
@click.option('--phmm-overlap-cutoff', '-ophmm', default=1, help='pHMM overlap significance cutoff (default=1)')
def single(fasta, outdir, prodigal_path, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff):
    """Run SmORFinder on a complete or draft genome assembly of a single species.
    
    This workflow uses Prodigal to predict genes from the genome, then filters for
    small proteins (≤50aa) and applies SmORFinder's deep learning models.
    
    FASTA: Genome assembly file (.fa, .fna, .fasta)
    
    Output files:
    • {outdir}.faa: Predicted small protein sequences
    • {outdir}.ffn: Predicted small nucleotide sequences  
    • {outdir}.gff: Gene predictions with SmORFinder annotations
    • {outdir}.tsv: Detailed results table
    """
    log_params(command='single', fasta=fasta, outdir=outdir, prodigal_path=prodigal_path, dsn1_model_path=dsn1_model_path,
               dsn2_model_path=dsn2_model_path, smorf_hmm_path=smorf_hmm_path, hmmsearch_path=hmmsearch_path, 
               force=force, dsn1_indiv_cutoff=dsn1_indiv_cutoff, dsn2_indiv_cutoff=dsn2_indiv_cutoff, 
               phmm_indiv_cutoff=phmm_indiv_cutoff, dsn1_overlap_cutoff=dsn1_overlap_cutoff, 
               dsn2_overlap_cutoff=dsn2_overlap_cutoff, phmm_overlap_cutoff=phmm_overlap_cutoff)

    _run(fasta, outdir, 1, prodigal_path, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, 
         force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, 
         dsn2_overlap_cutoff, phmm_overlap_cutoff, mode='single')


@cli.command(short_help='Run SmORFinder on a metagenomic assembly (multi-threaded)', help_priority=2)
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--outdir', '-o', default='smorf_output', help='Output directory for results')
@click.option('--threads', '-t', default=1, help='Number of threads for Prodigal')
@click.option('--prodigal-path', '-pp', default=PRODIGAL_PATH, type=click.Path(exists=True), help='Path to Prodigal executable')
@click.option('--dsn1-model-path', '-d1p', default=DSN1_MODEL_PATH, type=click.Path(exists=True), help='Path to DSN1 model file')
@click.option('--dsn2-model-path', '-d2p', default=DSN2_MODEL_PATH, type=click.Path(exists=True), help='Path to DSN2 model file')
@click.option('--smorf-hmm-path', '-shp', default=SMORFHMM_PATH, type=click.Path(exists=True), help='Path to SmORF HMM database')
@click.option('--hmmsearch-path', '-hp', default=HMMSEARCH_PATH, type=click.Path(exists=True), help='Path to HMMER hmmsearch executable')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory")
@click.option('--dsn1-indiv-cutoff', '-idsn1', default=0.9999, help='DSN1 individual significance cutoff (0-1, default=0.9999)')
@click.option('--dsn2-indiv-cutoff', '-idsn2', default=0.9999, help='DSN2 individual significance cutoff (0-1, default=0.9999)')
@click.option('--phmm-indiv-cutoff', '-iphmm', default=1e-6, help='pHMM individual significance cutoff (default=1e-6)')
@click.option('--dsn1-overlap-cutoff', '-odsn1', default=0.5, help='DSN1 overlap significance cutoff (0-1, default=0.5)')
@click.option('--dsn2-overlap-cutoff', '-odsn2', default=0.5, help='DSN2 overlap significance cutoff (0-1, default=0.5)')
@click.option('--phmm-overlap-cutoff', '-ophmm', default=1, help='pHMM overlap significance cutoff (default=1)')
def meta(fasta, outdir, threads, prodigal_path, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff):
    """Run SmORFinder on a metagenomic assembly using multi-threaded Prodigal.
    
    This workflow is optimized for metagenomic assemblies and uses Prodigal in
    metagenomic mode with multi-threading for faster gene prediction.
    
    FASTA: Metagenomic assembly file (.fa, .fna, .fasta)
    
    Output files:
    • {outdir}.faa: Predicted small protein sequences
    • {outdir}.ffn: Predicted small nucleotide sequences  
    • {outdir}.gff: Gene predictions with SmORFinder annotations
    • {outdir}.tsv: Detailed results table
    """
    log_params(command='meta', fasta=fasta, outdir=outdir, threads=threads, prodigal_path=prodigal_path,
               dsn1_model_path=dsn1_model_path, dsn2_model_path=dsn2_model_path, smorf_hmm_path=smorf_hmm_path, 
               hmmsearch_path=hmmsearch_path, force=force, dsn1_indiv_cutoff=dsn1_indiv_cutoff, 
               dsn2_indiv_cutoff=dsn2_indiv_cutoff, phmm_indiv_cutoff=phmm_indiv_cutoff, 
               dsn1_overlap_cutoff=dsn1_overlap_cutoff, dsn2_overlap_cutoff=dsn2_overlap_cutoff, 
               phmm_overlap_cutoff=phmm_overlap_cutoff)

    _run(fasta, outdir, threads, prodigal_path, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, 
         force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, 
         dsn2_overlap_cutoff, phmm_overlap_cutoff, mode='meta')


@cli.command(short_help='Run SmORFinder on pre-predicted proteins (≤50aa)', help_priority=3)
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('protein_faa', type=click.Path(exists=True))
@click.argument('gff', type=click.Path(exists=True))
@click.option('--outdir', '-o', default='smorf_output', help='Output directory for results')
@click.option('--dsn1-model-path', '-d1p', default=DSN1_MODEL_PATH, type=click.Path(exists=True), help='Path to DSN1 model file')
@click.option('--dsn2-model-path', '-d2p', default=DSN2_MODEL_PATH, type=click.Path(exists=True), help='Path to DSN2 model file')
@click.option('--smorf-hmm-path', '-shp', default=SMORFHMM_PATH, type=click.Path(exists=True), help='Path to SmORF HMM database')
@click.option('--hmmsearch-path', '-hp', default=HMMSEARCH_PATH, type=click.Path(exists=True), help='Path to HMMER hmmsearch executable')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory")
@click.option('--dsn1-indiv-cutoff', '-idsn1', default=0.9999, help='DSN1 individual significance cutoff (0-1, default=0.9999)')
@click.option('--dsn2-indiv-cutoff', '-idsn2', default=0.9999, help='DSN2 individual significance cutoff (0-1, default=0.9999)')
@click.option('--phmm-indiv-cutoff', '-iphmm', default=1e-6, help='pHMM individual significance cutoff (default=1e-6)')
@click.option('--dsn1-overlap-cutoff', '-odsn1', default=0.5, help='DSN1 overlap significance cutoff (0-1, default=0.5)')
@click.option('--dsn2-overlap-cutoff', '-odsn2', default=0.5, help='DSN2 overlap significance cutoff (0-1, default=0.5)')
@click.option('--phmm-overlap-cutoff', '-ophmm', default=1, help='pHMM overlap significance cutoff (default=1)')
def protein(fasta, protein_faa, gff, outdir, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff):
    """Run SmORFinder on pre-predicted proteins with custom locus tags.
    
    This workflow is for when you already have protein predictions and want to use
    your own locus tag formatting instead of Prodigal's default naming.
    
    Required files:
    • FASTA: Reference genome file (.fa, .fna, .fasta) for context extraction
    • PROTEIN_FAA: Protein sequences in FASTA format (≤50aa, headers must match GFF IDs)
    • GFF: GFF file with gene predictions (must have ID field in attributes)
    
    Note: Protein headers in the FASTA file must match the ID field in the GFF file
    for proper sequence extraction and annotation.
    
    Output files:
    • {outdir}.faa: Predicted small protein sequences
    • {outdir}.gff: Gene predictions with SmORFinder annotations
    • {outdir}.tsv: Detailed results table
    
    Alternative: Use 'single' or 'meta' commands with just a genome FASTA file
    to let SmORFinder run Prodigal for gene prediction automatically.
    """
    log_params(command='protein', fasta=fasta, protein_faa=protein_faa, gff=gff, 
               outdir=outdir, dsn1_model_path=dsn1_model_path, dsn2_model_path=dsn2_model_path,
               smorf_hmm_path=smorf_hmm_path, hmmsearch_path=hmmsearch_path, force=force, 
               dsn1_indiv_cutoff=dsn1_indiv_cutoff, dsn2_indiv_cutoff=dsn2_indiv_cutoff, 
               phmm_indiv_cutoff=phmm_indiv_cutoff, dsn1_overlap_cutoff=dsn1_overlap_cutoff, 
               dsn2_overlap_cutoff=dsn2_overlap_cutoff, phmm_overlap_cutoff=phmm_overlap_cutoff)

    _run_protein(fasta, protein_faa, gff, outdir, dsn1_model_path, dsn2_model_path, 
                smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, 
                phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff)


def log_params(**kwargs):
    click.echo("#### PARAMETERS ####")
    click.echo('\n'.join(list(map(lambda x: ': '.join(list(map(str, x))), kwargs.items()))))
    click.echo("####################")

if __name__ == '__main__':
    cli()
