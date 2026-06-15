#!/usr/bin/env python3

from pathlib import Path

import typer

from emcSP.core import SampleConfig, analyze_from_tsv, install_reference

app = typer.Typer(
    name="emcSP",
    help="EMC SigProfiler tool for mutational signature analysis",
    add_completion=False,
)


@app.command()
def install(
    reference: str = typer.Option(
        "GRCh37", "--reference", "-r", help="Reference genome (GRCh37 or GRCh38)"
    ),
    rsync: bool = typer.Option(
        True, "--rsync/--no-rsync", help="Use rsync for download"
    ),
    bash: bool = typer.Option(True, "--bash/--no-bash", help="Use bash for extraction"),
):
    """Install reference genome for SigProfiler."""
    typer.echo(f"Installing reference genome: {reference}")
    install_reference(reference, rsync=rsync, bash=bash)
    typer.echo(f"Reference genome {reference} installed successfully.")


@app.command()
def analyze(
    input_tsv: Path = typer.Argument(
        ..., exists=True, help="Input TSV file with CHROM, POS, REF, ALT columns"
    ),
    sample_name: str = typer.Option(
        ..., "--sample-name", "-s", help="Sample name for the analysis"
    ),
    output_dir: Path = typer.Option(
        Path.cwd(), "--output", "-o", help="Output directory"
    ),
    reference: str = typer.Option(
        "GRCh37", "--reference", "-r", help="Reference genome (GRCh37 or GRCh38)"
    ),
    plot: bool = typer.Option(
        False, "--plot/--no-plot", help="Generate plots for mutational matrices"
    ),
    context_type: str = typer.Option(
        "96", "--context-type", "-c", help="Mutation context type (e.g., 96, 192, 1536)"
    ),
    cosmic_version: float = typer.Option(
        3.4, "--cosmic-version", help="COSMIC signature version"
    ),
    exome: bool = typer.Option(False, "--exome/--no-exome", help="Use exome filtering"),
    export_probabilities: bool = typer.Option(
        True,
        "--export-probabilities/--no-export-probabilities",
        help="Export probability matrices",
    ),
    nnls_add_penalty: float = typer.Option(
        0.05, "--nnls-add-penalty", help="NNLS additional penalty (default: 0.05)"
    ),
    nnls_remove_penalty: float = typer.Option(
        0.01, "--nnls-remove-penalty", help="NNLS removal penalty (default: 0.01)"
    ),
    initial_remove_penalty: float = typer.Option(
        0.05, "--initial-remove-penalty", help="Initial removal penalty (default: 0.05)"
    ),
    make_plots: bool = typer.Option(
        False, "--make-plots/--no-make-plots", help="Generate signature plots"
    ),
    output_format: str = typer.Option(
        "tsv", "--format", "-f", help="Output format: 'tsv', 'html', or 'both'"
    ),
):
    """Analyze a single sample from TSV file to signature assignment."""
    if output_format not in ("tsv", "html", "both"):
        typer.echo(
            f"Error: Invalid output format '{output_format}'. Choose from 'tsv', 'html', or 'both'.",
            err=True,
        )
        raise typer.Exit(code=1)

    effective_plot = plot or (output_format in ("html", "both"))

    cfg = SampleConfig(
        sample_name=sample_name,
        output_dir=output_dir,
        reference=reference,
        plot=effective_plot,
        context_type=context_type,
        cosmic_version=cosmic_version,
        exome=exome,
        export_probabilities=export_probabilities,
        nnls_add_penalty=nnls_add_penalty,
        nnls_remove_penalty=nnls_remove_penalty,
        initial_remove_penalty=initial_remove_penalty,
        make_plots=make_plots,
        output_format=output_format,
    )
    analyze_from_tsv(input_tsv=input_tsv, cfg=cfg, logger=typer.echo)


if __name__ == "__main__":
    app()
