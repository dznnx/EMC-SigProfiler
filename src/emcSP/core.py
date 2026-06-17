import base64
import csv
import importlib.resources
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

import jinja2
from SigProfilerAssignment.Analyzer import cosmic_fit as cosmic_fit_func
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

from emcSP.signatures import SIGNATURE_DESCRIPTIONS


@dataclass
class SampleConfig:
    sample_name: str
    output_dir: Path
    vcf_path: Path | None = None
    reference: str = "GRCh37"
    context_type: str = "96"
    cosmic_version: float = 3.5
    exome: bool = False
    plot: bool = False
    make_plots: bool = False
    export_probabilities: bool = True
    nnls_add_penalty: float = 0.05
    nnls_remove_penalty: float = 0.01
    initial_remove_penalty: float = 0.05
    output_format: str = "tsv"


def analyze_sample(cfg: SampleConfig, logger=None):
    """Analyze a single sample from a VCF file to signature assignment"""

    def log(msg):
        (logger or print)(msg)

    vcf_dir = cfg.output_dir / "vcf_tmp"
    vcf_dir.mkdir(parents=True, exist_ok=True)

    vcf_path_in_dir = vcf_dir / f"{cfg.sample_name}.vcf"
    if cfg.vcf_path is None:
        return
    if cfg.vcf_path.resolve() != vcf_path_in_dir.resolve():
        shutil.copy(cfg.vcf_path, vcf_path_in_dir)

    log("Running SigProfilerMatrixGenerator...")
    matrix_output_dir = cfg.output_dir / "matrix"

    for attempt in ("first", "retry"):
        try:
            matGen.SigProfilerMatrixGeneratorFunc(
                project=cfg.sample_name,
                reference_genome=cfg.reference,
                path_to_input_files=str(vcf_dir),
                plot=cfg.plot,
                exome=cfg.exome,
                output_directory=str(matrix_output_dir),
            )
            break
        except Exception:
            if attempt == "first":
                log(
                    f"Reference genome {cfg.reference} not installed. Installing now..."
                )
                genInstall.install(cfg.reference, rsync=True, bash=True)
                log("Installation complete. Retrying SigProfilerMatrixGenerator...")
            else:
                raise

    matrix_input_dir = matrix_output_dir / "input"
    if matrix_input_dir.exists():
        shutil.rmtree(matrix_input_dir)

    final_vcf_path = cfg.output_dir / f"{cfg.sample_name}.vcf"
    if (
        vcf_path_in_dir.exists()
        and vcf_path_in_dir.resolve() != final_vcf_path.resolve()
    ):
        vcf_path_in_dir.rename(final_vcf_path)
    if vcf_dir.exists():
        shutil.rmtree(vcf_dir)

    log("Running SigProfilerAssignment.cosmic_fit...")
    output_path = cfg.output_dir / "output"
    matrix_file = (
        matrix_output_dir / "SBS" / f"{cfg.sample_name}.SBS{cfg.context_type}.all"
    )

    make_plots = cfg.make_plots
    sample_reconstruction_plots = "none"
    if cfg.output_format in ("html", "both"):
        make_plots = True
        if shutil.which("pdftoppm"):
            sample_reconstruction_plots = "both"
        else:
            log(
                "Warning: pdftoppm (poppler-utils) not found. Sample reconstruction plots will not be generated."
            )

    cosmic_fit_func(
        samples=str(matrix_file),
        output=str(output_path),
        input_type="matrix",
        context_type=cfg.context_type,
        genome_build=cfg.reference,
        cosmic_version=cfg.cosmic_version,
        exome=cfg.exome,
        export_probabilities=cfg.export_probabilities,
        nnls_add_penalty=cfg.nnls_add_penalty,
        nnls_remove_penalty=cfg.nnls_remove_penalty,
        initial_remove_penalty=cfg.initial_remove_penalty,
        make_plots=make_plots,
        sample_reconstruction_plots=sample_reconstruction_plots,
    )

    log(f"Analysis complete. Results in: {output_path}")
    log("Formatting signature outputs...")

    activities_file = (
        output_path
        / "Assignment_Solution"
        / "Activities"
        / "Assignment_Solution_Activities.txt"
    )
    if not activities_file.exists():
        return

    signatures = _parse_signatures(activities_file, cfg.sample_name)
    if not signatures:
        return

    headers = ["Signature", "Relative Contribution", "Mutation Count", "Etiology"]

    if cfg.output_format in ("tsv", "both"):
        _write_tsv(cfg, signatures, headers, log)

    if cfg.output_format in ("html", "both"):
        _write_html(cfg, output_path, signatures, log)


def _parse_signatures(activities_file, sample_name):
    """Parse signature activities for sample, returns list of (sig, desc, prob, count)"""
    with open(activities_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("Samples") != sample_name:
                continue

            total_counts = sum(
                float(v)
                for k, v in row.items()
                if k != "Samples" and v.replace(".", "", 1).isdigit()
            )
            if total_counts <= 0:
                return []

            signatures = []
            for k, v in row.items():
                if k == "Samples" or not v.replace(".", "", 1).isdigit():
                    continue
                count = float(v)
                if count > 0:
                    prob = count / total_counts
                    desc = SIGNATURE_DESCRIPTIONS.get(k, "Unknown")
                    signatures.append((k, desc, prob, count))

            signatures.sort(key=lambda x: x[3], reverse=True)
            return signatures

    return []


def _write_tsv(cfg, signatures, headers, log):
    tsv_out = cfg.output_dir / f"{cfg.sample_name}_signatures.tsv"
    with open(tsv_out, "w") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(headers)
        for sig, desc, prob, count in signatures:
            display_count = int(count) if count.is_integer() else count
            writer.writerow([sig, f"{prob:.3f}", display_count, desc])
    log(f"Signatures TSV created: {tsv_out}")


def _write_html(cfg, output_path, signatures, log):
    plot_b64 = None
    plot_path = _find_sbs_plot(cfg)
    if plot_path:
        plot_b64 = _pdf_to_base64_png(plot_path)

    stats = _parse_stats(output_path, cfg.sample_name, log)
    annotated_stats = _annotate_stats(stats)

    decomp_plot_b64 = _load_decomp_plot(output_path, cfg.sample_name, log)

    try:
        template_text = (
            importlib.resources.files("emcSP")
            .joinpath("templates/report.html")
            .read_text(encoding="utf-8")
        )
        template = jinja2.Template(template_text)
        html_content = template.render(
            sample_name=cfg.sample_name,
            signatures=[
                {
                    "name": sig,
                    "contribution": prob,
                    "count": int(count) if count.is_integer() else count,
                    "etiology": desc,
                }
                for sig, desc, prob, count in signatures
            ],
            plot_b64=plot_b64,
            stats=annotated_stats,
            decomp_plot_b64=decomp_plot_b64,
        )
    except Exception as e:
        log(f"Error rendering HTML template: {e}")
        html_content = (
            f"<html><body><h1>Mutational Signature Analysis - {cfg.sample_name}</h1>"
            f"<p>Error rendering report: {e}</p></body></html>"
        )

    html_out = cfg.output_dir / f"{cfg.sample_name}_signatures.html"
    with open(html_out, "w", encoding="utf-8") as out_f:
        out_f.write(html_content)
    log(f"Signatures HTML created: {html_out}")


def _parse_stats(output_path, sample_name, log):
    stats_file = (
        output_path
        / "Assignment_Solution"
        / "Solution_Stats"
        / "Assignment_Solution_Samples_Stats.txt"
    )
    if not stats_file.exists():
        return None

    try:
        with open(stats_file) as f:
            for s_row in csv.DictReader(f, delimiter="\t"):
                if s_row.get("Sample Names") == sample_name:
                    return {
                        "cosine_similarity": s_row.get("Cosine Similarity", "0.0"),
                        "correlation": s_row.get("Correlation", "0.0"),
                        "kl_divergence": s_row.get("KL Divergence", "0.0"),
                        "total_mutations": s_row.get("Total Mutations", "0"),
                    }
    except Exception as e:
        log(f"Error parsing stats file: {e}")
    return None


def _annotate_stats(stats_dict):
    """Add status field for each metric based on thresholds"""

    cs = float(stats_dict.get("cosine_similarity", 0))
    corr = float(stats_dict.get("correlation", 0))
    kl = float(stats_dict.get("kl_divergence", 0))
    mut = int(float(stats_dict.get("total_mutations", 0)))

    # cosine_similarity: ≥0.9 success, ≥0.8 warning, else error
    if cs >= 0.9:
        stats_dict["cosine_similarity_status"] = "success"
    elif cs >= 0.8:
        stats_dict["cosine_similarity_status"] = "warning"
    else:
        stats_dict["cosine_similarity_status"] = "error"

    # correlation: ≥0.9 success, ≥0.8 warning, else error
    if corr >= 0.9:
        stats_dict["correlation_status"] = "success"
    elif corr >= 0.8:
        stats_dict["correlation_status"] = "warning"
    else:
        stats_dict["correlation_status"] = "error"

    # kl_divergence: ≤0.1 success, ≤0.2 warning, else error
    if kl <= 0.1:
        stats_dict["kl_divergence_status"] = "success"
    elif kl <= 0.2:
        stats_dict["kl_divergence_status"] = "warning"
    else:
        stats_dict["kl_divergence_status"] = "error"

    # total_mutations: <50 low, 50-200 moderate, 200-500 high, >500 very_high
    if mut < 50:
        stats_dict["total_mutations_status"] = "error"
    elif mut < 200:
        stats_dict["total_mutations_status"] = "warning"
    elif mut < 500:
        stats_dict["total_mutations_status"] = "success"
    else:
        stats_dict["total_mutations_status"] = "info"

    return stats_dict


def _load_decomp_plot(output_path, sample_name, log):
    decomp_path = (
        output_path
        / "Assignment_Solution"
        / "Activities"
        / "SampleReconstruction"
        / "WebPNGs"
        / f"{sample_name}.png"
    )
    if not decomp_path.exists():
        return None

    try:
        with open(decomp_path, "rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")
    except Exception as e:
        log(f"Error reading decomposition plot: {e}")
        return None


def install_reference(reference: str, rsync: bool = True, bash: bool = True):
    """Install reference genome for SigProfiler."""
    genInstall.install(reference, rsync=rsync, bash=bash)


def analyze_from_tsv(input_tsv: Path, cfg: SampleConfig, logger=None):
    """Analyze a single sample from TSV file to signature assignment."""

    def log(msg):
        if logger:
            logger(msg)
        else:
            print(msg)

    log("Converting TSV to VCF...")
    vcf_path = _tsv_to_vcf(input_tsv, cfg.sample_name)
    log(f"VCF created: {vcf_path}")

    cfg = replace(cfg, vcf_path=vcf_path)
    analyze_sample(cfg, logger=logger)

    if (
        vcf_path.exists()
        and vcf_path.resolve() != (cfg.output_dir / f"{cfg.sample_name}.vcf").resolve()
    ):
        vcf_path.unlink()


def _tsv_to_vcf(tsv_path: Path, sample_name: str) -> Path:
    vcf_lines = [
        "##fileformat=VCFv4.2",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}",
    ]

    with open(tsv_path) as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
            vcf_lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t.\t.")

    vcf_path = tsv_path.with_suffix(".vcf")
    with open(vcf_path, "w") as f:
        f.write("\n".join(vcf_lines))
    return vcf_path


def _find_sbs_plot(cfg: SampleConfig) -> Path | None:
    plot_dir = cfg.output_dir / "matrix" / "plots"
    specific_plot = plot_dir / f"SBS_{cfg.context_type}_plots_{cfg.sample_name}.pdf"
    if specific_plot.exists():
        return specific_plot

    fallback_plot = plot_dir / f"SBS_96_plots_{cfg.sample_name}.pdf"
    if fallback_plot.exists():
        return fallback_plot

    return None


def _pdf_to_base64_png(pdf_path: Path) -> str | None:
    if not pdf_path.exists() or not shutil.which("pdftoppm"):
        return None
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "page"
            subprocess.run(
                ["pdftoppm", "-png", "-r", "150", str(pdf_path), str(prefix)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            png_files = sorted(Path(tmpdir).glob("page-*.png"))
            if png_files:
                with open(png_files[0], "rb") as f:
                    return base64.b64encode(f.read()).decode("utf-8")
    except Exception:
        pass
    return None
