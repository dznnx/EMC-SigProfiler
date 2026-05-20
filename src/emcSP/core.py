import csv
import shutil
from pathlib import Path

from SigProfilerAssignment.Analyzer import cosmic_fit as cosmic_fit_func
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

from emcSP.signatures import SIGNATURE_DESCRIPTIONS


def tsv_to_vcf(tsv_path: Path, sample_name: str) -> Path:
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


def analyze_sample(
    vcf_path: Path,
    sample_name: str,
    output_dir: Path,
    reference: str = "GRCh37",
    plot: bool = False,
    context_type: str = "96",
    cosmic_version: float = 3.4,
    exome: bool = False,
    export_probabilities: bool = True,
    nnls_add_penalty: float = 0.05,
    nnls_remove_penalty: float = 0.01,
    initial_remove_penalty: float = 0.05,
    make_plots: bool = False,
    logger=None,
):
    """Analyze a single sample from a VCF file to signature assignment."""

    def log(msg):
        if logger:
            logger(msg)
        else:
            print(msg)

    vcf_dir = output_dir / "vcf_tmp"
    vcf_dir.mkdir(parents=True, exist_ok=True)

    vcf_path_in_dir = vcf_dir / f"{sample_name}.vcf"
    if vcf_path.resolve() != vcf_path_in_dir.resolve():
        shutil.copy(vcf_path, vcf_path_in_dir)

    log("Running SigProfilerMatrixGenerator...")
    matrix_output_dir = output_dir / "matrix"
    try:
        log("matgen")
        matGen.SigProfilerMatrixGeneratorFunc(
            project=sample_name,
            reference_genome=reference,
            path_to_input_files=str(vcf_dir),
            plot=plot,
            exome=exome,
            output_directory=str(matrix_output_dir),
        )
    except Exception as e:
        log(f"Reference genome {reference} not installed. Installing now...")
        genInstall.install(reference, rsync=True, bash=True)
        log("Installation complete. Retrying SigProfilerMatrixGenerator...")
        matGen.SigProfilerMatrixGeneratorFunc(
            project=sample_name,
            reference_genome=reference,
            path_to_input_files=str(vcf_dir),
            plot=plot,
            exome=exome,
            output_directory=str(matrix_output_dir),
        )

    matrix_input_dir = matrix_output_dir / "input"
    if matrix_input_dir.exists():
        shutil.rmtree(matrix_input_dir)

    final_vcf_path = output_dir / f"{sample_name}.vcf"
    if (
        vcf_path_in_dir.exists()
        and vcf_path_in_dir.resolve() != final_vcf_path.resolve()
    ):
        vcf_path_in_dir.rename(final_vcf_path)
    if vcf_dir.exists():
        shutil.rmtree(vcf_dir)

    log("Running SigProfilerAssignment.cosmic_fit...")
    output_path = output_dir / "output"

    matrix_file = matrix_output_dir / "SBS" / f"{sample_name}.SBS{context_type}.all"

    cosmic_fit_func(
        samples=str(matrix_file),
        output=str(output_path),
        input_type="matrix",
        context_type=context_type,
        genome_build=reference,
        cosmic_version=cosmic_version,
        exome=exome,
        export_probabilities=export_probabilities,
        nnls_add_penalty=nnls_add_penalty,
        nnls_remove_penalty=nnls_remove_penalty,
        initial_remove_penalty=initial_remove_penalty,
        make_plots=make_plots,
    )

    log(f"Analysis complete. Results in: {output_path}")

    log("Formatting signature outputs...")
    activities_file = (
        output_path
        / "Assignment_Solution"
        / "Activities"
        / "Assignment_Solution_Activities.txt"
    )
    if activities_file.exists():
        with open(activities_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("Samples") == sample_name:
                    total_counts = sum(
                        float(v)
                        for k, v in row.items()
                        if k != "Samples" and v.replace(".", "", 1).isdigit()
                    )
                    if total_counts > 0:
                        signatures = []
                        for k, v in row.items():
                            if k != "Samples" and v.replace(".", "", 1).isdigit():
                                count = float(v)
                                if count > 0:
                                    prob = count / total_counts
                                    desc = SIGNATURE_DESCRIPTIONS.get(k, "Unknown")
                                    signatures.append((k, desc, prob, count))

                        signatures.sort(key=lambda x: x[3], reverse=True)

                        tsv_out = output_dir / f"{sample_name}_signatures.tsv"
                        with open(tsv_out, "w") as out_f:
                            writer = csv.writer(out_f, delimiter="\t")
                            writer.writerow(
                                ["Signature", "Probability", "Count", "Description"]
                            )
                            for sig, desc, prob, count in signatures:
                                display_count = (
                                    int(count) if count.is_integer() else count
                                )
                                writer.writerow(
                                    [sig, f"{prob:.3f}", display_count, desc]
                                )

                        log(f"Signatures TSV created: {tsv_out}")
                    break


def install_reference(reference: str, rsync: bool = True, bash: bool = True):
    """Install reference genome for SigProfiler."""
    genInstall.install(reference, rsync=rsync, bash=bash)


def analyze_from_tsv(
    input_tsv: Path,
    sample_name: str,
    output_dir: Path,
    reference: str = "GRCh37",
    plot: bool = False,
    context_type: str = "96",
    cosmic_version: float = 3.4,
    exome: bool = False,
    export_probabilities: bool = True,
    nnls_add_penalty: float = 0.05,
    nnls_remove_penalty: float = 0.01,
    initial_remove_penalty: float = 0.05,
    make_plots: bool = False,
    logger=None,
):
    """Analyze a single sample from TSV file to signature assignment."""

    def log(msg):
        if logger:
            logger(msg)
        else:
            print(msg)

    log("Converting TSV to VCF...")
    vcf_path = tsv_to_vcf(input_tsv, sample_name)
    log(f"VCF created: {vcf_path}")

    analyze_sample(
        vcf_path=vcf_path,
        sample_name=sample_name,
        output_dir=output_dir,
        reference=reference,
        plot=plot,
        context_type=context_type,
        cosmic_version=cosmic_version,
        exome=exome,
        export_probabilities=export_probabilities,
        nnls_add_penalty=nnls_add_penalty,
        nnls_remove_penalty=nnls_remove_penalty,
        initial_remove_penalty=initial_remove_penalty,
        make_plots=make_plots,
        logger=logger,
    )

    if (
        vcf_path.exists()
        and vcf_path.resolve() != (output_dir / f"{sample_name}.vcf").resolve()
    ):
        vcf_path.unlink()
