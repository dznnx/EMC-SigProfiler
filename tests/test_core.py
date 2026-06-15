from pathlib import Path

import pytest
from emcSP import core


@pytest.fixture
def sample_tsv(tmp_path: Path) -> Path:
    tsv_content = "CHROM\tPOS\tREF\tALT\nchr1\t123456\tA\tG\nchr2\t789012\tC\tT\nchrX\t456789\tG\tA\n"
    tsv_path = tmp_path / "sample.tsv"
    tsv_path.write_text(tsv_content)
    return tsv_path


class TestTsvToVcf:
    def test_tsv_to_vcf_creates_correct_vcf(self, sample_tsv: Path, tmp_path: Path):
        result = core._tsv_to_vcf(sample_tsv, "test_sample")

        assert result.exists()
        assert result.suffix == ".vcf"

        lines = result.read_text().split("\n")
        assert lines[0] == "##fileformat=VCFv4.2"
        assert (
            lines[1]
            == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest_sample"
        )
        assert lines[2] == "chr1\t123456\t.\tA\tG\t.\t.\t.\t.\t."
        assert lines[3] == "chr2\t789012\t.\tC\tT\t.\t.\t.\t.\t."
        assert lines[4] == "chrX\t456789\t.\tG\tA\t.\t.\t.\t.\t."

    def test_tsv_to_vcf_returns_vcf_path(self, sample_tsv: Path):
        vcf_path = core._tsv_to_vcf(sample_tsv, "sample1")
        expected = sample_tsv.with_suffix(".vcf")
        assert vcf_path == expected

    def test_tsv_to_vcf_skips_header(self, sample_tsv: Path, tmp_path: Path):
        vcf_path = core._tsv_to_vcf(sample_tsv, "sample1")
        lines = vcf_path.read_text().split("\n")
        assert "CHROM" not in lines[-1]


class TestAnalyzeFromTsv:
    def test_analyze_from_tsv_calls_external_funcs(
        self, sample_tsv: Path, tmp_path: Path, monkeypatch
    ):
        matGen_called = False
        cosmic_fit_called = False

        def mock_matGen(*args, **kwargs):
            nonlocal matGen_called
            matGen_called = True

        def mock_cosmic_fit(*args, **kwargs):
            nonlocal cosmic_fit_called
            cosmic_fit_called = True

        monkeypatch.setattr(core.matGen, "SigProfilerMatrixGeneratorFunc", mock_matGen)
        monkeypatch.setattr(core, "cosmic_fit_func", mock_cosmic_fit)

        cfg = core.SampleConfig(sample_name="test_sample", output_dir=tmp_path)
        core.analyze_from_tsv(input_tsv=sample_tsv, cfg=cfg)

        assert matGen_called is True
        assert cosmic_fit_called is True


class TestInstallReference:
    def test_install_reference_calls_gen_install(self, monkeypatch):
        installed_kwargs = {}

        def mock_install(ref, rsync=True, bash=True):
            installed_kwargs["ref"] = ref
            installed_kwargs["rsync"] = rsync
            installed_kwargs["bash"] = bash

        monkeypatch.setattr(core.genInstall, "install", mock_install)

        core.install_reference("GRCh38", rsync=False, bash=False)

        assert installed_kwargs["ref"] == "GRCh38"
        assert installed_kwargs["rsync"] is False
        assert installed_kwargs["bash"] is False


class TestOutputFormatting:
    def test_analyze_sample_generates_tsv_and_html_with_plot(
        self, sample_tsv: Path, tmp_path: Path, monkeypatch
    ):
        def mock_matGen(*args, **kwargs):
            plot_dir = tmp_path / "matrix" / "plots"
            plot_dir.mkdir(parents=True, exist_ok=True)
            src_pdf = Path(__file__).parent.parent / "src" / "emcSP" / "SBS_96_plots_BRCA_example.pdf"
            if src_pdf.exists():
                import shutil
                shutil.copy(src_pdf, plot_dir / "SBS_96_plots_test_sample.pdf")

        def mock_cosmic_fit(*args, **kwargs):
            activities_dir = tmp_path / "output" / "Assignment_Solution" / "Activities"
            activities_dir.mkdir(parents=True, exist_ok=True)
            activities_file = activities_dir / "Assignment_Solution_Activities.txt"
            activities_file.write_text("Samples\tSBS1\tSBS5\ntest_sample\t10.0\t20.0\n")

            # Mock stats file
            stats_dir = tmp_path / "output" / "Assignment_Solution" / "Solution_Stats"
            stats_dir.mkdir(parents=True, exist_ok=True)
            stats_file = stats_dir / "Assignment_Solution_Samples_Stats.txt"
            stats_file.write_text("Sample Names\tTotal Mutations\tCosine Similarity\tL1 Norm\tL1_Norm_%\tL2 Norm\tL2_Norm_%\tKL Divergence\tCorrelation\ntest_sample\t30\t0.985\t10.0\t100%\t2.0\t50%\t0.05\t0.990\n")

            # Mock decomposition plot file
            recon_dir = activities_dir / "SampleReconstruction" / "WebPNGs"
            recon_dir.mkdir(parents=True, exist_ok=True)
            decomp_plot = recon_dir / "test_sample.png"
            decomp_plot.write_bytes(b"dummy_png_data")

        monkeypatch.setattr(core.matGen, "SigProfilerMatrixGeneratorFunc", mock_matGen)
        monkeypatch.setattr(core, "cosmic_fit_func", mock_cosmic_fit)

        cfg = core.SampleConfig(
            sample_name="test_sample",
            output_dir=tmp_path,
            output_format="both",
            plot=True
        )

        core.analyze_from_tsv(input_tsv=sample_tsv, cfg=cfg)

        tsv_out = tmp_path / "test_sample_signatures.tsv"
        assert tsv_out.exists()
        tsv_content = tsv_out.read_text().splitlines()
        assert tsv_content[0] == "Signature\tRelative Contribution\tMutation Count\tEtiology"
        
        # Verify both rows exist with correct format
        # Sorted desc by mutation count: SBS5 first, then SBS1
        assert "SBS5\t0.667\t20\tUnknown (Clock-like)" in tsv_content[1]
        assert "SBS1\t0.333\t10\tDeamination of 5-methylcytosine" in tsv_content[2]

        html_out = tmp_path / "test_sample_signatures.html"
        assert html_out.exists()
        html_content = html_out.read_text()
        assert "Mutational Signature Analysis" in html_content
        assert "test_sample" in html_content
        assert "SBS1" in html_content
        assert "SBS5" in html_content
        assert "data:image/png;base64," in html_content

        # Verify stats and decomposition plot in HTML
        assert "Reconstruction Quality" in html_content
        assert "0.985" in html_content
        assert "0.990" in html_content
        assert "0.050" in html_content
        assert "30" in html_content
        assert "Signature Decomposition (Reconstruction) Plot" in html_content
        assert "data:image/png;base64,ZHVtbXlfcG5nX2RhdGE=" in html_content

