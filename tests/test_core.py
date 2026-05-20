import pytest
from pathlib import Path
from emcSP import core

@pytest.fixture
def sample_tsv(tmp_path: Path) -> Path:
    tsv_content = "CHROM\tPOS\tREF\tALT\nchr1\t123456\tA\tG\nchr2\t789012\tC\tT\nchrX\t456789\tG\tA\n"
    tsv_path = tmp_path / "sample.tsv"
    tsv_path.write_text(tsv_content)
    return tsv_path

class TestTsvToVcf:
    def test_tsv_to_vcf_creates_correct_vcf(self, sample_tsv: Path, tmp_path: Path):
        result = core.tsv_to_vcf(sample_tsv, "test_sample")

        assert result.exists()
        assert result.suffix == ".vcf"

        lines = result.read_text().split("\n")
        assert lines[0] == "##fileformat=VCFv4.2"
        assert lines[1] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest_sample"
        assert lines[2] == "chr1\t123456\t.\tA\tG\t.\t.\t.\t.\t."
        assert lines[3] == "chr2\t789012\t.\tC\tT\t.\t.\t.\t.\t."
        assert lines[4] == "chrX\t456789\t.\tG\tA\t.\t.\t.\t.\t."

    def test_tsv_to_vcf_returns_vcf_path(self, sample_tsv: Path):
        vcf_path = core.tsv_to_vcf(sample_tsv, "sample1")
        expected = sample_tsv.with_suffix(".vcf")
        assert vcf_path == expected

    def test_tsv_to_vcf_skips_header(self, sample_tsv: Path, tmp_path: Path):
        vcf_path = core.tsv_to_vcf(sample_tsv, "sample1")
        lines = vcf_path.read_text().split("\n")
        assert "CHROM" not in lines[-1]

class TestAnalyzeFromTsv:
    def test_analyze_from_tsv_calls_external_funcs(self, sample_tsv: Path, tmp_path: Path, monkeypatch):
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

        core.analyze_from_tsv(
            input_tsv=sample_tsv,
            sample_name="test_sample",
            output_dir=tmp_path
        )

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
