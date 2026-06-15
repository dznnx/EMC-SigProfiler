#!/usr/bin/env python3

from pathlib import Path

import pytest
from emcSP import core
from emcSP.cli import app
from typer.testing import CliRunner


@pytest.fixture
def sample_tsv(tmp_path: Path) -> Path:
    tsv_content = "CHROM\tPOS\tREF\tALT\nchr1\t123456\tA\tG\nchr2\t789012\tC\tT\nchrX\t456789\tG\tA\n"
    tsv_path = tmp_path / "sample.tsv"
    tsv_path.write_text(tsv_content)
    return tsv_path


@pytest.fixture
def runner():
    return CliRunner()


class TestAnalyzeCommand:
    def test_analyze_command_requires_input_file(self, runner: CliRunner):
        result = runner.invoke(app, ["analyze"])
        assert result.exit_code != 0
        assert "Missing argument" in result.output or "Error" in result.output

    def test_analyze_command_with_valid_args_calls_core_func(
        self, runner: CliRunner, sample_tsv: Path, tmp_path: Path, monkeypatch
    ):
        called_args = {}

        def mock_analyze(input_tsv, cfg, **kwargs):
            called_args["cfg"] = cfg
            called_args["input_tsv"] = input_tsv

        monkeypatch.setattr("emcSP.cli.analyze_from_tsv", mock_analyze)

        result = runner.invoke(
            app,
            [
                "analyze",
                str(sample_tsv),
                "--sample-name",
                "test_sample",
                "--output",
                str(tmp_path),
            ],
        )

        assert result.exit_code == 0
        assert called_args["input_tsv"] == sample_tsv
        assert called_args["cfg"].sample_name == "test_sample"

    def test_analyze_command_validates_input_exists(self, runner: CliRunner):
        result = runner.invoke(
            app, ["analyze", "nonexistent.tsv", "--sample-name", "test"]
        )
        assert result.exit_code != 0
        assert "does not exist" in result.output or "Error" in result.output

    def test_analyze_command_requires_sample_name(
        self, runner: CliRunner, sample_tsv: Path
    ):
        result = runner.invoke(app, ["analyze", str(sample_tsv)])
        assert result.exit_code != 0
        assert "Missing option" in result.output or "Error" in result.output

    def test_analyze_command_with_plot_option(
        self, runner: CliRunner, sample_tsv: Path, tmp_path: Path, monkeypatch
    ):
        called_cfg = {}

        def mock_analyze(input_tsv, cfg, **kwargs):
            called_cfg["cfg"] = cfg

        monkeypatch.setattr("emcSP.cli.analyze_from_tsv", mock_analyze)

        result = runner.invoke(
            app,
            [
                "analyze",
                str(sample_tsv),
                "--sample-name",
                "test_sample",
                "--output",
                str(tmp_path),
                "--plot",
            ],
        )

        assert result.exit_code == 0
        assert called_cfg["cfg"].plot is True

    def test_analyze_command_with_custom_context_type(
        self, runner: CliRunner, sample_tsv: Path, tmp_path: Path, monkeypatch
    ):
        called_cfg = {}

        def mock_analyze(input_tsv, cfg, **kwargs):
            called_cfg["cfg"] = cfg

        monkeypatch.setattr("emcSP.cli.analyze_from_tsv", mock_analyze)

        result = runner.invoke(
            app,
            [
                "analyze",
                str(sample_tsv),
                "--sample-name",
                "test_sample",
                "--output",
                str(tmp_path),
                "--context-type",
                "1536",
            ],
        )

        assert result.exit_code == 0
        assert called_cfg["cfg"].context_type == "1536"

    def test_analyze_command_with_format_html(
        self, runner: CliRunner, sample_tsv: Path, tmp_path: Path, monkeypatch
    ):
        called_cfg = {}

        def mock_analyze(input_tsv, cfg, **kwargs):
            called_cfg["cfg"] = cfg

        monkeypatch.setattr("emcSP.cli.analyze_from_tsv", mock_analyze)

        result = runner.invoke(
            app,
            [
                "analyze",
                str(sample_tsv),
                "--sample-name",
                "test_sample",
                "--output",
                str(tmp_path),
                "--format",
                "html",
            ],
        )

        assert result.exit_code == 0
        assert called_cfg["cfg"].output_format == "html"
        assert called_cfg["cfg"].plot is True

    def test_analyze_command_with_invalid_format(
        self, runner: CliRunner, sample_tsv: Path, tmp_path: Path
    ):
        result = runner.invoke(
            app,
            [
                "analyze",
                str(sample_tsv),
                "--sample-name",
                "test_sample",
                "--output",
                str(tmp_path),
                "--format",
                "invalid",
            ],
        )

        assert result.exit_code != 0
        assert "Invalid output format" in result.output


class TestInstallCommand:
    def test_install_command_exists(self, runner: CliRunner):
        result = runner.invoke(app, ["install", "--help"])
        assert result.exit_code == 0
        assert "Install reference genome" in result.output

    def test_install_command_calls_core_install(self, runner: CliRunner, monkeypatch):
        called_kwargs = {}

        def mock_install(ref, rsync=True, bash=True):
            called_kwargs["ref"] = ref

        monkeypatch.setattr("emcSP.cli.install_reference", mock_install)

        result = runner.invoke(app, ["install", "--reference", "GRCh38"])
        assert result.exit_code == 0
        assert called_kwargs.get("ref") == "GRCh38"

    def test_install_command_default_reference(self, runner: CliRunner, monkeypatch):
        called_kwargs = {}

        def mock_install(ref, rsync=True, bash=True):
            called_kwargs["ref"] = ref

        monkeypatch.setattr("emcSP.cli.install_reference", mock_install)

        result = runner.invoke(app, ["install"])
        assert result.exit_code == 0
        assert called_kwargs.get("ref") == "GRCh37"

    def test_install_command_rsync_option(self, runner: CliRunner, monkeypatch):
        called_kwargs = {}

        def mock_install(ref, rsync=True, bash=True):
            called_kwargs["rsync"] = rsync

        monkeypatch.setattr("emcSP.cli.install_reference", mock_install)

        result = runner.invoke(app, ["install", "--no-rsync"])
        assert result.exit_code == 0
        assert called_kwargs.get("rsync") is False

    def test_install_command_bash_option(self, runner: CliRunner, monkeypatch):
        called_kwargs = {}

        def mock_install(ref, rsync=True, bash=True):
            called_kwargs["bash"] = bash

        monkeypatch.setattr("emcSP.cli.install_reference", mock_install)

        result = runner.invoke(app, ["install", "--no-bash"])
        assert result.exit_code == 0
        assert called_kwargs.get("bash") is False
