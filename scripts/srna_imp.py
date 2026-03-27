#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import yaml

WORKFLOWS = ["preprocessing", "host_specific", "nonhost_specific", "assembly_based", "novel_ncrna", "summary"]
DOWNSTREAM_NEEDS_PREPROCESS = {"host_specific", "nonhost_specific", "assembly_based", "novel_ncrna"}


class PipelineContext:
    def __init__(self, config: dict[str, Any]):
        self.config = config
        self.sample_name = config["sample"]["name"]
        self.sample_type = config["sample"].get("type", "sRNA-Ext")
        self.root_dir = Path(config["project"]["root_dir"])
        self.results_dir = Path(config["project"]["results_dir"])
        self.logs_dir = Path(config["project"]["logs_dir"])
        self.sample_dir = self.results_dir / self.sample_type / self.sample_name
        self.analysis_root = self.sample_dir / "Analysis"
        self.analysis_type_dir = self.analysis_root / self.sample_type
        self.annotation_dir = self.analysis_root / "annotation"
        self.taxonomy_dir = self.analysis_root / "taxonomy"
        self.kraken_dir = self.taxonomy_dir / "kraken"
        self.preprocess_root = self.sample_dir / "Preprocessing"
        self.preprocess_dir = self.preprocess_root / self.sample_type
        self.stats_root = self.sample_dir / "Stats"
        self.stats_dir = self.stats_root / self.sample_type
        self.assembly_root = self.sample_dir / "Assembly"
        self.host_dir = self.analysis_type_dir / "human"
        self.nonhost_dir = self.analysis_type_dir / "nonhuman"
        self.assembly_dir = self.analysis_type_dir
        self.novel_dir = self.analysis_type_dir / "novel"
        self.summary_dir = self.stats_dir
        self.logs_cluster_jobs_root = self.sample_dir / "logs_cluster_jobs"
        self.report_logs_dir = self.logs_cluster_jobs_root / self.sample_type
        self.helper_scripts_dir = Path(config["project"].get("helper_scripts_dir", self.root_dir / "scripts"))
        self.raw_fastq = Path(config["input"]["fastq"]).resolve()
        self.assembly_fasta = Path(config["input"]["assembly_fasta"]).resolve() if config["input"].get("assembly_fasta") else None
        self.assembly_gff = Path(config["input"]["assembly_gff"]).resolve() if config["input"].get("assembly_gff") else None
        self.threads = int(config["settings"].get("threads", 8))
        self.trimmed_fastq = self.preprocess_dir / f"{self.sample_name}.trimmed.fastq.gz"
        self.host_fastq = self.preprocess_dir / f"{self.sample_name}.host.fastq.gz"
        self.nonhost_fastq = self.preprocess_dir / f"{self.sample_name}.nonhost.fastq.gz"
        self.preprocess_stats = self.preprocess_dir / f"{self.sample_name}.preprocessing.summary.tsv"
        self.metadata_json = self.sample_dir / "sample.json"

    def ensure_dirs(self) -> None:
        for path in [
            self.results_dir,
            self.logs_dir,
            self.sample_dir,
            self.analysis_root,
            self.analysis_type_dir,
            self.annotation_dir,
            self.taxonomy_dir,
            self.kraken_dir,
            self.preprocess_root,
            self.preprocess_dir,
            self.stats_root,
            self.stats_dir,
            self.assembly_root,
            self.host_dir,
            self.nonhost_dir,
            self.assembly_dir,
            self.novel_dir,
            self.summary_dir,
            self.logs_cluster_jobs_root,
            self.report_logs_dir,
            self.report_logs_dir / "human",
            self.report_logs_dir / "nonhuman",
            self.report_logs_dir / "novel",
            self.host_dir / "mirna",
            self.host_dir / "otherrna",
            self.host_dir / "rrna",
            self.host_dir / "trna",
            self.host_dir / "unknown",
            self.nonhost_dir / "collapsed",
            self.nonhost_dir / "mirna_like",
            self.nonhost_dir / "rfam",
            self.nonhost_dir / "rrna",
            self.nonhost_dir / "rrna_tax",
            self.nonhost_dir / "trna",
            self.nonhost_dir / "unknown",
        ]:
            path.mkdir(parents=True, exist_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the native sRNA-IMP pipeline.")
    parser.add_argument("--config", required=True, help="Path to YAML config file.")
    parser.add_argument("--workflow", choices=["all"] + WORKFLOWS, default="all", help="Run one workflow or the full pipeline.")
    parser.add_argument("--sample", action="append", default=[], help="Sample name to run. Repeat to select multiple samples from a multi-sample config.")
    parser.add_argument("--keep-going", action="store_true", help="Continue with remaining samples if one sample fails in multi-sample mode.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them.")
    return parser.parse_args()


def load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError("Config file must contain a YAML mapping.")
    return config



def load_samples_from_tsv(tsv_path: Path) -> list[dict[str, Any]]:
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        if reader.fieldnames is None:
            raise ValueError(f"Sample TSV is empty: {tsv_path}")
        required = {"name", "fastq"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"Sample TSV missing required column(s): {', '.join(missing)}")
        samples: list[dict[str, Any]] = []
        for row in reader:
            cleaned = {key: (value.strip() if isinstance(value, str) else value) for key, value in row.items()}
            if not any(cleaned.values()):
                continue
            sample = {
                "name": cleaned.get("name", ""),
                "type": cleaned.get("type", ""),
                "fastq": cleaned.get("fastq", ""),
                "assembly_fasta": cleaned.get("assembly_fasta", ""),
                "assembly_gff": cleaned.get("assembly_gff", ""),
            }
            if not sample["name"] or not sample["fastq"]:
                raise ValueError(f"Sample TSV row missing required values in {tsv_path}: {row}")
            samples.append(sample)
    if not samples:
        raise ValueError(f"Sample TSV contains no usable sample rows: {tsv_path}")
    return samples


def build_sample_config(base_config: dict[str, Any], sample_entry: dict[str, Any], sample_defaults: dict[str, Any], input_defaults: dict[str, Any]) -> dict[str, Any]:
    config = copy.deepcopy({
        key: value
        for key, value in base_config.items()
        if key not in {"sample", "input", "samples", "sample_defaults", "input_defaults"}
    })
    sample_block = sample_entry.get("sample", {})
    input_block = sample_entry.get("input", {})

    if "fastq" in sample_entry and "fastq" not in input_block:
        input_block["fastq"] = sample_entry["fastq"]
    if "assembly_fasta" in sample_entry and "assembly_fasta" not in input_block:
        input_block["assembly_fasta"] = sample_entry["assembly_fasta"]
    if "assembly_gff" in sample_entry and "assembly_gff" not in input_block:
        input_block["assembly_gff"] = sample_entry["assembly_gff"]

    sample_name = sample_entry.get("name", sample_block.get("name"))
    sample_type = sample_entry.get("type", sample_block.get("type", sample_defaults.get("type", "sRNA-Ext")))
    config["sample"] = {"name": sample_name, "type": sample_type}
    config["input"] = {
        "fastq": input_block.get("fastq", input_defaults.get("fastq")),
        "assembly_fasta": input_block.get("assembly_fasta", input_defaults.get("assembly_fasta")),
        "assembly_gff": input_block.get("assembly_gff", input_defaults.get("assembly_gff")),
    }
    return config


def expand_configs(config: dict[str, Any], selected_samples: list[str] | None = None) -> list[dict[str, Any]]:
    sample_defaults = config.get("sample_defaults", {})
    input_defaults = config.get("input_defaults", {})

    if "samples_tsv" in config:
        tsv_path = Path(str(config["samples_tsv"])).resolve()
        sample_entries = load_samples_from_tsv(tsv_path)
        configs = [build_sample_config(config, sample_entry, sample_defaults, input_defaults) for sample_entry in sample_entries]
    elif "samples" not in config:
        configs = [config]
    else:
        samples = config.get("samples")
        if not isinstance(samples, list) or not samples:
            raise ValueError("Multi-sample config must contain a non-empty 'samples' list.")
        configs = [build_sample_config(config, sample_entry, sample_defaults, input_defaults) for sample_entry in samples]

    sample_names = [cfg.get("sample", {}).get("name") for cfg in configs]
    if any(not name for name in sample_names):
        raise ValueError("Each sample config must define sample.name.")
    duplicates = sorted({name for name in sample_names if sample_names.count(name) > 1})
    if duplicates:
        raise ValueError(f"Duplicate sample names in config: {', '.join(duplicates)}")

    if selected_samples:
        wanted = set(selected_samples)
        configs = [cfg for cfg in configs if cfg["sample"]["name"] in wanted]
        found = {cfg["sample"]["name"] for cfg in configs}
        missing = sorted(wanted - found)
        if missing:
            raise ValueError(f"Requested sample(s) not found in config: {', '.join(missing)}")
    return configs


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def write_text(path: Path, content: str) -> None:
    ensure_parent(path)
    path.write_text(content, encoding="utf-8")


def append_text(path: Path, content: str) -> None:
    ensure_parent(path)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(content)


def run_command(command: list[str], log_path: Path, dry_run: bool, cwd: Path | None = None) -> None:
    ensure_parent(log_path)
    rendered = " ".join(command)
    prefix = f"(cd {cwd} && " if cwd else ""
    suffix = ")" if cwd else ""
    if dry_run:
        print(f"[dry-run] {prefix}{rendered}{suffix}")
        return
    with log_path.open("w", encoding="utf-8") as log_handle:
        process = subprocess.run(command, cwd=cwd, stdout=log_handle, stderr=subprocess.STDOUT)
    if process.returncode != 0:
        raise RuntimeError(f"Command failed ({process.returncode}): {rendered}")


def run_capture(command: list[str], dry_run: bool) -> str:
    if dry_run:
        return "0"
    process = subprocess.run(command, check=True, capture_output=True, text=True)
    return process.stdout.strip()


def run_command_stdout_to_file(command: list[str], output_path: Path, log_path: Path, dry_run: bool) -> None:
    ensure_parent(output_path)
    ensure_parent(log_path)
    rendered = " ".join(command)
    if dry_run:
        print(f"[dry-run] {rendered} > {output_path}")
        return
    with output_path.open("w", encoding="utf-8") as out_handle, log_path.open("w", encoding="utf-8") as log_handle:
        process = subprocess.run(command, stdout=out_handle, stderr=log_handle, text=True)
    if process.returncode != 0:
        raise RuntimeError(f"Command failed ({process.returncode}): {rendered}")


def validate_exists(path: Path | None, description: str) -> None:
    if path is None:
        raise ValueError(f"Missing required input: {description}")
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {path}")


def validate_executable(name_or_path: str) -> None:
    if "/" in name_or_path:
        if not Path(name_or_path).exists():
            raise FileNotFoundError(f"Required executable not found: {name_or_path}")
    elif shutil.which(name_or_path) is None:
        raise FileNotFoundError(f"Required executable not found in PATH: {name_or_path}")


def count_fastq_reads(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    line_count = 0
    with opener(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for _ in handle:
            line_count += 1
    return line_count // 4


def copy_or_symlink(src: Path, dst: Path, dry_run: bool) -> None:
    ensure_parent(dst)
    if dry_run:
        print(f"[dry-run] ln -s {src} {dst}")
        return
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src)


def gzip_empty_fastq(dst: Path, dry_run: bool) -> None:
    ensure_parent(dst)
    if dry_run:
        print(f"[dry-run] create empty fastq {dst}")
        return
    with gzip.open(dst, "wt", encoding="utf-8") as handle:
        handle.write("")


def stage_plain_copy(src: Path, dst: Path, dry_run: bool) -> None:
    ensure_parent(dst)
    if dry_run:
        print(f"[dry-run] stage plain {src} -> {dst}")
        return
    if src.suffix == ".gz" and dst.suffix != ".gz":
        with gzip.open(src, "rt", encoding="utf-8") as in_handle, dst.open("w", encoding="utf-8") as out_handle:
            shutil.copyfileobj(in_handle, out_handle)
    else:
        shutil.copy2(src, dst)


def stage_gzip_copy(src: Path, dst: Path, dry_run: bool) -> None:
    ensure_parent(dst)
    if dry_run:
        print(f"[dry-run] stage gzip {src} -> {dst}")
        return
    if src.suffix == ".gz":
        shutil.copy2(src, dst)
    else:
        with src.open("rb") as in_handle, gzip.open(dst, "wb") as out_handle:
            shutil.copyfileobj(in_handle, out_handle)


def sort_bed_in_place(path: Path, dry_run: bool) -> None:
    if dry_run:
        print(f"[dry-run] sort -k1,1 -k2,2n {path} -o {path}")
        return
    subprocess.run(["sort", "-k1,1", "-k2,2n", str(path), "-o", str(path)], check=True)


def format_normalized(value: float | int) -> str:
    return f"{float(value):.2f}"


def sum_count_table(path: Path) -> float:
    total = 0.0
    with path.open("r", encoding="utf-8") as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[-1]:
                total += float(parts[-1])
    return total


def get_tool(config: dict[str, Any], key: str, default: str | None = None) -> str:
    tools = config.get("tools", {})
    value = tools.get(key, default)
    if value is None:
        raise KeyError(f"Missing tool configuration: {key}")
    return str(value)


def get_workflow_config(config: dict[str, Any], workflow: str) -> dict[str, Any]:
    return config.get("workflows", {}).get(workflow, {})


def bowtie_index_exists(index_base: Path) -> bool:
    small = [index_base.with_suffix(suffix) for suffix in [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]]
    large = [index_base.with_suffix(suffix) for suffix in [".1.ebwtl", ".2.ebwtl", ".3.ebwtl", ".4.ebwtl", ".rev.1.ebwtl", ".rev.2.ebwtl"]]
    return all(path.exists() for path in small) or all(path.exists() for path in large)


def shared_bowtie_index_base(fasta: Path) -> Path:
    base = fasta
    if base.suffix == ".gz":
        base = base.with_suffix("")
    return base


def workflow_marker_path(ctx: PipelineContext, workflow: str) -> Path:
    return ctx.sample_dir / f".{workflow}.complete"


def workflow_required_outputs(ctx: PipelineContext, workflow: str) -> list[Path]:
    if workflow == "preprocessing":
        return [ctx.trimmed_fastq, ctx.host_fastq, ctx.nonhost_fastq, ctx.preprocess_stats]
    if workflow == "host_specific":
        return [ctx.host_dir / f"{ctx.sample_name}.human.summary.tsv"]
    if workflow == "nonhost_specific":
        return [
            ctx.nonhost_dir / "rrna" / f"{ctx.sample_name}.rRNA.fastq.gz",
            ctx.nonhost_dir / "trna" / f"{ctx.sample_name}.trna.counts.tsv",
            ctx.nonhost_dir / "mirna_like" / f"{ctx.sample_name}.mirna_like.counts.tsv",
            ctx.nonhost_dir / "rfam" / f"{ctx.sample_name}.rfam_class.counts.tsv",
        ]
    if workflow == "assembly_based":
        return [
            ctx.assembly_dir / f"{ctx.sample_name}.RawCounts.txt",
            ctx.assembly_dir / f"{ctx.sample_name}.WeightedCounts.txt",
        ]
    if workflow == "novel_ncrna":
        return [ctx.novel_dir / f"{ctx.sample_name}.novel.summary.tsv"]
    if workflow == "summary":
        return [ctx.summary_dir / f"{ctx.sample_name}.pipeline.summary.tsv"]
    if workflow == "multiqc":
        return [ctx.sample_dir / "multiqc" / "multiqc_report.html"]
    return []


def workflow_complete(ctx: PipelineContext, workflow: str) -> bool:
    marker = workflow_marker_path(ctx, workflow)
    required = workflow_required_outputs(ctx, workflow)
    return marker.exists() and all(path.exists() for path in required)


def mark_workflow_complete(ctx: PipelineContext, workflow: str, dry_run: bool) -> None:
    if dry_run:
        print(f"[dry-run] mark {workflow} complete")
        return
    write_text(workflow_marker_path(ctx, workflow), "complete\n")


def write_metadata(ctx: PipelineContext, dry_run: bool) -> None:
    data = {
        "sample": ctx.sample_name,
        "type": ctx.sample_type,
        "input_fastq": str(ctx.raw_fastq),
        "assembly_fasta": str(ctx.assembly_fasta) if ctx.assembly_fasta else None,
        "assembly_gff": str(ctx.assembly_gff) if ctx.assembly_gff else None,
        "sample_dir": str(ctx.sample_dir),
    }
    if dry_run:
        print(f"[dry-run] write metadata {ctx.metadata_json}")
        return
    write_text(ctx.metadata_json, json.dumps(data, indent=2))


def validate_config(config: dict[str, Any], workflow: str) -> None:
    for key in ["project", "sample", "input", "settings"]:
        if key not in config:
            raise KeyError(f"Missing config section: {key}")
    validate_exists(Path(config["input"]["fastq"]).resolve(), "input FASTQ")

    required_tools = {"python", "multiqc"}
    preprocessing_cfg = config.get("workflows", {}).get("preprocessing", {})
    host_filter_indexes = preprocessing_cfg.get("host_filter_indexes", [])
    host_reference_fastas = preprocessing_cfg.get("host_reference_fastas", [])
    needs_preprocessing = workflow in {"all", "preprocessing", "host_specific", "nonhost_specific", "assembly_based", "novel_ncrna"}
    if needs_preprocessing and not (host_filter_indexes or host_reference_fastas):
        raise ValueError(
            "preprocessing requires host filtering references: set workflows.preprocessing.host_reference_fastas "
            "or workflows.preprocessing.host_filter_indexes"
        )
    if needs_preprocessing:
        required_tools.update({"fastqc", "fastp", "seqkit", "samtools", "bowtie"})
    if host_reference_fastas:
        required_tools.add("bowtie-build")
    if workflow in {"all", "nonhost_specific", "assembly_based", "novel_ncrna"}:
        required_tools.add("bedtools")
    if workflow in {"all", "nonhost_specific"}:
        required_tools.add("ribodetector_cpu")
    if workflow in {"all", "assembly_based"}:
        required_tools.update({"bwa", "seqtk", "fastx_collapser", "gff2bed", "srnaMapper"})
    if workflow in {"all", "novel_ncrna"}:
        required_tools.update({"bowtie-build", "bowtie", "sort", "awk", "cut", "uniq"})
    for tool in sorted(required_tools):
        validate_executable(get_tool(config, tool, tool))


def run_preprocessing(ctx: PipelineContext, dry_run: bool) -> None:
    config = ctx.config
    workflow_cfg = get_workflow_config(config, "preprocessing")
    min_len = int(workflow_cfg.get("min_length", 15))
    max_len = int(workflow_cfg.get("max_length", 250))
    adapter = str(workflow_cfg.get("adapter_sequence", ""))
    host_indexes = list(workflow_cfg.get("host_filter_indexes", []))
    host_reference_fastas = workflow_cfg.get("host_reference_fastas", [])
    run_fastqc = bool(workflow_cfg.get("run_fastqc", True))

    write_metadata(ctx, dry_run)

    if run_fastqc:
        run_command(
            [get_tool(config, "fastqc", "fastqc"), "-o", str(ctx.preprocess_dir), str(ctx.raw_fastq)],
            ctx.report_logs_dir / f"{ctx.sample_name}.fastqc.raw.log",
            dry_run,
        )

    fastp_cmd = [
        get_tool(config, "fastp", "fastp"),
        "-i",
        str(ctx.raw_fastq),
        "-o",
        str(ctx.trimmed_fastq),
        "--length_required",
        str(min_len),
        "--length_limit",
        str(max_len),
        "--thread",
        str(ctx.threads),
        "--json",
        str(ctx.preprocess_dir / f"{ctx.sample_name}.fastp.json"),
        "--html",
        str(ctx.preprocess_dir / f"{ctx.sample_name}.fastp.html"),
    ]
    if adapter:
        fastp_cmd.extend(["--adapter_sequence", adapter])
    run_command(fastp_cmd, ctx.report_logs_dir / f"{ctx.sample_name}.fastp.log", dry_run)

    if run_fastqc:
        run_command(
            [get_tool(config, "fastqc", "fastqc"), "-o", str(ctx.preprocess_dir), str(ctx.trimmed_fastq)],
            ctx.report_logs_dir / f"{ctx.sample_name}.fastqc.trimmed.log",
            dry_run,
        )

    if host_reference_fastas:
        built_indexes = []
        for index_no, fasta_path in enumerate(host_reference_fastas, start=1):
            fasta = Path(fasta_path).resolve()
            validate_exists(fasta, f"host reference FASTA #{index_no}")
            index_base = shared_bowtie_index_base(fasta)
            if dry_run:
                print(f"[dry-run] reuse shared Bowtie index for {index_base} when available, otherwise build from {fasta}")
            elif not bowtie_index_exists(index_base):
                run_command(
                    [get_tool(config, "bowtie-build", "bowtie-build"), str(fasta), str(index_base)],
                    ctx.report_logs_dir / f"{ctx.sample_name}.host_index_build.{index_no}.log",
                    False,
                )
            built_indexes.append(index_base)
        host_indexes.extend(str(path) for path in built_indexes)

    if not host_indexes:
        raise ValueError(
            "preprocessing requires at least one host filtering reference from host_reference_fastas or host_filter_indexes"
        )
    else:
        current_input = ctx.trimmed_fastq
        host_accumulator = ctx.preprocess_dir / f"{ctx.sample_name}.host.accum.fastq"
        if not dry_run and host_accumulator.exists():
            host_accumulator.unlink()
        tmp_nonhost_files: list[Path] = []
        for index_no, index_path in enumerate(host_indexes, start=1):
            tmp_host = ctx.preprocess_dir / f"{ctx.sample_name}.host.part{index_no}.fastq"
            tmp_un = ctx.preprocess_dir / f"{ctx.sample_name}.nonhost.part{index_no}.fastq"
            cmd = [
                get_tool(config, "bowtie", "bowtie"),
                "-p",
                str(ctx.threads),
                "-v",
                "1",
                "-k",
                "1",
                "--best",
                "--al",
                str(tmp_host),
                "--un",
                str(tmp_un),
                str(index_path),
                str(current_input),
            ]
            run_command(cmd, ctx.report_logs_dir / f"{ctx.sample_name}.host_filter.{index_no}.log", dry_run)
            if not dry_run and tmp_host.exists():
                if host_accumulator.exists():
                    append_text(host_accumulator, tmp_host.read_text(encoding="utf-8"))
                    tmp_host.unlink()
                else:
                    tmp_host.rename(host_accumulator)
            tmp_nonhost_files.append(tmp_un)
            current_input = tmp_un

        if dry_run:
            print(f"[dry-run] gzip host/nonhost split to {ctx.host_fastq} and {ctx.nonhost_fastq}")
        else:
            if host_accumulator.exists():
                run_command(["gzip", "-f", str(host_accumulator)], ctx.report_logs_dir / f"{ctx.sample_name}.host.gzip.log", False)
                shutil.move(str(host_accumulator) + ".gz", ctx.host_fastq)
            else:
                gzip_empty_fastq(ctx.host_fastq, False)
            run_command(["gzip", "-f", str(current_input)], ctx.report_logs_dir / f"{ctx.sample_name}.nonhost.gzip.log", False)
            shutil.move(str(current_input) + ".gz", ctx.nonhost_fastq)
            for tmp_un in tmp_nonhost_files[:-1]:
                if tmp_un.exists():
                    tmp_un.unlink()

    if dry_run:
        print(f"[dry-run] write preprocessing summary {ctx.preprocess_stats}")
        return

    total = count_fastq_reads(ctx.raw_fastq)
    trimmed = count_fastq_reads(ctx.trimmed_fastq)
    host = count_fastq_reads(ctx.host_fastq)
    nonhost = count_fastq_reads(ctx.nonhost_fastq)
    write_text(
        ctx.preprocess_stats,
        "sample\tcategory\treads\n"
        f"{ctx.sample_name}\tinput\t{total}\n"
        f"{ctx.sample_name}\ttrimmed\t{trimmed}\n"
        f"{ctx.sample_name}\thost\t{host}\n"
        f"{ctx.sample_name}\tnon_host\t{nonhost}\n",
    )
    mark_workflow_complete(ctx, "preprocessing", False)


def run_host_specific(ctx: PipelineContext, dry_run: bool) -> None:
    cfg = get_workflow_config(ctx.config, "host_specific")
    required = ["mirna_index", "trna_index", "rrna_index", "other_ncrna_index"]
    missing = [key for key in required if not cfg.get(key)]
    if missing:
        raise ValueError(f"host_specific config missing: {', '.join(missing)}")

    host_input = ctx.host_fastq
    summary_path = ctx.host_dir / f"{ctx.sample_name}.human.summary.tsv"
    if not dry_run and host_input.exists() and host_input.stat().st_size == 0:
        write_text(
            summary_path,
            "sample\tcategory\treads\n"
            f"{ctx.sample_name}\ttrimmed\t0.00\n"
            f"{ctx.sample_name}\thuman_miRNA\t0.00\n"
            f"{ctx.sample_name}\thuman_tRNA\t0.00\n"
            f"{ctx.sample_name}\thuman_rRNA\t0.00\n"
            f"{ctx.sample_name}\thuman_otherRNA\t0.00\n"
            f"{ctx.sample_name}\tremaining_after_miRNA\t0.00\n"
            f"{ctx.sample_name}\tremaining_after_tRNA\t0.00\n"
            f"{ctx.sample_name}\tremaining_after_rRNA\t0.00\n"
            f"{ctx.sample_name}\tunknown_after_all\t0.00\n",
        )
        mark_workflow_complete(ctx, "host_specific", False)
        return

    helper_counts = ctx.helper_scripts_dir / "sam_fractional_counts.py"
    helper_summary = ctx.helper_scripts_dir / "summarize_human_subtree.py"
    validate_exists(helper_counts, "sam_fractional_counts.py")
    validate_exists(helper_summary, "summarize_human_subtree.py")

    categories = [
        ("mirna", cfg["mirna_index"], f"{ctx.sample_name}.human.mirna.sam", f"{ctx.sample_name}.human.mirna.counts.tsv", f"{ctx.sample_name}.human.non_mirna.fastq"),
        ("trna", cfg["trna_index"], f"{ctx.sample_name}.human.trna.sam", f"{ctx.sample_name}.human.trna.counts.tsv", f"{ctx.sample_name}.human.non_mirna.non_trna.fastq"),
        ("rrna", cfg["rrna_index"], f"{ctx.sample_name}.human.rrna.sam", f"{ctx.sample_name}.human.rrna.counts.tsv", f"{ctx.sample_name}.human.non_mirna.non_trna.non_rrna.fastq"),
        ("otherrna", cfg["other_ncrna_index"], f"{ctx.sample_name}.human.otherrna.sam", f"{ctx.sample_name}.human.otherrna.counts.tsv", f"{ctx.sample_name}.human.unknown.fastq"),
    ]

    current_input = host_input
    with tempfile.TemporaryDirectory(prefix=f"{ctx.sample_name}.human.") as tmpdir_name:
        tmpdir = Path(tmpdir_name)
        for category, index, sam_name, counts_name, unmapped_name in categories:
            out_dir = ctx.host_dir / category
            out_dir.mkdir(parents=True, exist_ok=True)
            sam_path = tmpdir / sam_name
            unmapped_plain = tmpdir / unmapped_name
            final_unmapped = out_dir / f"{unmapped_name}.gz"
            cmd = [
                get_tool(ctx.config, "bowtie", "bowtie"),
                "-p",
                str(ctx.threads),
                "-v",
                "1",
                "-a",
                "--best",
                "--sam",
                "--un",
                str(unmapped_plain),
                str(index),
                str(current_input),
            ]
            run_command_stdout_to_file(cmd, sam_path, ctx.report_logs_dir / f"{ctx.sample_name}.host.{category}.log", dry_run)
            counts_path = out_dir / counts_name
            run_command(
                [get_tool(ctx.config, "python", "python"), str(helper_counts), "-i", str(sam_path), "-o", str(counts_path)],
                ctx.report_logs_dir / f"{ctx.sample_name}.host.{category}.counts.log",
                dry_run,
            )
            if dry_run:
                print(f"[dry-run] gzip -f {unmapped_plain} && move to {final_unmapped}")
            else:
                run_command(["gzip", "-f", str(unmapped_plain)], ctx.report_logs_dir / f"{ctx.sample_name}.host.{category}.gzip.log", False)
                shutil.move(str(unmapped_plain) + ".gz", final_unmapped)
            current_input = final_unmapped

        run_command(
            [get_tool(ctx.config, "python", "python"), str(helper_summary), str(ctx.host_dir), "-o", str(summary_path)],
            ctx.report_logs_dir / f"{ctx.sample_name}.host.summary.log",
            dry_run,
        )

    if not dry_run:
        mark_workflow_complete(ctx, "host_specific", False)

def run_nonhost_specific(ctx: PipelineContext, dry_run: bool) -> None:
    cfg = get_workflow_config(ctx.config, "nonhost_specific")
    required = ["kraken_db", "bracken_db", "trna_index", "mirna_euk_index", "rfam_class_index", "rfamseq_tsv"]
    missing = [key for key in required if not cfg.get(key)]
    if missing:
        raise ValueError(f"nonhost_specific config missing: {', '.join(missing)}")

    helper_collapse = ctx.helper_scripts_dir / "collapse_fastq_to_counted_fasta.py"
    helper_mirna = ctx.helper_scripts_dir / "sum_collapsed_bowtie_hits.py"
    helper_trna = ctx.helper_scripts_dir / "sum_collapsed_trna_by_taxon.py"
    helper_rfam = ctx.helper_scripts_dir / "sum_collapsed_rfam_by_class_and_taxon.py"
    for helper in [helper_collapse, helper_mirna, helper_trna, helper_rfam]:
        validate_exists(helper, helper.name)

    nonhost_input = ctx.nonhost_fastq
    rrna_dir = ctx.nonhost_dir / "rrna"
    rrna_tax_dir = ctx.nonhost_dir / "rrna_tax"
    collapsed_dir = ctx.nonhost_dir / "collapsed"
    trna_dir = ctx.nonhost_dir / "trna"
    mirna_dir = ctx.nonhost_dir / "mirna_like"
    rfam_dir = ctx.nonhost_dir / "rfam"
    unknown_dir = ctx.nonhost_dir / "unknown"
    for path in [rrna_dir, rrna_tax_dir, collapsed_dir, trna_dir, mirna_dir, rfam_dir, unknown_dir]:
        path.mkdir(parents=True, exist_ok=True)

    non_rrna = rrna_dir / f"{ctx.sample_name}.nonrRNA.fastq"
    rrna_fastq = rrna_dir / f"{ctx.sample_name}.rRNA.fastq"
    if not dry_run:
        mean_len = run_capture([get_tool(ctx.config, "seqkit", "seqkit"), "stats", "-T", str(nonhost_input)], False)
        mean_field = mean_len.strip().splitlines()[-1].split("\t")[6]
        mean_read_length = str(int(float(mean_field)))
    else:
        mean_read_length = "50"

    run_command(
        [
            get_tool(ctx.config, "ribodetector_cpu", "ribodetector_cpu"),
            "--threads",
            str(ctx.threads),
            "--len",
            mean_read_length,
            "--input",
            str(nonhost_input),
            "--ensure",
            "none",
            "--chunk_size",
            "256",
            "--output",
            str(non_rrna),
            "--rrna",
            str(rrna_fastq),
        ],
        ctx.report_logs_dir / f"{ctx.sample_name}.ribodetector.log",
        dry_run,
    )

    rrna_fastq_gz = Path(str(rrna_fastq) + ".gz")
    non_rrna_fastq_gz = Path(str(non_rrna) + ".gz")
    run_command(["gzip", "-f", str(non_rrna), str(rrna_fastq)], ctx.report_logs_dir / f"{ctx.sample_name}.nonhost.gzip.log", dry_run)

    kraken_report = rrna_tax_dir / f"{ctx.sample_name}.rRNA.kraken.report"
    kraken_out = rrna_tax_dir / f"{ctx.sample_name}.rRNA.kraken.out"
    run_command(
        [
            get_tool(ctx.config, "kraken2", "kraken2"),
            "--db",
            str(cfg["kraken_db"]),
            "--threads",
            str(ctx.threads),
            "--gzip-compressed",
            "--memory-mapping",
            "--report",
            str(kraken_report),
            "--output",
            str(kraken_out),
            str(rrna_fastq_gz),
        ],
        ctx.report_logs_dir / f"{ctx.sample_name}.kraken2.log",
        dry_run,
    )

    collapsed_fa = collapsed_dir / f"{ctx.sample_name}.non_rRNA.collapsed.fa"
    run_command(
        [get_tool(ctx.config, "python", "python"), str(helper_collapse), "-i", str(non_rrna_fastq_gz), "-o", str(collapsed_fa)],
        ctx.report_logs_dir / f"{ctx.sample_name}.collapse.log",
        dry_run,
    )

    trna_hits = trna_dir / f"{ctx.sample_name}.trna.hits.txt"
    trna_un = trna_dir / f"{ctx.sample_name}.non_tRNA.fa"
    run_command_stdout_to_file(
        [
            get_tool(ctx.config, "bowtie", "bowtie"),
            "-p",
            str(ctx.threads),
            "-v",
            "1",
            "-k",
            "1",
            "--best",
            "-f",
            "--un",
            str(trna_un),
            str(cfg["trna_index"]),
            str(collapsed_fa),
        ],
        trna_hits,
        ctx.report_logs_dir / f"{ctx.sample_name}.nonhost.trna.log",
        dry_run,
    )

    run_command(
        [
            get_tool(ctx.config, "python", "python"),
            str(helper_trna),
            "-i",
            str(trna_hits),
            "--feature-out",
            str(trna_dir / f"{ctx.sample_name}.trna.counts.tsv"),
            "--phylum-out",
            str(trna_dir / f"{ctx.sample_name}.trna_phylum_class.counts.tsv"),
            "--species-out",
            str(trna_dir / f"{ctx.sample_name}.trna_species.counts.tsv"),
            "--aa-out",
            str(trna_dir / f"{ctx.sample_name}.trna_aa.counts.tsv"),
            "--anticodon-out",
            str(trna_dir / f"{ctx.sample_name}.trna_anticodon.counts.tsv"),
            "--phylum-species-out",
            str(trna_dir / f"{ctx.sample_name}.trna_phylum_species.counts.tsv"),
        ],
        ctx.report_logs_dir / f"{ctx.sample_name}.nonhost.trna.counts.log",
        dry_run,
    )

    mirna_in = mirna_dir / f"{ctx.sample_name}.non_tRNA.L18_26.fa"
    longer_in = mirna_dir / f"{ctx.sample_name}.non_tRNA.GT26.fa"
    run_command_stdout_to_file([get_tool(ctx.config, "seqkit", "seqkit"), "seq", "-g", "-m", "18", "-M", "26", str(trna_un)], mirna_in, ctx.report_logs_dir / f"{ctx.sample_name}.mirna.split.short.log", dry_run)
    run_command_stdout_to_file([get_tool(ctx.config, "seqkit", "seqkit"), "seq", "-g", "-m", "27", str(trna_un)], longer_in, ctx.report_logs_dir / f"{ctx.sample_name}.mirna.split.long.log", dry_run)

    mirna_hits = mirna_dir / f"{ctx.sample_name}.mirna_like.hits.txt"
    mirna_un = mirna_dir / f"{ctx.sample_name}.non_tRNA_non_miRNA.L18_26.fa"
    run_command_stdout_to_file(
        [
            get_tool(ctx.config, "bowtie", "bowtie"),
            "-p",
            str(ctx.threads),
            "-v",
            "0",
            "-k",
            "1",
            "--best",
            "-f",
            "--un",
            str(mirna_un),
            str(cfg["mirna_euk_index"]),
            str(mirna_in),
        ],
        mirna_hits,
        ctx.report_logs_dir / f"{ctx.sample_name}.nonhost.mirna.log",
        dry_run,
    )
    run_command(
        [get_tool(ctx.config, "python", "python"), str(helper_mirna), "-i", str(mirna_hits), "-o", str(mirna_dir / f"{ctx.sample_name}.mirna_like.counts.tsv")],
        ctx.report_logs_dir / f"{ctx.sample_name}.nonhost.mirna.counts.log",
        dry_run,
    )

    rfam_input = rfam_dir / f"{ctx.sample_name}.rfam_input.fa"
    rfam_hits = rfam_dir / f"{ctx.sample_name}.rfam.hits.txt"
    unknown_fa = unknown_dir / f"{ctx.sample_name}.unknown.fa"
    if dry_run:
        print(f"[dry-run] concatenate {mirna_un} and {longer_in} into {rfam_input}")
    else:
        write_text(rfam_input, mirna_un.read_text(encoding="utf-8") + longer_in.read_text(encoding="utf-8"))
    run_command_stdout_to_file(
        [
            get_tool(ctx.config, "bowtie", "bowtie"),
            "-p",
            str(ctx.threads),
            "-v",
            "1",
            "-k",
            "1",
            "--best",
            "-f",
            "--un",
            str(unknown_fa),
            str(cfg["rfam_class_index"]),
            str(rfam_input),
        ],
        rfam_hits,
        ctx.report_logs_dir / f"{ctx.sample_name}.nonhost.rfam.log",
        dry_run,
    )
    run_command(
        [
            get_tool(ctx.config, "python", "python"),
            str(helper_rfam),
            "-i",
            str(rfam_hits),
            "-s",
            str(cfg["rfamseq_tsv"]),
            "--class-out",
            str(rfam_dir / f"{ctx.sample_name}.rfam_class.counts.tsv"),
            "--taxon-out",
            str(rfam_dir / f"{ctx.sample_name}.rfam_taxon.counts.tsv"),
            "--class-taxon-out",
            str(rfam_dir / f"{ctx.sample_name}.rfam_class_taxon.counts.tsv"),
            "--family-out",
            str(rfam_dir / f"{ctx.sample_name}.rfam_family.counts.tsv"),
            "--accession-out",
            str(rfam_dir / f"{ctx.sample_name}.rfam_accession.counts.tsv"),
        ],
        ctx.report_logs_dir / f"{ctx.sample_name}.nonhost.rfam.counts.log",
        dry_run,
    )
    if not dry_run:
        rrna_reads = count_fastq_reads(rrna_fastq_gz)
        collapsed_sequences = sum(1 for line in collapsed_fa.read_text(encoding="utf-8").splitlines() if line.startswith(">"))
        unknown_sequences = sum(1 for line in unknown_fa.read_text(encoding="utf-8").splitlines() if line.startswith(">"))
        trna_total = sum_count_table(trna_dir / f"{ctx.sample_name}.trna.counts.tsv")
        mirna_total = sum_count_table(mirna_dir / f"{ctx.sample_name}.mirna_like.counts.tsv")
        rfam_total = sum_count_table(rfam_dir / f"{ctx.sample_name}.rfam_class.counts.tsv")
        write_text(
            ctx.nonhost_dir / f"{ctx.sample_name}.nonhuman.summary.tsv",
            "sample\tcategory\tcount\n"
            f"{ctx.sample_name}\trrna_reads\t{rrna_reads}\n"
            f"{ctx.sample_name}\tnon_rrna_collapsed_sequences\t{collapsed_sequences}\n"
            f"{ctx.sample_name}\ttrna_like_reads\t{format_normalized(trna_total)}\n"
            f"{ctx.sample_name}\tmirna_like_reads\t{format_normalized(mirna_total)}\n"
            f"{ctx.sample_name}\trfam_classified_reads\t{format_normalized(rfam_total)}\n"
            f"{ctx.sample_name}\tunknown_sequences\t{unknown_sequences}\n",
        )
        mark_workflow_complete(ctx, "nonhost_specific", False)


def run_assembly_based(ctx: PipelineContext, dry_run: bool) -> None:
    validate_exists(ctx.assembly_fasta, "assembly FASTA")
    validate_exists(ctx.assembly_gff, "assembly GFF")
    cfg = get_workflow_config(ctx.config, "assembly_based")
    min_len = int(cfg.get("min_length", 15))
    max_len = int(cfg.get("max_length", 51))
    tmp_dir = ctx.assembly_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    helper_add_xc = ctx.helper_scripts_dir / "add_xc_fast.py"
    helper_cmc = ctx.helper_scripts_dir / "cmc_py3_modern_final.py"
    helper_xw = ctx.helper_scripts_dir / "xw_sam_to_bed_external_sort.py"
    helper_clean = ctx.helper_scripts_dir / "clean_bed_two_pass.py"
    helper_len = ctx.helper_scripts_dir / "read_length_dist.py"
    helper_parse = ctx.helper_scripts_dir / "parse_annotation_counts_fast.py"
    helper_deseq = ctx.helper_scripts_dir / "prepareDESeq2_py3.py"
    helper_db = ctx.helper_scripts_dir / "make_db_matrices_from_cds.py"
    for helper in [helper_add_xc, helper_cmc, helper_xw, helper_clean, helper_len, helper_parse, helper_deseq, helper_db]:
        validate_exists(helper, helper.name)

    assembly_gz = ctx.assembly_dir / f"{ctx.sample_name}.assembly.fa.gz"
    annotation_gff = ctx.assembly_dir / f"{ctx.sample_name}.annotation.gff"
    annotation_bed = ctx.assembly_dir / f"{ctx.sample_name}.annotation.bed"
    if ctx.assembly_fasta is not None:
        if dry_run:
            print(f"[dry-run] stage assembly {ctx.assembly_fasta} -> {assembly_gz}")
        else:
            if ctx.assembly_fasta.suffix == ".gz":
                shutil.copy2(ctx.assembly_fasta, assembly_gz)
            else:
                with ctx.assembly_fasta.open("rb") as in_handle, gzip.open(assembly_gz, "wb") as out_handle:
                    shutil.copyfileobj(in_handle, out_handle)
    if ctx.assembly_gff is not None:
        stage_plain_copy(ctx.assembly_gff, annotation_gff, dry_run)
    if dry_run:
        print(f"[dry-run] gff2bed < {annotation_gff} > {annotation_bed}")
    else:
        with (ctx.report_logs_dir / f"{ctx.sample_name}.gff2bed.log").open("w", encoding="utf-8") as log_handle:
            proc = subprocess.run(
                [get_tool(ctx.config, "gff2bed", "gff2bed")],
                input=annotation_gff.read_text(encoding="utf-8"),
                text=True,
                capture_output=True,
                check=True,
            )
            log_handle.write(proc.stderr)
        write_text(annotation_bed, proc.stdout)

    uniq_fa = ctx.assembly_dir / f"{ctx.sample_name}.nonhost.uniq.fa"
    uniq_min_fa = ctx.assembly_dir / f"{ctx.sample_name}.nonhost.uniq.min{min_len}.fa"
    uniq_min_fastq = ctx.assembly_dir / f"{ctx.sample_name}.nonhost.uniq.min{min_len}.fastq"
    mapping_sam = ctx.assembly_dir / f"{ctx.sample_name}.mapping.sam"
    mapped_sam = ctx.assembly_dir / f"{ctx.sample_name}.mapped.sam"
    mapped_bam = ctx.assembly_dir / f"{ctx.sample_name}.mapped.bam"
    mapped_xs = ctx.assembly_dir / f"{ctx.sample_name}.mapped.XS.sam"
    mapped_xw = ctx.assembly_dir / f"{ctx.sample_name}.mapped.XW.sam"
    mapped_bed = ctx.assembly_dir / f"{ctx.sample_name}.mapped.XW.bed"
    mapped_annot = ctx.assembly_dir / f"{ctx.sample_name}.mapped.XW.bed.mapped"
    mapped_clean = ctx.assembly_dir / f"{ctx.sample_name}.mapped.XW.bed.mapped.cleaned"
    unmapped_sam = ctx.assembly_dir / f"{ctx.sample_name}.unmapped.sam"

    nonhost_plain = ctx.assembly_dir / f"{ctx.sample_name}.nonhost.fastq"
    if dry_run:
        print(f"[dry-run] gunzip -c {ctx.nonhost_fastq} > {nonhost_plain}")
    else:
        with gzip.open(ctx.nonhost_fastq, "rt", encoding="utf-8") as in_handle, nonhost_plain.open("w", encoding="utf-8") as out_handle:
            shutil.copyfileobj(in_handle, out_handle)

    run_command([get_tool(ctx.config, "fastx_collapser", "fastx_collapser"), "-v", "-i", str(nonhost_plain), "-o", str(uniq_fa)], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.fastx.log", dry_run)
    run_command([get_tool(ctx.config, "seqkit", "seqkit"), "seq", "-g", "-m", str(min_len), str(uniq_fa)], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.seqkit.filter.log", dry_run)
    if dry_run:
        print(f"[dry-run] seqkit seq -g -m {min_len} {uniq_fa} > {uniq_min_fa}")
    else:
        proc = subprocess.run([get_tool(ctx.config, "seqkit", "seqkit"), "seq", "-g", "-m", str(min_len), str(uniq_fa)], capture_output=True, text=True, check=True)
        write_text(uniq_min_fa, proc.stdout)
    run_command([get_tool(ctx.config, "seqtk", "seqtk"), "seq", "-F", "I", str(uniq_min_fa)], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.seqtk.log", dry_run)
    if dry_run:
        print(f"[dry-run] seqtk seq -F I {uniq_min_fa} > {uniq_min_fastq}")
    else:
        proc = subprocess.run([get_tool(ctx.config, "seqtk", "seqtk"), "seq", "-F", "I", str(uniq_min_fa)], capture_output=True, text=True, check=True)
        write_text(uniq_min_fastq, proc.stdout)

    run_command([get_tool(ctx.config, "bwa", "bwa"), "index", str(assembly_gz)], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.bwaindex.log", dry_run)
    run_command([get_tool(ctx.config, "srnaMapper", "srnaMapper"), "-t", str(ctx.threads), "-r", str(uniq_min_fastq), "-g", str(assembly_gz), "-o", str(mapping_sam)], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.srnaMapper.log", dry_run)
    run_command([get_tool(ctx.config, "samtools", "samtools"), "view", "-@", str(ctx.threads), "-f4", str(mapping_sam)], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.unmapped.log", dry_run)
    if dry_run:
        print(f"[dry-run] samtools view -@ {ctx.threads} -f4 {mapping_sam} > {unmapped_sam}")
    else:
        proc = subprocess.run([get_tool(ctx.config, "samtools", "samtools"), "view", "-@", str(ctx.threads), "-f4", str(mapping_sam)], capture_output=True, text=True, check=True)
        write_text(unmapped_sam, proc.stdout)
    run_command([get_tool(ctx.config, "samtools", "samtools"), "view", "-@", str(ctx.threads), "-F4", "-h", "-S", str(mapping_sam)], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.mappedsam.log", dry_run)
    if dry_run:
        print(f"[dry-run] samtools view -@ {ctx.threads} -F4 -h -S {mapping_sam} > {mapped_sam}")
    else:
        proc = subprocess.run([get_tool(ctx.config, "samtools", "samtools"), "view", "-@", str(ctx.threads), "-F4", "-h", "-S", str(mapping_sam)], capture_output=True, text=True, check=True)
        write_text(mapped_sam, proc.stdout)
    run_command([get_tool(ctx.config, "python", "python"), str(helper_add_xc), str(mapped_sam), str(mapped_xs)], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.addxc.log", dry_run)
    run_command([get_tool(ctx.config, "python", "python"), str(helper_cmc), "-i", str(mapped_xs), "-o", str(mapped_xw), "-n", "1", "-e", "3"], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.cmc.log", dry_run)
    run_command([get_tool(ctx.config, "python", "python"), str(helper_xw), "-i", str(mapped_xw), "-o", str(mapped_bed), "--chunk-records", "3000000", "--tmpdir", str(tmp_dir)], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.xw.log", dry_run)
    run_command([get_tool(ctx.config, "bedtools", "bedtools"), "intersect", "-a", str(mapped_bed), "-b", str(annotation_bed), "-wao"], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.intersect.log", dry_run)
    if dry_run:
        print(f"[dry-run] bedtools intersect -a {mapped_bed} -b {annotation_bed} -wao > {mapped_annot}")
    else:
        proc = subprocess.run([get_tool(ctx.config, "bedtools", "bedtools"), "intersect", "-a", str(mapped_bed), "-b", str(annotation_bed), "-wao"], capture_output=True, text=True, check=True)
        write_text(mapped_annot, proc.stdout)
    run_command([get_tool(ctx.config, "python", "python"), str(helper_clean), str(mapped_annot), str(mapped_clean)], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.clean.log", dry_run)
    run_command([get_tool(ctx.config, "python", "python"), str(helper_len), "--max-len", str(max_len), "--bed", str(mapped_clean), "--bed-use-weighted", "-o", str(ctx.assembly_dir / f"{ctx.sample_name}.ReadLengthDistribution.weighted.txt")], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.length.log", dry_run)
    run_command([get_tool(ctx.config, "python", "python"), str(helper_parse), str(mapped_clean), "--raw-out", str(ctx.assembly_dir / f"{ctx.sample_name}.RawCounts.txt"), "--weighted-out", str(ctx.assembly_dir / f"{ctx.sample_name}.WeightedCounts.txt")], ctx.report_logs_dir / f"{ctx.sample_name}.assembly.parse.log", dry_run)
    if not dry_run:
        mark_workflow_complete(ctx, "assembly_based", False)


def run_novel_ncrna(ctx: PipelineContext, dry_run: bool) -> None:
    validate_exists(ctx.assembly_fasta, "assembly FASTA")
    validate_exists(ctx.assembly_gff, "assembly GFF")
    cfg = get_workflow_config(ctx.config, "novel_ncrna")
    min_len = int(cfg.get("min_length", 18))
    max_len = int(cfg.get("max_length", 250))
    min_reads = int(cfg.get("min_reads", 10))
    merge_dist = int(cfg.get("merge_distance", 10))
    mismatches = int(cfg.get("mismatches", 1))
    run_aragorn = bool(cfg.get("run_aragorn", True))
    run_trnascan = bool(cfg.get("run_trnascan", False))
    run_cmscan = bool(cfg.get("run_cmscan", True))
    run_rnafold = bool(cfg.get("run_rnafold", True))
    rfam_cm = cfg.get("rfam_cm")

    ref_gz = ctx.novel_dir / f"{ctx.sample_name}.assembly.fa.gz"
    gff = ctx.novel_dir / f"{ctx.sample_name}.annotation.gff"
    stage_gzip_copy(ctx.assembly_fasta, ref_gz, dry_run)
    stage_plain_copy(ctx.assembly_gff, gff, dry_run)
    ref = ctx.novel_dir / f"{ctx.sample_name}.assembly.fa"
    if dry_run:
        print(f"[dry-run] gunzip -k -f {ref_gz}")
    else:
        run_command(["gunzip", "-k", "-f", str(ref_gz)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.refgunzip.log", False)

    sam_path = ctx.novel_dir / f"{ctx.sample_name}.sam"
    bam_path = ctx.novel_dir / f"{ctx.sample_name}.sorted.bam"
    plus_bed = ctx.novel_dir / f"{ctx.sample_name}.plus.bed"
    minus_bed = ctx.novel_dir / f"{ctx.sample_name}.minus.bed"
    clusters = ctx.novel_dir / f"{ctx.sample_name}.clusters.bed"
    clusters_len = ctx.novel_dir / f"{ctx.sample_name}.clusters.lenfiltered.bed"
    counts_bed = ctx.novel_dir / f"{ctx.sample_name}.clusters.counts.bed"
    supported = ctx.novel_dir / f"{ctx.sample_name}.clusters.supported.bed"
    annot_bed = ctx.novel_dir / f"{ctx.sample_name}.annot.bed"
    cds_bed = ctx.novel_dir / f"{ctx.sample_name}.cds.bed"
    intergenic = ctx.novel_dir / f"{ctx.sample_name}.candidates.intergenic.bed"
    antisense_full = ctx.novel_dir / f"{ctx.sample_name}.candidates.antisense_to_CDS.full.bed"
    antisense = ctx.novel_dir / f"{ctx.sample_name}.candidates.antisense_to_CDS.bed"

    run_command([get_tool(ctx.config, "bowtie-build", "bowtie-build"), str(ref_gz), str(ref)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.bowtiebuild.log", dry_run)
    run_command([get_tool(ctx.config, "bowtie", "bowtie"), "-v", str(mismatches), "-a", "--best", "--strata", "-S", "-p", str(ctx.threads), str(ref), str(ctx.nonhost_fastq)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.bowtie.log", dry_run)
    if dry_run:
        print(f"[dry-run] bowtie ... > {sam_path}")
    else:
        proc = subprocess.run([get_tool(ctx.config, "bowtie", "bowtie"), "-v", str(mismatches), "-a", "--best", "--strata", "-S", "-p", str(ctx.threads), str(ref), str(ctx.nonhost_fastq)], capture_output=True, text=True, check=True)
        write_text(sam_path, proc.stdout)
    run_command([get_tool(ctx.config, "samtools", "samtools"), "view", "-@", str(ctx.threads), "-bS", str(sam_path)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.samtools.log", dry_run)
    if dry_run:
        print(f"[dry-run] samtools view -@ {ctx.threads} -bS {sam_path} | samtools sort -@ {ctx.threads} -o {bam_path}")
    else:
        first = subprocess.Popen([get_tool(ctx.config, "samtools", "samtools"), "view", "-@", str(ctx.threads), "-bS", str(sam_path)], stdout=subprocess.PIPE)
        second = subprocess.run([get_tool(ctx.config, "samtools", "samtools"), "sort", "-@", str(ctx.threads), "-o", str(bam_path)], stdin=first.stdout, check=True)
        if first.stdout is not None:
            first.stdout.close()
        first.wait()
        _ = second
    run_command([get_tool(ctx.config, "samtools", "samtools"), "index", str(bam_path)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.index.log", dry_run)

    if dry_run:
        print(f"[dry-run] strand-aware BED extraction from {bam_path} into {plus_bed} and {minus_bed}")
    else:
        proc = subprocess.run([get_tool(ctx.config, "samtools", "samtools"), "view", "-@", str(ctx.threads), "-F", "4", str(bam_path)], capture_output=True, text=True, check=True)
        plus_lines: list[str] = []
        minus_lines: list[str] = []
        for line in proc.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 11:
                continue
            flag = int(parts[1])
            start = int(parts[3]) - 1
            end = start + len(parts[9])
            strand = "-" if flag & 16 else "+"
            record = f"{parts[2]}\t{start}\t{end}\t{parts[0]}\t0\t{strand}\n"
            if strand == "+":
                plus_lines.append(record)
            else:
                minus_lines.append(record)
        write_text(plus_bed, "".join(plus_lines))
        write_text(minus_bed, "".join(minus_lines))
        sort_bed_in_place(plus_bed, False)
        sort_bed_in_place(minus_bed, False)

    plus_clusters = ctx.novel_dir / f"{ctx.sample_name}.plus.clusters.bed"
    minus_clusters = ctx.novel_dir / f"{ctx.sample_name}.minus.clusters.bed"
    run_command([get_tool(ctx.config, "bedtools", "bedtools"), "merge", "-s", "-d", str(merge_dist), "-i", str(plus_bed), "-c", "6", "-o", "distinct"], ctx.report_logs_dir / f"{ctx.sample_name}.novel.merge.plus.log", dry_run)
    run_command([get_tool(ctx.config, "bedtools", "bedtools"), "merge", "-s", "-d", str(merge_dist), "-i", str(minus_bed), "-c", "6", "-o", "distinct"], ctx.report_logs_dir / f"{ctx.sample_name}.novel.merge.minus.log", dry_run)
    if dry_run:
        print(f"[dry-run] merge plus/minus clusters into {clusters}")
    else:
        proc_plus = subprocess.run([get_tool(ctx.config, "bedtools", "bedtools"), "merge", "-s", "-d", str(merge_dist), "-i", str(plus_bed), "-c", "6", "-o", "distinct"], capture_output=True, text=True, check=True)
        proc_minus = subprocess.run([get_tool(ctx.config, "bedtools", "bedtools"), "merge", "-s", "-d", str(merge_dist), "-i", str(minus_bed), "-c", "6", "-o", "distinct"], capture_output=True, text=True, check=True)
        lines = []
        counter = 0
        for chunk in [proc_plus.stdout, proc_minus.stdout]:
            for line in chunk.splitlines():
                parts = line.split("\t")
                if len(parts) >= 4:
                    counter += 1
                    lines.append(f"{parts[0]}\t{parts[1]}\t{parts[2]}\tcluster_{counter}\t0\t{parts[3]}\n")
        write_text(clusters, "".join(sorted(lines)))
        kept = []
        for line in lines:
            parts = line.rstrip("\n").split("\t")
            if (int(parts[2]) - int(parts[1])) >= min_len and (int(parts[2]) - int(parts[1])) <= max_len:
                kept.append(line)
        write_text(clusters_len, "".join(sorted(kept)))

    if dry_run:
        print(f"[dry-run] strand-specific read support for {clusters_len} into {counts_bed}")
    else:
        plus_len = ctx.novel_dir / f"{ctx.sample_name}.clusters.lenfiltered.plus.bed"
        minus_len = ctx.novel_dir / f"{ctx.sample_name}.clusters.lenfiltered.minus.bed"
        plus_counts = ctx.novel_dir / f"{ctx.sample_name}.clusters.counts.plus.bed"
        minus_counts = ctx.novel_dir / f"{ctx.sample_name}.clusters.counts.minus.bed"
        with clusters_len.open("r", encoding="utf-8") as handle:
            plus_lines = []
            minus_lines = []
            for line in handle:
                if line.rstrip("\n").endswith("\t+"):
                    plus_lines.append(line)
                elif line.rstrip("\n").endswith("\t-"):
                    minus_lines.append(line)
        write_text(plus_len, "".join(plus_lines))
        write_text(minus_len, "".join(minus_lines))
        sort_bed_in_place(plus_len, False)
        sort_bed_in_place(minus_len, False)
        if plus_lines:
            proc = subprocess.run([get_tool(ctx.config, "bedtools", "bedtools"), "intersect", "-sorted", "-c", "-a", str(plus_len), "-b", str(plus_bed)], capture_output=True, text=True, check=True)
            write_text(plus_counts, proc.stdout)
        else:
            write_text(plus_counts, "")
        if minus_lines:
            proc = subprocess.run([get_tool(ctx.config, "bedtools", "bedtools"), "intersect", "-sorted", "-c", "-a", str(minus_len), "-b", str(minus_bed)], capture_output=True, text=True, check=True)
            write_text(minus_counts, proc.stdout)
        else:
            write_text(minus_counts, "")
        write_text(counts_bed, plus_counts.read_text(encoding="utf-8") + minus_counts.read_text(encoding="utf-8"))
        supported_lines = []
        for line in counts_bed.read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            if len(parts) >= 7 and int(parts[6]) >= min_reads:
                supported_lines.append(line + "\n")
        write_text(supported, "".join(supported_lines))
        ann_lines = []
        for line in gff.read_text(encoding="utf-8").splitlines():
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 8 and parts[2] not in {"region", "source"}:
                start = max(int(parts[3]) - 1, 0)
                end = int(parts[4])
                if start > end:
                    start, end = end, start
                ann_lines.append(f"{parts[0]}\t{start}\t{end}\t{parts[2]}\t.\t{parts[6]}\n")
        write_text(annot_bed, "".join(sorted(ann_lines)))
        cds_lines = [line for line in ann_lines if "\tCDS\t" in line]
        write_text(cds_bed, "".join(sorted(cds_lines)))

    if run_aragorn:
        run_command([get_tool(ctx.config, "aragorn", "aragorn"), "-gcstd", "-w", "-o", str(ctx.novel_dir / f"{ctx.sample_name}.aragorn.txt"), str(ref)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.aragorn.log", dry_run)
    if run_trnascan:
        run_command([get_tool(ctx.config, "tRNAscan-SE", "tRNAscan-SE"), "-o", str(ctx.novel_dir / f"{ctx.sample_name}.trnascan.out"), str(ref)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.trnascan.log", dry_run)

    run_command([get_tool(ctx.config, "bedtools", "bedtools"), "intersect", "-s", "-v", "-a", str(supported), "-b", str(annot_bed)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.intergenic.log", dry_run)
    run_command([get_tool(ctx.config, "bedtools", "bedtools"), "intersect", "-S", "-wa", "-wb", "-a", str(supported), "-b", str(cds_bed)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.antisense.log", dry_run)
    if dry_run:
        print(f"[dry-run] bedtools intersect outputs into {intergenic} and {antisense_full}")
    else:
        proc = subprocess.run([get_tool(ctx.config, "bedtools", "bedtools"), "intersect", "-s", "-v", "-a", str(supported), "-b", str(annot_bed)], capture_output=True, text=True, check=True)
        write_text(intergenic, proc.stdout)
        proc = subprocess.run([get_tool(ctx.config, "bedtools", "bedtools"), "intersect", "-S", "-wa", "-wb", "-a", str(supported), "-b", str(cds_bed)], capture_output=True, text=True, check=True)
        write_text(antisense_full, proc.stdout)
        unique = []
        seen = set()
        for line in antisense_full.read_text(encoding="utf-8").splitlines():
            key = tuple(line.split("\t")[:6])
            if key not in seen:
                seen.add(key)
                unique.append("\t".join(key) + "\n")
        write_text(antisense, "".join(sorted(unique)))

    intergenic_fa = ctx.novel_dir / f"{ctx.sample_name}.candidates.intergenic.fa"
    antisense_fa = ctx.novel_dir / f"{ctx.sample_name}.candidates.antisense_to_CDS.fa"
    all_fa = ctx.novel_dir / f"{ctx.sample_name}.candidates.all.fa"
    run_command([get_tool(ctx.config, "bedtools", "bedtools"), "getfasta", "-s", "-fi", str(ref), "-bed", str(intergenic), "-fo", str(intergenic_fa)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.getfasta.intergenic.log", dry_run)
    run_command([get_tool(ctx.config, "bedtools", "bedtools"), "getfasta", "-s", "-fi", str(ref), "-bed", str(antisense), "-fo", str(antisense_fa)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.getfasta.antisense.log", dry_run)
    if dry_run:
        print(f"[dry-run] concatenate {intergenic_fa} and {antisense_fa} into {all_fa}")
    else:
        write_text(all_fa, intergenic_fa.read_text(encoding="utf-8") + antisense_fa.read_text(encoding="utf-8"))

    no_rfam_fa = ctx.novel_dir / f"{ctx.sample_name}.candidates.no_rfam.fa"
    if run_cmscan and rfam_cm:
        tbl = ctx.novel_dir / f"{ctx.sample_name}.candidates.rfam.tbl"
        run_command([get_tool(ctx.config, "cmscan", "cmscan"), "--cpu", str(ctx.threads), "--rfam", "--cut_ga", "--notrunc", "--tblout", str(tbl), str(rfam_cm), str(all_fa)], ctx.report_logs_dir / f"{ctx.sample_name}.novel.cmscan.log", dry_run)
        if dry_run:
            print(f"[dry-run] filter rfam hits from {all_fa} into {no_rfam_fa}")
        else:
            hit_ids = []
            for line in tbl.read_text(encoding="utf-8").splitlines():
                if line and not line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 3:
                        hit_ids.append(parts[2])
            if hit_ids:
                proc = subprocess.run([get_tool(ctx.config, "seqkit", "seqkit"), "grep", "-v", "-f", "/dev/stdin", str(all_fa)], input="\n".join(hit_ids) + "\n", capture_output=True, text=True, check=True)
                write_text(no_rfam_fa, proc.stdout)
            else:
                write_text(no_rfam_fa, all_fa.read_text(encoding="utf-8"))
    else:
        if dry_run:
            print(f"[dry-run] copy {all_fa} -> {no_rfam_fa}")
        else:
            write_text(no_rfam_fa, all_fa.read_text(encoding="utf-8"))

    if run_rnafold:
        run_command([get_tool(ctx.config, "RNAfold", "RNAfold"), "--noPS"], ctx.report_logs_dir / f"{ctx.sample_name}.novel.rnafold.log", dry_run)
        if dry_run:
            print(f"[dry-run] RNAfold --noPS < {no_rfam_fa} > {ctx.novel_dir / f'{ctx.sample_name}.candidates.no_rfam.RNAfold.txt'}")
        else:
            with no_rfam_fa.open("r", encoding="utf-8") as handle:
                proc = subprocess.run([get_tool(ctx.config, "RNAfold", "RNAfold"), "--noPS"], stdin=handle, capture_output=True, text=True, check=True)
            write_text(ctx.novel_dir / f"{ctx.sample_name}.candidates.no_rfam.RNAfold.txt", proc.stdout)

    if not dry_run:
        alignment_count = run_capture([get_tool(ctx.config, "samtools", "samtools"), "view", "-c", str(bam_path)], False)
        intergenic_count = sum(1 for line in intergenic.read_text(encoding="utf-8").splitlines() if line.strip())
        antisense_count = sum(1 for line in antisense.read_text(encoding="utf-8").splitlines() if line.strip())
        write_text(
            ctx.novel_dir / f"{ctx.sample_name}.novel.summary.tsv",
            "sample\tcategory\treads\n"
            f"{ctx.sample_name}\ttotal_alignments\t{alignment_count}\n"
            f"{ctx.sample_name}\tnovel_intergenic_candidates\t{intergenic_count}\n"
            f"{ctx.sample_name}\tnovel_antisense_candidates\t{antisense_count}\n",
        )
        mark_workflow_complete(ctx, "novel_ncrna", False)


def run_summary(ctx: PipelineContext, dry_run: bool) -> None:
    summary_files = [
        ctx.preprocess_stats,
        ctx.host_dir / f"{ctx.sample_name}.human.summary.tsv",
        ctx.nonhost_dir / f"{ctx.sample_name}.nonhuman.summary.tsv",
        ctx.novel_dir / f"{ctx.sample_name}.novel.summary.tsv",
    ]
    out_path = ctx.summary_dir / f"{ctx.sample_name}.pipeline.summary.tsv"
    if dry_run:
        print(f"[dry-run] aggregate summaries into {out_path}")
        return
    lines = ["source\tsample\tcategory\tvalue"]
    for summary in summary_files:
        if not summary.exists():
            continue
        with summary.open("r", encoding="utf-8") as handle:
            header = next(handle, None)
            _ = header
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    lines.append(f"{summary.name}\t{parts[0]}\t{parts[1]}\t{parts[2]}")
    write_text(out_path, "\n".join(lines) + "\n")
    mark_workflow_complete(ctx, "summary", False)


def run_multiqc(ctx: PipelineContext, dry_run: bool) -> None:
    if not ctx.config["settings"].get("run_multiqc", True):
        return
    run_command(
        [
            get_tool(ctx.config, "multiqc", "multiqc"),
            str(ctx.sample_dir),
            "-c",
            str(ctx.root_dir / "scripts" / "srna_imp_multiqc_config.yaml"),
            "-o",
            str(ctx.sample_dir / "multiqc"),
        ],
        ctx.report_logs_dir / f"{ctx.sample_name}.multiqc.log",
        dry_run,
    )
    if not dry_run:
        mark_workflow_complete(ctx, "multiqc", False)


def selected_workflows(requested: str) -> list[str]:
    if requested == "all":
        return WORKFLOWS
    workflows = [requested]
    if requested in DOWNSTREAM_NEEDS_PREPROCESS:
        workflows = ["preprocessing", requested]
    if requested == "summary":
        workflows = ["summary"]
    return workflows


def run_sample(config: dict[str, Any], workflow: str, dry_run: bool) -> None:
    validate_config(config, workflow)
    ctx = PipelineContext(config)
    ctx.ensure_dirs()

    for workflow_name in selected_workflows(workflow):
        if not dry_run and workflow_complete(ctx, workflow_name):
            print(f"[resume][{ctx.sample_name}] skipping completed workflow: {workflow_name}")
            continue
        if workflow_name == "preprocessing":
            run_preprocessing(ctx, dry_run)
        elif workflow_name == "host_specific":
            run_host_specific(ctx, dry_run)
        elif workflow_name == "nonhost_specific":
            run_nonhost_specific(ctx, dry_run)
        elif workflow_name == "assembly_based":
            run_assembly_based(ctx, dry_run)
        elif workflow_name == "novel_ncrna":
            run_novel_ncrna(ctx, dry_run)
        elif workflow_name == "summary":
            run_summary(ctx, dry_run)

    if workflow == "all":
        if not dry_run and workflow_complete(ctx, "multiqc"):
            print(f"[resume][{ctx.sample_name}] skipping completed workflow: multiqc")
        else:
            run_multiqc(ctx, dry_run)


def main() -> int:
    args = parse_args()
    base_config = load_config(Path(args.config).resolve())
    configs = expand_configs(base_config, args.sample)
    total = len(configs)
    exit_code = 0

    for index, config in enumerate(configs, start=1):
        sample_name = config.get("sample", {}).get("name", f"sample_{index}")
        prefix = f"[{index}/{total}][{sample_name}]" if total > 1 else f"[{sample_name}]"
        print(f"{prefix} starting workflow '{args.workflow}'")
        try:
            run_sample(config, args.workflow, args.dry_run)
        except Exception as exc:
            print(f"{prefix} failed: {exc}", file=sys.stderr)
            exit_code = 1
            if not args.keep_going:
                break
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
