#!/usr/bin/env python3
"""
Production-Ready Standalone Pipeline

This script provides a production-grade alternative to the Nextflow workflow with:
- Comprehensive logging (per-model failure tracking)
- Caching and resume capability
- HTML execution reports
- Error categorization and retry logic
- Configurable via CLI
- Real-time progress tracking

Replicates all 15 stages from end_to_end_with_validation_multibatch.nf
(includes ipSAE interface quality metrics calculation)
"""

import argparse
import csv
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any

import duckdb
import orjson

# Add repo to sys.path for direct imports
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_DIR = _SCRIPT_DIR.parent
if str(_REPO_DIR) not in sys.path:
    sys.path.insert(0, str(_REPO_DIR))


# ============================================================================
# PipelineLogger - Multi-file logging with per-model tracking
# ============================================================================

class PipelineLogger:
    """
    Comprehensive logging system that tracks:
    - General pipeline execution (logs/pipeline.log)
    - Per-model failures (logs/failed_models.txt, logs/errors.log)
    - Per-model successes (logs/successful_models.txt)
    - Per-stage failure details (logs/stage_XX_failures.json)
    """

    def __init__(self, output_dir: Path, log_level: str = "INFO"):
        self.output_dir = output_dir
        self.log_dir = output_dir / "logs"
        self.log_dir.mkdir(exist_ok=True, parents=True)

        # Log files
        self.pipeline_log = self.log_dir / "pipeline.log"
        self.errors_log = self.log_dir / "errors.log"
        self.failed_models_txt = self.log_dir / "failed_models.txt"
        self.successful_models_txt = self.log_dir / "successful_models.txt"

        # In-memory tracking
        self.failed_models: Dict[str, List[Dict[str, str]]] = {}  # stage -> [model failures]
        self.successful_models: Dict[str, List[str]] = {}  # stage -> [model ids]

        # Setup Python logging
        self.logger = logging.getLogger("production_pipeline")
        self.logger.setLevel(getattr(logging, log_level.upper()))

        # Console handler (INFO level)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S')
        console_handler.setFormatter(console_formatter)

        # File handler (DEBUG level)
        file_handler = logging.FileHandler(self.pipeline_log, mode='w')
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(
            '[%(asctime)s] [%(levelname)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_formatter)

        # Error file handler
        error_handler = logging.FileHandler(self.errors_log, mode='w')
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(file_formatter)

        self.logger.addHandler(console_handler)
        self.logger.addHandler(file_handler)
        self.logger.addHandler(error_handler)

        # Initialize empty files
        self.failed_models_txt.write_text("")
        self.successful_models_txt.write_text("")

    def info(self, message: str, indent: int = 0):
        """Log info message with optional indentation"""
        prefix = "  " * indent
        self.logger.info(f"{prefix}{message}")

    def debug(self, message: str, indent: int = 0):
        """Log debug message"""
        prefix = "  " * indent
        self.logger.debug(f"{prefix}{message}")

    def error(self, message: str, indent: int = 0):
        """Log error message"""
        prefix = "  " * indent
        self.logger.error(f"{prefix}{message}")

    def warning(self, message: str, indent: int = 0):
        """Log warning message"""
        prefix = "  " * indent
        self.logger.warning(f"{prefix}{message}")

    def log_model_failure(self, stage: str, model_id: str, error_message: str, error_category: str = "UNKNOWN"):
        """
        Track a model failure

        Args:
            stage: Stage identifier (e.g., "stage_10")
            model_id: Model identifier
            error_message: Error description
            error_category: Error category (FILE_NOT_FOUND, VALIDATION_ERROR, etc.)
        """
        if stage not in self.failed_models:
            self.failed_models[stage] = []

        self.failed_models[stage].append({
            "model_id": model_id,
            "error_message": error_message,
            "error_category": error_category,
            "timestamp": datetime.now(timezone.utc).isoformat()
        })

        # Append to failed_models.txt
        with open(self.failed_models_txt, 'a') as f:
            f.write(f"{model_id}\n")

        # Log to errors.log
        self.error(f"[{stage}] Model {model_id} failed: {error_message}")

    def log_model_success(self, stage: str, model_id: str):
        """Track a model success"""
        if stage not in self.successful_models:
            self.successful_models[stage] = []

        self.successful_models[stage].append(model_id)

        # Append to successful_models.txt
        with open(self.successful_models_txt, 'a') as f:
            f.write(f"{model_id}\n")

    def save_stage_failures(self, stage: str):
        """Save per-stage failure details to JSON"""
        if stage in self.failed_models and self.failed_models[stage]:
            stage_failure_file = self.log_dir / f"{stage}_failures.json"
            with open(stage_failure_file, 'w') as f:
                json.dump(self.failed_models[stage], f, indent=2)
            self.debug(f"Saved {len(self.failed_models[stage])} failures to {stage_failure_file}")

    def get_failure_summary(self) -> Dict[str, int]:
        """Get summary of failures by stage"""
        return {stage: len(failures) for stage, failures in self.failed_models.items()}

    def get_success_summary(self) -> Dict[str, int]:
        """Get summary of successes by stage"""
        return {stage: len(models) for stage, models in self.successful_models.items()}


# ============================================================================
# PipelineCache - Resume capability with state tracking
# ============================================================================

class PipelineCache:
    """
    Pipeline state caching for resume capability

    Tracks completed stages and allows skipping them on resume.
    Cache is invalidated if configuration changes.
    """

    def __init__(self, output_dir: Path, config_hash: str):
        self.output_dir = output_dir
        self.cache_file = output_dir / ".pipeline_cache.json"
        self.config_hash = config_hash
        self.cache_data = self._load_cache()

    def _load_cache(self) -> Dict[str, Any]:
        """Load cache from file"""
        if not self.cache_file.exists():
            return self._create_empty_cache()

        try:
            with open(self.cache_file) as f:
                data = json.load(f)

            # Validate config hash
            if data.get("config_hash") != self.config_hash:
                print(f"[CACHE] Configuration changed - invalidating cache")
                return self._create_empty_cache()

            return data
        except Exception as e:
            print(f"[CACHE] Failed to load cache: {e}")
            return self._create_empty_cache()

    def _create_empty_cache(self) -> Dict[str, Any]:
        """Create empty cache structure"""
        return {
            "run_id": datetime.now().strftime("%Y%m%d_%H%M%S"),
            "config_hash": self.config_hash,
            "stages": {},
            "failed_models": {}
        }

    def is_stage_completed(self, stage_name: str) -> bool:
        """Check if a stage is marked as completed"""
        stage_data = self.cache_data["stages"].get(stage_name, {})
        return stage_data.get("completed", False)

    def mark_stage_completed(
        self,
        stage_name: str,
        duration: float,
        models_succeeded: int = 0,
        models_failed: int = 0,
        success: bool = True
    ):
        """Mark a stage as completed"""
        self.cache_data["stages"][stage_name] = {
            "completed": success,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "duration": duration,
            "models_succeeded": models_succeeded,
            "models_failed": models_failed
        }
        self._save_cache()

    def mark_model_failed(self, stage_name: str, model_id: str):
        """Track a failed model"""
        if stage_name not in self.cache_data["failed_models"]:
            self.cache_data["failed_models"][stage_name] = []
        self.cache_data["failed_models"][stage_name].append(model_id)

    def get_failed_models(self, stage_name: str) -> List[str]:
        """Get list of failed models for a stage"""
        return self.cache_data["failed_models"].get(stage_name, [])

    def _save_cache(self):
        """Save cache to file"""
        try:
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache_data, f, indent=2)
        except Exception as e:
            print(f"[CACHE] Failed to save cache: {e}")

    def clear_cache(self):
        """Clear the cache file"""
        if self.cache_file.exists():
            self.cache_file.unlink()
        self.cache_data = self._create_empty_cache()


# ============================================================================
# ErrorTracker - Categorize and track errors
# ============================================================================

class ErrorTracker:
    """
    Analyze and categorize errors from command output

    Error categories:
    - FILE_NOT_FOUND: Missing input files
    - VALIDATION_ERROR: Schema or data validation failures
    - MEMORY_ERROR: Out of memory
    - TIMEOUT: Process timeout
    - DEPENDENCY_ERROR: Missing packages
    - UNKNOWN: Uncategorized
    """

    ERROR_PATTERNS = {
        "FILE_NOT_FOUND": [
            r"No such file or directory",
            r"FileNotFoundError",
            r"does not exist",
            r"cannot find"
        ],
        "VALIDATION_ERROR": [
            r"ValidationError",
            r"Invalid.*format",
            r"Schema validation failed",
            r"Unexpected.*value"
        ],
        "MEMORY_ERROR": [
            r"MemoryError",
            r"Out of memory",
            r"Cannot allocate memory",
            r"OOM"
        ],
        "TIMEOUT": [
            r"TimeoutError",
            r"timed out",
            r"Timeout exceeded"
        ],
        "DEPENDENCY_ERROR": [
            r"ModuleNotFoundError",
            r"ImportError",
            r"No module named",
            r"command not found"
        ]
    }

    SUGGESTIONS = {
        "FILE_NOT_FOUND": "Check that input files exist and paths are correct",
        "VALIDATION_ERROR": "Verify input data format matches expected schema",
        "MEMORY_ERROR": "Increase memory allocation or reduce batch size/workers",
        "TIMEOUT": "Increase timeout limit or check for hung processes",
        "DEPENDENCY_ERROR": "Install missing dependencies (check requirements.txt)",
        "UNKNOWN": "Review error logs for details"
    }

    def __init__(self):
        self.error_counts: Dict[str, int] = {}

    def categorize_error(self, error_text: str) -> str:
        """
        Categorize an error based on its text

        Returns:
            Error category string
        """
        for category, patterns in self.ERROR_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, error_text, re.IGNORECASE):
                    self.error_counts[category] = self.error_counts.get(category, 0) + 1
                    return category

        self.error_counts["UNKNOWN"] = self.error_counts.get("UNKNOWN", 0) + 1
        return "UNKNOWN"

    def get_suggestion(self, error_category: str) -> str:
        """Get actionable suggestion for error category"""
        return self.SUGGESTIONS.get(error_category, self.SUGGESTIONS["UNKNOWN"])

    def get_error_summary(self) -> Dict[str, int]:
        """Get summary of errors by category"""
        return self.error_counts.copy()


# ============================================================================
# Configuration - CLI arguments and defaults
# ============================================================================

@dataclass
class Config:
    """Pipeline configuration"""
    # Paths
    repo_dir: Path
    input_dir: Path
    output_dir: Path
    mapping_file: Optional[Path] = None
    chain_mapping: Optional[Path] = None
    dataset_config: Optional[Path] = None
    provider_json: Optional[Path] = None
    uniprot_db: Optional[Path] = None
    modelcif_template: Optional[Path] = None

    # Execution
    python_cmd: List[str] = None  # Will default to ["python"]
    workers: int = 30
    model_version: str = "v1"
    batch_size: int = 1000

    # ipSAE parameters
    pae_cutoff: float = 10.0
    dist_cutoff: float = 15.0

    # Clash/Interface analysis parameters
    clash_cutoff: float = 0.4
    interface_cutoff: float = 8.0
    analysis_batch_size: int = 128
    clash_device: str = "cuda"  # "cpu", "cuda", or "auto" (resolved in __post_init__)
    analyses: Set[str] = None  # None → all three; resolved in __post_init__

    # Enrichment
    enrichment_metrics: List[str] = None  # Will default in __post_init__
    cif_qa_metrics: Optional[str] = None  # Comma-separated or "auto"

    # Features
    resume: bool = False
    no_cache: bool = False
    retry: int = 2
    skip_stages: List[str] = None
    dry_run: bool = False
    dssp_algorithm: str = "psea"  # "psea", "pydssp", or "tmalign"
    parallel_stages: bool = False  # wave-based parallel execution of independent stages

    # Dataset metadata
    tool_used: str = "ColabFold v1.6.0 / AlphaFold-Multimer"

    # Heterodimer mode
    heterodimers: bool = False
    provider_id: str = "NVDA"
    provider_name: str = "NVIDIA"

    # Reporting
    log_level: str = "INFO"
    report_file: Optional[Path] = None

    def __post_init__(self):
        """Convert strings to Paths and resolve to absolute paths"""
        for field in ['repo_dir', 'input_dir', 'output_dir', 'mapping_file',
                      'chain_mapping', 'dataset_config', 'provider_json',
                      'uniprot_db', 'modelcif_template']:
            value = getattr(self, field)
            if value is not None:
                if not isinstance(value, Path):
                    value = Path(value)
                # Resolve to absolute path to avoid cwd issues
                setattr(self, field, value.resolve())

        if self.skip_stages is None:
            self.skip_stages = []

        if self.analyses is None:
            self.analyses = {"interface", "backbone_clashes", "heavy_atom_clashes"}

        if self.clash_device == "auto":
            try:
                import torch
                self.clash_device = "cuda" if torch.cuda.is_available() else "cpu"
            except ImportError:
                self.clash_device = "cpu"

        if self.enrichment_metrics is None:
            self.enrichment_metrics = [
                "pae_cutoff", "dist_cutoff", "iptm_af",
                "ipsae", "ipsae_d0chn", "ipsae_d0dom",
                "iptm_d0chn", "pDockQ2", "LIS",
                "n0res", "n0dom", "d0res", "d0dom",
                "nres1", "nres2", "dist_nres1", "dist_nres2",
                "pDockQ", "n0chn", "N_clash_backbone", "N_clash_heavyAtom",
            ]

        # Set default python command if not provided
        if self.python_cmd is None:
            self.python_cmd = ["python"]

        # Select modelcif template based on tool_used if not explicitly provided
        if self.modelcif_template is None:
            is_colabfold = (self.tool_used or "").lower().startswith("colabfold")
            if is_colabfold:
                template_name = "colabfold_modelcif_metadata.json"
            else:
                template_name = "openfold_modelcif_metadata.json"
            self.modelcif_template = self.repo_dir / "uniprot/templates" / template_name

    def get_hash(self) -> str:
        """Generate hash of configuration for cache validation"""
        # Hash relevant fields that affect output
        config_str = (
            f"{self.repo_dir}:{self.input_dir}:{self.mapping_file}:"
            f"{self.chain_mapping}:{self.workers}:{self.model_version}:"
            f"{self.batch_size}:{self.heterodimers}:{self.provider_id}"
        )
        return hashlib.sha256(config_str.encode()).hexdigest()[:16]


def parse_arguments() -> Config:
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description="Production-ready standalone pipeline for AFDB integration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic run with required arguments
  %(prog)s --output-dir /path/to/output \\
    --input-dir /path/to/input \\
    --mapping-file /path/to/mapping.tsv \\
    --chain-mapping /path/to/manifest.csv \\
    --dataset-config /path/to/config.json \\
    --provider-json /path/to/provider.json \\
    --uniprot-db /path/to/uniprot.duckdb

  # Resume from previous run
  %(prog)s --output-dir /path/to/output [...] --resume

  # Use uv for Python environment
  %(prog)s --output-dir /path/to/output [...] --python-cmd "uv run python"

  # Skip specific stages
  %(prog)s --output-dir /path/to/output [...] --skip-stages stage_12,stage_14

  # Dry run (show what would be executed)
  %(prog)s --output-dir /path/to/output [...] --dry-run

  # Heterodimer mode (pre-built manifest + DuckDB)
  %(prog)s --output-dir /path/to/output \\
    --input-dir /path/to/raw_colabfold \\
    --heterodimers \\
    --chain-mapping /path/to/manifest.csv \\
    --uniprot-db /path/to/uniprot.duckdb
        """
    )

    # Required arguments
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for pipeline results"
    )

    # Path arguments - repo-dir auto-detected, others required
    parser.add_argument(
        "--repo-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="Repository root directory (default: auto-detected from script location)"
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Input directory containing PDB and JSON files "
             "(raw ColabFold outputs when --heterodimers is set)"
    )

    parser.add_argument(
        "--mapping-file",
        type=Path,
        default=None,
        help="Model ID mapping file (TSV format). "
             "Required in homodimer mode. Not needed in heterodimer mode "
             "(model IDs are derived from --chain-mapping)."
    )

    parser.add_argument(
        "--chain-mapping",
        type=Path,
        default=None,
        help="Chain-to-accession mapping CSV (model_entity_id, chain_id, uniprot_ac). "
             "Required in homodimer mode. In heterodimer mode, pass with "
             "--uniprot-db for production (skips API fetching)."
    )

    parser.add_argument(
        "--dataset-config",
        type=Path,
        default=None,
        help="Dataset configuration JSON file. "
             "Auto-generated if not provided (stoichiometry and modelCreatedDate "
             "are inferred from input files)."
    )

    parser.add_argument(
        "--provider-json",
        type=Path,
        default=None,
        help="Provider configuration JSON file. "
             "Required in homodimer mode. Auto-generated in heterodimer mode "
             "if not provided."
    )

    parser.add_argument(
        "--uniprot-db",
        type=Path,
        default=None,
        help="UniProt DuckDB database file. "
             "Required in homodimer mode. In heterodimer mode, pass with "
             "--chain-mapping for production (skips API fetching)."
    )

    # Dataset metadata
    parser.add_argument(
        "--tool-used",
        type=str,
        choices=[
            "ColabFold v1.6.0 / AlphaFold-Multimer",
            "OpenFold-TRT / AlphaFold-Multimer",
        ],
        default="ColabFold v1.6.0 / AlphaFold-Multimer",
        help="Prediction tool recorded in dataset_config.json (default: %(default)s)."
    )

    # Heterodimer mode
    parser.add_argument(
        "--heterodimers",
        action="store_true",
        help="Enable heterodimer mode. --input-dir may contain raw ColabFold "
             "outputs (long suffixes are detected automatically). "
             "Requires --chain-mapping + --uniprot-db. Model IDs are derived "
             "from --chain-mapping; --dataset-config and --provider-json are "
             "auto-generated if not provided."
    )

    parser.add_argument(
        "--provider-id",
        type=str,
        default="NVDA",
        help="(Heterodimer mode) Dataset provider ID (default: %(default)s)."
    )

    parser.add_argument(
        "--provider-name",
        type=str,
        default="NVIDIA",
        help="(Heterodimer mode) Display name for the provider (default: %(default)s)."
    )

    # Execution parameters
    parser.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count() or 8,
        help="Number of parallel workers (default: all available CPUs)"
    )

    parser.add_argument(
        "--model-version",
        default="v1",
        help="Model version string (default: %(default)s)"
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=1000,
        help="Batch size for metadata combining (default: %(default)s)"
    )

    parser.add_argument(
        "--python-cmd",
        default="python",
        help="Python command to use (default: %(default)s). Use 'uv run python' for uv environments."
    )

    # Features
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from previous run (skip completed stages)"
    )

    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable caching (force re-run all stages)"
    )

    parser.add_argument(
        "--retry",
        type=int,
        default=2,
        help="Number of retries for failed stages (default: %(default)s)"
    )

    parser.add_argument(
        "--skip-stages",
        type=str,
        default="",
        help="Comma-separated list of stages to skip (e.g., stage_12,stage_14)"
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be executed without running"
    )

    parser.add_argument(
        "--dssp-algorithm",
        choices=["psea", "pydssp", "tmalign"],
        default="pydssp",
        help="Algorithm for secondary structure: 'psea' (geometry), 'pydssp' (H-bond), or 'tmalign' (CA-CA distance). Default: %(default)s"
    )

    parser.add_argument(
        "--parallel-stages",
        action="store_true",
        default=False,
        help="Run independent pipeline stages in parallel waves (splits worker budget among concurrent stages)"
    )

    # Reporting
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: %(default)s)"
    )

    parser.add_argument(
        "--report-file",
        type=Path,
        help="Path to HTML report file (default: output_dir/report.html)"
    )

    # ipSAE parameters
    parser.add_argument(
        "--pae-cutoff",
        type=float,
        default=10.0,
        help="PAE cutoff for ipSAE calculation (default: %(default)s)"
    )

    parser.add_argument(
        "--dist-cutoff",
        type=float,
        default=15.0,
        help="Distance cutoff for ipSAE calculation (default: %(default)s)"
    )

    # Clash/Interface analysis parameters
    parser.add_argument(
        "--clash-cutoff",
        type=float,
        default=0.4,
        help="VDW overlap fraction for clash detection (default: %(default)s, ~0.4 A overlap)"
    )

    parser.add_argument(
        "--interface-cutoff",
        type=float,
        default=8.0,
        help="CA-CA distance cutoff for interface residues (default: %(default)s Angstrom)"
    )

    parser.add_argument(
        "--analysis-batch-size",
        type=int,
        default=4,
        help="Batch size for clash/interface analysis (default: %(default)s)"
    )

    parser.add_argument(
        "--interface-clash-analysis",
        nargs="+",
        choices=["interface", "backbone_clashes", "heavy_atom_clashes"],
        default=["interface", "backbone_clashes", "heavy_atom_clashes"],
        help="Which analyses to run in stage 13 (default: all three)"
    )

    parser.add_argument(
        "--clash-device",
        choices=["cpu", "cuda", "auto"],
        default="cuda",
        help="Device for stage 13 clash/interface analysis (default: cuda). Use 'cpu' to force CPU.",
    )

    parser.add_argument(
        "--enrichment-metrics",
        nargs="+",
        default=[
            "pae_cutoff", "dist_cutoff", "iptm_af",
            "ipsae", "ipsae_d0chn", "ipsae_d0dom",
            "iptm_d0chn", "pDockQ2", "LIS",
            "n0res", "n0dom", "d0res", "d0dom",
            "nres1", "nres2", "dist_nres1", "dist_nres2",
            "pDockQ", "n0chn", "N_clash_backbone", "N_clash_heavyAtom",
        ],
        help="iPSAE/clash metric base names to include in model/chain metadata JSONs "
             "(default: all known metrics)"
    )

    parser.add_argument(
        "--cif-qa-metrics",
        type=str,
        default="auto",
        help="Comma-separated CIF QA metric short names to embed in mmCIF files "
             "(e.g. 'ipsae_AB,iptm_af,N_clash_backbone'). Use 'auto' for all (default). "
             "Requires model JSONs enriched with complexPredictionAccuracy_* keys."
    )

    args = parser.parse_args()

    # ---- Argument validation for heterodimer vs homodimer mode ----
    if args.heterodimers:
        # Heterodimer mode: require chain-mapping + uniprot-db.
        # mapping-file, dataset-config, and provider-json are optional
        # (derived / auto-generated).
        heterodimer_required = {
            "--chain-mapping": args.chain_mapping,
            "--uniprot-db": args.uniprot_db,
        }
        missing = [name for name, val in heterodimer_required.items() if val is None]
        if missing:
            parser.error(
                f"--heterodimers requires: {', '.join(missing)}"
            )
    else:
        # Homodimer mode: all config paths are required
        homodimer_required = {
            "--mapping-file": args.mapping_file,
            "--chain-mapping": args.chain_mapping,
            "--provider-json": args.provider_json,
            "--uniprot-db": args.uniprot_db,
        }
        missing = [name for name, val in homodimer_required.items() if val is None]
        if missing:
            parser.error(
                f"The following arguments are required (unless --heterodimers is set): "
                f"{', '.join(missing)}"
            )

    # Process skip_stages
    skip_stages = [s.strip() for s in args.skip_stages.split(',') if s.strip()]

    # Process python_cmd (split by space to support "uv run python")
    python_cmd = args.python_cmd.split()

    # Create Config object
    config = Config(
        repo_dir=args.repo_dir,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        mapping_file=args.mapping_file,
        chain_mapping=args.chain_mapping,
        dataset_config=args.dataset_config,
        provider_json=args.provider_json,
        uniprot_db=args.uniprot_db,
        python_cmd=python_cmd,
        workers=args.workers,
        model_version=args.model_version,
        batch_size=args.batch_size,
        pae_cutoff=args.pae_cutoff,
        dist_cutoff=args.dist_cutoff,
        clash_cutoff=args.clash_cutoff,
        interface_cutoff=args.interface_cutoff,
        analysis_batch_size=args.analysis_batch_size,
        clash_device=args.clash_device,
        analyses=set(args.interface_clash_analysis),
        enrichment_metrics=args.enrichment_metrics,
        cif_qa_metrics=args.cif_qa_metrics,
        resume=args.resume,
        no_cache=args.no_cache,
        retry=args.retry,
        skip_stages=skip_stages,
        dry_run=args.dry_run,
        dssp_algorithm=args.dssp_algorithm,
        parallel_stages=args.parallel_stages,
        tool_used=args.tool_used,
        heterodimers=args.heterodimers,
        provider_id=args.provider_id,
        provider_name=args.provider_name,
        log_level=args.log_level,
        report_file=args.report_file or (args.output_dir / "report.html")
    )

    return config


# ============================================================================
# Helper Functions
# ============================================================================

def run_command(
    cmd: List[str],
    stage: str,
    logger: PipelineLogger,
    error_tracker: ErrorTracker,
    cwd: Optional[Path] = None
) -> Tuple[float, bool, str, str]:
    """
    Run a command and return results

    Returns:
        (duration, success, stdout, stderr)
    """
    start = time.time()
    logger.debug(f"Running command: {' '.join(str(c) for c in cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            cwd=str(cwd) if cwd else None
        )
        duration = time.time() - start
        logger.info(f"✓ Completed in {duration:.2f}s", indent=1)

        # Show first few lines of output
        if result.stdout:
            for line in result.stdout.strip().split('\n')[:3]:
                if line.strip():
                    logger.debug(line.strip(), indent=2)

        return duration, True, result.stdout, ""

    except subprocess.CalledProcessError as e:
        duration = time.time() - start
        logger.error(f"✗ Failed after {duration:.2f}s", indent=1)

        # Categorize error
        error_category = error_tracker.categorize_error(e.stderr)
        suggestion = error_tracker.get_suggestion(error_category)

        logger.error(f"Error category: {error_category}", indent=2)
        logger.error(f"Suggestion: {suggestion}", indent=2)

        if e.stderr:
            lines = e.stderr.strip().split('\n')
            for line in lines[-30:]:
                if line.strip():
                    logger.error(line.strip(), indent=2)

        return duration, False, "", e.stderr


def load_model_ids(mapping_file: Path, logger: PipelineLogger) -> List[str]:
    """Load model IDs from mapping file (homodimer mode)."""
    logger.info(f"Loading model IDs from {mapping_file}")
    with open(mapping_file) as f:
        model_ids = [line.strip().split('\t')[0] for line in f if line.strip()]
    logger.info(f"Loaded {len(model_ids)} model IDs")
    return model_ids


def load_model_ids_from_manifest(chain_mapping: Path, logger: PipelineLogger) -> List[str]:
    """Derive unique model IDs from chain-mapping CSV using DuckDB.

    Uses DuckDB's native C++ CSV reader for speed at scale (~1-2s for 30M rows
    vs ~30s for csv.DictReader).
    """
    logger.info(f"Loading model IDs from manifest {chain_mapping}")
    query = """
        SELECT DISTINCT model_entity_id
        FROM read_csv_auto(?)
        ORDER BY model_entity_id
    """
    rows = duckdb.execute(query, [str(chain_mapping)]).fetchall()
    model_ids = [row[0] for row in rows]
    logger.info(f"Loaded {len(model_ids)} unique model IDs from manifest")
    return model_ids


# ---------------------------------------------------------------------------
# Enrichment helpers -- iPSAE + clash data → model/chain metadata JSONs
# ---------------------------------------------------------------------------

# Regex: strip _{X}{Y} chain-pair suffix (e.g. "_AB", "_BA") from column name
_CHAIN_PAIR_SUFFIX_RE = re.compile(r"_([A-Z])([A-Z])$")


def _model_id_from_ipsae_pdb_path(pdb_path: str) -> Optional[str]:
    """Extract model ID from ipSAE CSV pdb_path (supports homodimer and heterodimer naming).

    Paths are like .../AF_0000000067800110_AF_0000000067800429-model_v1.pdb or
    .../AF-1234567890123456-model_v1.pdb. We derive ID by stripping the trailing
    -model_v1.pdb (or -model_v1) suffix from the filename.
    """
    if not pdb_path or not pdb_path.strip():
        return None
    name = Path(pdb_path).stem
    for suffix in ("-model_v1", "-model-v1"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return None


def _parse_ipsae_csv(csv_path: Path) -> Dict[str, Dict[str, Any]]:
    """Read iPSAE summary CSV and return {model_id: {col: value, ...}}.

    Non-metric columns (``pdb_path``, ``processing_time_ms``) are skipped.
    The model ID is extracted from the ``pdb_path`` filename (supports both
    homodimer and heterodimer naming, e.g. AF_..._AF_... or AF-...).
    """
    result: Dict[str, Dict[str, Any]] = {}
    if not csv_path.exists():
        return result

    skip_cols = {"pdb_path", "processing_time_ms"}

    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pdb_path = row.get("pdb_path", "")
            model_id = _model_id_from_ipsae_pdb_path(pdb_path)
            if not model_id:
                continue
            metrics: Dict[str, Any] = {}
            for col, val in row.items():
                if col in skip_cols:
                    continue
                try:
                    metrics[col] = float(val)
                except (ValueError, TypeError):
                    metrics[col] = val
            result[model_id] = metrics
    return result


def _model_id_from_clash_filename(stem: str) -> Optional[str]:
    """Extract model ID from clash file stem (e.g. AF_..._AF_...-model_v1_clashes)."""
    suffix = "-model_v1_clashes"
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return None


def _parse_clash_jsons(clash_dir: Path) -> Dict[str, Dict[str, int]]:
    """Read ``*_clashes.json`` files and return {model_id: {"N_clash_backbone": n, "N_clash_heavyAtom": n}}.

    Looks for site labels ``backbone_clashes`` and ``heavy_atom_clashes``
    inside each JSON's ``sites`` array. Model ID is derived from filename (supports heterodimer naming).
    """
    result: Dict[str, Dict[str, int]] = {}
    if not clash_dir.exists():
        return result

    label_map = {
        "backbone_clashes": "N_clash_backbone",
        "heavy_atom_clashes": "N_clash_heavyAtom",
        # Backwards compat with old runs that still say "side_chain_clashes"
        "side_chain_clashes": "N_clash_heavyAtom",
    }

    for json_file in clash_dir.glob("*_clashes.json"):
        model_id = _model_id_from_clash_filename(json_file.stem)
        if not model_id:
            continue
        try:
            data = json.loads(json_file.read_bytes())
        except Exception:
            continue
        counts: Dict[str, int] = {}
        for site in data.get("sites", []):
            key = label_map.get(site.get("label"))
            if key:
                counts[key] = site.get("additional_site_annotations", {}).get("n_clashes", 0)
        # Ensure both keys are always present (default 0 if a site was omitted)
        counts.setdefault("N_clash_backbone", 0)
        counts.setdefault("N_clash_heavyAtom", 0)
        result[model_id] = counts
    return result


def _parse_interface_jsons(clash_dir: Path) -> Dict[str, int]:
    """Read ``*_interface.json`` files and return {model_id: number_of_interactions}."""
    result: Dict[str, int] = {}
    if not clash_dir.exists():
        return result

    interface_suffix = "-model_v1_interface"
    for json_file in clash_dir.glob("*_interface.json"):
        stem = json_file.stem
        if not stem.endswith(interface_suffix):
            continue
        model_id = stem[: -len(interface_suffix)]
        try:
            data = json.loads(json_file.read_bytes())
        except Exception:
            continue
        total = sum(
            len(site.get("additional_site_annotations", {}).get("interactions", []))
            for site in data.get("sites", [])
        )
        result[model_id] = total
    return result


def _base_name(col: str) -> str:
    """Strip ``_{X}{Y}`` chain-pair suffix to get the base metric name.

    ``ipsae_AB`` → ``ipsae``, ``pDockQ2_BA`` → ``pDockQ2``,
    ``iptm_af`` → ``iptm_af`` (no two-char suffix → returned as-is).
    """
    m = _CHAIN_PAIR_SUFFIX_RE.search(col)
    return col[: m.start()] if m else col


def _ipsae_json_key(col: str) -> str:
    """Map an iPSAE CSV column name to its JSON key.

    Metrics that are NOT iPSAE-specific get ``complexPredictionAccuracy_{col}``,
    iPSAE-derived metrics get ``complexPredictionAccuracy_ipsae_{col}``.
    """
    non_ipsae_bases = {"iptm_af", "pDockQ2", "LIS", "pDockQ", "N_clash_backbone", "N_clash_heavyAtom"}
    if _base_name(col) in non_ipsae_bases:
        return f"complexPredictionAccuracy_{col}"
    # Columns already starting with "ipsae_" don't need the prefix doubled
    if col.startswith("ipsae_"):
        return f"complexPredictionAccuracy_{col}"
    return f"complexPredictionAccuracy_ipsae_{col}"


def _build_model_enrichment(
    ipsae_row: Dict[str, Any],
    clash_row: Dict[str, int],
    metrics_filter: List[str],
) -> Dict[str, Any]:
    """Build ``complexPredictionAccuracy_*`` dict for one model."""
    filter_set = set(metrics_filter)
    out: Dict[str, Any] = {}

    for col, val in ipsae_row.items():
        if _base_name(col) in filter_set:
            out[_ipsae_json_key(col)] = val

    for key, val in clash_row.items():
        if key in filter_set:
            out[f"complexPredictionAccuracy_{key}"] = val

    return out


def _build_chain_enrichment(
    ipsae_row: Dict[str, Any],
    chain_id: str,
    metrics_filter: List[str],
) -> Dict[str, Any]:
    """Build ``complexPredictionAccuracy_*`` dict for one chain.

    General metrics (no chain-pair suffix) go to all chains.
    Chain-specific metrics (e.g. ``ipsae_AB``) go only to the chain whose
    ID matches the **first** letter of the suffix (A for ``_AB``).
    """
    filter_set = set(metrics_filter)
    out: Dict[str, Any] = {}

    for col, val in ipsae_row.items():
        base = _base_name(col)
        if base not in filter_set:
            continue

        m = _CHAIN_PAIR_SUFFIX_RE.search(col)
        if m:
            # Chain-specific: include only if first letter matches chain_id
            if m.group(1) == chain_id:
                out[_ipsae_json_key(col)] = val
        else:
            # General metric: all chains get it
            out[_ipsae_json_key(col)] = val

    return out


# ---------------------------------------------------------------------------
# Input file discovery -- supports canonical and ColabFold naming conventions
# ---------------------------------------------------------------------------

CANONICAL_PDB_SUFFIX = "-model_v1.pdb"
CANONICAL_META_SUFFIX = "-meta_v1.json"
COLABFOLD_PDB_SUFFIX = (
    ".merged_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb"
)
COLABFOLD_SCORE_SUFFIX = (
    ".merged_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json"
)
COLABFOLD_PDB_SUFFIX_ALT = (
    "_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb"
)
COLABFOLD_SCORE_SUFFIX_ALT = (
    "_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json"
)

PDB_SUFFIXES = [CANONICAL_PDB_SUFFIX, COLABFOLD_PDB_SUFFIX, COLABFOLD_PDB_SUFFIX_ALT]
META_SUFFIXES = [CANONICAL_META_SUFFIX, COLABFOLD_SCORE_SUFFIX, COLABFOLD_SCORE_SUFFIX_ALT]

_AFDB_PREFIX = "AFDB_"

MAX_IO_WORKERS = 16


def _strip_afdb_prefix(model_id: str) -> str:
    """Strip the ``AFDB_`` prefix that ColabFold batch outputs prepend."""
    if model_id.startswith(_AFDB_PREFIX):
        return model_id[len(_AFDB_PREFIX):]
    return model_id


def _canonicalize_model_id(model_id: str) -> str:
    """Normalize AF_ to AF- for simple IDs. Compound heterodimer IDs (AF_X_AF_Y) are returned unchanged so they match the underscore form stored in batch_ids."""
    if model_id.startswith("AF_") and "_AF_" not in model_id:
        return "AF-" + model_id[3:]
    return model_id


def scan_input_dir(input_dir: Path) -> Dict[str, Tuple[Path, Path]]:
    """Single os.scandir() pass to build {model_id: (pdb_path, meta_path)} lookup.

    Tries each known suffix convention (canonical ``-model_v1.pdb`` first, then
    ColabFold long suffix).  One directory read instead of N * stat() calls per
    model -- critical on networked filesystems like Lustre.
    """
    pdbs: Dict[str, Path] = {}
    metas: Dict[str, Path] = {}

    for entry in os.scandir(input_dir):
        if not entry.is_file():
            continue
        name = entry.name
        for sfx in PDB_SUFFIXES:
            if name.endswith(sfx):
                raw_id = _strip_afdb_prefix(name[: -len(sfx)])
                pdbs[_canonicalize_model_id(raw_id)] = Path(entry.path)
                break
        for sfx in META_SUFFIXES:
            if name.endswith(sfx):
                raw_id = _strip_afdb_prefix(name[: -len(sfx)])
                metas[_canonicalize_model_id(raw_id)] = Path(entry.path)
                break

    matched = pdbs.keys() & metas.keys()
    return {mid: (pdbs[mid], metas[mid]) for mid in matched}


def _symlink_model(
    model_id: str,
    pdb_src: Path,
    meta_src: Path,
    staging_dir: Path,
) -> Optional[str]:
    """Symlink one model's files to staging with canonical names.

    Returns the model_id on failure, ``None`` on success.
    """
    try:
        pdb_dst = staging_dir / f"{model_id}{CANONICAL_PDB_SUFFIX}"
        meta_dst = staging_dir / f"{model_id}{CANONICAL_META_SUFFIX}"
        for dst in (pdb_dst, meta_dst):
            if dst.is_symlink() or dst.exists():
                dst.unlink()
        os.symlink(pdb_src, pdb_dst)
        os.symlink(meta_src, meta_dst)
        return None
    except Exception:
        return model_id


# ---------------------------------------------------------------------------
# Inline config generation (heterodimer mode)
# ---------------------------------------------------------------------------

def _infer_model_created_date(input_dir: Path) -> str:
    """Infer model creation date from the earliest PDB file modification time."""
    earliest_mtime = None
    for entry in os.scandir(input_dir):
        if not entry.is_file() or not entry.name.endswith(".pdb"):
            continue
        mtime = entry.stat().st_mtime
        if earliest_mtime is None or mtime < earliest_mtime:
            earliest_mtime = mtime
    if earliest_mtime is None:
        return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return datetime.fromtimestamp(earliest_mtime, tz=timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )


def generate_dataset_config(
    path: Path,
    provider_id: str,
    input_dir: Path,
    tool_used: str,
) -> None:
    """Write a dataset configuration JSON with dynamically inferred fields.

    modelCreatedDate is derived from the earliest PDB file mtime in input_dir.
    """
    model_created_date = _infer_model_created_date(input_dir)
    config = {
        "toolUsed": tool_used,
        "providerId": provider_id,
        "entityType": "protein",
        "isUniProt": True,
        "modelCreatedDate": model_created_date,
        "versionTag": "v1",
        "ordinal": 1,
        "latestVersion": 1,
        "allVersions": [1],
    }
    path.write_bytes(orjson.dumps(config, option=orjson.OPT_INDENT_2))


def generate_provider_json(
    path: Path, provider_id: str, provider_name: str
) -> None:
    """Write a provider configuration JSON for heterodimer mode."""
    provider = {
        "providerId": provider_id,
        "providerName": provider_name,
        "providerUrl": "https://alphafold.ebi.ac.uk",
        "copyrights": [
            f"Copyright 2026 {provider_name}. All rights reserved."
        ],
    }
    path.write_bytes(orjson.dumps(provider, option=orjson.OPT_INDENT_2))


def preflight_checks(config: "Config", logger: PipelineLogger) -> bool:
    """Run pre-flight checks before pipeline execution.

    Checks for required binaries and attempts to compile them if missing.

    Returns:
        True if all checks pass, False otherwise.
    """
    logger.info("Running pre-flight checks...")
    all_passed = True

    # Check ipSAE binary
    ipsae_binary = config.repo_dir / "afdb_integration_kit" / "ipsae" / "ipsae_cpp"
    ipsae_cpp_dir = ipsae_binary.parent

    if not ipsae_binary.exists():
        logger.warning("ipSAE binary not found, attempting to compile...")
        makefile = ipsae_cpp_dir / "Makefile"

        if not makefile.exists():
            logger.error("=" * 60)
            logger.error("PREFLIGHT CHECK FAILED: ipSAE Makefile missing")
            logger.error("=" * 60)
            logger.error(f"  Makefile expected at: {makefile}")
            logger.error("=" * 60)
            all_passed = False
        else:
            result = subprocess.run(
                ["make"],
                cwd=ipsae_cpp_dir,
                capture_output=True,
                text=True
            )
            if result.returncode != 0 or not ipsae_binary.exists():
                logger.error("=" * 60)
                logger.error("PREFLIGHT CHECK FAILED: ipSAE binary compilation failed")
                logger.error("=" * 60)
                logger.error(f"  Binary expected at: {ipsae_binary}")
                logger.error(f"  To compile manually, run:")
                logger.error(f"    cd {ipsae_cpp_dir} && make")
                if result.stderr:
                    logger.error(f"  Compiler error: {result.stderr[:200]}")
                logger.error("=" * 60)
                all_passed = False
            else:
                logger.info("  Successfully compiled ipSAE binary")
    else:
        logger.info("  ipSAE binary: OK")

    if all_passed:
        logger.info("Pre-flight checks passed")

    return all_passed


# ============================================================================
# Stage Functions
# ============================================================================

def stage_01_prepare_assets(
    model_ids: List[str],
    config: Config,
    logger: PipelineLogger
) -> Dict[str, Any]:
    """STAGE 1: BATCH_PREPARE_AF_ASSETS - Symlink input files to staging.

    Supports both canonical (``-model_v1.pdb``) and ColabFold long-suffix
    naming conventions.  Uses a single ``os.scandir()`` pass for discovery
    and ``ThreadPoolExecutor`` for parallel symlinking.
    """
    logger.info(f"STAGE 1: PREPARE_ASSETS ({len(model_ids)} models)")

    staging_dir = config.output_dir / "staging"
    staging_dir.mkdir(exist_ok=True, parents=True)

    start = time.time()

    # Single directory scan -- O(1) regardless of naming convention
    file_map = scan_input_dir(config.input_dir)
    logger.debug(f"Directory scan found {len(file_map)} matched pairs", indent=1)

    # Split into found / missing
    model_id_set = set(model_ids)
    found_ids = [mid for mid in model_ids if mid in file_map]
    missing_ids = [mid for mid in model_ids if mid not in file_map]

    # Parallel symlinking
    symlink_failures: List[str] = []
    with ThreadPoolExecutor(max_workers=MAX_IO_WORKERS) as executor:
        futures = {
            executor.submit(
                _symlink_model, mid, file_map[mid][0], file_map[mid][1], staging_dir
            ): mid
            for mid in found_ids
        }
        for future in as_completed(futures):
            failed_id = future.result()
            if failed_id is not None:
                symlink_failures.append(failed_id)

    # Batch logging -- single write instead of per-model file I/O
    succeeded = [mid for mid in found_ids if mid not in set(symlink_failures)]
    all_failed = missing_ids + symlink_failures

    if succeeded:
        logger.successful_models.setdefault("stage_01", []).extend(succeeded)
        with open(logger.successful_models_txt, "a") as f:
            f.write("\n".join(succeeded) + "\n")

    for mid in missing_ids:
        logger.log_model_failure("stage_01", mid, "Missing PDB or JSON file", "FILE_NOT_FOUND")
    for mid in symlink_failures:
        logger.log_model_failure("stage_01", mid, "Symlink creation failed", "UNKNOWN")

    duration = time.time() - start
    success_count = len(succeeded)
    logger.info(
        f"✓ Completed in {duration:.2f}s ({success_count}/{len(model_ids)} succeeded)",
        indent=1,
    )

    return {
        "duration": duration,
        "success": success_count > 0,
        "models_succeeded": success_count,
        "models_failed": len(all_failed),
    }


def stage_02_validate_assets(
    model_ids: List[str],
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 2: BATCH_VALIDATE_ASSETS - Validate PDB/JSON consistency"""
    logger.info(f"STAGE 2: VALIDATE_ASSETS ({len(model_ids)} models, {config.workers} workers)")

    # Write model IDs file
    ids_file = config.output_dir / "model_ids.txt"
    ids_file.write_text("\n".join(model_ids))

    staging_dir = config.output_dir / "staging"

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "uniprot/scripts/batch_validate_assets.py"),
        "--ids", str(ids_file),
        "--input-dir", str(staging_dir),
        "--output", str(config.output_dir / "validation_results.tsv"),
        "--workers", str(config.workers)
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_02", logger, error_tracker, config.repo_dir)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": len(model_ids) if success else 0,
        "models_failed": 0 if success else len(model_ids)
    }


def stage_03_convert_colabfold(
    model_ids: List[str],
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 3: BATCH_CONVERT_COLABFOLD - Convert ColabFold to AFDB format"""
    logger.info(f"STAGE 3: CONVERT_COLABFOLD ({len(model_ids)} models, {config.workers} workers)")

    ids_file = config.output_dir / "model_ids.txt"
    staging_dir = config.output_dir / "staging"
    scores_dir = config.output_dir / "scores"
    chain_manifest_dir = config.output_dir / "chain_manifests"
    model_manifest_dir = config.output_dir / "model_manifests"

    for d in [scores_dir, chain_manifest_dir, model_manifest_dir]:
        d.mkdir(exist_ok=True, parents=True)

    failed_ids_file = config.output_dir / "failed_models.tsv"

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "uniprot/scripts/batch_convert_colabfold.py"),
        "--model-ids-file", str(ids_file),
        "--input-dir", str(staging_dir),
        "--manifest", str(config.chain_mapping),
        "--duckdb", str(config.uniprot_db),
        "--output-dir", str(scores_dir),
        "--chain-manifest-dir", str(chain_manifest_dir),
        "--model-manifest-dir", str(model_manifest_dir),
        "--workers", str(config.workers),
        "--failed-ids-file", str(failed_ids_file),
        "--stage-name", "stage_03_convert_colabfold",
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_03", logger, error_tracker, config.repo_dir)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": len(model_ids) if success else 0,
        "models_failed": 0 if success else len(model_ids)
    }


def stage_04_merge_manifests(
    config: Config,
    logger: PipelineLogger
) -> Dict[str, Any]:
    """STAGE 4: MERGE_MANIFESTS - Merge chain and model manifests"""
    logger.info("STAGE 4: MERGE_MANIFESTS")
    start = time.time()

    chain_manifest_dir = config.output_dir / "chain_manifests"
    model_manifest_dir = config.output_dir / "model_manifests"
    merged_manifest_dir = config.output_dir / "merged_manifests"
    merged_manifest_dir.mkdir(exist_ok=True, parents=True)

    # Merge chain manifests
    chain_files = list(chain_manifest_dir.glob("*_afid_mapping.csv"))
    chain_output = merged_manifest_dir / "uniprot_afid_mapping.csv"

    with open(chain_output, 'w') as out:
        header_written = False
        for f in sorted(chain_files):
            with open(f) as inp:
                lines = inp.readlines()
                if not header_written:
                    out.writelines(lines)
                    header_written = True
                else:
                    out.writelines(lines[1:])  # Skip header

    # Merge model manifests
    model_files = list(model_manifest_dir.glob("*_model_metadata.csv"))
    model_output = merged_manifest_dir / "uniprot_model_metadata.csv"

    with open(model_output, 'w') as out:
        header_written = False
        for f in sorted(model_files):
            with open(f) as inp:
                lines = inp.readlines()
                if not header_written:
                    out.writelines(lines)
                    header_written = True
                else:
                    out.writelines(lines[1:])

    duration = time.time() - start
    logger.info(f"✓ Completed in {duration:.2f}s", indent=1)
    logger.debug(f"Chain manifest: {len(chain_files)} files merged", indent=2)
    logger.debug(f"Model manifest: {len(model_files)} files merged", indent=2)

    return {
        "duration": duration,
        "success": True,
        "models_succeeded": 0,
        "models_failed": 0
    }


def stage_05_export_model_metadata(
    model_ids: List[str],
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 5: BATCH_EXPORT_MODEL_METADATA"""
    logger.info(f"STAGE 5: EXPORT_MODEL_METADATA ({len(model_ids)} models, {config.workers} workers)")

    ids_file = config.output_dir / "model_ids.txt"
    chain_manifest = config.output_dir / "merged_manifests" / "uniprot_afid_mapping.csv"
    model_manifest = config.output_dir / "merged_manifests" / "uniprot_model_metadata.csv"
    model_json_dir = config.output_dir / "model_jsons"
    model_json_dir.mkdir(exist_ok=True, parents=True)

    failed_ids_file = config.output_dir / "failed_models.tsv"

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "uniprot/scripts/batch_export_metadata.py"),
        "--model-ids", str(ids_file),
        "--db", str(config.uniprot_db),
        "--config", str(config.dataset_config),
        "--mapping", str(chain_manifest),
        "--model-manifest", str(model_manifest),
        "--output-dir", str(model_json_dir),
        "--export-type", "model",
        "--workers", str(config.workers),
        "--failed-ids-file", str(failed_ids_file),
        "--stage-name", "stage_05_export_model_metadata",
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_05", logger, error_tracker, config.repo_dir)

    # ---- Enrich model JSONs with iPSAE + clash + interface metrics ----
    if success:
        clash_dir = config.output_dir / "clash_interface_analysis"
        ipsae_data = _parse_ipsae_csv(config.output_dir / "ipsae" / "ipsae_summary.csv")
        clash_data = _parse_clash_jsons(clash_dir)
        interface_data = _parse_interface_jsons(clash_dir)
        enriched = 0
        for json_file in model_json_dir.glob("*.json"):
            model_id = json_file.stem
            enrichment = _build_model_enrichment(
                ipsae_data.get(model_id, {}),
                clash_data.get(model_id, {}),
                config.enrichment_metrics,
            )
            n_interactions = interface_data.get(model_id)
            if enrichment or n_interactions is not None:
                data = json.loads(json_file.read_bytes())
                data.update(enrichment)
                if n_interactions is not None:
                    data["numberOfInteractions"] = n_interactions
                json_file.write_bytes(orjson.dumps(data, option=orjson.OPT_INDENT_2) + b"\n")
                enriched += 1
        logger.info(f"Enriched {enriched} model JSONs with iPSAE/clash/interface metrics", indent=1)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": len(model_ids) if success else 0,
        "models_failed": 0 if success else len(model_ids)
    }


def stage_06_export_chain_metadata(
    model_ids: List[str],
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 6: BATCH_EXPORT_CHAIN_METADATA"""
    logger.info(f"STAGE 6: EXPORT_CHAIN_METADATA ({len(model_ids)} models, {config.workers} workers)")

    ids_file = config.output_dir / "model_ids.txt"
    chain_manifest = config.output_dir / "merged_manifests" / "uniprot_afid_mapping.csv"
    model_manifest = config.output_dir / "merged_manifests" / "uniprot_model_metadata.csv"
    chain_json_dir = config.output_dir / "chain_jsons"
    chain_json_dir.mkdir(exist_ok=True, parents=True)

    failed_ids_file = config.output_dir / "failed_models.tsv"

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "uniprot/scripts/batch_export_metadata.py"),
        "--model-ids", str(ids_file),
        "--db", str(config.uniprot_db),
        "--config", str(config.dataset_config),
        "--mapping", str(chain_manifest),
        "--model-manifest", str(model_manifest),
        "--output-dir", str(chain_json_dir),
        "--export-type", "chain",
        "--workers", str(config.workers),
        "--failed-ids-file", str(failed_ids_file),
        "--stage-name", "stage_06_export_chain_metadata",
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_06", logger, error_tracker, config.repo_dir)

    # ---- Enrich chain JSONs with iPSAE metrics ----
    if success:
        ipsae_data = _parse_ipsae_csv(config.output_dir / "ipsae" / "ipsae_summary.csv")
        enriched = 0
        for json_file in chain_json_dir.glob("*.json"):
            model_id = json_file.stem
            ipsae_row = ipsae_data.get(model_id, {})
            if not ipsae_row:
                continue
            chain_records = json.loads(json_file.read_bytes())
            modified = False
            for record in chain_records:
                chain_id = record.get("uniqueId", "").rsplit("_", 1)[-1]
                enrichment = _build_chain_enrichment(
                    ipsae_row,
                    chain_id,
                    config.enrichment_metrics,
                )
                if enrichment:
                    record.update(enrichment)
                    modified = True
            if modified:
                json_file.write_bytes(orjson.dumps(chain_records, option=orjson.OPT_INDENT_2) + b"\n")
                enriched += 1
        logger.info(f"Enriched {enriched} chain JSONs with iPSAE metrics", indent=1)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": len(model_ids) if success else 0,
        "models_failed": 0 if success else len(model_ids)
    }


def stage_07_combine_model_metadata(
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 7: COMBINE_MODEL_METADATA"""
    logger.info("STAGE 7: COMBINE_MODEL_METADATA")

    model_json_dir = config.output_dir / "model_jsons"
    model_batch_dir = config.output_dir / "model_batches"
    model_batch_dir.mkdir(exist_ok=True, parents=True)

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "uniprot/scripts/combine_metadata.py"),
        "--input-dir", str(model_json_dir),
        "--output-dir", str(model_batch_dir),
        "--output-prefix", "AF-metadata",
        "--chunk-size", str(config.batch_size)
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_07", logger, error_tracker, config.repo_dir)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": 0,
        "models_failed": 0
    }


def stage_08_combine_chain_metadata(
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 8: COMBINE_CHAIN_METADATA"""
    logger.info("STAGE 8: COMBINE_CHAIN_METADATA")

    chain_json_dir = config.output_dir / "chain_jsons"
    chain_batch_dir = config.output_dir / "chain_batches"
    chain_batch_dir.mkdir(exist_ok=True, parents=True)

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "uniprot/scripts/combine_metadata.py"),
        "--input-dir", str(chain_json_dir),
        "--output-dir", str(chain_batch_dir),
        "--output-prefix", "AF-chain-metadata",
        "--chunk-size", str(config.batch_size)
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_08", logger, error_tracker, config.repo_dir)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": 0,
        "models_failed": 0
    }


def stage_09_export_modelcif_input(
    model_ids: List[str],
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 9: BATCH_EXPORT_MODELCIF_INPUT"""
    logger.info(f"STAGE 9: EXPORT_MODELCIF_INPUT ({len(model_ids)} models, {config.workers} workers)")

    ids_file = config.output_dir / "model_ids.txt"
    chain_manifest = config.output_dir / "merged_manifests" / "uniprot_afid_mapping.csv"
    modelcif_input_dir = config.output_dir / "modelcif_input"
    modelcif_input_dir.mkdir(exist_ok=True, parents=True)

    failed_ids_file = config.output_dir / "failed_models.tsv"

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "uniprot/scripts/batch_export_modelcif_input.py"),
        "--model-ids", str(ids_file),
        "--manifest", str(chain_manifest),
        "--db", str(config.uniprot_db),
        "--template", str(config.modelcif_template),
        "--output-dir", str(modelcif_input_dir),
        "--workers", str(config.workers),
        "--failed-ids-file", str(failed_ids_file),
        "--stage-name", "stage_09_export_modelcif_input",
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_09", logger, error_tracker, config.repo_dir)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": len(model_ids) if success else 0,
        "models_failed": 0 if success else len(model_ids)
    }


def stage_10_batch_modelcif_gen(
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 10: BATCH_RUN_MODELCIF_GEN"""
    logger.info(f"STAGE 10: BATCH_MODELCIF_GEN ({config.workers} workers)")

    staging_dir = config.output_dir / "staging"
    modelcif_input_dir = config.output_dir / "modelcif_input"
    modelcif_dir = config.output_dir / "modelcif"
    modelcif_dir.mkdir(exist_ok=True, parents=True)

    # Create symlinks for PDB and metadata files
    pdb_dir = config.output_dir / "pdb_staging"
    metadata_dir = config.output_dir / "metadata_staging"
    pdb_dir.mkdir(exist_ok=True, parents=True)
    metadata_dir.mkdir(exist_ok=True, parents=True)

    for pdb in staging_dir.glob("*-model_v1.pdb"):
        link_path = pdb_dir / pdb.name
        if link_path.exists():
            link_path.unlink()
        link_path.symlink_to(pdb)

    for json_file in modelcif_input_dir.glob("*.json"):
        link_path = metadata_dir / json_file.name
        if link_path.exists():
            link_path.unlink()
        link_path.symlink_to(json_file)

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "main.py"),
        "run-batch-modelcif-gen",
        "--input-dir", str(pdb_dir),
        "--metadata-dir", str(metadata_dir),
        "--output-dir", str(modelcif_dir),
        "--model-version", config.model_version,
        "--skip-validation",
        "--skip-alignment",
        "--workers", str(config.workers)
    ]

    # Optionally inject QA metrics from model JSONs into CIF files
    model_json_dir = config.output_dir / "model_jsons"
    if config.cif_qa_metrics and model_json_dir.exists():
        cmd.extend(["--model-json-dir", str(model_json_dir)])
        cmd.extend(["--cif-qa-metrics", config.cif_qa_metrics])

    duration, success, stdout, stderr = run_command(cmd, "stage_10", logger, error_tracker, config.repo_dir)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": len(list(modelcif_dir.glob("*.cif"))) if success else 0,
        "models_failed": 0
    }


def stage_11_batch_dssp(
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 11: BATCH_RUN_DSSP"""
    algo_names = {"psea": "P-SEA", "pydssp": "PyDSSP", "tmalign": "TM-align"}
    algo_name = algo_names.get(config.dssp_algorithm, config.dssp_algorithm)
    device = config.clash_device
    logger.info(f"STAGE 11: BATCH_DSSP ({config.workers} workers, {algo_name}, device={device})")

    modelcif_dir = config.output_dir / "modelcif"

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "main.py"),
        "batch-dssp",
        "--input-dir", str(modelcif_dir),
        "--output-dir", str(modelcif_dir),
        "--workers", str(config.workers),
        "--algorithm", config.dssp_algorithm,
        "--device", device,
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_11", logger, error_tracker, config.repo_dir)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": len(list(modelcif_dir.glob("*.cif"))) if success else 0,
        "models_failed": 0
    }

def stage_12_batch_ipsae(
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 12: BATCH_IPSAE - Calculate interface scores (ipSAE, pDockQ, LIS)"""
    logger.info(f"STAGE 12: BATCH_IPSAE (pae_cutoff={config.pae_cutoff}, dist_cutoff={config.dist_cutoff}, {config.workers} workers)")

    # Directories
    staging_dir = config.output_dir / "staging"
    ipsae_dir = config.output_dir / "ipsae"
    ipsae_input_dir = ipsae_dir / "input"
    ipsae_dir.mkdir(exist_ok=True, parents=True)
    ipsae_input_dir.mkdir(exist_ok=True, parents=True)

    # Create symlinks for PDB files and meta JSON files in the same directory
    # ipsae_cpp expects *-model_v1.pdb and *-meta_v1.json pairs
    start = time.time()
    pdb_count = 0
    json_count = 0

    # Link PDB files
    for pdb_file in staging_dir.glob("*-model_v1.pdb"):
        link_path = ipsae_input_dir / pdb_file.name
        if link_path.exists() or link_path.is_symlink():
            link_path.unlink()
        link_path.symlink_to(pdb_file)
        pdb_count += 1

    # Link meta JSON files (contain PAE data)
    for json_file in staging_dir.glob("*-meta_v1.json"):
        link_path = ipsae_input_dir / json_file.name
        if link_path.exists() or link_path.is_symlink():
            link_path.unlink()
        link_path.symlink_to(json_file)
        json_count += 1

    logger.debug(f"Linked {pdb_count} PDB files and {json_count} meta JSON files", indent=1)

    # Run ipsae_cpp
    ipsae_binary = config.repo_dir / "afdb_integration_kit" / "ipsae" / "ipsae_cpp"
    summary_csv = ipsae_dir / "ipsae_summary.csv"

    cmd = [
        str(ipsae_binary),
        "--batch", str(ipsae_input_dir),
        str(config.pae_cutoff),
        str(config.dist_cutoff),
        "--summary", str(summary_csv),
        "--workers", str(config.workers),
        "--quiet",
        "--no-individual",
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_12", logger, error_tracker, config.repo_dir)

    # Count output lines (excluding header)
    models_succeeded = 0
    if success and summary_csv.exists():
        with open(summary_csv) as f:
            models_succeeded = max(0, sum(1 for _ in f) - 1)  # Subtract header line

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": models_succeeded,
        "models_failed": 0
    }


def _analyze_chunk(
    args: Tuple[List[str], str, float, float, int, Optional[Set[str]], str]
) -> Tuple[int, int, List[str]]:
    """
    Helper function to analyze a chunk of PDB files using direct import.

    Args:
        args: Tuple of (pdb_paths, output_dir, clash_cutoff, interface_cutoff,
              batch_size, analyses, device)

    Returns:
        Tuple of (success_count, fail_count, failed_model_ids)
    """
    pdb_paths, output_dir, clash_cutoff, interface_cutoff, batch_size, analyses, device = args

    try:
        repo_root = Path(__file__).resolve().parent.parent
        if str(repo_root) not in sys.path:
            sys.path.insert(0, str(repo_root))
        from afdb_integration_kit.gpu.analyze import analyze_pdb_files

        results = analyze_pdb_files(
            paths=pdb_paths,
            output_dir=output_dir,
            batch_size=batch_size,
            n_workers=1,  # Single worker per chunk
            device=device,
            clash_cutoff=clash_cutoff,
            interface_cutoff=interface_cutoff,
            progress=False,
            analyses=analyses,
        )
        return (len(results), 0, [])
    except Exception as e:
        # If the whole chunk fails, report all as failed; log the real error for debugging
        logging.exception("Stage 13 chunk failed: %s", e)
        failed_ids = [Path(p).stem for p in pdb_paths]
        return (0, len(pdb_paths), failed_ids)


def stage_13_batch_clash_and_interface(
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 13: BATCH_CLASH_AND_INTERFACE - Compute clashes and interface residues

    Uses analyze_pdb_files_pipelined which overlaps PDB parsing (ProcessPool),
    GPU computation, and JSON writing for maximum throughput.
    """
    n_workers = config.workers
    logger.info(
        f"STAGE 13: BATCH_CLASH_AND_INTERFACE (clash_cutoff={config.clash_cutoff}, "
        f"interface_cutoff={config.interface_cutoff}, batch_size={config.analysis_batch_size}, "
        f"device={config.clash_device}, {n_workers} workers)"
    )

    staging_dir = config.output_dir / "staging"
    analysis_dir = config.output_dir / "clash_interface_analysis"
    analysis_dir.mkdir(exist_ok=True, parents=True)

    pdb_files = list(staging_dir.glob("*-model_v1.pdb"))
    logger.debug(f"Found {len(pdb_files)} PDB files to analyze", indent=1)

    if not pdb_files:
        logger.warning("No PDB files found for analysis", indent=1)
        return {
            "duration": 0.0,
            "success": True,
            "models_succeeded": 0,
            "models_failed": 0
        }

    pdb_paths = [str(p) for p in pdb_files]

    start = time.time()
    total_success = 0
    total_fail = 0

    try:
        repo_root = Path(__file__).resolve().parent.parent
        if str(repo_root) not in sys.path:
            sys.path.insert(0, str(repo_root))
        from afdb_integration_kit.gpu.analyze import analyze_pdb_files_pipelined

        results = analyze_pdb_files_pipelined(
            paths=pdb_paths,
            output_dir=str(analysis_dir),
            batch_size=config.analysis_batch_size,
            n_workers=n_workers,
            device=config.clash_device,
            clash_cutoff=config.clash_cutoff,
            interface_cutoff=config.interface_cutoff,
            progress=False,
            analyses=config.analyses,
        )
        total_success = len(results)
    except Exception as e:
        logging.exception("Stage 13 pipelined analysis failed: %s", e)
        total_fail = len(pdb_paths)
        for pdb_path in pdb_paths:
            model_id = Path(pdb_path).stem
            logger.log_model_failure("stage_13", model_id, str(e), "UNKNOWN")

    duration = time.time() - start
    overall_success = total_fail == 0

    logger.info(
        f"{'✓' if overall_success else '✗'} Completed in {duration:.2f}s ({total_success}/{len(pdb_files)} succeeded)",
        indent=1
    )

    return {
        "duration": duration,
        "success": overall_success,
        "models_succeeded": total_success,
        "models_failed": total_fail
    }


def stage_14_batch_modelpdb(
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 14: BATCH_RUN_MODELPDB - Enrich PDB files"""
    logger.info(f"STAGE 14: BATCH_MODELPDB ({config.workers} workers)")

    staging_dir = config.output_dir / "staging"
    modelcif_dir = config.output_dir / "modelcif"
    modelpdb_dir = config.output_dir / "modelpdb"
    modelpdb_dir.mkdir(exist_ok=True, parents=True)

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "main.py"),
        "run-batch-modelpdb-gen",
        "--cif-dir", str(modelcif_dir),
        "--pdb-dir", str(staging_dir),
        "--provider", str(config.provider_json),
        "--output-dir", str(modelpdb_dir),
        "--model-version", config.model_version,
        "--workers", str(config.workers)
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_14", logger, error_tracker, config.repo_dir)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": len(list(modelpdb_dir.glob("*.pdb"))) if success else 0,
        "models_failed": 0
    }


def stage_15_batch_cif2bcif(
    config: Config,
    logger: PipelineLogger,
    error_tracker: ErrorTracker
) -> Dict[str, Any]:
    """STAGE 15: BATCH_RUN_CIF2BCIF - Convert to BinaryCIF"""

    logger.info(f"STAGE 15: BATCH_CIF2BCIF ({config.workers} workers)")

    modelcif_dir = config.output_dir / "modelcif"
    bcif_dir = config.output_dir / "bcif"
    bcif_dir.mkdir(exist_ok=True, parents=True)

    cmd = [
        *config.python_cmd,
        str(config.repo_dir / "main.py"),
        "batch-cif2bcif",
        "--input-dir", str(modelcif_dir),
        "--output-dir", str(bcif_dir),
        "--workers", str(config.workers)
    ]

    duration, success, stdout, stderr = run_command(cmd, "stage_15", logger, error_tracker, config.repo_dir)

    return {
        "duration": duration,
        "success": success,
        "models_succeeded": len(list(bcif_dir.glob("*.bcif"))) if success else 0,
        "models_failed": 0
    }


CLEANUP_KEEP_DIRS = {
    "modelcif",
    "modelpdb",
    "bcif",
    "chain_batches",
    "model_batches",
    "model_jsons",
    "chain_jsons",
    "scores",
    "ipsae",
    "clash_interface_analysis",
    "input",
    "logs",
}


def _rm_rf(path: Path) -> tuple:
    """Remove a directory tree using system rm for Lustre performance."""
    result = subprocess.run(
        ["rm", "-rf", str(path)],
        capture_output=True, text=True,
    )
    return path.name, result.returncode == 0


def stage_16_cleanup(
    config: Config,
    logger: PipelineLogger
) -> Dict[str, Any]:
    """STAGE 16: BATCH_CLEANUP_FILES - Remove intermediate subdirectories."""
    start = time.time()
    logger.info("STAGE 16: CLEANUP")

    dirs_to_remove = [
        entry for entry in sorted(config.output_dir.iterdir())
        if entry.is_dir() and entry.name not in CLEANUP_KEEP_DIRS
    ]

    if not dirs_to_remove:
        logger.info("  Nothing to clean up")
        return {"duration": 0.0, "success": True, "models_succeeded": 0, "models_failed": 0}

    logger.info(f"  Removing {len(dirs_to_remove)} intermediate directories in parallel")
    removed = []
    failed = []
    with ThreadPoolExecutor(max_workers=min(len(dirs_to_remove), 8)) as pool:
        futures = {pool.submit(_rm_rf, d): d for d in dirs_to_remove}
        for future in as_completed(futures):
            name, success = future.result()
            if success:
                removed.append(name)
            else:
                failed.append(name)
                logger.warning(f"  Failed to remove: {name}/")

    duration = time.time() - start
    logger.info(
        f"  Cleanup complete: removed {len(removed)} directories in {duration:.1f}s"
    )
    return {
        "duration": duration,
        "success": len(failed) == 0,
        "models_succeeded": len(removed),
        "models_failed": len(failed)
    }

# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Main pipeline execution"""
    # Parse configuration
    config = parse_arguments()

    # Create output directory
    config.output_dir.mkdir(exist_ok=True, parents=True)

    # Initialize components
    config_hash = config.get_hash()
    logger = PipelineLogger(config.output_dir, config.log_level)
    cache = PipelineCache(config.output_dir, config_hash) if not config.no_cache else None
    error_tracker = ErrorTracker()

    # Print banner
    logger.info("=" * 80)
    logger.info("PRODUCTION-READY STANDALONE PIPELINE")
    logger.info("=" * 80)
    logger.info("")
    mode_label = "HETERODIMER" if config.heterodimers else "HOMODIMER (default)"
    logger.info(f"Mode: {mode_label}")
    logger.info("Configuration:")
    logger.info(f"  Output directory: {config.output_dir}", indent=1)
    logger.info(f"  Workers: {config.workers}", indent=1)
    logger.info(f"  Model version: {config.model_version}", indent=1)
    logger.info(f"  Resume: {config.resume}", indent=1)
    logger.info(f"  Retry: {config.retry}", indent=1)
    logger.info("")

    if config.dry_run:
        logger.info("DRY RUN MODE - No commands will be executed")
        logger.info("")

    # ---- Load model IDs and auto-generate configs if needed ----
    if config.heterodimers:
        logger.info("Heterodimer mode: deriving model IDs from chain-mapping")
        model_ids = load_model_ids_from_manifest(config.chain_mapping, logger)

        # Auto-generate dataset_config and provider_json if not provided
        config_dir = config.output_dir / "config"
        config_dir.mkdir(exist_ok=True, parents=True)

        if config.dataset_config is None:
            config.dataset_config = config_dir / "dataset_config.json"
            generate_dataset_config(
                config.dataset_config, config.provider_id,
                config.input_dir, config.tool_used,
            )
            logger.info(f"  Generated {config.dataset_config}", indent=1)

        if config.provider_json is None:
            config.provider_json = config_dir / "provider.json"
            generate_provider_json(
                config.provider_json, config.provider_id, config.provider_name
            )
            logger.info(f"  Generated {config.provider_json}", indent=1)

        # Write model_ids.txt (needed by subprocess stages)
        mapping_file = config_dir / "af_mapping.tsv"
        mapping_file.write_text("\n".join(model_ids) + "\n")
        config.mapping_file = mapping_file.resolve()

        logger.info("")
    else:
        model_ids = load_model_ids(config.mapping_file, logger)

        if config.dataset_config is None:
            config_dir = config.output_dir / "config"
            config_dir.mkdir(exist_ok=True, parents=True)
            config.dataset_config = config_dir / "dataset_config.json"
            generate_dataset_config(
                config.dataset_config, config.provider_id,
                config.input_dir, config.tool_used,
            )
            logger.info(f"  Auto-generated {config.dataset_config}", indent=1)

        logger.info("")

    # Run pre-flight checks
    if not config.dry_run:
        if not preflight_checks(config, logger):
            logger.error("Pre-flight checks failed. Exiting.")
            sys.exit(1)
        logger.info("")

    # Stage registry: id -> (name, callable)
    stage_registry = {
        "stage_01": ("Prepare assets", lambda: stage_01_prepare_assets(model_ids, config, logger)),
        "stage_02": ("Validate assets", lambda: stage_02_validate_assets(model_ids, config, logger, error_tracker)),
        "stage_03": ("Convert ColabFold", lambda: stage_03_convert_colabfold(model_ids, config, logger, error_tracker)),
        "stage_04": ("Merge manifests", lambda: stage_04_merge_manifests(config, logger)),
        "stage_12": ("Calculate ipSAE scores", lambda: stage_12_batch_ipsae(config, logger, error_tracker)),
        "stage_13": ("Analyze clashes/interfaces", lambda: stage_13_batch_clash_and_interface(config, logger, error_tracker)),
        "stage_05": ("Export model metadata", lambda: stage_05_export_model_metadata(model_ids, config, logger, error_tracker)),
        "stage_06": ("Export chain metadata", lambda: stage_06_export_chain_metadata(model_ids, config, logger, error_tracker)),
        "stage_07": ("Combine model metadata", lambda: stage_07_combine_model_metadata(config, logger, error_tracker)),
        "stage_08": ("Combine chain metadata", lambda: stage_08_combine_chain_metadata(config, logger, error_tracker)),
        "stage_09": ("Export ModelCIF input", lambda: stage_09_export_modelcif_input(model_ids, config, logger, error_tracker)),
        "stage_10": ("Generate ModelCIF", lambda: stage_10_batch_modelcif_gen(config, logger, error_tracker)),
        "stage_11": ("Add DSSP annotations", lambda: stage_11_batch_dssp(config, logger, error_tracker)),
        "stage_14": ("Enrich PDB files", lambda: stage_14_batch_modelpdb(config, logger, error_tracker)),
        "stage_15": ("Convert to BinaryCIF", lambda: stage_15_batch_cif2bcif(config, logger, error_tracker)),
        "stage_16": ("Cleanup", lambda: stage_16_cleanup(config, logger)),
    }

    # Sequential order: used when parallel_stages is off
    sequential_stage_order = [
        "stage_01", "stage_02", "stage_03", "stage_12", "stage_13",
        "stage_04", "stage_05", "stage_06", "stage_09",
        "stage_07", "stage_08", "stage_10",
        "stage_11", "stage_14", "stage_15", "stage_16",
    ]

    # Parallel waves: groups of stages that run concurrently.
    # Within a wave, heavy stages share the CPU worker budget evenly.
    # Waves execute sequentially; all stages in a wave must finish
    # before the next wave starts.
    execution_waves = [
        ["stage_01"],
        ["stage_02"],
        ["stage_03", "stage_12", "stage_13"],
        ["stage_04"],
        ["stage_05", "stage_06", "stage_09"],
        ["stage_07", "stage_08", "stage_10"],
        ["stage_11"],
        ["stage_14", "stage_15"],
        ["stage_16"],
    ]

    HEAVY_STAGES = {
        "stage_03", "stage_05", "stage_06", "stage_09", "stage_10",
        "stage_11", "stage_12", "stage_13", "stage_14", "stage_15",
    }

    overall_start = time.time()
    results = {}
    total_workers = config.workers

    def _exec_stage(stage_id: str) -> None:
        stage_name, stage_func = stage_registry[stage_id]

        if stage_id in config.skip_stages:
            logger.info(f"⊘ {stage_id} ({stage_name}): Skipped")
            results[stage_id] = {"duration": 0.0, "success": True, "skipped": True}
            return

        if config.resume and cache and cache.is_stage_completed(stage_id):
            logger.info(f"✓ {stage_id} ({stage_name}): Resumed")
            results[stage_id] = {"duration": 0.0, "success": True, "resumed": True}
            return

        if config.dry_run:
            logger.info(f"→ {stage_id} ({stage_name}): Dry run")
            results[stage_id] = {"duration": 0.0, "success": True, "dry_run": True}
            return

        logger.info(f"▶ {stage_id}: {stage_name.upper()}")
        result = stage_func()
        results[stage_id] = result

        if cache and result["success"]:
            cache.mark_stage_completed(
                stage_id, result["duration"],
                result.get("models_succeeded", 0),
                result.get("models_failed", 0),
                result["success"],
            )

        logger.save_stage_failures(stage_id)

        if not result["success"]:
            raise RuntimeError(f"Stage {stage_id} failed")

    execution_mode = "parallel (wave-based DAG)" if config.parallel_stages else "sequential"
    logger.info(f"Execution mode: {execution_mode}")

    stage_failed = False
    try:
        if not config.parallel_stages:
            for stage_id in sequential_stage_order:
                _exec_stage(stage_id)
        else:
            for wave in execution_waves:
                heavy_in_wave = sum(1 for s in wave if s in HEAVY_STAGES)
                if heavy_in_wave > 1:
                    config.workers = max(1, total_workers // heavy_in_wave)
                else:
                    config.workers = total_workers

                if len(wave) == 1:
                    _exec_stage(wave[0])
                else:
                    logger.info(f"  [parallel wave: {', '.join(wave)} | {config.workers} workers each]")
                    with ThreadPoolExecutor(max_workers=len(wave)) as executor:
                        futures = {executor.submit(_exec_stage, sid): sid for sid in wave}
                        for future in as_completed(futures):
                            future.result()

                config.workers = total_workers

    except KeyboardInterrupt:
        logger.info("")
        logger.warning("Pipeline interrupted by user")
        sys.exit(1)
    except RuntimeError as e:
        logger.error(str(e))
        stage_failed = True
    except Exception as e:
        logger.error(f"Pipeline failed with exception: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)

    overall_duration = time.time() - overall_start

    # Summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("PIPELINE SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total execution time: {overall_duration:.2f}s ({overall_duration/60:.2f} min)")
    logger.info(f"Models processed: {len(model_ids)}")
    logger.info(f"Workers used: {config.workers}")
    logger.info("")

    # Stage durations
    total_processing = sum(r.get("duration", 0) for r in results.values())
    logger.info("Stage durations:")
    for stage_id, result in results.items():
        if result.get("skipped"):
            logger.info(f"  {stage_id}: SKIPPED", indent=1)
        elif result.get("resumed"):
            logger.info(f"  {stage_id}: RESUMED", indent=1)
        elif result.get("dry_run"):
            logger.info(f"  {stage_id}: DRY RUN", indent=1)
        else:
            status = "✓" if result.get("success") else "✗"
            logger.info(f"  {stage_id}: {status} {result.get('duration', 0):.2f}s", indent=1)

    logger.info(f"\nTotal processing: {total_processing:.2f}s")
    logger.info("")

    # Error summary
    error_summary = error_tracker.get_error_summary()
    if error_summary:
        logger.info("Error summary:")
        for category, count in error_summary.items():
            logger.info(f"  {category}: {count}", indent=1)
        logger.info("")

    # Failure summary
    failure_summary = logger.get_failure_summary()
    if failure_summary:
        logger.info("Failed models by stage:")
        for stage, count in failure_summary.items():
            logger.info(f"  {stage}: {count} models", indent=1)
        logger.info("")

    # Count output files
    logger.info("Output files generated:")
    output_counts = {
        "Confidence JSONs": len(list((config.output_dir / "scores").glob("*-confidence_v1.json"))) if (config.output_dir / "scores").exists() else 0,
        "PAE JSONs": len(list((config.output_dir / "scores").glob("*-predicted_aligned_error_v1.json"))) if (config.output_dir / "scores").exists() else 0,
        "ModelCIF files": len(list((config.output_dir / "modelcif").glob("*.cif"))) if (config.output_dir / "modelcif").exists() else 0,
        "DSSP-annotated CIFs": len(list((config.output_dir / "modelcif").glob("*.cif"))) if (config.output_dir / "modelcif").exists() else 0,
        "Enriched PDBs": len(list((config.output_dir / "modelpdb").glob("*.pdb"))) if (config.output_dir / "modelpdb").exists() else 0,
        "BCIF files": len(list((config.output_dir / "bcif").glob("*.bcif"))) if (config.output_dir / "bcif").exists() else 0,
    }
    for name, count in output_counts.items():
        logger.info(f"  {name}: {count}", indent=1)
    logger.info("")

    # Save results
    results_file = config.output_dir / "pipeline_results.json"
    results["total_duration"] = overall_duration
    results["model_count"] = len(model_ids)
    results["workers"] = config.workers
    results["error_summary"] = error_summary
    results["failure_summary"] = failure_summary
    results["output_counts"] = output_counts

    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"Results saved to {results_file}")
    logger.info(f"Logs saved to {logger.log_dir}")
    logger.info("")
    logger.info("=" * 80)

    if stage_failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
