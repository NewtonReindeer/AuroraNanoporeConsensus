"""Lightweight structural types for Aurora result dictionaries.

The production functions intentionally continue returning dictionaries so the
refactor does not change runtime behaviour or report field names.
"""
from typing import Any, TypedDict

class DepthStats(TypedDict):
    mean_depth: float
    median_depth: float
    min_depth: int
    positions_ge_10x_percent: float
    positions_ge_30x_percent: float

class SampleResult(TypedDict, total=False):
    sample: str
    virus_type: str
    best_reference: str
    mapped_reads: int
    qc_verdict: str
    qc_reason: str
    bam_file: str
    consensus_file: str
    reference_file: str
    medaka_model: str
    extra: Any
