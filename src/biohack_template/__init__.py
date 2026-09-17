"""BioHackathon spatial-omics notebook helpers for KIDS26-Team18."""

from biohack_template.celltyping_compare import compare_celltyping_methods
from biohack_template.celltyping_rule_based import assign_celltypes_rule_based, rule_based_summary
from biohack_template.celltyping_scanpy import (
    apply_scanpy_cache,
    assign_celltypes_scanpy,
    build_scanpy_cache,
    run_expensive_scanpy_analysis,
    run_scanpy_clustering,
)
from biohack_template.panels import load_demo_marker_panels
from biohack_template.paths import build_analyst_output_dirs
from biohack_template.plotting import (
    CELL_TYPE_COLORS,
    plot_agreement_bars,
    plot_composition,
    plot_confusion_matrix,
    plot_qc_bar,
    plot_spatial_celltyping_comparison,
    plot_spatial_core,
)
from biohack_template.qc import basic_qc_filter, cells_per_core_summary
from biohack_template.results_io import load_result, save_result
from biohack_template.synthetic_tma import load_data, make_synthetic_tma_adata

__all__ = [
    "CELL_TYPE_COLORS",
    "apply_scanpy_cache",
    "assign_celltypes_rule_based",
    "assign_celltypes_scanpy",
    "basic_qc_filter",
    "build_analyst_output_dirs",
    "build_scanpy_cache",
    "cells_per_core_summary",
    "compare_celltyping_methods",
    "load_data",
    "load_demo_marker_panels",
    "load_result",
    "make_synthetic_tma_adata",
    "plot_agreement_bars",
    "plot_composition",
    "plot_confusion_matrix",
    "plot_qc_bar",
    "plot_spatial_celltyping_comparison",
    "plot_spatial_core",
    "rule_based_summary",
    "run_expensive_scanpy_analysis",
    "run_scanpy_clustering",
    "save_result",
]
