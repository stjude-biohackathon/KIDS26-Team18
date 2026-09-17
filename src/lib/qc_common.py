"""QC helpers vendored from cell_typing_v3.lib.clustering_common."""

from __future__ import annotations

import logging
import os
from typing import Any

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp


def ensure_output_dir(output_path: str) -> None:
    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)


def ensure_csr_float32(x):
    if sp.isspmatrix_csr(x):
        return x.astype(np.float32, copy=False)
    if sp.issparse(x):
        return x.tocsr().astype(np.float32, copy=False)
    return np.asarray(x, dtype=np.float32)


def make_ribo_mask(adata: ad.AnnData, cfg: dict[str, Any]) -> None:
    ribo_genes = cfg.get("ribo_genes", None)

    if ribo_genes:
        ribo_set = {str(g).upper() for g in ribo_genes}
        gene_names = adata.var_names.astype(str).str.upper()
        adata.var["ribo"] = gene_names.isin(ribo_set)
        logging.info("Ribosomal whitelist enabled: %d genes matched", int(adata.var["ribo"].sum()))
        return

    ribo_regex = cfg.get("ribo_regex", None)
    if ribo_regex:
        gene_names = adata.var_names.astype(str).str.upper()
        adata.var["ribo"] = gene_names.str.match(str(ribo_regex))
        logging.info("Ribosomal regex enabled: %d genes matched", int(adata.var["ribo"].sum()))
        return

    logging.info("No ribosomal gene rule provided; pct_counts_ribo will not be computed.")


def compute_missing_qc_metrics(adata: ad.AnnData, cfg: dict[str, Any]) -> None:
    four_qc_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]

    if all(col in adata.obs.columns for col in four_qc_cols):
        logging.info(
            "Reusing existing QC columns to call QC_flag: %s",
            ", ".join(four_qc_cols),
        )
        return

    has_basic = all(col in adata.obs.columns for col in ["total_counts", "n_genes_by_counts"])
    has_mt = "pct_counts_mt" in adata.obs.columns
    has_ribo = "pct_counts_ribo" in adata.obs.columns

    qc_vars = []

    want_mt = cfg.get("pct_counts_mt_max", cfg.get("pct_mt_max", None)) is not None or bool(
        cfg.get("compute_mt", True)
    )
    mt_prefix = cfg.get("mt_prefix", None)

    if want_mt and not has_mt and mt_prefix is not None:
        gene_names = adata.var_names.astype(str).str.upper()
        adata.var["mt"] = gene_names.str.startswith(str(mt_prefix).upper())
        qc_vars.append("mt")
        logging.info("MT genes matched: %d", int(adata.var["mt"].sum()))
    elif has_mt:
        logging.info("Reusing existing pct_counts_mt column.")

    want_ribo = cfg.get("pct_counts_ribo_max", None) is not None or bool(cfg.get("compute_ribo", True))

    if want_ribo and not has_ribo:
        make_ribo_mask(adata, cfg)
        if "ribo" in adata.var.columns and int(adata.var["ribo"].sum()) > 0:
            qc_vars.append("ribo")
    elif has_ribo:
        logging.info("Reusing existing pct_counts_ribo column.")

    if has_basic and not qc_vars:
        logging.info(
            "Reusing existing total_counts and n_genes_by_counts; no missing QC metrics to compute."
        )
        return

    logging.info("Computing missing QC metrics with qc_vars=%s", qc_vars)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_vars,
        percent_top=cfg.get("percent_top", None),
        log1p=bool(cfg.get("qc_log1p", False)),
        inplace=True,
    )


def call_qc_flag(adata: ad.AnnData, cfg: dict[str, Any]) -> ad.AnnData:
    qc_flag_col = str(cfg.get("qc_flag_col", "QC_flag"))

    total_counts_min = cfg.get("total_counts_min", None)
    n_genes_min = cfg.get("n_genes_min", None)
    pct_mt_max = cfg.get("pct_counts_mt_max", cfg.get("pct_mt_max", None))
    pct_ribo_max = cfg.get("pct_counts_ribo_max", None)

    fail_mask = np.zeros(adata.n_obs, dtype=bool)

    if total_counts_min is not None:
        fail_mask |= adata.obs["total_counts"].to_numpy() < float(total_counts_min)

    if n_genes_min is not None:
        fail_mask |= adata.obs["n_genes_by_counts"].to_numpy() < float(n_genes_min)

    if pct_mt_max is not None:
        if "pct_counts_mt" not in adata.obs.columns:
            logging.warning("pct_counts_mt_max was set, but pct_counts_mt is missing. Ignoring MT threshold.")
        else:
            fail_mask |= adata.obs["pct_counts_mt"].to_numpy() > float(pct_mt_max)

    if pct_ribo_max is not None:
        if "pct_counts_ribo" not in adata.obs.columns:
            logging.warning(
                "pct_counts_ribo_max was set, but pct_counts_ribo is missing. Ignoring ribo threshold."
            )
        else:
            fail_mask |= adata.obs["pct_counts_ribo"].to_numpy() > float(pct_ribo_max)

    adata.obs[qc_flag_col] = np.where(fail_mask, "fail", "pass")

    n_pass = int((adata.obs[qc_flag_col] == "pass").sum())
    n_fail = int((adata.obs[qc_flag_col] == "fail").sum())

    logging.info(
        "QC_flag assigned using thresholds: total_counts_min=%s, n_genes_min=%s, "
        "pct_counts_mt_max=%s, pct_counts_ribo_max=%s",
        total_counts_min,
        n_genes_min,
        pct_mt_max,
        pct_ribo_max,
    )
    logging.info("QC_flag counts: pass=%d, fail=%d, total=%d", n_pass, n_fail, adata.n_obs)

    return adata


def subset_raw_adata_copy(adata: ad.AnnData, mask, reason: str = "subset") -> ad.AnnData:
    logging.info("%s: keeping %d / %d cells", reason, int(np.sum(mask)), adata.n_obs)
    adata_sub = adata[mask].copy()
    adata_sub.X = ensure_csr_float32(adata_sub.X)
    return adata_sub


def run_qc_flagging_and_filtering(adata: ad.AnnData, cfg: dict[str, Any]) -> ad.AnnData:
    if not bool(cfg.get("do_qc", True)):
        logging.info("QC disabled")
        return adata

    compute_missing_qc_metrics(adata, cfg)
    adata = call_qc_flag(adata, cfg)

    if bool(cfg.get("filter_qc_pass", True)):
        qc_flag_col = str(cfg.get("qc_flag_col", "QC_flag"))
        keep_mask = adata.obs[qc_flag_col].to_numpy() == "pass"
        return subset_raw_adata_copy(
            adata=adata,
            mask=keep_mask,
            reason="QC filtering",
        )

    logging.info("filter_qc_pass=false; keeping all cells with QC_flag annotation.")
    return adata
