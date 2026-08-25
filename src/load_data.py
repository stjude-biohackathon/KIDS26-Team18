"""Data loading helpers extracted from Utils.py."""

from __future__ import annotations

import json
import re
import shutil
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

def _debug_log(location, message, data, hypothesis_id, run_id="pre-fix"):
    payload = {
        "sessionId": "92bf40",
        "runId": run_id,
        "hypothesisId": hypothesis_id,
        "location": location,
        "message": message,
        "data": data,
        "timestamp": int(time.time() * 1000),
    }
    with _DEBUG_LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(payload) + "\n")
# #endregion



def _load_metadata(
    merged_parquet=None,
    persample_parquets=None,
    per_sample_column="celltype_from_leiden_unintegrated_r0.7",
    verbose=True,
):
    """
    Load merged metadata and optionally merge per-sample annotation columns.

    per_sample_column controls which columns are pulled from persample_parquets:
    - str (default): merge that single column name
    - list/tuple: merge only columns present in each per-sample file
    - "all": merge all data columns (excluding cell_id used as index)
    """

    def _select_columns(df):
        if per_sample_column == "all":
            return [c for c in df.columns if c != "cell_id"]
        if isinstance(per_sample_column, (list, tuple)):
            return [c for c in per_sample_column if c in df.columns]
        if per_sample_column in df.columns:
            return [per_sample_column]
        return []

    # ------------------------------------------------------------
    # Load merged metadata
    # ------------------------------------------------------------
    if merged_parquet is None:
        raise ValueError("merged_parquet must be provided.")

    merged_parquet = Path(merged_parquet)

    if not merged_parquet.is_file():
        raise FileNotFoundError(merged_parquet)

    meta = pd.read_parquet(merged_parquet)

    if "cell_id" in meta.columns:
        meta = meta.set_index("cell_id")

    meta.index = meta.index.astype(str)

    # ------------------------------------------------------------
    # Load per-sample annotations (optional)
    # ------------------------------------------------------------
    if persample_parquets is None:
        return meta

    if isinstance(persample_parquets, (str, Path)):
        persample_parquets = [persample_parquets]

    parquet_files = []

    for p in persample_parquets:

        p = Path(p)

        if p.is_dir():
            parquet_files.extend(sorted(p.glob("*.parquet")))
        else:
            parquet_files.append(p)

    if len(parquet_files) == 0:
        return meta

    annotation_list = []

    for p in parquet_files:

        df = pd.read_parquet(p)

        if "cell_id" in df.columns:
            df = df.set_index("cell_id")

        df.index = df.index.astype(str)

        cols = _select_columns(df)
        if not cols:
            continue

        annotation_list.append(df[cols])

    if len(annotation_list):

        annotations = pd.concat(annotation_list)
        annotations = annotations[~annotations.index.duplicated(keep="first")]

        new_cols = [c for c in annotations.columns if c not in meta.columns]
        if new_cols:
            meta = meta.join(annotations[new_cols], how="left")

    return meta

# def _load_metadata(
#     metadata_parquet_path,
#     per_sample_column="celltype_from_leiden_unintegrated_r0.7",
#     verbose=True,
# ):

#     # allow single path or list
#     if isinstance(metadata_parquet_path, (str, Path)):
#         metadata_parquet_path = [metadata_parquet_path]

#     integrated = None

#     for p in metadata_parquet_path:
#         p = Path(p)
#         df = pd.read_parquet(p)

#         if "cell_id" in df.columns:
#             df = df.set_index("cell_id")

#         df.index = df.index.astype(str)

#         # integrated metadata
#         if integrated is None:
#             integrated = df.copy()
#             continue

#         # per-sample metadata
#         if per_sample_column in df.columns:
#             integrated = integrated.join(
#                 df[[per_sample_column]],
#                 how="left",
#                 rsuffix="_tmp",
#             )

#     return integrated


def _h5ad_obs_column_names(h5ad_path):
    with h5py.File(h5ad_path, "r") as f:
        return [key for key in f["obs"].keys() if key != "_index"]


def _restore_umaps(adata, verbose=True):
    for prefix in [
        "X_umap_unintegrated",
        "X_umap_harmony_rsc",
        "X_umap_scVI",
    ]:
        c1 = f"{prefix}_1"
        c2 = f"{prefix}_2"
        if c1 in adata.obs.columns and c2 in adata.obs.columns:
            adata.obsm[prefix] = (
                adata.obs[[c1, c2]]
                .astype(float)
                .to_numpy()
            )
            if verbose:
                print(f"Restored {prefix}")


def _resolve_gene_columns(genes, var_names, verbose=True):
    """
    Resolve ``genes`` to present gene names and h5ad column indices.

    Returns
    -------
    present_genes : list[str]
        Gene names to keep, in request order.
    col_indices : np.ndarray or None
        Integer column positions for h5ad ``X``, or None for all genes.
    n_var_total : int
        Total number of genes in the source object.
    """
    var_names = pd.Index(var_names).astype(str)
    n_var_total = len(var_names)

    if genes is None or genes == "all_genes":
        return list(var_names), None, n_var_total

    if isinstance(genes, str):
        raise ValueError(
            "genes must be None, 'all_genes', or a list of gene names; "
            f"got {genes!r}."
        )

    genes = [str(g) for g in genes]
    if len(genes) == 0:
        raise ValueError("genes list is empty.")

    col_pos = var_names.get_indexer(genes)
    present = []
    col_indices = []
    missing = []

    for gene, pos in zip(genes, col_pos):
        if pos >= 0:
            present.append(gene)
            col_indices.append(int(pos))
        else:
            missing.append(gene)

    if verbose and missing:
        print(f"Requested genes missing from adata: {missing}")

    if len(present) == 0:
        raise ValueError(
            "None of the requested genes were found in adata.var_names."
        )

    if verbose and len(present) < n_var_total:
        print(f"Loading {len(present):,} / {n_var_total:,} genes.")

    return present, np.asarray(col_indices, dtype=np.int64), n_var_total


def _build_filtered_csr(data, indices, indptr, n_rows, n_var, col_indices=None):
    """Build a CSR matrix from raw arrays, optionally subsetting columns."""
    if col_indices is not None:
        col_indices = np.asarray(col_indices, dtype=np.int64)
        n_cols = len(col_indices)
        col_map = np.full(n_var, -1, dtype=np.int32)
        col_map[col_indices] = np.arange(n_cols, dtype=np.int32)

        mapped = col_map[indices]
        keep = mapped >= 0
        f_data = data[keep]
        f_indices = mapped[keep]

        new_indptr = np.empty(n_rows + 1, dtype=np.int64)
        new_indptr[0] = 0
        pos = 0
        for r in range(n_rows):
            start, end_r = indptr[r], indptr[r + 1]
            pos += int(keep[start:end_r].sum())
            new_indptr[r + 1] = pos

        return sp.csr_matrix(
            (f_data, f_indices, new_indptr),
            shape=(n_rows, n_cols),
        )

    return sp.csr_matrix((data, indices, indptr), shape=(n_rows, n_var))


def _read_csr_row_slice(h5ad_path, row_start, row_end, col_indices=None, verbose=False):
    """
    Read CSR rows ``[row_start, row_end)`` from an h5ad ``X`` group.

    Unlike :func:`_read_csr_rows`, this reads only the requested row block
    from disk (no prefix re-read from row 0).
    """
    row_start = int(row_start)
    row_end = int(row_end)
    if row_start < 0 or row_end <= row_start:
        raise ValueError(
            f"Invalid row slice: row_start={row_start}, row_end={row_end}"
        )

    with h5py.File(h5ad_path, "r") as f:
        Xg = f["X"]
        attrs = dict(Xg.attrs)
        encoding = attrs.get("encoding-type", "csr_matrix")
        if encoding not in ("csr_matrix", "csc_matrix"):
            raise ValueError(f"Unsupported X encoding: {encoding}")
        if encoding == "csc_matrix":
            raise ValueError(
                "X is stored as CSC; this fast loader expects CSR. "
                "Fall back to backed loading for CSC."
            )

        shape = tuple(int(s) for s in attrs["shape"])
        n_var = shape[1]
        n_rows = row_end - row_start

        indptr_full = Xg["indptr"]
        data_start = int(indptr_full[row_start])
        data_end = int(indptr_full[row_end])

        if verbose:
            gene_msg = (
                f", {len(col_indices):,} genes"
                if col_indices is not None
                else ""
            )
            print(
                f"Reading X rows[{row_start:,}:{row_end:,}] "
                f"(nnz={data_end - data_start:,}{gene_msg}) ..."
            )

        indptr = (
            indptr_full[row_start : row_end + 1].astype(np.int64) - data_start
        )
        data = Xg["data"].astype("float32")[data_start:data_end]
        indices = Xg["indices"][data_start:data_end]

    mat = _build_filtered_csr(
        data,
        indices,
        indptr,
        n_rows,
        n_var,
        col_indices=col_indices,
    )

    if not mat.has_sorted_indices:
        mat.sort_indices()

    return mat


def _read_csr_rows(h5ad_path, row_pos, col_indices=None, verbose=True):
    """
    Bulk-read the CSR expression matrix from an h5ad ``X`` group and return a
    ``scipy.sparse.csr_matrix`` restricted to ``row_pos`` (integer positions),
    with data cast to ``float32``.

    This reads the ``data``/``indices``/``indptr`` arrays sequentially (fast,
    I/O-bound) and slices the requested rows in memory. It is dramatically
    faster than backed fancy row-indexing (``adata[cells, :].to_memory()``),
    which reads one row-group at a time and is the main load bottleneck for
    very large matrices.

    The read is bounded by the largest requested row position, so when the
    requested rows are a small prefix (e.g. a quick test subset) only that
    portion of the file is read.

    When ``col_indices`` is provided, only those columns are materialized.
    """
    row_pos = np.asarray(row_pos, dtype=np.int64)
    hi = int(row_pos.max()) + 1  # exclusive upper bound of rows to read

    with h5py.File(h5ad_path, "r") as f:
        Xg = f["X"]
        attrs = dict(Xg.attrs)
        encoding = attrs.get("encoding-type", "csr_matrix")
        if encoding not in ("csr_matrix", "csc_matrix"):
            raise ValueError(f"Unsupported X encoding: {encoding}")
        if encoding == "csc_matrix":
            raise ValueError(
                "X is stored as CSC; this fast loader expects CSR. "
                "Fall back to backed loading for CSC."
            )

        shape = tuple(int(s) for s in attrs["shape"])
        n_var = shape[1]

        indptr = Xg["indptr"][: hi + 1].astype(np.int64)
        end = int(indptr[-1])

        if verbose:
            gene_msg = (
                f", {len(col_indices):,} genes"
                if col_indices is not None
                else ""
            )
            print(
                f"Bulk-reading X rows[0:{hi:,}] of {shape[0]:,} "
                f"(nnz={end:,}{gene_msg}) ..."
            )

        # Chunked read + cast to float32 (avoids a full float64 temporary).
        data = Xg["data"].astype("float32")[:end]
        indices = Xg["indices"][:end]

    full = _build_filtered_csr(
        data,
        indices,
        indptr,
        hi,
        n_var,
        col_indices=col_indices,
    )
    del data, indices, indptr

    sub = full[row_pos, :]
    del full

    if not sub.has_sorted_indices:
        sub.sort_indices()

    return sub


def load_expression_matrix(
    h5ad_path,
    cells=None,
    metadata_parquet_path=None,
    genes=None,
    max_cells=None,
    seed=66,
    verbose=True,
):
    """
    Load the expression matrix once as a ``float32`` CSR AnnData.

    The returned AnnData has ``X`` + ``var`` populated and a minimal ``obs``
    (index only). Use :func:`attach_metadata` to build annotated objects
    (e.g. per parquet) that share/subset this matrix without re-reading the
    huge h5ad file.

    Parameters
    ----------
    cells
        Optional iterable of cell ids to load. If given, only those cells are
        materialized. Ignored when ``metadata_parquet_path`` is provided.
    metadata_parquet_path
        Optional parquet path/dir whose index defines the cells to load.
    genes
        Gene subset to load. ``None`` or ``"all_genes"`` loads every gene;
        otherwise pass a list/tuple/set of gene names (loaded in that order).
        Missing genes are skipped with a warning.
    max_cells
        Optional cap on the number of cells for quick end-to-end testing.
        To keep the read fast this keeps the first ``max_cells`` cells in
        on-disk order (a contiguous prefix), so only that slice of the file is
        read. Defaults to all cells. NOTE: this is a smoke-test convenience,
        not a representative random sample.
    """
    h5ad_path = str(h5ad_path)

    adata_backed = sc.read_h5ad(h5ad_path, backed="r")
    adata_backed.obs_names = adata_backed.obs_names.astype(str)
    all_names = adata_backed.obs_names
    var = adata_backed.var.copy()
    var.index = var.index.astype(str)
    adata_backed.file.close()

    if metadata_parquet_path is not None:
        meta_index = _load_metadata(metadata_parquet_path, verbose=verbose).index
        cells = meta_index.intersection(all_names)
    elif cells is not None:
        cells = pd.Index(pd.Series(cells).astype(str)).intersection(all_names)
    else:
        cells = all_names

    if len(cells) == 0:
        raise ValueError("No overlapping cells found for expression load.")

    row_pos = all_names.get_indexer(cells)
    if (row_pos < 0).any():
        raise ValueError("Some requested cells are absent from the h5ad.")

    if max_cells is not None and max_cells < len(cells):
        # Keep the cells with the smallest on-disk positions so the bulk read
        # only touches a prefix of the file (fast smoke test).
        keep = np.argsort(row_pos, kind="stable")[: int(max_cells)]
        keep = np.sort(keep)
        cells = cells[keep]
        row_pos = row_pos[keep]
        if verbose:
            print(
                f"[testing] Using the first {len(cells):,} cells "
                f"(max_cells; contiguous prefix)."
            )

    present_genes, col_indices, _ = _resolve_gene_columns(
        genes,
        var.index,
        verbose=verbose,
    )
    var = var.loc[present_genes]

    X = _read_csr_rows(
        h5ad_path,
        row_pos,
        col_indices=col_indices,
        verbose=verbose,
    )

    obs = pd.DataFrame(index=pd.Index(cells, name=None).astype(str))
    adata = ad.AnnData(X=X, obs=obs, var=var)

    if verbose:
        print(f"\nLoaded expression matrix: {adata.n_obs:,} x {adata.n_vars:,}")
        print(f"X dtype: {adata.X.dtype}, nnz: {adata.X.nnz:,}")

    return adata


def attach_metadata(
    adata_expr,
    metadata_parquet_path,
    subset_to_metadata=True,
    restore_umaps=True,
    verbose=True,
):
    """
    Build a new annotated AnnData from an already-loaded expression object
    (see :func:`load_expression_matrix`) by setting ``.obs`` to the parquet
    metadata. ``X`` is subset/copied to the overlapping cells; ``adata_expr``
    is left unchanged. Avoids re-reading the h5ad file.
    """
    # #region agent log
    _debug_log(
        "read_adata_parquet_mode.py:attach_metadata:entry",
        "attach_metadata called",
        {
            "module_file": __file__,
            "input_obs_n_cols": int(adata_expr.obs.shape[1]),
            "input_obs_cols_sample": list(adata_expr.obs.columns[:8]),
            "metadata_parquet_path": str(metadata_parquet_path),
        },
        "H1",
    )
    # #endregion

    meta = _load_metadata(metadata_parquet_path, verbose=verbose)

    # #region agent log
    _debug_log(
        "read_adata_parquet_mode.py:attach_metadata:meta_loaded",
        "metadata loaded from parquet",
        {
            "meta_n_cols": int(meta.shape[1]),
            "meta_cols": list(meta.columns),
        },
        "H3",
    )
    # #endregion

    if subset_to_metadata:
        cells = meta.index.intersection(adata_expr.obs_names)
    else:
        cells = adata_expr.obs_names

    if len(cells) == 0:
        raise ValueError(
            "No overlapping cells between expression matrix and metadata."
        )

    if verbose:
        print(f"\nCells in metadata parquet(s): {len(meta):,}")
        print(f"Cells selected from expression matrix: {len(cells):,}")

    subset = adata_expr[cells, :]
    obs_new = meta.reindex(cells).copy()
    obs_new.index = obs_new.index.astype(str)

    X = subset.X
    if hasattr(X, "copy"):
        X = X.copy()

    adata = ad.AnnData(
        X=X,
        obs=obs_new,
        var=subset.var.copy(),
    )

    # #region agent log
    _debug_log(
        "read_adata_parquet_mode.py:attach_metadata:after_fresh_anndata",
        "obs columns on freshly built AnnData",
        {
            "obs_n_cols": int(adata.obs.shape[1]),
            "obs_cols": list(adata.obs.columns),
            "matches_meta_cols": list(adata.obs.columns) == list(meta.columns),
        },
        "H2",
    )
    # #endregion

    if list(adata.obs.columns) != list(meta.columns):
        raise RuntimeError(
            f"attach_metadata obs has {len(adata.obs.columns)} columns but parquet "
            f"metadata has {len(meta.columns)}. Use "
            "importlib.reload(read_adata_parquet_mode) and call "
            "read_adata_parquet_mode.attach_metadata(...) instead of a stale "
            "`from ... import attach_metadata` binding."
        )

    if restore_umaps:
        _restore_umaps(adata, verbose=verbose)

    # #region agent log
    _debug_log(
        "read_adata_parquet_mode.py:attach_metadata:return",
        "final obs columns before return",
        {
            "obs_n_cols": int(adata.obs.shape[1]),
            "obs_cols": list(adata.obs.columns),
            "obsm_keys": list(adata.obsm.keys()),
        },
        "H5",
    )
    # #endregion

    return adata


def _extract_gene_panel_script(h5ad_path):
    """Locate extract_gene_panel.py relative to the project h5ad layout."""
    candidate = (
        Path(h5ad_path).resolve().parents[2]
        / "scripts/prepare_obj/extract_gene_panel.py"
    )
    if candidate.is_file():
        return candidate
    fallback = Path(
        "/mnt/scratch2/Maycon/Illumina_SpatialData/after_May282026/"
        "scripts/prepare_obj/extract_gene_panel.py"
    )
    if fallback.is_file():
        return fallback
    raise FileNotFoundError(
        "Could not locate extract_gene_panel.py. "
        "Run scripts/prepare_obj/extract_gene_panel.py manually."
    )


def _genes_for_panel(gene_panel_h5ad, genes):
    """Resolve gene list from argument or sibling *_genes.txt file."""
    if genes is not None and genes != "all_genes":
        if isinstance(genes, str):
            return [genes]
        return list(genes)

    genes_path = Path(gene_panel_h5ad).with_name(
        f"{Path(gene_panel_h5ad).stem}_genes.txt"
    )
    if genes_path.is_file():
        return [
            line.strip()
            for line in genes_path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        ]
    return None


def _build_gene_panel_h5ad(
    h5ad_path,
    gene_panel_h5ad,
    genes,
    metadata_parquet=None,
    verbose=True,
):
    """Run extract_gene_panel.py to create a missing gene panel."""
    import subprocess
    import sys

    panel_path = Path(gene_panel_h5ad)
    gene_list = _genes_for_panel(gene_panel_h5ad, genes)
    if not gene_list:
        raise FileNotFoundError(
            f"Gene panel not found: {gene_panel_h5ad}\n"
            "Provide genes=... or create "
            f"{panel_path.with_name(panel_path.stem + '_genes.txt')}, "
            "then run extract_gene_panel.py."
        )

    extract_script = _extract_gene_panel_script(h5ad_path)
    cmd = [
        sys.executable,
        str(extract_script),
        "--h5ad",
        str(h5ad_path),
        "--output-dir",
        str(panel_path.parent),
        "--panel-name",
        panel_path.stem,
        "--genes",
        *gene_list,
    ]
    if metadata_parquet is not None:
        cmd.extend(["--metadata-parquet", str(metadata_parquet)])

    if verbose:
        print(
            f"Gene panel missing; building {panel_path} "
            f"({len(gene_list)} genes). This may take hours."
        )
        print("Command:", " ".join(cmd))

    subprocess.run(cmd, check=True)

    if not panel_path.is_file():
        raise FileNotFoundError(
            f"Gene panel build finished but file not found: {panel_path}"
        )


def _check_gene_panel_manifest(gene_panel_h5ad, h5ad_path, verbose=True):
    """Warn if the gene panel was built from an older source h5ad."""
    panel_path = Path(gene_panel_h5ad)
    manifest_path = panel_path.with_name(f"{panel_path.stem}_manifest.json")
    if not manifest_path.is_file():
        return

    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        if verbose:
            print(f"Warning: could not read gene panel manifest: {manifest_path}")
        return

    source_h5ad = manifest.get("source_h5ad")
    source_mtime = manifest.get("source_h5ad_mtime")
    if source_h5ad and Path(source_h5ad).resolve() != Path(h5ad_path).resolve():
        if verbose:
            print(
                f"Warning: gene panel source h5ad differs from h5ad_path.\n"
                f"  panel manifest: {source_h5ad}\n"
                f"  h5ad_path:      {h5ad_path}"
            )
        return

    if source_mtime is not None and Path(h5ad_path).is_file():
        current_mtime = Path(h5ad_path).stat().st_mtime
        if current_mtime > float(source_mtime) + 1:
            if verbose:
                print(
                    "Warning: source h5ad is newer than the gene panel manifest. "
                    "Re-run extract_gene_panel.py."
                )


def _join_metadata_obs(adata, meta, h5ad_obs=None, verbose=True):
    """Attach parquet metadata (and optional h5ad-only obs cols) to adata."""
    if h5ad_obs is not None:
        adata.obs = adata.obs.join(h5ad_obs, how="left")

    meta_aligned = meta.reindex(adata.obs_names)
    overlap = adata.obs.columns.intersection(meta_aligned.columns)

    if len(overlap) > 0:
        if verbose:
            print(
                f"\nDropping overlapping columns from adata.obs: "
                f"{list(overlap)}"
            )
        adata.obs = adata.obs.drop(columns=overlap)

    adata.obs = adata.obs.join(meta_aligned, how="left")
    return adata


def _read_h5ad_only_obs(h5ad_path, cells, meta, verbose=True):
    """Load obs columns present in h5ad but absent from parquet metadata."""
    adata_backed = sc.read_h5ad(h5ad_path, backed="r")
    adata_backed.obs_names = adata_backed.obs_names.astype(str)

    h5ad_only_cols = [
        col for col in adata_backed.obs.columns
        if col not in meta.columns
    ]

    if h5ad_only_cols:
        h5ad_obs = adata_backed.obs.loc[cells, h5ad_only_cols].copy()
        h5ad_obs.index = h5ad_obs.index.astype(str)
    else:
        h5ad_obs = None

    adata_backed.file.close()
    return h5ad_obs


def read_adata(
    h5ad_path,
    metadata_parquet_path=None,
    merged_parquet=None,
    persample_parquets=None,
    per_sample_column="celltype_from_leiden_unintegrated_r0.7",
    load_expression=False,
    subset_to_metadata=True,
    keep_var=True,
    restore_umaps=True,
    genes=None,
    gene_panel_h5ad=None,
    build_gene_panel_if_missing=False,
    max_cells=None,
    verbose=True,
):
    """
    Load h5ad metadata and update `.obs` from parquet.

    Metadata can be supplied in one of two ways:

    1. metadata_parquet_path
       Original API. Can be a parquet file or directory.

    2. merged_parquet + persample_parquets
       merged_parquet:
           Single merged metadata parquet used as the base metadata.

       persample_parquets:
           Optional parquet, list of parquets, or directory containing
           per-sample metadata. The specified annotation column(s) are merged
           onto the base metadata.

       per_sample_column:
           Which column(s) to merge from persample_parquets. Can be a single
           column name (default: ``celltype_from_leiden_unintegrated_r0.7``),
           a list of column names, or ``"all"`` to merge every data column
           from each per-sample parquet.

    By default this returns a lightweight in-memory AnnData object with
    updated `.obs` and `.var`, without loading the expression matrix.

    Set ``load_expression=True`` to materialize the expression matrix.

    genes
        Applied only when ``load_expression=True`` and ``gene_panel_h5ad`` is
        not set. ``None`` or ``"all_genes"`` loads every gene; otherwise pass
        a list of gene names to load only those columns (in request order).
        Missing genes are skipped with a warning.

    gene_panel_h5ad
        Optional path to a pre-extracted gene-panel h5ad (see
        ``extract_gene_panel.py``). When set with ``load_expression=True``,
        expression is loaded from this small file instead of the source h5ad.

    build_gene_panel_if_missing
        If ``gene_panel_h5ad`` is missing, run ``extract_gene_panel.py``
        once (requires ``genes=`` or a sibling ``*_genes.txt`` file).
        Uses ``merged_parquet`` / ``metadata_parquet_path`` to keep only
        metadata cells when building the panel.
    """

    h5ad_path = str(h5ad_path)
    if gene_panel_h5ad is not None:
        gene_panel_h5ad = str(gene_panel_h5ad)

    # ------------------------------------------------------------
    # Load metadata
    # ------------------------------------------------------------
    if merged_parquet is not None:

        meta = _load_metadata(
            merged_parquet=merged_parquet,
            persample_parquets=persample_parquets,
            per_sample_column=per_sample_column,
            verbose=verbose,
        )

    elif metadata_parquet_path is not None:

        meta = _load_metadata(
            merged_parquet=metadata_parquet_path,
            verbose=verbose,
        )

    else:
        raise ValueError(
            "Provide either `metadata_parquet_path` or "
            "`merged_parquet`."
        )

    if not load_expression:
        if genes is not None and genes != "all_genes" and verbose:
            print("Note: genes ignored because load_expression=False.")
        if gene_panel_h5ad is not None and verbose:
            print("Note: gene_panel_h5ad ignored because load_expression=False.")

    use_gene_panel = load_expression and gene_panel_h5ad is not None

    if use_gene_panel:
        panel_path = Path(gene_panel_h5ad)
        if not panel_path.is_file():
            if build_gene_panel_if_missing:
                meta_parquet = merged_parquet or metadata_parquet_path
                _build_gene_panel_h5ad(
                    h5ad_path,
                    gene_panel_h5ad,
                    genes,
                    metadata_parquet=meta_parquet,
                    verbose=verbose,
                )
            else:
                raise FileNotFoundError(
                    f"Gene panel not found: {gene_panel_h5ad}\n"
                    "Run scripts/prepare_obj/extract_gene_panel.py first, e.g.:\n"
                    "  python extract_gene_panel.py \\\n"
                    "    --genes EPCAM KLK3 AR MSMB CD3D CD79A PECAM1 COL1A1 \\\n"
                    "    --panel-name fig1_markers_8genes \\\n"
                    "    --metadata-parquet <merged__cellmetadata.parquet>\n"
                    "Or pass build_gene_panel_if_missing=True with genes=..."
                )

        _check_gene_panel_manifest(gene_panel_h5ad, h5ad_path, verbose=verbose)

        panel_backed = sc.read_h5ad(gene_panel_h5ad, backed="r")
        panel_backed.obs_names = panel_backed.obs_names.astype(str)

        if subset_to_metadata:
            cells = meta.index.intersection(panel_backed.obs_names)
        else:
            cells = panel_backed.obs_names

        if len(cells) == 0:
            panel_backed.file.close()
            raise ValueError(
                "No overlapping cells between gene panel and metadata index."
            )

        if max_cells is not None and max_cells < len(cells):
            panel_positions = panel_backed.obs_names.get_indexer(cells)
            keep = np.sort(np.argsort(panel_positions, kind="stable")[: int(max_cells)])
            cells = cells[keep]

            if verbose:
                print(
                    f"[testing] Using the first {len(cells):,} cells "
                    f"(max_cells; contiguous prefix)."
                )

        if verbose:
            print(f"\nCells in metadata parquet(s): {len(meta):,}")
            print(f"Cells selected from gene panel: {len(cells):,}")
            print(f"Loading expression from gene panel: {gene_panel_h5ad}")

        if genes is not None and genes != "all_genes" and verbose:
            print(
                "Note: subsetting gene panel further with genes= "
                "(panel genes are loaded first)."
            )

        h5ad_obs = _read_h5ad_only_obs(h5ad_path, cells, meta, verbose=verbose)

        adata = panel_backed[list(cells)].to_memory()
        panel_backed.file.close()

        if genes is not None and genes != "all_genes":
            present_genes, _, _ = _resolve_gene_columns(
                genes,
                adata.var_names,
                verbose=verbose,
            )
            adata = adata[:, present_genes].copy()

        adata = _join_metadata_obs(adata, meta, h5ad_obs=h5ad_obs, verbose=verbose)
        adata.uns["__gene_panel_h5ad__"] = gene_panel_h5ad
        adata.uns["__h5ad_path__"] = h5ad_path

    elif load_expression:
        adata_backed = sc.read_h5ad(h5ad_path, backed="r")
        adata_backed.obs_names = adata_backed.obs_names.astype(str)

        if subset_to_metadata:
            cells = meta.index.intersection(adata_backed.obs_names)
        else:
            cells = adata_backed.obs_names

        if len(cells) == 0:
            adata_backed.file.close()
            raise ValueError(
                "No overlapping cells found between h5ad obs_names and metadata index."
            )

        if max_cells is not None and max_cells < len(cells):
            pos = adata_backed.obs_names.get_indexer(cells)
            keep = np.sort(np.argsort(pos, kind="stable")[: int(max_cells)])
            cells = cells[keep]

            if verbose:
                print(
                    f"[testing] Using the first {len(cells):,} cells "
                    f"(max_cells; contiguous prefix)."
                )

        if verbose:
            print(f"\nCells in metadata parquet(s): {len(meta):,}")
            print(f"Cells selected from h5ad: {len(cells):,}")

        all_names = adata_backed.obs_names

        var = adata_backed.var.copy()
        var.index = var.index.astype(str)

        present_genes, col_indices, _ = _resolve_gene_columns(
            genes,
            var.index,
            verbose=verbose,
        )
        var = var.loc[present_genes]

        h5ad_only_cols = [
            col for col in adata_backed.obs.columns
            if col not in meta.columns
        ]

        if h5ad_only_cols:
            h5ad_obs = adata_backed.obs.loc[cells, h5ad_only_cols].copy()
            h5ad_obs.index = h5ad_obs.index.astype(str)
        else:
            h5ad_obs = None

        row_pos = all_names.get_indexer(cells)

        adata_backed.file.close()

        if (row_pos < 0).any():
            raise ValueError(
                "Some requested cells are absent from the h5ad."
            )

        X = _read_csr_rows(
            h5ad_path,
            row_pos,
            col_indices=col_indices,
            verbose=verbose,
        )

        obs_index = pd.Index(cells).astype(str)

        adata = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=obs_index),
            var=var,
        )

        adata = _join_metadata_obs(adata, meta, h5ad_obs=h5ad_obs, verbose=verbose)

    else:
        adata_backed = sc.read_h5ad(h5ad_path, backed="r")
        adata_backed.obs_names = adata_backed.obs_names.astype(str)

        if subset_to_metadata:
            cells = meta.index.intersection(adata_backed.obs_names)
        else:
            cells = adata_backed.obs_names

        if len(cells) == 0:
            adata_backed.file.close()
            raise ValueError(
                "No overlapping cells found between h5ad obs_names and metadata index."
            )

        if max_cells is not None and max_cells < len(cells):
            pos = adata_backed.obs_names.get_indexer(cells)
            keep = np.sort(np.argsort(pos, kind="stable")[: int(max_cells)])
            cells = cells[keep]

            if verbose:
                print(
                    f"[testing] Using the first {len(cells):,} cells "
                    f"(max_cells; contiguous prefix)."
                )

        if verbose:
            print(f"\nCells in metadata parquet(s): {len(meta):,}")
            print(f"Cells selected from h5ad: {len(cells):,}")

        h5ad_only_cols = [
            col for col in _h5ad_obs_column_names(h5ad_path)
            if col not in meta.columns
        ]

        if h5ad_only_cols:

            h5ad_obs = adata_backed.obs.loc[
                cells,
                h5ad_only_cols,
            ].copy()

            h5ad_obs.index = h5ad_obs.index.astype(str)

            obs = h5ad_obs.join(
                meta.loc[cells],
                how="left",
            )

        else:

            obs = meta.loc[cells].copy()

        obs.index = obs.index.astype(str)

        if keep_var:

            var = adata_backed.var.copy()
            var.index = var.index.astype(str)

        else:
            var = None

        adata_backed.file.close()

        if keep_var:
            adata = ad.AnnData(
                X=None,
                obs=obs,
                var=var,
            )
        else:
            adata = ad.AnnData(
                X=None,
                obs=obs,
            )

        adata.uns["__h5ad_path__"] = h5ad_path

    if restore_umaps:
        _restore_umaps(
            adata,
            verbose=verbose,
        )

    if verbose:

        missing = adata.obs[meta.columns].isna().sum()
        missing = missing[missing > 0].sort_values(ascending=False)

        print("\nMetadata columns with missing values:")
        print(missing if len(missing) > 0 else "None")

        print("\nFinal adata:")
        print(adata)

        print("\nIs backed:")
        print(adata.isbacked)

        print("\nHas X:")
        print(adata.X is not None)

        print("\nAvailable UMAPs:")
        print(list(adata.obsm.keys()))

        print("\nAvailable cluster columns:")
        print([
            c for c in adata.obs.columns
            if c.startswith("leiden")
        ])

        print("\nAvailable annotation columns:")
        print([
            c for c in adata.obs.columns
            if "cell_type" in c.lower()
            or "anno" in c.lower()
        ])

    return adata


DEFAULT_GENE_ANNO_PATH = (
    "/mnt/scratch2/Maycon/Bissler_extvis_bRNAseq/Round_1/Objects/gene_anno.csv"
)


def _strip_ensembl_version(gene_id):
    return re.sub(r"\.\d+$", "", str(gene_id))


def _load_gene_biotype_lookup(gene_anno_path=DEFAULT_GENE_ANNO_PATH):
    anno = pd.read_csv(gene_anno_path)
    if "geneSymbol" not in anno.columns or "bioType" not in anno.columns:
        raise ValueError(
            f"gene_anno.csv must contain geneSymbol and bioType columns: {gene_anno_path}"
        )

    anno = anno.copy()
    anno["geneSymbol"] = anno["geneSymbol"].astype(str)
    anno["bioType"] = anno["bioType"].astype(str)

    sym_df = anno.drop_duplicates("geneSymbol", keep="first")
    sym_lookup = dict(zip(sym_df["geneSymbol"], sym_df["bioType"]))

    id_lookup = {}
    if "geneID" in anno.columns:
        anno["geneID_base"] = anno["geneID"].astype(str).map(_strip_ensembl_version)
        id_df = anno.drop_duplicates("geneID_base", keep="first")
        id_lookup = dict(zip(id_df["geneID_base"], id_df["bioType"]))

    return sym_lookup, id_lookup


def _classify_gene_biotype(gene_name, sym_lookup, id_lookup):
    if gene_name in sym_lookup:
        return sym_lookup[gene_name]

    gene_base = _strip_ensembl_version(gene_name)
    if gene_base in id_lookup:
        return id_lookup[gene_base]

    return None


def read_var_names(h5ad_path):
    """
    Cheaply read the ``var`` index (gene names) from an h5ad without loading X.
    """
    with h5py.File(str(h5ad_path), "r") as f:
        var_group = f["var"]
        index_name = var_group.attrs.get("_index", "_index")
        if isinstance(index_name, bytes):
            index_name = index_name.decode()
        raw = var_group[index_name][:]

    return [x.decode() if isinstance(x, bytes) else str(x) for x in raw]


def get_noncoding_genes(
    var_names,
    gene_anno_path=DEFAULT_GENE_ANNO_PATH,
    biotype_exclude=("protein_coding",),
    include_unclassified=False,
    verbose=True,
):
    """
    Return non-coding genes present in ``var_names``, preserving input order.

    A gene is kept when its biotype is matched in ``gene_anno.csv`` and is not
    in ``biotype_exclude``. Set ``include_unclassified=True`` to also keep genes
    with no matching biotype in the annotation table.

    Classification matches each gene by symbol first, then by version-stripped
    Ensembl ID (logic ported from the noncoding_clusters pipeline).
    """
    var_names = [str(g) for g in var_names]
    exclude = {str(x) for x in biotype_exclude}
    sym_lookup, id_lookup = _load_gene_biotype_lookup(gene_anno_path)

    biotypes = [
        _classify_gene_biotype(g, sym_lookup, id_lookup) for g in var_names
    ]

    selected = []
    for gene, biotype in zip(var_names, biotypes):
        if biotype is None:
            if include_unclassified:
                selected.append(gene)
        elif biotype not in exclude:
            selected.append(gene)

    if verbose:
        n_total = len(var_names)
        n_matched = sum(b is not None for b in biotypes)
        n_unclassified = n_total - n_matched
        n_noncoding = len(selected)
        print(
            "Non-coding gene selection: "
            f"total={n_total} matched={n_matched} noncoding={n_noncoding} "
            f"unclassified={n_unclassified} excluded={n_total - n_noncoding} "
            f"include_unclassified={include_unclassified}"
        )

    return selected


# # Usage
# H5AD_PATH = (
#     "/mnt/scratch2/Maycon/Illumina_SpatialData/"
#     "after_May282026/results/prepare_obj/adata_16smp.h5ad"
# )

# # Base metadata (required)
# merged_parquet = (
#     "/mnt/scratch2/Maycon/Illumina_SpatialData/"
#     "after_May282026/results/cell_typing_v3/"
#     "rapids_merged_harmony/merged_harmony/clustering/rapids/"
#     "merged__cellmetadata.parquet"
# )

# # Optional per-sample annotations
# persample_parquets = (
#     "/mnt/scratch2/Maycon/Illumina_SpatialData/"
#     "after_May282026/results/cell_typing_v3/"
#     "scanpy_per_sample_res07/per_sample/annotations"
# )

# adata = Utils.read_adata(
#     h5ad_path=H5AD_PATH,
#     per_sample_column="celltype_from_leiden_unintegrated_r0.7",
#     merged_parquet=merged_parquet,
#     persample_parquets=persample_parquets,
#     load_expression=True,
#     genes=["EPCAM", "KLK3", "AR"],  # or genes="all_genes"
#     gene_panel_h5ad=".../gene_panels/fig1_markers_8genes.h5ad",
# )
