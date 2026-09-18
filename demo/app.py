"""
Streamlit demo: cell typing method comparison.

Run from repo root:
    conda activate spatialdata
    streamlit run demo/app.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from demo.cluster_plots import (
    plot_ncount_violin,
    plot_rank_genes_dotplot,
    plot_umap_leiden_vs_scanpy,
    unassigned_ncount_qc,
)
from demo.datasets import (
    DATASETS,
    DatasetSpec,
    assign_rule_labels,
    cluster_annotations_table,
    get_dataset,
    load_celltyping_config,
    load_demo_adata,
    resolve_rule_panels,
    unassigned_column,
    uses_rule_engine,
)
from demo.labels import format_label_display, normalize_cell_label
from demo.metrics import compare_methods
from demo.plots import (
    counts_bar,
    crosstab_heatmap,
    disagreement_scatter,
    spatial_scatter,
    spatial_triptych,
    subsample_df,
)


def _h5ad_sig(path: Path) -> str:
    stat = path.stat()
    return f"{path}:{stat.st_mtime}:{stat.st_size}"


@st.cache_data(show_spinner="Loading cell data…")
def load_demo_adata_cached(dataset_id: str) -> ad.AnnData:
    return load_demo_adata(dataset_id)


@st.cache_data(show_spinner="Preparing marker panels…")
def load_rule_panel_bundle(dataset_id: str, _h5ad_sig: str) -> tuple[dict, dict, dict]:
    spec = get_dataset(dataset_id)
    adata = load_demo_adata_cached(dataset_id)
    markers, panel_min_hits, rule_cfg = resolve_rule_panels(adata, spec)
    return markers, panel_min_hits, rule_cfg


@st.cache_data(show_spinner="Rendering UMAP…")
def cached_umap_figure(dataset_id: str, _h5ad_sig: str) -> bytes:
    adata = load_demo_adata_cached(dataset_id)
    fig = plot_umap_leiden_vs_scanpy(adata)
    return _fig_to_png(fig)


@st.cache_data(show_spinner="Rendering marker dotplot…")
def cached_dotplot_figure(dataset_id: str, _h5ad_sig: str, n_genes: int) -> bytes:
    adata = load_demo_adata_cached(dataset_id)
    fig = plot_rank_genes_dotplot(adata, n_genes=n_genes)
    return _fig_to_png(fig)


def _fig_to_png(fig: plt.Figure) -> bytes:
    import io

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def obs_to_frame(adata: ad.AnnData, spec: DatasetSpec) -> pd.DataFrame:
    cols = ["coord_x", "coord_y", spec.method_a_col, spec.method_b_col]
    if "final_CT" in adata.obs.columns:
        cols.append("final_CT")

    df = adata.obs[[c for c in cols if c in adata.obs.columns or c in ("coord_x", "coord_y")]].copy()
    if "coord_x" not in df.columns and "spatial" in adata.obsm:
        df["coord_x"] = adata.obsm["spatial"][:, 0]
        df["coord_y"] = adata.obsm["spatial"][:, 1]

    return df.reset_index(names="cell_id")


def default_spatial_ranges(
    adata: ad.AnnData,
    spec: DatasetSpec,
) -> tuple[tuple[float, float], tuple[float, float]]:
    x_min = float(adata.obs["coord_x"].min())
    x_max = float(adata.obs["coord_x"].max())
    y_min = float(adata.obs["coord_y"].min())
    y_max = float(adata.obs["coord_y"].max())

    if spec.config_path and spec.config_path.is_file():
        cfg = load_celltyping_config(spec)
        crop = cfg.get("spatial_crop", {})
        if crop.get("enabled"):
            x_min = max(x_min, float(crop["x_min"]))
            x_max = min(x_max, float(crop["x_max"]))
            y_min = max(y_min, float(crop["y_min"]))
            y_max = min(y_max, float(crop["y_max"]))

    return (x_min, x_max), (y_min, y_max)


def crop_adata(adata: ad.AnnData, x_range: tuple[float, float], y_range: tuple[float, float]) -> ad.AnnData:
    mask = adata.obs["coord_x"].between(x_range[0], x_range[1]) & adata.obs["coord_y"].between(
        y_range[0], y_range[1]
    )
    return adata[mask].copy()


def agreement_series(df: pd.DataFrame, col_a: str, col_b: str) -> pd.Series:
    return pd.Series(
        np.where(df[col_a].astype(str) == df[col_b].astype(str), "Agree", "Disagree"),
        index=df.index,
    )


def prepare_comparison_columns(df: pd.DataFrame, spec: DatasetSpec) -> tuple[pd.DataFrame, str, str]:
    """Add canonical comparison columns when label normalization is enabled."""
    if not spec.label_normalization:
        return df, spec.method_a_col, spec.method_b_col

    out = df.copy()
    cmp_a = f"{spec.method_a_col}__cmp"
    cmp_b = f"{spec.method_b_col}__cmp"
    out[cmp_a] = out[spec.method_a_col].astype(str).map(normalize_cell_label)
    out[cmp_b] = out[spec.method_b_col].astype(str).map(normalize_cell_label)
    return out, cmp_a, cmp_b


def apply_display_labels(df: pd.DataFrame, spec: DatasetSpec) -> pd.DataFrame:
    """Apply canonical display labels so matched types share the same legend name."""
    if not spec.label_normalization:
        return df

    cmp_a = f"{spec.method_a_col}__cmp"
    cmp_b = f"{spec.method_b_col}__cmp"
    out = df.copy()
    out[spec.method_a_col] = [
        format_label_display(raw, canonical=norm)
        for raw, norm in zip(out[spec.method_a_col].astype(str), out[cmp_a].astype(str))
    ]
    out[spec.method_b_col] = [
        format_label_display(raw, canonical=norm)
        for raw, norm in zip(out[spec.method_b_col].astype(str), out[cmp_b].astype(str))
    ]
    return out


def format_crosstab_labels(ct: pd.DataFrame) -> pd.DataFrame:
    out = ct.copy()
    out.index = [format_label_display(str(x), canonical=str(x)) for x in out.index]
    out.columns = [format_label_display(str(x), canonical=str(x)) for x in out.columns]
    out = out.groupby(out.index).sum()
    return out.T.groupby(out.columns).sum().T


AGREEMENT_HELP = (
    "Agreement is the percentage of cells in the current view where both methods assign "
    "the same cell type after label normalization (case, plurals, and cell/cells suffixes). "
    "It is a strict label match, not a biological similarity score — 100% means every cell "
    "received identical type names from both methods; lower values mean more naming disagreements."
)


def count_unassigned(df: pd.DataFrame, col: str, unassigned_labels: tuple[str, ...]) -> tuple[int, float]:
    labels = df[col].astype(str)
    n_unassigned = int(labels.isin(unassigned_labels).sum())
    pct_unassigned = 100.0 * n_unassigned / len(df) if len(df) else 0.0
    return n_unassigned, pct_unassigned


def spatial_panels(spec: DatasetSpec) -> list[tuple[str, str]]:
    return [
        (spec.method_a_col, spec.method_a_label),
        (spec.method_b_col, spec.method_b_label),
    ]


def main() -> None:
    st.set_page_config(
        page_title="Cell typing comparison",
        page_icon="🧬",
        layout="wide",
    )

    st.title("Cell typing comparison")

    dataset_options = list(DATASETS.keys())
    default_dataset = "dataset01" if "dataset01" in DATASETS else dataset_options[0]

    with st.sidebar:
        st.header("Dataset")
        dataset_id = st.selectbox(
            "Select dataset",
            dataset_options,
            format_func=lambda key: DATASETS[key].label,
            index=dataset_options.index(default_dataset),
        )

    try:
        spec = get_dataset(dataset_id)
        adata = load_demo_adata_cached(dataset_id)
    except FileNotFoundError as exc:
        st.error(str(exc))
        st.stop()

    h5ad_sig = _h5ad_sig(spec.clustered_h5ad)
    has_rules = uses_rule_engine(spec)
    marker_panels: dict[str, list[str]] = {}
    panel_min_hits: dict[str, int] = {}
    rule_cfg: dict = {"min_hits": 1, "count_threshold": 0.0, "max_hits_slider": 4}
    if has_rules:
        marker_panels, panel_min_hits, rule_cfg = load_rule_panel_bundle(dataset_id, h5ad_sig)

    has_leiden = has_rules and "leiden" in adata.obs.columns and "X_umap" in adata.obsm
    has_published = "final_CT" in adata.obs.columns
    comparison_label = f"{spec.method_a_label} vs {spec.method_b_label}"

    st.caption(spec.caption)

    (def_x_min, def_x_max), (def_y_min, def_y_max) = default_spatial_ranges(adata, spec)
    full_x_min, full_x_max = float(adata.obs["coord_x"].min()), float(adata.obs["coord_x"].max())
    full_y_min, full_y_max = float(adata.obs["coord_y"].min()), float(adata.obs["coord_y"].max())
    x_pad = (full_x_max - full_x_min) * 0.02
    y_pad = (full_y_max - full_y_min) * 0.02

    color_options = ["Side-by-side", spec.method_a_label, spec.method_b_label, "Agreement"]
    if has_published:
        color_options.insert(3, "Published")

    max_hits = int(rule_cfg.get("max_hits_slider", 4))
    default_min_hits = int(rule_cfg.get("min_hits", 1))
    default_count_threshold = int(rule_cfg.get("count_threshold", 0))

    with st.sidebar:
        st.header("Controls")

        view_mode = st.radio("Color by", color_options, index=0)

        st.subheader("Spatial crop")
        x_range = st.slider(
            "X range",
            full_x_min - x_pad,
            full_x_max + x_pad,
            (def_x_min, def_x_max),
            step=10.0,
        )
        y_range = st.slider(
            "Y range",
            full_y_min - y_pad,
            full_y_max + y_pad,
            (def_y_min, def_y_max),
            step=10.0,
        )

        st.subheader("Display")
        max_points = st.slider("Max points plotted", 5_000, 50_000, 20_000, step=1_000)
        point_size = st.slider("Point size", 1.0, 6.0, 2.0, step=0.5)
        opacity = st.slider("Opacity", 0.2, 1.0, 0.65, step=0.05)
        show_final_ct = False
        if has_published:
            show_final_ct = st.checkbox(
                "Show published (final_CT)",
                value=False,
                help="Include the published reference panel in the side-by-side spatial view.",
            )

        if has_rules:
            st.subheader("Rule-based tuning")
            use_tuned_rules = st.checkbox("Recompute rule labels live", value=False)
            min_hits = st.slider("Min marker hits", 1, max_hits, default_min_hits)
            count_threshold = st.slider("Count threshold", 0, 10, default_count_threshold)
        else:
            use_tuned_rules = False
            min_hits = default_min_hits
            count_threshold = default_count_threshold

        filter_disagree = st.checkbox("Show disagreements only", value=False)
        cell_type_filter = st.multiselect(
            f"Filter {spec.method_a_label} types",
            sorted(adata.obs[spec.method_a_col].astype(str).unique()),
        )

        st.divider()
        st.caption(f"Data: `{spec.clustered_h5ad.name}`")
        if spec.reference_h5ad and spec.reference_h5ad.is_file() and has_published:
            st.caption(f"Reference: `{spec.reference_h5ad.name}` (final_CT)")
        if spec.config_path and spec.config_path.is_file():
            st.caption(f"Config: `{spec.config_path.name}`")
        st.caption(f"{adata.n_obs:,} cells × {adata.n_vars:,} genes")
        if has_leiden:
            st.caption(f"{adata.obs['leiden'].nunique()} Leiden clusters — see **Clustering** tab")

    adata_view = crop_adata(adata, x_range, y_range)
    df = obs_to_frame(adata_view, spec)

    if has_rules:
        if use_tuned_rules:
            rule_labels = assign_rule_labels(
                adata_view,
                spec,
                marker_panels,
                panel_min_hits,
                min_hits=min_hits,
                count_threshold=float(count_threshold),
                use_config_panels=False,
            )
            adata_view.obs[spec.method_b_col] = rule_labels
            df[spec.method_b_col] = rule_labels
        elif spec.method_b_col not in df.columns:
            rule_labels = assign_rule_labels(
                adata_view,
                spec,
                marker_panels,
                panel_min_hits,
                min_hits=min_hits,
                count_threshold=float(count_threshold),
                use_config_panels=True,
            )
            adata_view.obs[spec.method_b_col] = rule_labels
            df[spec.method_b_col] = rule_labels

    df, cmp_a_col, cmp_b_col = prepare_comparison_columns(df, spec)
    df["agreement"] = agreement_series(df, cmp_a_col, cmp_b_col)

    if filter_disagree:
        df = df[df["agreement"] == "Disagree"].copy()
    if cell_type_filter:
        df = df[df[spec.method_a_col].astype(str).isin(cell_type_filter)].copy()

    if df.empty:
        st.warning("No cells match the current filters. Widen the spatial crop or clear filters.")
        st.stop()

    df = apply_display_labels(df, spec)

    df_plot = subsample_df(
        df,
        max_points=max_points,
        stratify_col=spec.method_a_col,
        secondary_stratify_col=spec.method_b_col,
    )
    metrics = compare_methods(
        df,
        method_a_key=cmp_a_col,
        method_b_key=cmp_b_col,
        method_a_label=spec.method_a_label,
        method_b_label=spec.method_b_label,
    )

    unassigned_col = unassigned_column(spec)
    unassigned_label_name = spec.method_a_label if spec.unassigned_on == "method_a" else spec.method_b_label
    unassigned_cmp_col = (
        f"{unassigned_col}__cmp"
        if spec.label_normalization and f"{unassigned_col}__cmp" in df.columns
        else unassigned_col
    )
    unassigned_values = tuple(normalize_cell_label(x) for x in spec.unassigned_labels)
    n_unassigned, pct_unassigned = count_unassigned(df, unassigned_cmp_col, unassigned_values)

    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Cells in view", f"{len(df):,}")
    m2.metric("Agreement", f"{metrics['pct_agree']:.1f}%", help=AGREEMENT_HELP)
    m3.metric("Disagreements", f"{metrics['n_disagree']:,}")
    summary_row = metrics["summary"].iloc[0]
    m4.metric(f"Cohen's κ ({comparison_label})", f"{summary_row['cohen_kappa']:.3f}")
    m5.metric(
        f"Unassigned ({unassigned_label_name})",
        f"{pct_unassigned:.1f}%",
        help=f"{n_unassigned:,} cells in view",
    )

    panels = spatial_panels(spec)
    plot_title = f"Same tissue — {comparison_label}"

    if view_mode == "Side-by-side":
        if show_final_ct and "final_CT" not in df_plot.columns:
            st.warning("Published `final_CT` not available — showing the two main panels only.")
        st.plotly_chart(
            spatial_triptych(
                df_plot,
                panels=panels,
                point_size=point_size,
                opacity=opacity,
                include_final_ct=show_final_ct,
                plot_title=plot_title,
            ),
            width="stretch",
        )
    elif view_mode == "Agreement":
        st.plotly_chart(
            disagreement_scatter(df_plot, point_size=point_size + 0.5, opacity=opacity),
            width="stretch",
        )
    else:
        col_map = {
            spec.method_a_label: spec.method_a_col,
            spec.method_b_label: spec.method_b_col,
            "Published": "final_CT",
        }
        color_col = col_map[view_mode]
        if color_col not in df_plot.columns:
            st.warning(f"Column `{color_col}` not available in this file.")
        else:
            st.plotly_chart(
                spatial_scatter(
                    df_plot,
                    color_col,
                    title=f"Spatial — {view_mode}",
                    point_size=point_size,
                    opacity=opacity,
                ),
                width="stretch",
            )

    tab_names = ["Metrics", "Crosstab", "Type counts"]
    if has_leiden:
        tab_names = ["Clustering", *tab_names, "Marker panels"]
    elif spec.dataset_id == "dataset05":
        tab_names = [*tab_names, "Hierarchy"]
    else:
        if has_rules:
            tab_names = [*tab_names, "Marker panels"]

    tabs = st.tabs(tab_names)
    tab_idx = 0

    if has_leiden:
        with tabs[tab_idx]:
            st.markdown(
                "Notebook-style cluster views on the **full sample** (not spatial crop). "
                "Leiden clusters are annotated manually via cluster annotations in the notebook/config."
            )
            try:
                st.image(
                    cached_umap_figure(dataset_id, h5ad_sig),
                    caption="UMAP — Leiden clusters vs manual Scanpy labels",
                )
            except ValueError as exc:
                st.error(str(exc))

            n_genes = st.slider("Genes per cluster (dotplot)", 3, 10, 5, key="cluster_n_genes")
            try:
                st.image(
                    cached_dotplot_figure(dataset_id, h5ad_sig, n_genes),
                    caption=f"Top {n_genes} DE genes per Leiden cluster",
                )
            except ValueError as exc:
                st.error(str(exc))

            st.markdown("**Manual cluster annotations** (from notebook/config):")
            st.dataframe(cluster_annotations_table(spec), width="stretch", hide_index=True)
        tab_idx += 1

    with tabs[tab_idx]:
        st.dataframe(metrics["summary"], width="stretch", hide_index=True)
        if spec.label_normalization:
            st.info(
                "Agreement metrics normalize labels before comparison (case, plurals, and "
                "'cell'/'cells' suffixes) so formatting differences are not counted as disagreements."
            )
        if has_published:
            st.info(
                "Published `final_CT` uses finer-grained labels than either demo method, "
                "so accuracy vs truth is a loose reference — focus on spatial patterns and method disagreement."
            )
    tab_idx += 1

    with tabs[tab_idx]:
        crosstab = format_crosstab_labels(metrics["crosstab"])
        st.plotly_chart(crosstab_heatmap(crosstab), width="stretch")
    tab_idx += 1

    with tabs[tab_idx]:
        c1, c2 = st.columns(2)
        with c1:
            st.plotly_chart(
                counts_bar(
                    df[spec.method_a_col].astype(str).value_counts(),
                    title=f"{spec.method_a_label} labels",
                ),
                width="stretch",
            )
        with c2:
            st.plotly_chart(
                counts_bar(
                    df[spec.method_b_col].astype(str).value_counts(),
                    title=f"{spec.method_b_label} labels",
                ),
                width="stretch",
            )

        count_col = spec.count_col
        st.markdown(f"**Transcript counts (`{count_col}`) by cell type** — current spatial view")
        violin_options = [spec.method_b_col, spec.method_a_col]
        violin_groupby = st.radio(
            "Violin groupby",
            violin_options,
            format_func=lambda x: spec.method_b_label if x == spec.method_b_col else spec.method_a_label,
            horizontal=True,
            key=f"violin_groupby_{dataset_id}",
        )
        if count_col not in adata_view.obs.columns:
            st.warning(f"`{count_col}` not found in cell metadata.")
        else:
            adata_violin = adata_view[df["cell_id"].astype(str).tolist()].copy()
            label_lookup = df.set_index("cell_id")
            for col in violin_options:
                if col in label_lookup.columns:
                    adata_violin.obs[col] = (
                        adata_violin.obs_names.astype(str).map(label_lookup[col].astype(str))
                    )
            if violin_groupby not in adata_violin.obs.columns:
                st.warning(f"Column `{violin_groupby}` not available — pick another grouping.")
            else:
                adata_violin.obs[violin_groupby] = adata_violin.obs[violin_groupby].astype(str)
                try:
                    fig_violin = plot_ncount_violin(adata_violin, key=count_col, groupby=violin_groupby)
                    st.pyplot(fig_violin, clear_figure=True)
                except ValueError as exc:
                    st.error(str(exc))

                st.markdown(f"**Unassigned cells ({unassigned_label_name})** — low `{count_col}` QC")
                ncount_threshold = st.slider(
                    f"{count_col} threshold",
                    min_value=1,
                    max_value=20,
                    value=5,
                    key=f"unassigned_ncount_threshold_{dataset_id}",
                    help=f"Flag unassigned cells with transcript count below this value.",
                )
                try:
                    unassigned_labels_for_qc = (
                        tuple(format_label_display(x) for x in spec.unassigned_labels)
                        if spec.label_normalization
                        else spec.unassigned_labels
                    )
                    fig_unassigned, n_unassigned_qc, pct_below = unassigned_ncount_qc(
                        adata_violin,
                        key=count_col,
                        rule_key=unassigned_col,
                        unassigned_labels=unassigned_labels_for_qc,
                        threshold=float(ncount_threshold),
                    )
                    if n_unassigned_qc == 0:
                        st.info("No unassigned cells in the current view.")
                    else:
                        m_u, m_p = st.columns(2)
                        m_u.metric("Unassigned cells", f"{n_unassigned_qc:,}")
                        m_p.metric(
                            f"With {count_col} < {ncount_threshold}",
                            f"{pct_below:.2f}%",
                        )
                        if fig_unassigned is not None:
                            st.pyplot(fig_unassigned, clear_figure=True)
                except ValueError as exc:
                    st.error(str(exc))
    tab_idx += 1

    if has_leiden:
        with tabs[tab_idx]:
            st.markdown("Marker panels used for rule assignment:")
            panel_rows = []
            for cell_type, genes in marker_panels.items():
                row = {"cell_type": cell_type, "markers": ", ".join(genes)}
                if spec.rule_engine == "dataset04":
                    row["min_hits"] = panel_min_hits.get(cell_type, rule_cfg.get("min_hits", 1))
                panel_rows.append(row)
            st.dataframe(pd.DataFrame(panel_rows), width="stretch", hide_index=True)
    elif spec.dataset_id == "dataset05":
        with tabs[tab_idx]:
            st.markdown("Hierarchical annotation breakdown for the current spatial view.")
            hier_cols = [c for c in ["hier_level1", "hier_level2", "hier_level3"] if c in adata_view.obs.columns]
            for level_col in hier_cols:
                st.markdown(f"**{level_col}**")
                st.dataframe(
                    adata_view.obs.loc[df["cell_id"].astype(str), level_col]
                    .astype(str)
                    .value_counts()
                    .rename("cells")
                    .reset_index()
                    .rename(columns={"index": level_col}),
                    width="stretch",
                    hide_index=True,
                )

            if "hier_score" in adata_view.obs.columns:
                st.markdown("**hier_label summary**")
                summary = (
                    adata_view.obs.loc[df["cell_id"].astype(str), ["hier_label", "hier_depth", "hier_score", "hier_margin"]]
                    .astype({"hier_label": str, "hier_depth": float, "hier_score": float, "hier_margin": float})
                    .groupby("hier_label", observed=True)
                    .agg(cells=("hier_label", "size"), mean_depth=("hier_depth", "mean"), mean_score=("hier_score", "mean"), mean_margin=("hier_margin", "mean"))
                    .sort_values("cells", ascending=False)
                    .reset_index()
                )
                st.dataframe(summary, width="stretch", hide_index=True)
    elif has_rules:
        with tabs[tab_idx]:
            st.markdown("Marker panels used for rule assignment:")
            panel_rows = [{"cell_type": ct, "markers": ", ".join(genes)} for ct, genes in marker_panels.items()]
            st.dataframe(pd.DataFrame(panel_rows), width="stretch", hide_index=True)


if __name__ == "__main__":
    main()
