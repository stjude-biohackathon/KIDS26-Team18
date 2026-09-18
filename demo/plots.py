"""Plotly helpers for the cell-typing comparison demo."""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

UNASSIGNED_LABEL = "Unassigned"
AMBIGUOUS_LABEL = "Ambiguous"
DEFAULT_MARKER_COLOR = "#7a7a7a"
UNASSIGNED_COLOR = DEFAULT_MARKER_COLOR  # kept for backward compatibility
SPECIAL_LABEL_COLORS = {
    UNASSIGNED_LABEL: DEFAULT_MARKER_COLOR,
    AMBIGUOUS_LABEL: "#c0c0c0",
    "unassigned": DEFAULT_MARKER_COLOR,
    "low_quality": "#a8a8a8",
    "Low quality": "#a8a8a8",
}

# Colorblind-friendly palette for categorical labels
PALETTE = px.colors.qualitative.Safe
PLOTLY_TEMPLATE = "plotly_white"
LEGEND_STYLE = {
    "itemsizing": "constant",
    "itemwidth": 48,
    "font": {"size": 13},
}


def _finalize(fig: go.Figure) -> go.Figure:
    fig.update_layout(template=PLOTLY_TEMPLATE, paper_bgcolor="white", plot_bgcolor="white")
    return fig


def add_legend_toggle_buttons(fig: go.Figure) -> go.Figure:
    """Add Show all / Hide all buttons for Plotly legend traces."""
    if not fig.data:
        return fig
    n_traces = len(fig.data)
    fig.update_layout(
        updatemenus=[
            {
                "type": "buttons",
                "direction": "right",
                "active": -1,
                "x": 1.0,
                "xanchor": "right",
                "y": 1.01,
                "yanchor": "bottom",
                "showactive": False,
                "buttons": [
                    {
                        "label": "Show all",
                        "method": "restyle",
                        "args": [{"visible": [True] * n_traces}],
                    },
                    {
                        "label": "Hide all",
                        "method": "restyle",
                        "args": [{"visible": ["legendonly"] * n_traces}],
                    },
                ],
            }
        ],
    )
    return fig


def _sort_labels(labels: Iterable[str]) -> list[str]:
    uniq = sorted(set(str(x) for x in labels))
    for special in ("low_quality", AMBIGUOUS_LABEL, "unassigned", UNASSIGNED_LABEL):
        if special in uniq:
            uniq.remove(special)
            uniq.append(special)
    return uniq


def _color_map(categories: Iterable[str]) -> dict[str, str]:
    cats = _sort_labels(categories)
    cmap: dict[str, str] = {}
    palette_idx = 0
    for cat in cats:
        if cat in SPECIAL_LABEL_COLORS:
            cmap[cat] = SPECIAL_LABEL_COLORS[cat]
        else:
            cmap[cat] = PALETTE[palette_idx % len(PALETTE)]
            palette_idx += 1
    return cmap


def _y_axis_name(col_idx: int) -> str:
    return "y" if col_idx == 1 else f"y{col_idx}"


def spatial_scatter(
    df: pd.DataFrame,
    color_col: str,
    *,
    title: str,
    point_size: float = 2.0,
    opacity: float = 0.65,
    height: int = 520,
) -> go.Figure:
    plot_df = df.copy()
    plot_df[color_col] = plot_df[color_col].astype(str)
    cmap = _color_map(plot_df[color_col].unique())
    hover_cols = [c for c in ["cell_type_scanpy", "cell_type_rule", "final_CT"] if c in df.columns]
    fig = px.scatter(
        plot_df,
        x="coord_x",
        y="coord_y",
        color=color_col,
        color_discrete_map=cmap,
        category_orders={color_col: _sort_labels(plot_df[color_col].unique())},
        hover_data=hover_cols or None,
        title=title,
        height=height,
        render_mode="webgl",
    )
    fig.update_traces(marker={"size": point_size, "opacity": opacity})
    fig.update_layout(
        xaxis={"title": "x", "scaleanchor": "y", "scaleratio": 1, "showgrid": False},
        yaxis={"title": "y", "showgrid": False, "autorange": "reversed"},
        legend={
            "title": color_col,
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            **LEGEND_STYLE,
        },
        margin={"l": 20, "r": 20, "t": 60, "b": 20},
    )
    return add_legend_toggle_buttons(_finalize(fig))


def spatial_triptych(
    df: pd.DataFrame,
    *,
    panels: list[tuple[str, str]] | None = None,
    point_size: float = 1.5,
    opacity: float = 0.6,
    height: int = 480,
    include_final_ct: bool = False,
    plot_title: str | None = None,
) -> go.Figure:
    if panels is None:
        panels = [
            ("cell_type_scanpy", "Cluster + manual (Scanpy)"),
            ("cell_type_rule", "Rule-based"),
        ]
    if include_final_ct and "final_CT" in df.columns:
        panels = [*panels, ("final_CT", "Published (final_CT)")]

    active_panels = [(col, title) for col, title in panels if col in df.columns]
    ncols = len(active_panels)
    fig = make_subplots(
        rows=1,
        cols=ncols,
        subplot_titles=[title for _, title in active_panels],
        horizontal_spacing=0.03,
    )

    all_labels = pd.concat([df[col].astype(str) for col, _ in active_panels])
    cmap = _color_map(all_labels.unique())
    legend_shown: set[str] = set()

    for col_idx, (col, _) in enumerate(active_panels, start=1):
        for label in _sort_labels(df[col].astype(str).unique()):
            sub = df[df[col].astype(str) == label]
            showlegend = label not in legend_shown
            if showlegend:
                legend_shown.add(label)

            fig.add_trace(
                go.Scattergl(
                    x=sub["coord_x"],
                    y=sub["coord_y"],
                    mode="markers",
                    name=str(label),
                    legendgroup=str(label),
                    showlegend=showlegend,
                    marker={
                        "size": point_size,
                        "opacity": opacity,
                        "color": cmap.get(str(label), DEFAULT_MARKER_COLOR),
                    },
                    hovertemplate=(
                        f"{col}: %{{customdata[0]}}<br>"
                        "x: %{x:.0f}<br>y: %{y:.0f}<extra></extra>"
                    ),
                    customdata=sub[[col]].astype(str).to_numpy(),
                ),
                row=1,
                col=col_idx,
            )

        y_axis = _y_axis_name(col_idx)
        fig.update_yaxes(autorange="reversed", showgrid=False, row=1, col=col_idx)
        fig.update_xaxes(showgrid=False, scaleanchor=y_axis, scaleratio=1, row=1, col=col_idx)

    title = plot_title or (
        "Same tissue — Scanpy vs rule-based"
        if not include_final_ct
        else "Same tissue — three annotation views"
    )
    fig.update_layout(
        height=height,
        title=title,
        legend={
            "orientation": "h",
            "yanchor": "top",
            "y": -0.08,
            "x": 0.5,
            "xanchor": "center",
            **LEGEND_STYLE,
        },
        margin={"l": 10, "r": 10, "t": 70, "b": 90},
    )
    return add_legend_toggle_buttons(_finalize(fig))


def disagreement_scatter(
    df: pd.DataFrame,
    *,
    point_size: float = 2.5,
    opacity: float = 0.75,
    height: int = 520,
) -> go.Figure:
    colors = {"Agree": "#4daf4a", "Disagree": "#e41a1c"}
    hover_cols = [c for c in ["cell_type_scanpy", "cell_type_rule", "final_CT"] if c in df.columns]
    fig = px.scatter(
        df,
        x="coord_x",
        y="coord_y",
        color="agreement",
        color_discrete_map=colors,
        hover_data=hover_cols or ["cell_type_scanpy", "cell_type_rule"],
        title="Method agreement map",
        height=height,
        render_mode="webgl",
    )
    fig.update_traces(marker={"size": point_size, "opacity": opacity})
    fig.update_layout(
        xaxis={"scaleanchor": "y", "scaleratio": 1, "showgrid": False},
        yaxis={"autorange": "reversed", "showgrid": False},
        legend={"title": "Status", **LEGEND_STYLE},
    )
    return add_legend_toggle_buttons(_finalize(fig))


def _aggregate_duplicate_axis_labels(ct: pd.DataFrame) -> pd.DataFrame:
    out = ct.groupby(ct.index).sum()
    return out.T.groupby(out.columns).sum().T


def crosstab_heatmap(ct: pd.DataFrame, *, title: str = "Scanpy × Rule") -> go.Figure:
    ct = _aggregate_duplicate_axis_labels(ct)
    fig = px.imshow(
        ct,
        labels={"x": "Rule-based", "y": "Scanpy", "color": "Cells"},
        x=[str(c) for c in ct.columns],
        y=[str(r) for r in ct.index],
        color_continuous_scale="Blues",
        text_auto=True,
        aspect="auto",
        title=title,
    )
    fig.update_xaxes(tickangle=45, side="bottom")
    fig.update_layout(
        height=max(420, 28 * len(ct.index) + 120),
        margin={"l": 120, "b": 120},
    )
    return _finalize(fig)


def counts_bar(counts: pd.Series, *, title: str) -> go.Figure:
    s = counts.copy()
    s.index = s.index.astype(str)
    s = s.sort_values(ascending=True)

    top = s.tail(20)
    for special in (UNASSIGNED_LABEL, AMBIGUOUS_LABEL):
        if special in s.index and special not in top.index:
            top = pd.concat([s.loc[[special]], top])

    fig = px.bar(
        x=top.values,
        y=top.index.astype(str),
        orientation="h",
        title=title,
        labels={"x": "Cells", "y": "Cell type"},
        height=max(360, 22 * len(top) + 80),
        color=top.index.astype(str),
        color_discrete_map=_color_map(top.index),
    )
    fig.update_layout(showlegend=False, margin={"l": 140})
    return _finalize(fig)


def top_types_with_unassigned(counts: pd.Series, n: int = 12) -> pd.Index:
    counts = counts.astype(str).value_counts()
    top = counts.head(n).index.tolist()
    for special in ("low_quality", UNASSIGNED_LABEL, AMBIGUOUS_LABEL, "unassigned"):
        if special in counts.index and special not in top:
            top.append(special)
    return pd.Index(top)


def subsample_df(
    df: pd.DataFrame,
    *,
    max_points: int,
    stratify_col: str | None = "cell_type_scanpy",
    secondary_stratify_col: str | None = None,
    seed: int = 42,
) -> pd.DataFrame:
    if len(df) <= max_points:
        return df.copy()

    work = df.copy()
    if secondary_stratify_col and stratify_col and stratify_col in work.columns and secondary_stratify_col in work.columns:
        work["_stratify"] = (
            work[stratify_col].astype(str) + " | " + work[secondary_stratify_col].astype(str)
        )
        stratify_col = "_stratify"
    elif "cell_type_rule" in work.columns and "cell_type_scanpy" in work.columns:
        work["_stratify"] = (
            work["cell_type_scanpy"].astype(str) + " | " + work["cell_type_rule"].astype(str)
        )
        stratify_col = "_stratify"

    rng = np.random.default_rng(seed)
    if stratify_col and stratify_col in work.columns:
        parts = []
        for _, group in work.groupby(stratify_col, observed=True):
            n = max(1, int(round(max_points * len(group) / len(work))))
            n = min(n, len(group))
            parts.append(group.sample(n=n, random_state=int(rng.integers(0, 2**31 - 1))))
        out = pd.concat(parts, ignore_index=True)
        if stratify_col == "_stratify":
            out = out.drop(columns="_stratify", errors="ignore")
        if len(out) > max_points:
            out = out.sample(n=max_points, random_state=seed)
        return out.reset_index(drop=True)

    return work.sample(n=max_points, random_state=seed).reset_index(drop=True)
