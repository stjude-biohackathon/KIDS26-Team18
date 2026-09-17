"""Shared spatialdata-plot helpers for Xenium and CosMx QC notebooks."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from anndata import AnnData
    from spatialdata import SpatialData

# spatialdata_plot auto-switches to datashader above this threshold; datashader
# misaligns large Xenium transcript overlays (see Dataset_03 DEBUGGING.md).
DATASHADER_AUTO_THRESHOLD = 10_000
DEFAULT_POINTS_RENDER_METHOD = "matplotlib"


def attach_morphology_mip(
    sdata: SpatialData,
    mip_path: str | Path,
    *,
    image_key: str = "morphology_mip",
) -> bool:
    """Attach a 2D morphology MIP to ``sdata`` if not already present.

    Returns True when an image was added, False if it was already present.
    """
    if image_key in sdata.images:
        return False

    import dask.array as da
    import tifffile
    from spatialdata.models import Image2DModel
    from spatialdata.transformations import Identity

    path = Path(mip_path)
    if not path.exists():
        raise FileNotFoundError(f"Morphology MIP not found: {path}")

    mip = tifffile.imread(path)
    arr = mip[None, ...] if mip.ndim == 2 else mip.max(axis=0)[None, ...]
    sdata.images[image_key] = Image2DModel.parse(
        da.from_array(arr),
        dims=("c", "y", "x"),
        transformations={"global": Identity()},
        c_coords=["DAPI"],
    )
    return True


def _gene_expression(adata: AnnData, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


_BASAL_CO_MARKERS = ("KRT5", "KRT14")


def _roi_cell_mask(
    coords: np.ndarray, cx: float, cy: float, half: float
) -> np.ndarray:
    return (
        (coords[:, 0] >= cx - half)
        & (coords[:, 0] <= cx + half)
        & (coords[:, 1] >= cy - half)
        & (coords[:, 1] <= cy + half)
    )


def _candidate_centers(adata: AnnData, gene: str, coords: np.ndarray) -> list[tuple[float, float]]:
    centers: list[tuple[float, float]] = []
    seen: set[tuple[float, float]] = set()

    def _add(g: str) -> None:
        if g not in adata.var_names:
            return
        expr = _gene_expression(adata, g)
        if expr.max() <= 0:
            return
        idx = int(expr.argmax())
        center = (float(coords[idx, 0]), float(coords[idx, 1]))
        if center not in seen:
            seen.add(center)
            centers.append(center)

    _add(gene)
    if gene in _BASAL_CO_MARKERS:
        for marker in _BASAL_CO_MARKERS:
            if marker != gene:
                _add(marker)
    return centers


def _score_crop(
    sdata: SpatialData,
    adata: AnnData,
    gene: str,
    cx: float,
    cy: float,
    half: float,
    obsm_key: str,
) -> tuple[tuple[int, int, float], int, int]:
    """Return sortable score, circle count, and transcript count for a crop."""
    coords = adata.obsm[obsm_key]
    expr = _gene_expression(adata, gene)
    mask = _roi_cell_mask(coords, cx, cy, half)
    crop = sdata.query.bounding_box(
        axes=("x", "y"),
        min_coordinate=[cx - half, cy - half],
        max_coordinate=[cx + half, cy + half],
        target_coordinate_system="global",
        filter_table=False,
    )
    region = adata.uns["spatialdata_attrs"]["region"]
    circles = len(crop.shapes.get(region, []))
    transcripts = len(crop.points.get("transcripts", []))
    score = (circles, transcripts, float(expr[mask].sum()))
    return score, circles, transcripts


def _get_xy_scale(sdata: SpatialData, element_key: str) -> float:
    """Return the isotropic x/y scale factor (raw → global) for a SpatialData element.

    Returns 1.0 if no Scale transform is found.
    """
    from spatialdata.transformations import get_transformation

    try:
        t = get_transformation(sdata[element_key], "global")
    except Exception:
        return 1.0

    t_str = str(t)
    # Parse first numeric value from Scale representation
    for line in t_str.splitlines():
        line = line.strip()
        if line and not line.startswith("Scale") and not line.startswith("["):
            try:
                return float(line.split()[0])
            except (ValueError, IndexError):
                pass
        if line.startswith("["):
            try:
                return float(line.strip("[]").split()[0])
            except (ValueError, IndexError):
                pass
    return 1.0


def _infer_transcript_feature_col(sdata: SpatialData, points_key: str) -> str:
    cols = set(sdata[points_key].columns)
    for candidate in ("feature_name", "target", "gene"):
        if candidate in cols:
            return candidate
    raise ValueError(f"Cannot infer transcript feature column from {sorted(cols)}")


def _infer_instance_key(sdata: SpatialData, points_key: str = "transcripts") -> str:
    """Infer the cell-instance column from table metadata or points columns."""
    table = sdata.get("table")
    if table is not None:
        attrs = table.uns.get("spatialdata_attrs", {})
        if "instance_key" in attrs:
            return str(attrs["instance_key"])

    cols = set(sdata[points_key].columns)
    for candidate in ("cell_id", "cell_uid", "cell_ID"):
        if candidate in cols:
            return candidate
    return "cell_id"


def _points_coordinate_map(columns: list[str] | set[str]) -> dict[str, str]:
    """Build PointsModel coordinate map (x/y required; z optional)."""
    cols = set(columns)
    coords = {"x": "x", "y": "y"}
    if "z" in cols:
        coords["z"] = "z"
    return coords


def transcript_roi_center(
    sdata: SpatialData,
    gene: str,
    *,
    points_key: str = "transcripts",
    transcript_feature_col: str | None = None,
    half: float = 500,
    grid: int = 12,
) -> tuple[float, float]:
    """Return (cx, cy) in **global pixel coordinates** maximizing transcript
    molecule density for ``gene`` inside a half×2 window.

    Transcript x/y are stored in Xenium internal units (pre-Scale). This function
    converts them to global coordinates before grid-searching.
    """
    if transcript_feature_col is None:
        transcript_feature_col = _infer_transcript_feature_col(sdata, points_key)

    scale = _get_xy_scale(sdata, points_key)

    gene_tx = (
        sdata[points_key][sdata[points_key][transcript_feature_col] == gene][["x", "y"]]
        .compute()
    )
    if gene_tx.empty:
        raise ValueError(f"No transcript molecules found for {gene!r}")

    # convert stored coords to global pixel coordinates
    gx = gene_tx["x"] * scale
    gy = gene_tx["y"] * scale

    x_lo, x_hi = float(gx.quantile(0.15)), float(gx.quantile(0.85))
    y_lo, y_hi = float(gy.quantile(0.15)), float(gy.quantile(0.85))
    xs = np.linspace(x_lo, x_hi, grid)
    ys = np.linspace(y_lo, y_hi, grid)

    best_n = -1
    best_center = (float(gx.median()), float(gy.median()))
    for cx in xs:
        for cy in ys:
            n = int(
                ((gx >= cx - half) & (gx <= cx + half)
                 & (gy >= cy - half) & (gy <= cy + half)).sum()
            )
            if n > best_n:
                best_n = n
                best_center = (float(cx), float(cy))
    return best_center


def cell_density_roi_center(
    coords: np.ndarray,
    half: float = 2000,
    grid: int = 20,
) -> tuple[float, float]:
    """Return (cx, cy) maximizing cell count inside a half×2 window.

    Grid-searches the central 70–90% of the x/y extent of ``coords``.
    """
    x = coords[:, 0]
    y = coords[:, 1]
    x_lo, x_hi = float(np.quantile(x, 0.15)), float(np.quantile(x, 0.85))
    y_lo, y_hi = float(np.quantile(y, 0.15)), float(np.quantile(y, 0.85))
    xs = np.linspace(x_lo, x_hi, grid)
    ys = np.linspace(y_lo, y_hi, grid)

    best_n = -1
    best_center = (float(np.median(x)), float(np.median(y)))
    for cx in xs:
        for cy in ys:
            n = int(
                ((x >= cx - half) & (x <= cx + half)
                 & (y >= cy - half) & (y <= cy + half)).sum()
            )
            if n > best_n:
                best_n = n
                best_center = (float(cx), float(cy))
    return best_center


def build_plot_crop(
    sdata: SpatialData,
    cx: float,
    cy: float,
    half: float,
    *,
    gene: str | None = None,
    points_key: str = "transcripts",
    transcript_feature_col: str | None = None,
    instance_key: str | None = None,
    filter_table: bool = False,
) -> SpatialData:
    """BBox-crop ``sdata`` and replace points with a proper spatial transcript filter.

    ``sdata.query.bounding_box`` under-represents Xenium transcript points; this helper
    keeps shapes/images from the bbox crop but re-filters transcripts on x/y (and
    optionally on ``gene``).
    """
    if transcript_feature_col is None:
        transcript_feature_col = _infer_transcript_feature_col(sdata, points_key)
    if instance_key is None:
        instance_key = _infer_instance_key(sdata, points_key)

    crop = sdata.query.bounding_box(
        axes=("x", "y"),
        min_coordinate=[cx - half, cy - half],
        max_coordinate=[cx + half, cy + half],
        target_coordinate_system="global",
        filter_table=filter_table,
    )
    xmin, xmax = cx - half, cx + half
    ymin, ymax = cy - half, cy + half
    tx = sdata[points_key]

    # Transcript x/y are in raw (pre-Scale) units; convert bbox to raw space before filtering
    scale = _get_xy_scale(sdata, points_key)
    raw_xmin, raw_xmax = xmin / scale, xmax / scale
    raw_ymin, raw_ymax = ymin / scale, ymax / scale

    mask = (
        (tx["x"] >= raw_xmin)
        & (tx["x"] <= raw_xmax)
        & (tx["y"] >= raw_ymin)
        & (tx["y"] <= raw_ymax)
    )
    if gene is not None:
        mask = mask & (tx[transcript_feature_col] == gene)
    filtered = tx[mask]

    # Boolean filtering leaves sparse dask-expr partitions (npartitions != non-empty
    # chunks), which breaks spatialdata transform during plotting on large datasets.
    from spatialdata.models import PointsModel
    from spatialdata.transformations import get_transformation

    # reset_index(drop=True): boolean masking preserves original (non-monotonic,
    # duplicate-across-partition) indices, which break dask repartition/reindex on
    # write and on len(). A clean RangeIndex avoids "cannot reindex on an axis with
    # duplicate labels".
    pdf = filtered.compute().reset_index(drop=True)
    transform = get_transformation(sdata[points_key], "global")
    crop.points[points_key] = PointsModel.parse(
        pdf,
        coordinates=_points_coordinate_map(pdf.columns),
        feature_key=transcript_feature_col,
        instance_key=instance_key,
        transformations={"global": transform},
    )
    return crop


def expression_roi_center(
    adata: AnnData,
    gene: str,
    *,
    sdata: SpatialData | None = None,
    half: float = 500,
    obsm_key: str = "spatial",
    method: str = "argmax",
    percentile: float = 90,
) -> tuple[float, float]:
    """Return (cx, cy) for an ROI centered on high expression of ``gene``.

    When ``sdata`` is given, pick among candidate peak centers (the gene itself and,
    for basal keratins, co-markers KRT5/KRT14) using crop coverage and expression.
    """
    expr = _gene_expression(adata, gene)
    coords = adata.obsm[obsm_key]

    if sdata is not None:
        candidates = _candidate_centers(adata, gene, coords)
        if candidates:
            best_center = candidates[0]
            best_score: tuple[int, int, float] | None = None
            for cx, cy in candidates:
                score, _, _ = _score_crop(sdata, adata, gene, cx, cy, half, obsm_key)
                if best_score is None or score > best_score:
                    best_score = score
                    best_center = (cx, cy)
            return best_center

    if method == "argmax" and expr.max() > 0:
        idx = int(expr.argmax())
        return float(coords[idx, 0]), float(coords[idx, 1])

    expressing = expr > 0
    if not expressing.any():
        center = np.median(coords, axis=0)
        return float(center[0]), float(center[1])

    thr = np.percentile(expr[expressing], percentile)
    high = expr >= thr
    weights = expr[high].astype(float)
    cx = float(np.average(coords[high, 0], weights=weights))
    cy = float(np.average(coords[high, 1], weights=weights))
    return cx, cy


def plot_gene_seg_transcripts(
    sdata,
    gene,
    *,
    shapes_key=None,
    points_key="transcripts",
    transcript_feature_col=None,
    coordinate_system="global",
    figsize=(8, 8),
    point_size=1.0,
    point_alpha=0.6,
    palette="orange",
    points_render_method="matplotlib",
):
    """
    Plot cell-level gene expression on segmentation polygons together
    with individual transcript molecules for the same gene.

    Parameters
    ----------
    sdata
        SpatialData object.

    gene
        Gene to visualize.

    shapes_key
        Shape element corresponding to segmented cells.
        If None, infer it from table spatialdata metadata.

    points_key
        SpatialData point element containing transcript molecules.

    transcript_feature_col
        Column identifying transcript target/gene.
        If None, infer from feature_name / target / gene.

    coordinate_system
        Coordinate system used for plotting.

    figsize
        Figure size.

    point_size
        Transcript point size.

    point_alpha
        Transcript point transparency.

    palette
        Color used for transcript molecules.

    points_render_method
        spatialdata-plot point rendering backend.

    Returns
    -------
    Plot object returned by spatialdata-plot.
    """

    import spatialdata_plot  # noqa: F401

    # --------------------------------------------------------
    # Table
    # --------------------------------------------------------

    table = sdata["table"]

    # --------------------------------------------------------
    # Infer shape element
    # --------------------------------------------------------

    if shapes_key is None:

        attrs = table.uns.get(
            "spatialdata_attrs",
            {},
        )

        region = attrs.get("region")

        if isinstance(region, (list, tuple)):
            if len(region) != 1:
                raise ValueError(
                    "Could not uniquely infer shapes element "
                    f"from table region: {region}"
                )
            shapes_key = region[0]

        else:
            shapes_key = region

    if shapes_key is None:
        raise ValueError(
            "Could not infer shapes_key from "
            "table.uns['spatialdata_attrs']."
        )

    if shapes_key not in sdata.shapes:
        raise KeyError(
            f"Shape element {shapes_key!r} "
            "not found in sdata.shapes."
        )

    # --------------------------------------------------------
    # Validate gene
    # --------------------------------------------------------

    if gene not in table.var_names:
        raise KeyError(
            f"Gene {gene!r} not found in "
            "sdata['table'].var_names."
        )

    # --------------------------------------------------------
    # Validate points
    # --------------------------------------------------------

    if points_key not in sdata.points:
        raise KeyError(
            f"Point element {points_key!r} "
            "not found in sdata.points."
        )

    # --------------------------------------------------------
    # Infer transcript feature column
    # --------------------------------------------------------

    if transcript_feature_col is None:

        point_columns = set(
            sdata[points_key].columns
        )

        for candidate in (
            "feature_name",
            "target",
            "gene",
        ):

            if candidate in point_columns:
                transcript_feature_col = candidate
                break

        if transcript_feature_col is None:
            raise ValueError(
                "Could not infer transcript feature column. "
                f"Available columns: {sorted(point_columns)}"
            )

    # --------------------------------------------------------
    # Report what is being plotted
    # --------------------------------------------------------

    print(
        f"Plotting {gene}: "
        f"shapes={shapes_key!r}, "
        f"points={points_key!r}, "
        f"feature={transcript_feature_col!r}"
    )

    # --------------------------------------------------------
    # IMPORTANT:
    #
    # spatialdata-plot 0.3.3 expects the element to be supplied
    # explicitly as the `element=` keyword.
    # --------------------------------------------------------

    plot = (
        sdata.pl.render_shapes(
            element=shapes_key,
            color=gene,
            fill_alpha=0.45,
            outline_alpha=0.8,
            outline_color="white",
        )
        .pl.render_points(
            element=points_key,
            color=transcript_feature_col,
            groups=gene,
            palette=palette,
            size=point_size,
            alpha=point_alpha,
            method=points_render_method,
        )
        .pl.show(
            coordinate_systems=coordinate_system,
            figsize=figsize,
            title=(
                f"{gene}: expression + "
                "segmentation + transcripts"
            ),
        )
    )

    return plot
