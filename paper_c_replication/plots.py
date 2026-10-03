"""Paper plotting functions, parameterized by supplied data only."""
from __future__ import annotations
from pathlib import Path
from typing import Any, Mapping
import math
import tempfile
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from .statistics import official_deciles, spearman

def separate_labels(fig):
    """Keep close endpoint annotations readable for arbitrary supplied data."""
    from matplotlib.text import Annotation
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for axis in fig.axes:
        groups={}
        for item in axis.texts:
            if isinstance(item,Annotation) and item.anncoords=='offset points':
                groups.setdefault(float(item.xy[0]),[]).append(item)
        for items in groups.values():
            top=-float('inf')
            for item in sorted(items,key=lambda x:x.get_window_extent(renderer).y0):
                box=item.get_window_extent(renderer)
                shift=max(0,top+4-box.y0)
                if shift:
                    x,y=item.get_position();item.set_position((x,y+shift*72/fig.dpi))
                top=box.y1+shift

UNIVERSE_ORDER = ("J", "N", "A")


TREATMENT_ORDER = ("RAW", "REF")


INPUT_CELLS = ("Raw_J", "Filtered_J", "Raw_N", "Filtered_N", "Raw_OA", "Filtered_OA")


CELL_COLUMNS = {
    "Raw_J": ("J", "RAW", "ais_raw_j"),
    "Filtered_J": ("J", "REF", "ais_filtered_j"),
    "Raw_N": ("N", "RAW", "ais_raw_n"),
    "Filtered_N": ("N", "REF", "ais_filtered_n"),
    "Raw_OA": ("A", "RAW", "ais_raw_oa"),
    "Filtered_OA": ("A", "REF", "ais_filtered_oa"),
}


GROUP_COLUMNS = (
    "year", "support_scope", "support_id", "cell_id", "universe", "treatment",
    "official_ais_decile", "oa_field", "stable_n_ge_300", "n",
    "official_ais_mean", "mae", "scaled_mae", "rmse", "nrmse",
    "signed_mean_error", "absolute_error_sum", "squared_error_sum",
)


DISPLAY_CELL_ORDER = INPUT_CELLS


DISPLAY_CELL_LABELS = {
    "Raw_J": r"$\mathit{J}$-RAW",
    "Filtered_J": r"$\mathit{J}$-REF",
    "Raw_N": r"$\mathit{N}$-RAW",
    "Filtered_N": r"$\mathit{N}$-REF",
    "Raw_OA": r"$\mathit{A}$-RAW",
    "Filtered_OA": r"$\mathit{A}$-REF",
}


def _require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    if missing := sorted(columns - set(frame.columns)):
        raise AssertionError(f"{label} misses required columns: {missing}")


def _save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="wb", prefix=f".{path.stem}.", suffix=path.suffix,
        dir=path.parent, delete=False,
    )
    temporary = Path(handle.name)
    handle.close()
    try:
        fig.savefig(temporary, dpi=240, bbox_inches="tight", facecolor="white")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
        plt.close(fig)


def _nice_axis_ceiling(maximum: float) -> float:
    if not math.isfinite(maximum) or maximum <= 0:
        raise AssertionError("AIS axis bound must be positive and finite")
    step = 10.0 ** math.floor(math.log10(maximum / 8.0))
    return float(math.ceil(maximum / step * 1.02) * step)


def _error_group_record(
    group: pd.DataFrame, score_column: str, *, year: int, cell_id: str,
    official_ais_decile: str | None = None, oa_field: str | None = None,
    stable_n_ge_300: bool | None = None,
) -> dict[str, Any]:
    official = pd.to_numeric(group["official_ais"], errors="raise").to_numpy(float)
    opened = pd.to_numeric(group[score_column], errors="raise").to_numpy(float)
    if not np.isfinite(official).all() or not np.isfinite(opened).all():
        raise AssertionError("Error-profile group contains a non-finite score")
    mean_official = float(official.mean())
    if mean_official <= 0:
        raise AssertionError("Scaled error requires positive mean official AIS")
    error = opened - official
    absolute = np.abs(error)
    squared = np.square(error)
    universe, treatment, _ = CELL_COLUMNS[cell_id]
    return {
        "year": year,
        "support_scope": "annual_s6",
        "support_id": "s6",
        "cell_id": cell_id,
        "universe": universe,
        "treatment": treatment,
        "official_ais_decile": official_ais_decile,
        "oa_field": oa_field,
        "stable_n_ge_300": stable_n_ge_300,
        "n": int(len(group)),
        "official_ais_mean": mean_official,
        "mae": float(absolute.mean()),
        "scaled_mae": float(absolute.mean() / mean_official),
        "rmse": float(np.sqrt(squared.mean())),
        "nrmse": float(np.sqrt(squared.mean()) / mean_official),
        "signed_mean_error": float(error.mean()),
        "absolute_error_sum": float(absolute.sum()),
        "squared_error_sum": float(squared.sum()),
    }


def prepare_group_error_tables(
    annual_panels: Mapping[int, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return aggregate decile and field error profiles for 2023 and 2024."""

    decile_rows: list[dict[str, Any]] = []
    field_rows: list[dict[str, Any]] = []
    for year, panel in sorted(annual_panels.items()):
        work = panel.copy()
        work["official_ais_decile"] = official_deciles(work["official_ais"])
        work["oa_field"] = work["oa_field"].astype("string").fillna("Unknown")
        for cell_id, (_, _, score_column) in CELL_COLUMNS.items():
            for number in range(1, 11):
                label = f"D{number}"
                group = work.loc[work["official_ais_decile"].eq(label)]
                if group.empty:
                    raise AssertionError(f"{year} {cell_id} {label} is empty")
                decile_rows.append(
                    _error_group_record(
                        group, score_column, year=year, cell_id=cell_id,
                        official_ais_decile=label,
                    )
                )
            for field, group in work.groupby("oa_field", sort=True, dropna=False):
                field_rows.append(
                    _error_group_record(
                        group, score_column, year=year, cell_id=cell_id,
                        oa_field=str(field), stable_n_ge_300=bool(len(group) >= 300),
                    )
                )
    deciles = pd.DataFrame(decile_rows).loc[:, GROUP_COLUMNS]
    fields = pd.DataFrame(field_rows).loc[:, GROUP_COLUMNS]
    if deciles.duplicated(["year", "cell_id", "official_ais_decile"]).any():
        raise AssertionError("Decile aggregate key is duplicated")
    if fields.duplicated(["year", "cell_id", "oa_field"]).any():
        raise AssertionError("Field aggregate key is duplicated")
    for year, panel in annual_panels.items():
        expected_n = int(len(panel))
        decile_n = deciles.loc[deciles["year"].eq(year)].groupby("cell_id")["n"].sum()
        if not decile_n.eq(expected_n).all():
            raise AssertionError(f"{year} decile groups do not conserve S6 rows")
        field_n = fields.loc[fields["year"].eq(year)].groupby("cell_id")["n"].sum()
        if not field_n.eq(expected_n).all():
            raise AssertionError(f"{year} field groups do not conserve S6 rows")
    return deciles, fields


def plot_headline_connected(frame: pd.DataFrame, path: Path, *, year: int = 2024) -> dict[str, Any]:
    """Draw four connected RAW-to-REF aggregate metrics."""

    data = frame.loc[frame["year"].eq(year)].copy()
    if len(data) != 6:
        raise AssertionError("Connected metric figure requires six annual cells")
    colors = {"J": "#1769aa", "N": "#2a9d8f", "A": "#7a5195"}
    specs = (
        ("mae", "Mean absolute error", "{:.3f}"),
        ("rmse", "Root mean squared error", "{:.3f}"),
        ("spearman", r"Spearman $\rho$", "{:.3f}"),
        ("signed_mean_error", "Signed mean error", "{:+.3f}"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(9.8, 7.4))
    for axis, (metric, label, number_format) in zip(axes.flat, specs):
        for universe in UNIVERSE_ORDER:
            group = data.loc[data["universe"].eq(universe)].copy()
            group["x"] = group["treatment"].map({"RAW": 0, "REF": 1})
            group = group.sort_values("x")
            axis.plot(
                group["x"], group[metric], marker="o", markersize=6.5,
                linewidth=1.9, color=colors[universe], label=rf"$\mathit{{{universe}}}$",
            )
            for row in group.itertuples(index=False):
                x = 0 if row.treatment == "RAW" else 1
                value = float(getattr(row, metric))
                offset = {"J": -12, "N": 7, "A": 7}[universe]
                axis.annotate(
                    number_format.format(value), (x, value), xytext=(0, offset),
                    textcoords="offset points", ha="center", fontsize=7.4,
                )
        axis.set_xlim(-0.18, 1.18)
        axis.set_xticks([0, 1], ["RAW", "REF"])
        axis.set_ylabel(label)
        axis.grid(axis="y", color="#dddddd", linewidth=0.55)
        if metric == "signed_mean_error":
            axis.axhline(0, color="#555555", linewidth=0.8)
            minimum = float(data[metric].min())
            maximum = float(data[metric].max())
            padding = max(0.035, (maximum - minimum) * 0.24)
            axis.set_ylim(minimum - padding, maximum + padding)
        if metric == "spearman":
            lower = max(0.0, float(data[metric].min()) - 0.025)
            upper = min(1.0, float(data[metric].max()) + 0.012)
            axis.set_ylim(lower, upper)
        elif metric in {"mae", "rmse"}:
            axis.set_ylim(0, float(data[metric].max()) * 1.15)
    handles = [
        Line2D([0], [0], color=colors[u], marker="o", label=rf"$\mathit{{{u}}}$")
        for u in UNIVERSE_ORDER
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, title="Universe")
    fig.suptitle(f"{year} aggregate agreement with AIS")
    fig.tight_layout(rect=(0, 0.08, 1, 0.96))
    separate_labels(fig)
    _save_figure(fig, path)
    return {"year": year, "panels": [item[0] for item in specs], "n_points": 24}


def plot_six_panel_scatter(panel: pd.DataFrame, path: Path, *, year: int) -> dict[str, Any]:
    """Draw zoomed J/N/A by RAW/REF scatters with full-range insets.

    The six main panels share a data-derived pooled 99th-percentile ceiling.
    Each inset shows the complete range for its universe and marks that cell's
    ten largest absolute gaps.  This preserves cross-panel comparability in the
    dense region without hiding the observations that determine the full range.
    """

    official = pd.to_numeric(panel["official_ais"], errors="raise").to_numpy(float)
    score_columns = [value[2] for value in CELL_COLUMNS.values()]
    scores = panel[score_columns].apply(pd.to_numeric, errors="raise")
    all_values = np.concatenate([official, scores.to_numpy(float).ravel()])
    if not np.isfinite(all_values).all() or (all_values < 0).any():
        raise AssertionError("Scatter panel contains invalid AIS")
    cell_zoom_quantiles = [
        float(np.quantile(np.concatenate([official, scores[column].to_numpy(float)]), 0.99))
        for column in score_columns
    ]
    zoom_axis_max = _nice_axis_ceiling(max(cell_zoom_quantiles))
    inset_axis_max = {
        universe: _nice_axis_ceiling(
            max(
                float(official.max()),
                *[
                    float(scores[score_column].max())
                    for cell, (cell_universe, _, score_column) in CELL_COLUMNS.items()
                    if cell_universe == universe
                ],
            )
        )
        for universe in UNIVERSE_ORDER
    }
    fig, axes = plt.subplots(3, 2, figsize=(10.2, 13.2), sharex=True, sharey=True)
    metadata: dict[str, Any] = {}
    for row_index, universe in enumerate(UNIVERSE_ORDER):
        for column_index, treatment in enumerate(TREATMENT_ORDER):
            cell_id = next(
                cell for cell, parts in CELL_COLUMNS.items()
                if parts[0] == universe and parts[1] == treatment
            )
            opened = scores[CELL_COLUMNS[cell_id][2]].to_numpy(float)
            axis = axes[row_index, column_index]
            gap = np.abs(opened - official)
            top_gap = np.argsort(-gap, kind="stable")[:10]
            official_zero = official == 0
            open_zero = opened == 0

            axis.plot(
                [0, zoom_axis_max], [0, zoom_axis_max], color="#222222",
                linewidth=0.9, zorder=1,
            )
            axis.scatter(
                official, opened, s=5.5, alpha=0.16, color="#2166ac",
                edgecolors="none", rasterized=True, zorder=2,
            )
            # Axis rugs show numerical zeros without imposing a log pseudocount.
            axis.scatter(
                np.zeros(int(official_zero.sum())), opened[official_zero], marker="_",
                s=26, color="#f28e2b", linewidths=0.85, rasterized=True,
                clip_on=False, zorder=4,
            )
            axis.scatter(
                official[open_zero], np.zeros(int(open_zero.sum())), marker="|",
                s=26, color="#f28e2b", linewidths=0.85, rasterized=True,
                clip_on=False, zorder=4,
            )
            axis.set_xlim(0, zoom_axis_max)
            axis.set_ylim(0, zoom_axis_max)
            axis.set_aspect("equal", adjustable="box")
            axis.grid(True, color="#e3e3e3", linewidth=0.4)
            axis.set_title(rf"$\mathit{{{universe}}}$-{treatment}")
            axis.text(
                0.03, 0.96, f"Zoom: 0-{zoom_axis_max:g}", transform=axis.transAxes,
                ha="left", va="top", fontsize=7.2,
            )
            axis.text(
                0.97, 0.04, f"n = {len(panel):,}", transform=axis.transAxes,
                ha="right", va="bottom", fontsize=7.5,
            )

            full_max = inset_axis_max[universe]
            high_value = np.argsort(
                -np.maximum(opened, official), kind="stable"
            )[:10]
            inset = axis.inset_axes([0.57, 0.55, 0.39, 0.39])
            inset.plot([0, full_max], [0, full_max], color="#222222", linewidth=0.6)
            inset.scatter(
                official, opened, s=1.9, alpha=0.13, color="#2166ac",
                edgecolors="none", rasterized=True,
            )
            inset.scatter(
                official[high_value], opened[high_value], s=20, marker="+",
                color="#2166ac", linewidths=0.9, rasterized=True, zorder=3,
            )
            inset.scatter(
                official[top_gap], opened[top_gap], s=18, marker="x", color="#d73027",
                linewidths=0.85, rasterized=True, zorder=4,
            )
            inset.set_xlim(0, full_max)
            inset.set_ylim(0, full_max)
            inset.set_aspect("equal", adjustable="box")
            inset.set_xticks([0, full_max], ["0", f"{full_max:g}"])
            inset.set_yticks([0, full_max], ["0", f"{full_max:g}"])
            inset.tick_params(labelsize=5.7, length=2, pad=1)
            inset.set_title(f"Full range: 0-{full_max:g}", fontsize=6.2, pad=1.5)
            metadata[cell_id] = {
                "official_zero_n": int(official_zero.sum()),
                "open_zero_n": int(open_zero.sum()),
                "top_absolute_gap_n": int(len(top_gap)),
                "top_plotted_value_n": int(len(high_value)),
                "marker_overlap_n": int(len(set(top_gap) & set(high_value))),
                "zoom_axis_max": zoom_axis_max,
                "full_range_axis_max": full_max,
            }
    for axis in axes[-1, :]:
        axis.set_xlabel("AIS")
    for axis in axes[:, 0]:
        axis.set_ylabel("ANS")
    handles = [
        Line2D([0], [0], color="#222222", linewidth=0.9, label="45-degree line"),
        Line2D([0], [0], marker="x", color="#d73027", linestyle="none",
               label="Ten largest absolute gaps (insets)"),
        Line2D([0], [0], marker="+", color="#2166ac", linestyle="none",
               label="Ten highest plotted AIS values (insets)"),
        Line2D([0], [0], marker="|", color="#f28e2b", linestyle="none",
               label="Zero on either axis (main-panel rug)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.suptitle(f"{year} AIS and ANS on the six-cell comparison sample")
    fig.tight_layout(rect=(0, 0.065, 1, 0.97))
    _save_figure(fig, path)
    return {
        "year": year,
        "n": int(len(panel)),
        "zoom_quantile": 0.99,
        "zoom_axis_min": 0.0,
        "zoom_axis_max": zoom_axis_max,
        "full_range_axis_max_by_universe": inset_axis_max,
        "panels": metadata,
        "full_range_insets_present": True,
        "binned_median_present": False,
    }


def _stable_field_mask(frame: pd.DataFrame) -> pd.Series:
    values = frame["stable_n_ge_300"]
    if pd.api.types.is_bool_dtype(values):
        return values.astype(bool)
    parsed = values.astype("string").str.strip().str.casefold().map(
        {"true": True, "false": False}
    )
    if parsed.isna().any():
        raise AssertionError("stable_n_ge_300 contains a non-Boolean value")
    return parsed.astype(bool)


def _nice_color_ceiling(maximum: float) -> float:
    if not math.isfinite(maximum) or maximum <= 0:
        raise AssertionError("Heatmap maximum must be positive and finite")
    return float(math.ceil(maximum * 10.0 - 1e-12) / 10.0)


def shared_heatmap_ceiling(
    deciles: pd.DataFrame, fields: pd.DataFrame, *, metric: str, year: int,
) -> float:
    """Return one non-clipping scale ceiling for the decile and field panels."""

    if metric not in {"scaled_mae", "nrmse"}:
        raise ValueError(f"Unsupported heatmap metric: {metric}")
    decile_values = pd.to_numeric(
        deciles.loc[
            deciles["year"].eq(year) & deciles["cell_id"].isin(DISPLAY_CELL_ORDER),
            metric,
        ],
        errors="raise",
    )
    field_rows = fields.loc[
        fields["year"].eq(year)
        & fields["cell_id"].isin(DISPLAY_CELL_ORDER)
        & _stable_field_mask(fields)
    ]
    field_values = pd.to_numeric(field_rows[metric], errors="raise")
    combined = np.concatenate(
        [decile_values.to_numpy(dtype=float), field_values.to_numpy(dtype=float)]
    )
    if len(combined) == 0 or not np.isfinite(combined).all() or (combined < 0).any():
        raise AssertionError(f"{metric} heatmap values are invalid")
    return _nice_color_ceiling(float(combined.max()))


def _plot_annotated_heatmap(
    matrix: pd.DataFrame, path: Path, *, title: str, y_label: str,
    colorbar_label: str, color_scale_max: float, figure_size: tuple[float, float],
    wrap_y_labels: bool, annotation_fontsize: float,
) -> None:
    values = matrix.to_numpy(dtype=float)
    if values.size == 0 or not np.isfinite(values).all() or (values < 0).any():
        raise AssertionError("Heatmap matrix is empty or contains invalid values")
    if float(values.max()) > color_scale_max + 1e-12:
        raise AssertionError("Heatmap color scale clips a reported value")
    fig, axis = plt.subplots(figsize=figure_size, constrained_layout=True)
    image = axis.imshow(
        values, aspect="auto", cmap="YlOrRd", vmin=0.0, vmax=color_scale_max
    )
    axis.set_xticks(
        range(len(matrix.columns)),
        [DISPLAY_CELL_LABELS[str(value)] for value in matrix.columns],
    )
    y_labels = [
        textwrap.fill(str(value), width=34) if wrap_y_labels else str(value)
        for value in matrix.index
    ]
    axis.set_yticks(range(len(matrix.index)), y_labels)
    axis.set_ylabel(y_label)
    axis.set_title(title)
    for y in range(values.shape[0]):
        for x in range(values.shape[1]):
            value = float(values[y, x])
            axis.text(
                x, y, f"{value:.2f}", ha="center", va="center",
                fontsize=annotation_fontsize,
                color="white" if value / color_scale_max > 0.58 else "black",
            )
    colorbar = fig.colorbar(image, ax=axis, shrink=0.88, pad=0.018)
    colorbar.set_label(colorbar_label)
    _save_figure(fig, path)


def plot_decile_heatmap(
    frame: pd.DataFrame, path: Path, *, metric: str, year: int,
    color_scale_max: float,
) -> dict[str, Any]:
    """Draw a six-cell error heatmap with D1 at the top and D10 at the bottom."""

    rows = [f"D{i}" for i in range(1, 11)]
    data = frame.loc[
        frame["year"].eq(year) & frame["cell_id"].isin(DISPLAY_CELL_ORDER)
    ].copy()
    matrix = data.pivot(
        index="official_ais_decile", columns="cell_id", values=metric
    ).reindex(index=rows, columns=list(DISPLAY_CELL_ORDER))
    if matrix.isna().any().any():
        raise AssertionError(f"{year} {metric} decile heatmap is incomplete")
    label = "Scaled MAE" if metric == "scaled_mae" else "Normalized RMSE"
    _plot_annotated_heatmap(
        matrix,
        path,
        title=f"{year} {label} by official AIS decile",
        y_label="Official AIS decile (D1 highest)",
        colorbar_label=label,
        color_scale_max=color_scale_max,
        figure_size=(9.6, 6.0),
        wrap_y_labels=False,
        annotation_fontsize=8.5,
    )
    return {
        "year": year,
        "metric": metric,
        "rows": 10,
        "columns": 6,
        "d1_highest": True,
        "color_scale_min": 0.0,
        "color_scale_max": color_scale_max,
    }


def plot_field_heatmap(
    frame: pd.DataFrame, path: Path, *, metric: str, year: int,
    color_scale_max: float, expected_fields: int | None = None,
) -> dict[str, Any]:
    """Draw the six-cell error heatmap for OA fields with at least 300 journals."""

    stable = frame.loc[
        frame["year"].eq(year)
        & frame["cell_id"].isin(DISPLAY_CELL_ORDER)
        & _stable_field_mask(frame)
    ].copy()
    field_order = sorted(stable["oa_field"].astype(str).unique())
    if expected_fields is not None and len(field_order) != expected_fields:
        raise AssertionError(
            f"{year} stable-field inventory changed: {len(field_order)} != {expected_fields}"
        )
    matrix = stable.pivot(
        index="oa_field", columns="cell_id", values=metric
    ).reindex(index=field_order, columns=list(DISPLAY_CELL_ORDER))
    if matrix.isna().any().any():
        raise AssertionError(f"{year} {metric} field heatmap is incomplete")
    label = "Scaled MAE" if metric == "scaled_mae" else "Normalized RMSE"
    _plot_annotated_heatmap(
        matrix,
        path,
        title=f"{year} {label} by OpenAlex field (n >= 300)",
        y_label="",
        colorbar_label=label,
        color_scale_max=color_scale_max,
        figure_size=(11.4, 9.2),
        wrap_y_labels=True,
        annotation_fontsize=10.0,
    )
    return {
        "year": year,
        "metric": metric,
        "rows": len(field_order),
        "columns": 6,
        "minimum_field_n": 300,
        "fields": field_order,
        "color_scale_min": 0.0,
        "color_scale_max": color_scale_max,
    }
