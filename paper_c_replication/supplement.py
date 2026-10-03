"""Supplement plotting and profiles; no input discovery or data constants."""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from .statistics import official_deciles
from .plots import separate_labels
SHORT = {"J": "j", "N": "n", "A": "oa", "L": "l"}


def specs(universes: tuple[str, ...], supports: dict[str, str]) -> list[dict[str, str]]:
    out = []
    for u in universes:
        for treatment, token in (("RAW", "raw"), ("REF", "filtered")):
            out.append({"id": f"{token.title()}_{u}", "universe": u, "treatment": treatment,
                        "score": f"ais_{token}_{SHORT[u]}", "support": supports[u]})
    return out


def subset(frame: pd.DataFrame, spec: dict[str, str]) -> pd.DataFrame:
    column = "in_support_s6" if spec["support"] == "s6" else f"in_support_{spec['support']}"
    return frame.loc[frame[column].astype(bool), ["official_ais", "oa_field", spec["score"]]].dropna().copy()


def summary(frame: pd.DataFrame, spec_list: list[dict[str, str]]) -> pd.DataFrame:
    rows = []
    for spec in spec_list:
        data = subset(frame, spec); error = data[spec["score"]] - data["official_ais"]
        rows.append({**spec, "n": len(data), "mae": error.abs().mean(), "rmse": np.sqrt((error ** 2).mean()),
                     "spearman": data["official_ais"].corr(data[spec["score"]], method="spearman"), "signed": error.mean()})
    return pd.DataFrame(rows)


def profile(frame: pd.DataFrame, spec_list: list[dict[str, str]], field: bool) -> pd.DataFrame:
    rows = []
    for spec in spec_list:
        data = subset(frame, spec)
        key = "field" if field else "decile"
        data[key] = data["oa_field"].astype("string").fillna("Unknown") if field else official_deciles(data["official_ais"])
        for group, block in data.groupby(key, sort=field):
            if field and len(block) < 300: continue
            error = block[spec["score"]] - block["official_ais"]
            rows.append({"id": spec["id"], "group": str(group), "value": error.abs().mean() / block["official_ais"].mean()})
    return pd.DataFrame(rows)


def headline(data: pd.DataFrame, path: Path, year: int) -> None:
    metrics = (("mae", "Mean absolute error", "{:.3f}"), ("rmse", "Root mean squared error", "{:.3f}"),
               ("spearman", r"Spearman $\rho$", "{:.3f}"), ("signed", "Signed mean error", "{:+.3f}"))
    colours = {"J": "#1769aa", "N": "#2a9d8f", "A": "#7a5195", "L": "#c65d7b"}
    fig, axes = plt.subplots(2, 2, figsize=(9.8, 7.4))
    for axis, (metric, label, fmt) in zip(axes.flat, metrics):
        for u, group in data.groupby("universe", sort=False):
            group = group.assign(x=group.treatment.map({"RAW": 0, "REF": 1})).sort_values("x")
            axis.plot(group.x, group[metric], marker="o", color=colours[u])
            for row in group.itertuples(): axis.annotate(fmt.format(getattr(row, metric)), (row.x, getattr(row, metric)), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=7)
        axis.set(xticks=(0, 1), xticklabels=("RAW", "REF"), ylabel=label); axis.grid(axis="y", color="#dddddd", linewidth=.6)
        if metric == "signed": axis.axhline(0, color="#555555", linewidth=.8)
        if metric in {"mae", "rmse"}: axis.set_ylim(bottom=0)
    handles = [Line2D([0], [0], marker="o", color=colours[u], label=u) for u in data.universe.drop_duplicates()]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, title="Universe")
    fig.suptitle(f"{year} aggregate agreement with AIS"); fig.tight_layout(rect=(0, .07, 1, .96))
    separate_labels(fig)
    fig.savefig(path, dpi=240, bbox_inches="tight"); plt.close(fig)


def scatter(frame: pd.DataFrame, spec_list: list[dict[str, str]], path: Path, year: int) -> None:
    nrows = len(spec_list) // 2; fig, axes = plt.subplots(nrows, 2, figsize=(10.2, 4.35 * nrows), squeeze=False)
    samples = [subset(frame, spec) for spec in spec_list]
    values = np.concatenate([np.concatenate((x.official_ais.to_numpy(), x[s["score"]].to_numpy())) for x, s in zip(samples, spec_list)])
    limit = float(np.ceil(np.quantile(values, .99) * 10) / 10)
    full_limits = {}
    for universe in dict.fromkeys(spec["universe"] for spec in spec_list):
        universe_values = np.concatenate([
            np.concatenate((data.official_ais.to_numpy(), data[spec["score"]].to_numpy()))
            for data, spec in zip(samples, spec_list) if spec["universe"] == universe
        ])
        maximum = float(universe_values.max())
        step = 10.0 ** np.floor(np.log10(maximum / 8.0))
        full_limits[universe] = float(np.ceil(maximum / step * 1.02) * step)
    for axis, data, spec in zip(axes.flat, samples, spec_list):
        official, opened = data.official_ais.to_numpy(), data[spec["score"]].to_numpy()
        gap = np.abs(opened - official)
        high_value = np.argsort(-np.maximum(official, opened), kind="stable")[:10]
        top_gap = np.argsort(-gap, kind="stable")[:10]
        axis.plot((0, limit), (0, limit), color="#222222", linewidth=.8)
        axis.scatter(official, opened, s=4.5, alpha=.15, color="#2166ac", edgecolors="none", rasterized=True)
        axis.set(xlim=(0, limit), ylim=(0, limit), aspect="equal", title=f"{spec['universe']}-{spec['treatment']}  (n = {len(data):,})")
        axis.grid(color="#e3e3e3", linewidth=.4)
        full_limit = full_limits[spec["universe"]]
        inset = axis.inset_axes([.57, .55, .39, .39])
        inset.plot((0, full_limit), (0, full_limit), color="#222222", linewidth=.6)
        inset.scatter(official, opened, s=1.8, alpha=.13, color="#2166ac", edgecolors="none", rasterized=True)
        inset.scatter(official[high_value], opened[high_value], s=20, marker="+", color="#2166ac", linewidths=.9, rasterized=True, zorder=3)
        inset.scatter(official[top_gap], opened[top_gap], s=18, marker="x", color="#d73027", linewidths=.85, rasterized=True, zorder=4)
        inset.set(xlim=(0, full_limit), ylim=(0, full_limit), aspect="equal", xticks=(0, full_limit), yticks=(0, full_limit))
        inset.tick_params(labelsize=5.7, length=2, pad=1)
        inset.set_title(f"Full range: 0-{full_limit:g}", fontsize=6.2, pad=1.5)
    for axis in axes[-1, :]: axis.set_xlabel("AIS")
    for axis in axes[:, 0]: axis.set_ylabel("ANS")
    handles = [Line2D([0], [0], color="#222222", linewidth=.8, label="45-degree line"),
               Line2D([0], [0], marker="x", color="#d73027", linestyle="none", label="Ten largest absolute gaps (insets)"),
               Line2D([0], [0], marker="+", color="#2166ac", linestyle="none", label="Ten highest AIS values (insets)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.suptitle(f"{year} AIS and ANS on stated matched samples"); fig.tight_layout(rect=(0, .055, 1, .97))
    fig.savefig(path, dpi=240, bbox_inches="tight"); plt.close(fig)


def heatmap(data: pd.DataFrame, spec_list: list[dict[str, str]], path: Path, year: int, kind: str) -> None:
    groups = [f"D{i}" for i in range(1, 11)] if kind == "decile" else sorted(set.intersection(*[set(data.loc[data.id.eq(s["id"]), "group"]) for s in spec_list]))
    matrix = data.pivot(index="group", columns="id", values="value").reindex(index=groups, columns=[s["id"] for s in spec_list])
    ceiling = float(np.ceil(matrix.to_numpy(float).max() * 10) / 10)
    fig, axis = plt.subplots(figsize=(9.6, 6.0) if kind == "decile" else (11.4, 9.2), constrained_layout=True)
    image = axis.imshow(matrix.to_numpy(float), aspect="auto", cmap="YlOrRd", vmin=0, vmax=ceiling)
    axis.set_xticks(range(len(matrix.columns)), [x.replace("Filtered", "REF").replace("Raw", "RAW") for x in matrix.columns]); axis.set_yticks(range(len(groups)), groups)
    axis.set_title(f"{year} scaled MAE by {'official AIS decile' if kind == 'decile' else 'OpenAlex field (n >= 300)'}")
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]):
            value = float(matrix.iloc[y, x]); axis.text(x, y, f"{value:.2f}", ha="center", va="center", fontsize=8, color="white" if value / ceiling > .58 else "black")
    fig.colorbar(image, ax=axis, shrink=.88, pad=.018, label="Scaled MAE"); fig.savefig(path, dpi=240, bbox_inches="tight"); plt.close(fig)
