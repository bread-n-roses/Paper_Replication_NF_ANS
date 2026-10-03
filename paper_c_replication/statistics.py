"""Registered aggregate statistics for Paper C.

This module is deliberately free of repository I/O.  It receives already-frozen,
row-aligned private panels and returns aggregate records only.  In particular, it
does not choose a support, join official data, or write row-level values.  Those
responsibilities belong to the support freezer and the analysis driver.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
import hashlib
from typing import Any

import numpy as np
import pandas as pd


BASE_SEED = 20260821
BOOTSTRAP_DRAWS = 4_999
PRIMARY_CELLS = (
    "Raw_J",
    "Filtered_J",
    "Raw_N",
    "Filtered_N",
    "Raw_OA",
    "Filtered_OA",
)
CUBE_CELLS = tuple(f"E{e}T{t}A{a}" for e in (0, 1) for t in (0, 1) for a in (0, 1))


def child_seed(label: str, *, base_seed: int = BASE_SEED) -> int:
    """Return the registered 32-bit SHA256-derived child seed."""

    if not isinstance(label, str) or not label.strip():
        raise ValueError("Bootstrap label must be a nonblank string")
    digest = hashlib.sha256(f"{base_seed}|{label}".encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def _finite_numeric(values: Sequence[float] | pd.Series, *, label: str) -> np.ndarray:
    out = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=np.float64)
    if out.size == 0:
        raise ValueError(f"{label} is empty")
    if not np.isfinite(out).all():
        raise ValueError(f"{label} contains missing or nonfinite values")
    return out


def _nonnegative(values: Sequence[float] | pd.Series, *, label: str) -> np.ndarray:
    out = _finite_numeric(values, label=label)
    if (out < 0).any():
        raise ValueError(f"{label} contains negative values")
    return out


def spearman(x: Sequence[float] | pd.Series, y: Sequence[float] | pd.Series) -> float:
    """Spearman correlation with average ranks and explicit constant handling."""

    left = _finite_numeric(x, label="Spearman x")
    right = _finite_numeric(y, label="Spearman y")
    if len(left) != len(right):
        raise ValueError("Spearman vectors have different lengths")
    if len(left) < 2:
        return float("nan")
    left_rank = pd.Series(left).rank(method="average").to_numpy(dtype=np.float64)
    right_rank = pd.Series(right).rank(method="average").to_numpy(dtype=np.float64)
    if np.ptp(left_rank) == 0 or np.ptp(right_rank) == 0:
        return float("nan")
    return float(np.corrcoef(left_rank, right_rank)[0, 1])


def max_rank_percentile(values: Sequence[float] | pd.Series) -> np.ndarray:
    """Return pandas max-rank percentiles, preserving ties at their upper rank."""

    numeric = _finite_numeric(values, label="Rank values")
    return pd.Series(numeric).rank(method="max", pct=True).to_numpy(dtype=np.float64)


def official_deciles(values: Sequence[float] | pd.Series) -> np.ndarray:
    """Label D1 (highest) through D10 using registered max-rank percentiles.

    Exact boundaries remain in the lower-ranked decile.  Thus the registered
    top decile is exactly ``percentile > .90`` rather than ``>= .90``.
    """

    numeric = _finite_numeric(values, label="Rank values")
    # Integer ceiling avoids floating-point ambiguity at exact boundaries:
    # D = 11 - ceil(10 * max_rank / N).  A max rank of .90N is D2, while
    # any rank above .90N is D1, exactly matching Top10(X) = rank_pct > .90.
    ranks = pd.Series(numeric).rank(method="max").to_numpy(dtype=np.int64)
    quantile_number = (10 * ranks + len(ranks) - 1) // len(ranks)
    decile = 11 - quantile_number
    return np.asarray([f"D{int(value)}" for value in decile], dtype=object)


def _zero_cross_tab(official: np.ndarray, opened: np.ndarray) -> dict[str, int]:
    official_zero = official == 0
    open_zero = opened == 0
    return {
        "official_zero_open_zero_n": int((official_zero & open_zero).sum()),
        "official_zero_open_positive_n": int((official_zero & ~open_zero).sum()),
        "official_positive_open_zero_n": int((~official_zero & open_zero).sum()),
        "official_positive_open_positive_n": int((~official_zero & ~open_zero).sum()),
    }


def cell_statistics(
    official: Sequence[float] | pd.Series,
    opened: Sequence[float] | pd.Series,
) -> dict[str, Any]:
    """Compute the registered aggregate statistics for one AIS cell."""

    c = _nonnegative(official, label="Official AIS")
    o = _nonnegative(opened, label="Open AIS")
    if len(c) != len(o):
        raise ValueError("Official and open AIS vectors have different lengths")
    error = o - c
    absolute = np.abs(error)
    c_pct = max_rank_percentile(c)
    o_pct = max_rank_percentile(o)
    official_top = c_pct > 0.90
    open_top = o_pct > 0.90
    official_top_n = int(official_top.sum())
    if official_top_n == 0:
        top_retention = float("nan")
    else:
        top_retention = float((official_top & open_top).sum() / official_top_n)
    result: dict[str, Any] = {
        "n": int(len(c)),
        "mae": float(absolute.mean()),
        "mse": float(np.square(error).mean()),
        "spearman": spearman(c, o),
        "signed_mean_error": float(error.mean()),
        "median_absolute_error": float(np.median(absolute)),
        "official_top_decile_n": official_top_n,
        "open_top_decile_n": int(open_top.sum()),
        "official_top_decile_retained_n": int((official_top & open_top).sum()),
        "official_top_decile_retention": top_retention,
    }
    result.update(_zero_cross_tab(c, o))
    return result


def decile_error_profile(
    official: Sequence[float] | pd.Series,
    opened: Sequence[float] | pd.Series,
) -> pd.DataFrame:
    """Return D1-highest error summaries on one exact support."""

    c = _nonnegative(official, label="Official AIS")
    o = _nonnegative(opened, label="Open AIS")
    if len(c) != len(o):
        raise ValueError("Official and open AIS vectors have different lengths")
    work = pd.DataFrame({"official": c, "open": o, "decile": official_deciles(c)})
    work["error"] = work["open"] - work["official"]
    work["absolute_error"] = work["error"].abs()
    work["squared_error"] = work["error"].pow(2)
    rows: list[dict[str, Any]] = []
    for number in range(1, 11):
        label = f"D{number}"
        group = work.loc[work["decile"].eq(label)]
        if group.empty:
            continue
        rows.append(
            {
                "official_ais_decile": label,
                "n": int(len(group)),
                "mae": float(group["absolute_error"].mean()),
                "mse": float(group["squared_error"].mean()),
                "signed_mean_error": float(group["error"].mean()),
            }
        )
    if sum(row["n"] for row in rows) != len(work):
        raise AssertionError("Official-AIS decile accounting lost rows")
    return pd.DataFrame(rows)


def positive_log_calibration(
    official: Sequence[float] | pd.Series,
    opened: Sequence[float] | pd.Series,
) -> dict[str, Any]:
    """Fit the registered positive-pair log calibration without a pseudocount."""

    c = _nonnegative(official, label="Official AIS")
    o = _nonnegative(opened, label="Open AIS")
    if len(c) != len(o):
        raise ValueError("Official and open AIS vectors have different lengths")
    keep = (c > 0) & (o > 0)
    n = int(keep.sum())
    if n < 2 or np.ptp(np.log(c[keep])) == 0:
        alpha = beta = float("nan")
    else:
        design = np.column_stack([np.ones(n), np.log(c[keep])])
        alpha, beta = np.linalg.lstsq(design, np.log(o[keep]), rcond=None)[0]
    return {
        "n_all": int(len(c)),
        "n_positive_pairs": n,
        "excluded_zero_pairs": int((~keep).sum()),
        "alpha": float(alpha),
        "beta": float(beta),
    }


def named_contrast_statistics(
    official: Sequence[float] | pd.Series,
    score_a: Sequence[float] | pd.Series,
    score_b: Sequence[float] | pd.Series,
) -> dict[str, float | int]:
    """Return benchmark-error changes and movement for ``b minus a``."""

    c = _nonnegative(official, label="Official AIS")
    a = _nonnegative(score_a, label="Construction a AIS")
    b = _nonnegative(score_b, label="Construction b AIS")
    if not (len(c) == len(a) == len(b)):
        raise ValueError("Named-contrast vectors have different lengths")
    error_a = a - c
    error_b = b - c
    return {
        "n": int(len(c)),
        "mean_delta_ae": float((np.abs(error_b) - np.abs(error_a)).mean()),
        "mean_delta_se": float((np.square(error_b) - np.square(error_a)).mean()),
        "mean_abs_movement": float(np.abs(b - a).mean()),
        "delta_spearman": float(spearman(c, b) - spearman(c, a)),
    }


def interaction_statistics(
    official: Sequence[float] | pd.Series,
    raw_j: Sequence[float] | pd.Series,
    filtered_j: Sequence[float] | pd.Series,
    raw_u: Sequence[float] | pd.Series,
    filtered_u: Sequence[float] | pd.Series,
) -> dict[str, float | int]:
    """Return the registered Raw/Filtered-by-universe interaction summaries."""

    c = _nonnegative(official, label="Official AIS")
    rj = _nonnegative(raw_j, label="Raw J AIS")
    fj = _nonnegative(filtered_j, label="Filtered J AIS")
    ru = _nonnegative(raw_u, label="Raw comparison-universe AIS")
    fu = _nonnegative(filtered_u, label="Filtered comparison-universe AIS")
    if len({len(c), len(rj), len(fj), len(ru), len(fu)}) != 1:
        raise ValueError("Interaction vectors have different lengths")
    interaction = (fu - ru) - (fj - rj)
    did_ae = (np.abs(fu - c) - np.abs(ru - c)) - (
        np.abs(fj - c) - np.abs(rj - c)
    )
    did_se = (np.square(fu - c) - np.square(ru - c)) - (
        np.square(fj - c) - np.square(rj - c)
    )
    rank_interaction = (
        spearman(c, fu)
        - spearman(c, ru)
        - spearman(c, fj)
        + spearman(c, rj)
    )
    return {
        "n": int(len(c)),
        "mean_interaction": float(interaction.mean()),
        "mean_abs_interaction": float(np.abs(interaction).mean()),
        "did_mean_ae": float(did_ae.mean()),
        "did_mean_se": float(did_se.mean()),
        "rank_interaction": float(rank_interaction),
    }


def stratified_paired_bootstrap(
    frame: pd.DataFrame,
    statistic: Callable[[pd.DataFrame], Mapping[str, float]],
    *,
    label: str,
    strata_column: str = "oa_field",
    unit_column: str = "issn_l",
    draws: int = BOOTSTRAP_DRAWS,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Return fixed-seed percentile intervals from paired field-stratified draws.

    Whole rows are sampled, so every score vector remains paired.  Each stratum
    contributes exactly its observed number of rows to every draw.  Inputs are
    sorted by stratum and unit before sampling, making results independent of
    incoming row order.
    """

    if draws < 1:
        raise ValueError("Bootstrap draws must be positive")
    required = {strata_column, unit_column}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Bootstrap frame misses columns {sorted(missing)}")
    if frame.empty:
        raise ValueError("Bootstrap frame is empty")
    if frame[unit_column].isna().any() or frame[unit_column].duplicated().any():
        raise ValueError("Bootstrap units must be nonmissing and unique")
    work = frame.copy()
    work[strata_column] = work[strata_column].astype("string").fillna("Unknown")
    work[unit_column] = work[unit_column].astype("string")
    work = work.sort_values([strata_column, unit_column], kind="mergesort").reset_index(drop=True)
    group_indices = [
        group.index.to_numpy(dtype=np.int64)
        for _, group in work.groupby(strata_column, sort=True, dropna=False)
    ]
    estimate = {key: float(value) for key, value in statistic(work).items()}
    if len(estimate) != 1:
        raise ValueError(
            "Each registered scalar statistic needs its own bootstrap label and child seed"
        )
    if not np.isfinite(list(estimate.values())).all():
        raise ValueError("Bootstrap point estimate is nonfinite")
    seed = child_seed(label)
    generator = np.random.default_rng(seed)
    samples = {key: np.empty(draws, dtype=np.float64) for key in estimate}
    for draw in range(draws):
        selected = np.concatenate(
            [generator.choice(index, size=len(index), replace=True) for index in group_indices]
        )
        values = statistic(work.iloc[selected])
        if set(values) != set(estimate):
            raise AssertionError("Bootstrap statistic keys changed between draws")
        for key, value in values.items():
            samples[key][draw] = float(value)
    rows: list[dict[str, Any]] = []
    for key in estimate:
        values = samples[key]
        finite = np.isfinite(values)
        if not finite.all():
            raise AssertionError(f"Bootstrap statistic {key!r} produced nonfinite draws")
        lower, upper = np.quantile(values, [0.025, 0.975])
        rows.append(
            {
                "statistic": key,
                "estimate": estimate[key],
                "ci_lower_95": float(lower),
                "ci_upper_95": float(upper),
                "draws": int(draws),
                "bootstrap_label": label,
                "child_seed": seed,
            }
        )
    receipt = {
        "label": label,
        "child_seed": seed,
        "base_seed": BASE_SEED,
        "draws": int(draws),
        "n": int(len(work)),
        "strata_column": strata_column,
        "strata_count": int(len(group_indices)),
        "paired_rows": True,
    }
    return pd.DataFrame(rows), receipt


class PreparedStratifiedBootstrap:
    """Chunked exact paired bootstrap for registered scalar estimands.

    The frame and strata layout are sorted once and reused across scalar
    statistics on the same support.  Mean estimands use multinomial row counts
    and matrix products.  Rank estimands reconstruct the exact average ranks of
    the resampled observations from those counts, avoiding repeated sorting of
    the same value vectors.  Every scalar call nevertheless has its own label,
    seed, and independent resampling stream, as registered in the plan.
    """

    def __init__(
        self,
        frame: pd.DataFrame,
        *,
        strata_column: str = "oa_field",
        unit_column: str = "issn_l",
    ) -> None:
        required = {strata_column, unit_column}
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"Bootstrap frame misses columns {sorted(missing)}")
        if frame.empty:
            raise ValueError("Bootstrap frame is empty")
        if frame[unit_column].isna().any() or frame[unit_column].duplicated().any():
            raise ValueError("Bootstrap units must be nonmissing and unique")
        work = frame.copy()
        work[strata_column] = work[strata_column].astype("string").fillna("Unknown")
        work[unit_column] = work[unit_column].astype("string")
        self.frame = work.sort_values(
            [strata_column, unit_column], kind="mergesort"
        ).reset_index(drop=True)
        self.strata_column = strata_column
        self.unit_column = unit_column
        sizes = self.frame.groupby(strata_column, sort=True, dropna=False).size().to_numpy(int)
        endpoints = np.cumsum(sizes)
        starts = np.r_[0, endpoints[:-1]]
        self._strata = tuple(
            (
                int(start),
                int(stop),
                np.full(int(stop - start), 1.0 / int(stop - start), dtype=np.float64),
            )
            for start, stop in zip(starts, endpoints)
        )
        self._rank_plans: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    @property
    def n(self) -> int:
        return len(self.frame)

    @property
    def strata_count(self) -> int:
        return len(self._strata)

    def values(self, column: str) -> np.ndarray:
        if column not in self.frame:
            raise ValueError(f"Bootstrap frame misses column {column!r}")
        return _finite_numeric(self.frame[column], label=column)

    def _draw_counts(self, generator: np.random.Generator, batch: int) -> np.ndarray:
        counts = np.zeros((batch, self.n), dtype=np.int32)
        for start, stop, probabilities in self._strata:
            size = stop - start
            counts[:, start:stop] = generator.multinomial(
                size, probabilities, size=batch
            ).astype(np.int32, copy=False)
        if not np.equal(counts.sum(axis=1), self.n).all():
            raise AssertionError("Stratified bootstrap draw changed the sample size")
        return counts

    def _receipt(self, label: str, draws: int) -> dict[str, Any]:
        return {
            "label": label,
            "child_seed": child_seed(label),
            "base_seed": BASE_SEED,
            "draws": int(draws),
            "n": self.n,
            "strata_column": self.strata_column,
            "strata_count": self.strata_count,
            "paired_rows": True,
        }

    @staticmethod
    def _interval_record(
        statistic: str,
        estimate: float,
        draws_values: np.ndarray,
        label: str,
        draws: int,
    ) -> dict[str, Any]:
        if not np.isfinite(estimate) or not np.isfinite(draws_values).all():
            raise AssertionError(f"Bootstrap statistic {statistic!r} is nonfinite")
        lower, upper = np.quantile(draws_values, [0.025, 0.975])
        return {
            "statistic": statistic,
            "estimate": float(estimate),
            "ci_lower_95": float(lower),
            "ci_upper_95": float(upper),
            "draws": int(draws),
            "bootstrap_label": label,
            "child_seed": child_seed(label),
        }

    def mean_interval(
        self,
        values: Sequence[float] | pd.Series | np.ndarray,
        *,
        statistic: str,
        label: str,
        draws: int = BOOTSTRAP_DRAWS,
        batch_size: int = 64,
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        """Bootstrap one row-level mean estimand with its own seed."""

        vector = _finite_numeric(values, label=statistic)
        if len(vector) != self.n:
            raise ValueError("Mean-estimand vector does not align with bootstrap frame")
        if draws < 1 or batch_size < 1:
            raise ValueError("Bootstrap draws and batch size must be positive")
        generator = np.random.default_rng(child_seed(label))
        sampled = np.empty(draws, dtype=np.float64)
        cursor = 0
        while cursor < draws:
            batch = min(batch_size, draws - cursor)
            counts = self._draw_counts(generator, batch)
            sampled[cursor : cursor + batch] = counts @ vector / self.n
            cursor += batch
        return (
            self._interval_record(statistic, float(vector.mean()), sampled, label, draws),
            self._receipt(label, draws),
        )

    def _rank_plan(self, column: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if column in self._rank_plans:
            return self._rank_plans[column]
        values = self.values(column)
        order = np.argsort(values, kind="mergesort")
        ordered = values[order]
        starts = np.r_[0, np.flatnonzero(ordered[1:] != ordered[:-1]) + 1].astype(np.int64)
        sorted_group = np.cumsum(
            np.r_[True, ordered[1:] != ordered[:-1]]
        ).astype(np.int64) - 1
        group_by_original = np.empty(self.n, dtype=np.int64)
        group_by_original[order] = sorted_group
        plan = (order, starts, group_by_original)
        self._rank_plans[column] = plan
        return plan

    def _rank_matrix(self, counts: np.ndarray, column: str) -> np.ndarray:
        order, starts, group_by_original = self._rank_plan(column)
        group_counts = np.add.reduceat(counts[:, order], starts, axis=1)
        before = np.cumsum(group_counts, axis=1) - group_counts
        group_ranks = before + (group_counts + 1.0) / 2.0
        return group_ranks[:, group_by_original]

    @staticmethod
    def _weighted_correlations(
        counts: np.ndarray, rank_x: np.ndarray, rank_y: np.ndarray
    ) -> np.ndarray:
        total = counts.sum(axis=1).astype(np.float64)
        mean_x = (counts * rank_x).sum(axis=1) / total
        mean_y = (counts * rank_y).sum(axis=1) / total
        mean_x2 = (counts * np.square(rank_x)).sum(axis=1) / total
        mean_y2 = (counts * np.square(rank_y)).sum(axis=1) / total
        mean_xy = (counts * rank_x * rank_y).sum(axis=1) / total
        variance_x = np.maximum(mean_x2 - np.square(mean_x), 0.0)
        variance_y = np.maximum(mean_y2 - np.square(mean_y), 0.0)
        denominator = np.sqrt(variance_x * variance_y)
        result = np.full(len(total), np.nan, dtype=np.float64)
        valid = denominator > 0
        result[valid] = (mean_xy[valid] - mean_x[valid] * mean_y[valid]) / denominator[valid]
        result[valid] = np.clip(result[valid], -1.0, 1.0)
        return result

    def rank_expression_interval(
        self,
        terms: Sequence[tuple[float, str, str]],
        *,
        statistic: str,
        label: str,
        draws: int = BOOTSTRAP_DRAWS,
        batch_size: int = 16,
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        """Bootstrap a linear expression of exact resampled Spearman correlations."""

        if not terms:
            raise ValueError("Rank expression has no terms")
        columns = sorted({column for _, x, y in terms for column in (x, y)})
        for column in columns:
            self.values(column)
        estimate = float(
            sum(coefficient * spearman(self.frame[x], self.frame[y]) for coefficient, x, y in terms)
        )
        generator = np.random.default_rng(child_seed(label))
        sampled = np.empty(draws, dtype=np.float64)
        cursor = 0
        while cursor < draws:
            batch = min(batch_size, draws - cursor)
            counts = self._draw_counts(generator, batch)
            ranks = {column: self._rank_matrix(counts, column) for column in columns}
            value = np.zeros(batch, dtype=np.float64)
            for coefficient, x, y in terms:
                value += coefficient * self._weighted_correlations(
                    counts, ranks[x], ranks[y]
                )
            sampled[cursor : cursor + batch] = value
            cursor += batch
        return (
            self._interval_record(statistic, estimate, sampled, label, draws),
            self._receipt(label, draws),
        )


def field_cell_statistics(
    frame: pd.DataFrame,
    *,
    official_column: str,
    cell_columns: Mapping[str, str],
    field_column: str = "oa_field",
    stable_threshold: int = 300,
) -> pd.DataFrame:
    """Compute per-field six-cell summaries without creating a macro score."""

    required = {official_column, field_column, *cell_columns.values()}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Field panel misses columns {sorted(missing)}")
    if stable_threshold < 1:
        raise ValueError("Field stability threshold must be positive")
    work = frame.copy()
    work[field_column] = work[field_column].astype("string").fillna("Unknown")
    rows: list[dict[str, Any]] = []
    for field, group in work.groupby(field_column, sort=True, dropna=False):
        for cell_id, column in cell_columns.items():
            stats = cell_statistics(group[official_column], group[column])
            rows.append(
                {
                    "oa_field": str(field),
                    "cell_id": cell_id,
                    "stable_n_ge_300": bool(len(group) >= stable_threshold),
                    **stats,
                }
            )
    return pd.DataFrame(rows)


def ef_bridge_statistics(
    official_ef: Sequence[float] | pd.Series,
    open_ef: Sequence[float] | pd.Series,
    *,
    official_observed_reported_total: float,
) -> dict[str, Any]:
    """Compute the registered compact EF bridge on a shared support."""

    c = _nonnegative(official_ef, label="Official EF")
    o = _nonnegative(open_ef, label="Open EF")
    if len(c) != len(o):
        raise ValueError("Official and open EF vectors have different lengths")
    c_sum = float(c.sum())
    o_sum = float(o.sum())
    if c_sum <= 0 or o_sum <= 0:
        raise ValueError("EF bridge requires positive mass on shared support")
    if not np.isfinite(official_observed_reported_total) or official_observed_reported_total <= 0:
        raise ValueError("Official observed reported EF total must be finite and positive")
    p_c = c / c_sum
    p_o = o / o_sum
    return {
        "n": int(len(c)),
        "spearman": spearman(c, o),
        "total_variation": float(0.5 * np.abs(p_o - p_c).sum()),
        "open_native_ef_mass_on_support": o_sum,
        "official_reported_ef_mass_on_support": c_sum,
        "official_share_of_observed_reported_ef_mass": float(
            c_sum / official_observed_reported_total
        ),
    }


def component_cube_effects(cube: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute registered journal-level and aggregate 2024 component-cube effects."""

    required = {"issn_l", "cell_id", "ais"}
    missing = required - set(cube.columns)
    if missing:
        raise ValueError(f"Component cube misses columns {sorted(missing)}")
    if cube.duplicated(["issn_l", "cell_id"]).any():
        raise ValueError("Component cube has duplicate journal-cell rows")
    observed_cells = set(cube["cell_id"].astype(str).unique())
    if observed_cells != set(CUBE_CELLS):
        raise ValueError("Component cube does not contain exactly the eight registered cells")
    wide = cube.pivot(index="issn_l", columns="cell_id", values="ais")
    wide = wide.reindex(columns=CUBE_CELLS)
    if wide.isna().any().any():
        raise ValueError("Component-cube support contains undefined AIS")
    y = {cell: wide[cell].to_numpy(dtype=np.float64) for cell in CUBE_CELLS}
    effects = pd.DataFrame(index=wide.index)
    effects["tau_E"] = sum(
        y[f"E1T{t}A{a}"] - y[f"E0T{t}A{a}"] for t in (0, 1) for a in (0, 1)
    ) / 4.0
    effects["tau_T"] = sum(
        y[f"E{e}T1A{a}"] - y[f"E{e}T0A{a}"] for e in (0, 1) for a in (0, 1)
    ) / 4.0
    effects["tau_A"] = sum(
        y[f"E{e}T{t}A1"] - y[f"E{e}T{t}A0"] for e in (0, 1) for t in (0, 1)
    ) / 4.0
    effects["tau_ET"] = sum(
        (y[f"E1T1A{a}"] - y[f"E0T1A{a}"])
        - (y[f"E1T0A{a}"] - y[f"E0T0A{a}"])
        for a in (0, 1)
    ) / 2.0
    effects["tau_EA"] = sum(
        (y[f"E1T{t}A1"] - y[f"E0T{t}A1"])
        - (y[f"E1T{t}A0"] - y[f"E0T{t}A0"])
        for t in (0, 1)
    ) / 2.0
    effects["tau_TA"] = sum(
        (y[f"E{e}T1A1"] - y[f"E{e}T0A1"])
        - (y[f"E{e}T1A0"] - y[f"E{e}T0A0"])
        for e in (0, 1)
    ) / 2.0
    effects["tau_ETA"] = (
        (y["E1T1A1"] - y["E0T1A1"])
        - (y["E1T0A1"] - y["E0T0A1"])
        - (y["E1T1A0"] - y["E0T1A0"])
        + (y["E1T0A0"] - y["E0T0A0"])
    )
    effects = effects.reset_index()
    rows = []
    for column in ("tau_E", "tau_T", "tau_A", "tau_ET", "tau_EA", "tau_TA", "tau_ETA"):
        rows.append(
            {
                "effect": column,
                "n": int(len(effects)),
                "signed_mean": float(effects[column].mean()),
                "mean_absolute_movement": float(effects[column].abs().mean()),
            }
        )
    return effects, pd.DataFrame(rows)


__all__ = [
    "BASE_SEED",
    "BOOTSTRAP_DRAWS",
    "CUBE_CELLS",
    "PRIMARY_CELLS",
    "PreparedStratifiedBootstrap",
    "cell_statistics",
    "child_seed",
    "component_cube_effects",
    "decile_error_profile",
    "ef_bridge_statistics",
    "field_cell_statistics",
    "interaction_statistics",
    "max_rank_percentile",
    "named_contrast_statistics",
    "official_deciles",
    "positive_log_calibration",
    "spearman",
    "stratified_paired_bootstrap",
]
