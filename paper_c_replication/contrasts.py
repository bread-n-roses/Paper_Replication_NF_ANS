"""Portable contrast definitions and component summaries; no repository access."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np
import pandas as pd
from .statistics import *
@dataclass(frozen=True)
class ContrastSpec:
    contrast_id: str
    kind: str
    cells: tuple[str, ...]
    maximum_support_id: str
    appendix: bool = False


CONTRAST_SPECS = (
    ContrastSpec("Filtered_minus_Raw_J", "pair", ("Raw_J", "Filtered_J"), "max_treatment_J"),
    ContrastSpec("Filtered_minus_Raw_N", "pair", ("Raw_N", "Filtered_N"), "max_treatment_N"),
    ContrastSpec("Filtered_minus_Raw_OA", "pair", ("Raw_OA", "Filtered_OA"), "max_treatment_OA"),
    ContrastSpec("N_minus_J_Raw", "pair", ("Raw_J", "Raw_N"), "max_universe_Raw_J_N"),
    ContrastSpec(
        "N_minus_J_Filtered", "pair", ("Filtered_J", "Filtered_N"),
        "max_universe_Filtered_J_N",
    ),
    ContrastSpec("OA_minus_J_Raw", "pair", ("Raw_J", "Raw_OA"), "max_universe_Raw_J_OA"),
    ContrastSpec(
        "OA_minus_J_Filtered", "pair", ("Filtered_J", "Filtered_OA"),
        "max_universe_Filtered_J_OA",
    ),
    ContrastSpec(
        "RawFiltered_by_N_minus_J", "interaction",
        ("Raw_J", "Filtered_J", "Raw_N", "Filtered_N"), "max_rectangle_J_N",
    ),
    ContrastSpec(
        "RawFiltered_by_OA_minus_J", "interaction",
        ("Raw_J", "Filtered_J", "Raw_OA", "Filtered_OA"), "max_rectangle_J_OA",
    ),
    ContrastSpec(
        "Filtered_minus_Raw_L", "pair", ("Raw_L", "Filtered_L"),
        "max_treatment_L", appendix=True,
    ),
    ContrastSpec(
        "L_minus_J_Raw", "pair", ("Raw_J", "Raw_L"),
        "max_universe_Raw_J_L", appendix=True,
    ),
    ContrastSpec(
        "L_minus_J_Filtered", "pair", ("Filtered_J", "Filtered_L"),
        "max_universe_Filtered_J_L", appendix=True,
    ),
    ContrastSpec(
        "RawFiltered_by_L_minus_J", "interaction",
        ("Raw_J", "Filtered_J", "Raw_L", "Filtered_L"),
        "max_rectangle_J_L", appendix=True,
    ),
)


def _score_key(metric: str, cell: str) -> str:
    treatment, universe = cell.split("_", maxsplit=1)
    return f"{metric}_{treatment.lower()}_{universe.lower()}"


def _contrast_statistic(spec: ContrastSpec, sample: pd.DataFrame) -> dict[str, float]:
    official = sample["official_ais"]
    if spec.kind == "pair":
        a, b = (_score_key("ais", cell) for cell in spec.cells)
        result = named_contrast_statistics(official, sample[a], sample[b])
    elif spec.kind == "interaction":
        rj, fj, ru, fu = (_score_key("ais", cell) for cell in spec.cells)
        result = interaction_statistics(
            official, sample[rj], sample[fj], sample[ru], sample[fu]
        )
    else:  # pragma: no cover - frozen module constant
        raise ValueError(f"Unknown contrast kind {spec.kind!r}")
    return {key: float(value) for key, value in result.items() if key != "n"}


def _contrast_supports(spec: ContrastSpec) -> tuple[tuple[str, str], ...]:
    if spec.appendix:
        return (("maximum_support", spec.maximum_support_id),)
    return (("annual_s6", "s6"), ("maximum_support", spec.maximum_support_id))


def _support_mask(panel: pd.DataFrame, support_id: str) -> pd.Series:
    column = "in_support_s6" if support_id == "s6" else f"in_support_{support_id}"
    if column not in panel:
        raise AssertionError(f"Frozen support column is missing: {column}")
    mask = panel[column].astype(bool)
    if not mask.any():
        raise AssertionError(f"Frozen support is empty: {support_id}")
    return mask


def _cube_outputs(
    cube: pd.DataFrame, panel_2024: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    s6_ids = set(panel_2024.loc[panel_2024["in_support_s6"], "issn_l"].astype(str))
    counts = cube.groupby("issn_l", sort=False)["ais"].agg(
        rows="size", defined=lambda values: int(values.notna().sum())
    )
    complete_ids = set(counts.loc[(counts["rows"] == 8) & (counts["defined"] == 8)].index.astype(str))
    if not s6_ids <= complete_ids:
        raise AssertionError("2024 S6 is not mechanism-complete in the component cube")
    summaries: list[pd.DataFrame] = []
    for scope, identities in (
        ("annual_s6", s6_ids),
        ("mechanism_complete", complete_ids),
    ):
        _, aggregate = component_cube_effects(cube.loc[cube["issn_l"].isin(identities)])
        aggregate.insert(
            0, "support_id", "s6" if scope == "annual_s6" else "mechanism_complete"
        )
        aggregate.insert(0, "support_scope", scope)
        aggregate.insert(0, "year", 2024)
        summaries.append(aggregate)

    benchmark_rows = []
    official = panel_2024.set_index("issn_l")["official_ais"]
    for cell_id, group in cube.loc[cube["issn_l"].isin(s6_ids)].groupby("cell_id", sort=True):
        ordered = group.sort_values("issn_l")
        c = official.reindex(ordered["issn_l"])
        if c.isna().any() or ordered["ais"].isna().any():
            raise AssertionError(f"2024 cube benchmark cell {cell_id} is incomplete")
        benchmark_rows.append(
            {
                "year": 2024,
                "support_scope": "annual_s6",
                "support_id": "s6",
                "cell_id": cell_id,
                **cell_statistics(c, ordered["ais"]),
            }
        )
    return pd.concat(summaries, ignore_index=True), pd.DataFrame(benchmark_rows)
