#!/usr/bin/env python3
"""Simulation-calibrated global null for a frozen family of low-l statistics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
Alternative = str


@dataclass(frozen=True)
class GlobalNullResult:
    observed_local_p: FloatArray
    observed_max_score: float
    simulation_max_scores: FloatArray
    global_p: float


def _tail_p(value: float, reference: FloatArray, alternative: Alternative) -> float:
    denominator = len(reference) + 1.0
    upper = (1.0 + float(np.sum(reference >= value))) / denominator
    lower = (1.0 + float(np.sum(reference <= value))) / denominator
    if alternative == "greater":
        return upper
    if alternative == "less":
        return lower
    if alternative == "two-sided":
        return min(1.0, 2.0 * min(lower, upper))
    raise ValueError(f"Unknown alternative: {alternative}")


def calibrate_global_null(
    observed_statistics: FloatArray,
    simulated_statistics: FloatArray,
    alternatives: Sequence[Alternative],
) -> GlobalNullResult:
    """Return local ranks and a leave-one-out simulation maximum-statistic p-value."""
    observed = np.asarray(observed_statistics, dtype=float)
    simulations = np.asarray(simulated_statistics, dtype=float)
    if observed.ndim != 1 or simulations.ndim != 2 or simulations.shape[1] != len(observed):
        raise ValueError("Observed and simulated statistic arrays have incompatible shapes")
    if simulations.shape[0] < 20:
        raise ValueError("At least 20 simulations are required for a global null")
    if len(alternatives) != len(observed):
        raise ValueError("Every statistic requires a frozen tail alternative")
    if not np.all(np.isfinite(observed)) or not np.all(np.isfinite(simulations)):
        raise ValueError("Global-null inputs must be finite")

    observed_p = np.array(
        [
            _tail_p(observed[index], simulations[:, index], alternatives[index])
            for index in range(len(observed))
        ],
        dtype=float,
    )
    observed_score = float(np.max(-np.log10(observed_p)))

    simulation_local_p = np.empty_like(simulations)
    simulation_count = simulations.shape[0]
    for index, alternative in enumerate(alternatives):
        column = simulations[:, index]
        ordered = np.sort(column)
        upper = (simulation_count - np.searchsorted(ordered, column, side="left")) / simulation_count
        lower = np.searchsorted(ordered, column, side="right") / simulation_count
        if alternative == "greater":
            simulation_local_p[:, index] = upper
        elif alternative == "less":
            simulation_local_p[:, index] = lower
        elif alternative == "two-sided":
            simulation_local_p[:, index] = np.minimum(1.0, 2.0 * np.minimum(lower, upper))
        else:
            raise ValueError(f"Unknown alternative: {alternative}")
    simulation_scores = np.max(-np.log10(simulation_local_p), axis=1)
    global_p = (1.0 + float(np.sum(simulation_scores >= observed_score))) / (
        len(simulation_scores) + 1.0
    )
    return GlobalNullResult(
        observed_local_p=observed_p,
        observed_max_score=observed_score,
        simulation_max_scores=simulation_scores,
        global_p=global_p,
    )
