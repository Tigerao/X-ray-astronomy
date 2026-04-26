#!/usr/bin/env python3
"""Utilities for loading and preparing Chandra timing-search TXT products."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Any, Optional

import numpy as np


def get_default_paths(base, srcid, suffix="_p90"):
    """Return default source/background/epoch/output paths for a source id."""
    base = Path(base)
    srcid_str = str(srcid)
    txt_dir = base / "merge_data" / "timing" / "txt" / f"txt_all_obs{suffix}"
    out_dir = base / "merge_data" / "timing" / "period_search" / f"src_{srcid_str}"
    return {
        "src": txt_dir / f"{srcid_str}.txt",
        "bkg": txt_dir / f"{srcid_str}_bkg.txt",
        "epoch": txt_dir / f"epoch_src_{srcid_str}.txt",
        "outdir": out_dir,
    }


def _load_numeric_txt(filename, expected_cols, names):
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"Required input file not found: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"Input path is not a regular file: {path}")

    try:
        arr = np.genfromtxt(path, comments="#", dtype=float)
    except Exception as exc:
        raise ValueError(f"Failed reading numeric data from {path}: {exc}") from exc

    if arr.size == 0:
        return {k: np.array([], dtype=float) for k in names}

    arr = np.atleast_2d(arr)
    if arr.shape[1] != expected_cols:
        raise ValueError(
            f"File {path} has {arr.shape[1]} columns, expected {expected_cols}."
        )

    arr = arr[~np.any(~np.isfinite(arr), axis=1)]
    if arr.size == 0:
        return {k: np.array([], dtype=float) for k in names}

    return {name: arr[:, i].astype(float) for i, name in enumerate(names)}


def load_photon_txt(filename):
    """Load source/background photon TXT with columns: time, energy, obsid."""
    return _load_numeric_txt(filename, expected_cols=3, names=("time", "energy", "obsid"))


def load_epoch_txt(filename):
    """Load epoch TXT with columns: tstart, tstop, obsid, exposure."""
    epochs = _load_numeric_txt(
        filename, expected_cols=4, names=("tstart", "tstop", "obsid", "exposure")
    )
    if len(epochs["tstart"]) > 0:
        good = epochs["tstop"] > epochs["tstart"]
        for k in epochs:
            epochs[k] = epochs[k][good]
    return epochs


def filter_energy(data, emin=500, emax=8000):
    """Filter photon dictionary by energy band [emin, emax]."""
    if len(data["time"]) == 0:
        return {k: np.array([], dtype=float) for k in data}
    mask = (data["energy"] >= emin) & (data["energy"] <= emax)
    return {k: np.asarray(v)[mask] for k, v in data.items()}


def make_gti_mask(times, tstart, tstop):
    """Mask events that fall within any GTI/epoch interval [tstart, tstop)."""
    times = np.asarray(times, dtype=float)
    tstart = np.asarray(tstart, dtype=float)
    tstop = np.asarray(tstop, dtype=float)
    mask = np.zeros(times.shape, dtype=bool)
    if times.size == 0 or tstart.size == 0:
        return mask
    for ts, te in zip(tstart, tstop):
        if te <= ts:
            continue
        mask |= (times >= ts) & (times < te)
    return mask


def filter_time_by_epochs(data, epochs):
    """Filter photons by epoch windows."""
    if len(data["time"]) == 0:
        return {k: np.array([], dtype=float) for k in data}
    mask = make_gti_mask(data["time"], epochs["tstart"], epochs["tstop"])
    return {k: np.asarray(v)[mask] for k, v in data.items()}


def summarize_timing_data(src, bkg, epochs):
    """Produce a compact summary dictionary for timing inputs."""
    nsrc = int(len(src["time"]))
    nbkg = int(len(bkg["time"]))
    nepoch = int(len(epochs["tstart"]))

    if nepoch > 0:
        baseline = float(np.max(epochs["tstop"]) - np.min(epochs["tstart"]))
        total_exposure = float(np.sum(np.maximum(0.0, epochs["tstop"] - epochs["tstart"])))
    else:
        baseline = 0.0
        total_exposure = 0.0

    return {
        "n_source_events": nsrc,
        "n_background_events": nbkg,
        "n_epochs": nepoch,
        "baseline_s": baseline,
        "total_exposure_s": total_exposure,
    }


def save_summary_json(summary, outfile):
    """Write summary dictionary to JSON file."""
    out = Path(outfile)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)


def make_even_lightcurve_from_events(times, epochs, dt, t0: Optional[float] = None):
    """
    Build an evenly sampled light curve and exact per-bin exposure from event times + epochs.

    Important behavior:
    - Uses bins over the full baseline so gaps remain explicit.
    - Counts are from events only (no artificial filling).
    - exposure_per_bin is exact overlap of each bin with all epoch windows.
    """
    times = np.asarray(times, dtype=float)
    tstart = np.asarray(epochs["tstart"], dtype=float)
    tstop = np.asarray(epochs["tstop"], dtype=float)

    if dt <= 0:
        raise ValueError("dt must be > 0")
    if tstart.size == 0:
        return {
            "t_left": np.array([], dtype=float),
            "t_right": np.array([], dtype=float),
            "t_mid": np.array([], dtype=float),
            "counts": np.array([], dtype=float),
            "exposure_per_bin": np.array([], dtype=float),
            "dt": float(dt),
        }

    start = np.min(tstart) if t0 is None else float(t0)
    stop = np.max(tstop)
    edges = np.arange(start, stop + dt, dt, dtype=float)
    if len(edges) < 2:
        edges = np.array([start, start + dt], dtype=float)

    counts, _ = np.histogram(times, bins=edges)
    exposure = np.zeros(len(edges) - 1, dtype=float)

    left = edges[:-1]
    right = edges[1:]

    for ts, te in zip(tstart, tstop):
        if te <= ts:
            continue
        overlap_left = np.maximum(left, ts)
        overlap_right = np.minimum(right, te)
        exposure += np.clip(overlap_right - overlap_left, 0.0, None)

    return {
        "t_left": left,
        "t_right": right,
        "t_mid": 0.5 * (left + right),
        "counts": counts.astype(float),
        "exposure_per_bin": exposure,
        "dt": float(dt),
    }
