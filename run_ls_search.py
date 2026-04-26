#!/usr/bin/env python3
"""Run Lomb-Scargle period search from Chandra TXT timing products."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle

from timing_io import (
    filter_energy,
    filter_time_by_epochs,
    get_default_paths,
    load_epoch_txt,
    load_photon_txt,
    make_even_lightcurve_from_events,
    save_summary_json,
    summarize_timing_data,
)

DITHER_PERIODS = (706.96, 999.96)


def _get_frequency_grid(baseline, pmin, pmax, oversample):
    if baseline <= 0:
        raise ValueError("Cannot build frequency grid: baseline <= 0")
    if pmin <= 0 or pmax <= 0 or pmax <= pmin:
        raise ValueError("Require 0 < pmin < pmax")
    fmin = 1.0 / pmax
    fmax = 1.0 / pmin
    df = 1.0 / (baseline * oversample)
    n = int(np.floor((fmax - fmin) / df)) + 1
    if n < 2:
        raise ValueError("Frequency grid too small; adjust period bounds or oversample")
    return fmin + np.arange(n) * df


def _dither_warning(best_period):
    msgs = []
    for pd in DITHER_PERIODS:
        for k in range(1, 7):
            for rel, label in ((k * pd, "harmonic"), (pd / k, "subharmonic")):
                frac = abs(best_period - rel) / rel
                if frac <= 0.02:
                    msgs.append(
                        f"Best period {best_period:.3f}s is within {frac*100:.2f}% of dither {label} {rel:.3f}s (base {pd:.2f}s)."
                    )
    return sorted(set(msgs))


def _fold_profile_from_lc(times, counts, exposure, period, nphase):
    phase = (times / period) % 1.0
    bins = np.linspace(0, 1, nphase + 1)
    idx = np.clip(np.digitize(phase, bins) - 1, 0, nphase - 1)

    csum = np.bincount(idx, weights=counts, minlength=nphase).astype(float)
    esum = np.bincount(idx, weights=exposure, minlength=nphase).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.divide(csum, esum, out=np.full(nphase, np.nan), where=esum > 0)
        rate_err = np.divide(np.sqrt(csum), esum, out=np.full(nphase, np.nan), where=esum > 0)
    phase_center = 0.5 * (bins[:-1] + bins[1:])
    return phase_center, csum, esum, rate, rate_err


def _save_two_cycle_folded_txt(outfile, phase, rate, rate_err, counts, exposure):
    phase2 = np.concatenate([phase, phase + 1.0])
    rate2 = np.concatenate([rate, rate])
    err2 = np.concatenate([rate_err, rate_err])
    c2 = np.concatenate([counts, counts])
    e2 = np.concatenate([exposure, exposure])
    arr = np.column_stack([phase2, rate2, err2, c2, e2])
    np.savetxt(outfile, arr, header="phase rate rate_err counts exposure")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True)
    parser.add_argument("--srcid", required=True)
    parser.add_argument("--suffix", default="_p90")
    parser.add_argument("--emin", type=float, default=500)
    parser.add_argument("--emax", type=float, default=8000)
    parser.add_argument("--dt", type=float, default=50)
    parser.add_argument("--pmin", type=float, default=50)
    parser.add_argument("--pmax", type=float, default=50000)
    parser.add_argument("--oversample", type=float, default=5)
    parser.add_argument("--nphase", type=int, default=16)
    args = parser.parse_args()

    paths = get_default_paths(args.base, args.srcid, suffix=args.suffix)
    outdir = paths["outdir"]
    outdir.mkdir(parents=True, exist_ok=True)

    src = load_photon_txt(paths["src"])
    bkg = load_photon_txt(paths["bkg"])
    epochs = load_epoch_txt(paths["epoch"])

    src = filter_time_by_epochs(filter_energy(src, args.emin, args.emax), epochs)
    bkg = filter_time_by_epochs(filter_energy(bkg, args.emin, args.emax), epochs)

    summary = summarize_timing_data(src, bkg, epochs)
    summary.update({
        "method": "LombScargle",
        "srcid": str(args.srcid),
        "emin": args.emin,
        "emax": args.emax,
        "dt": args.dt,
        "pmin": args.pmin,
        "pmax": args.pmax,
        "oversample": args.oversample,
        "nphase": args.nphase,
        "warnings": [],
    })

    if summary["n_source_events"] == 0:
        raise RuntimeError("No source events remain after filtering; cannot run LS search.")

    baseline = summary["baseline_s"]
    if args.pmax > baseline / 2.0:
        warnings.warn("pmax > baseline/2. Sensitivity and interpretation at long periods may be weak.")
        summary["warnings"].append("pmax > baseline/2")

    lc = make_even_lightcurve_from_events(src["time"], epochs, dt=args.dt)
    valid = lc["exposure_per_bin"] > 0
    if not np.any(valid):
        raise RuntimeError("All bins have zero exposure; check epoch file or dt.")

    t = lc["t_mid"][valid]
    y = lc["counts"][valid] / lc["exposure_per_bin"][valid]

    freq = _get_frequency_grid(baseline, args.pmin, args.pmax, args.oversample)
    ls = LombScargle(t, y)
    power = ls.power(freq)

    ibest = int(np.nanargmax(power))
    best_freq = float(freq[ibest])
    best_period = 1.0 / best_freq
    best_power = float(power[ibest])

    fap = None
    try:
        fap = float(
            ls.false_alarm_probability(
                best_power, minimum_frequency=freq[0], maximum_frequency=freq[-1], method="baluev"
            )
        )
    except Exception:
        pass

    total_exposure = summary["total_exposure_s"]
    if total_exposure / best_period < 3:
        msg = "total exposure covers fewer than 3 cycles at best period"
        warnings.warn(msg)
        summary["warnings"].append(msg)
    if baseline / best_period < 3:
        msg = "baseline covers fewer than 3 cycles at best period"
        warnings.warn(msg)
        summary["warnings"].append(msg)

    for msg in _dither_warning(best_period):
        warnings.warn(msg)
        summary["warnings"].append(msg)

    per = 1.0 / freq
    np.savetxt(
        outdir / f"src_{args.srcid}_ls_periodogram.txt",
        np.column_stack([freq, per, power]),
        header="frequency period power",
    )

    ph, csum, esum, rate, rate_err = _fold_profile_from_lc(
        lc["t_mid"], lc["counts"], lc["exposure_per_bin"], best_period, args.nphase
    )
    _save_two_cycle_folded_txt(
        outdir / f"src_{args.srcid}_ls_folded.txt", ph, rate, rate_err, csum, esum
    )

    summary.update(
        {
            "best_frequency_hz": best_freq,
            "best_period_s": best_period,
            "best_power": best_power,
            "false_alarm_probability": fap,
            "n_freq": int(len(freq)),
        }
    )
    save_summary_json(summary, outdir / f"src_{args.srcid}_ls_summary.json")

    fig = plt.figure(figsize=(8, 5))
    plt.plot(per, power, lw=1)
    plt.xlabel("Period (s)")
    plt.ylabel("LS power")
    if args.pmax / args.pmin > 20:
        plt.xscale("log")
    plt.tight_layout()
    fig.savefig(outdir / f"src_{args.srcid}_ls_periodogram.png", dpi=150)
    plt.close(fig)

    fig = plt.figure(figsize=(8, 4))
    plt.step(lc["t_mid"], lc["counts"], where="mid", lw=0.8)
    plt.xlabel("Time (s)")
    plt.ylabel(f"Counts per {args.dt:g}s bin")
    plt.tight_layout()
    fig.savefig(outdir / f"src_{args.srcid}_lightcurve.png", dpi=150)
    plt.close(fig)

    ph2 = np.concatenate([ph, ph + 1.0])
    rate2 = np.concatenate([rate, rate])
    err2 = np.concatenate([rate_err, rate_err])
    fig = plt.figure(figsize=(8, 5))
    plt.errorbar(ph2, rate2, yerr=err2, fmt="o", ms=4, capsize=2)
    plt.xlim(0, 2)
    plt.xlabel("Phase")
    plt.ylabel("Rate (counts/s)")
    plt.tight_layout()
    fig.savefig(outdir / f"src_{args.srcid}_ls_folded.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
