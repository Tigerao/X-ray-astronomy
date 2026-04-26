#!/usr/bin/env python3
"""Run independent unbinned GL-like period search from Chandra TXT timing products."""

from __future__ import annotations

import argparse
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np

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


def _gl_like_stat_for_hist(counts):
    """Poisson deviance vs constant-rate model with BIC-like penalty."""
    n = counts.astype(float)
    ntot = np.sum(n)
    m = len(n)
    if ntot <= 0 or m < 2:
        return -np.inf
    mu = ntot / m
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(n > 0, n * np.log(n / mu), 0.0)
    llr = 2.0 * np.sum(term)
    bic_penalty = (m - 1) * np.log(max(ntot, 2.0))
    return float(llr - bic_penalty)


def _search_gl_like(times, freq, mmax, tref):
    stats = np.full(len(freq), np.nan)
    best_m = np.full(len(freq), -1, dtype=int)

    for i, f in enumerate(freq):
        phase = ((times - tref) * f) % 1.0
        smax = -np.inf
        mopt = -1
        for m in range(2, mmax + 1):
            hist, _ = np.histogram(phase, bins=np.linspace(0, 1, m + 1))
            s = _gl_like_stat_for_hist(hist)
            if s > smax:
                smax = s
                mopt = m
        stats[i] = smax
        best_m[i] = mopt
    return stats, best_m


def _simulate_uniform_in_epochs(n, tstart, tstop, rng):
    durations = np.clip(tstop - tstart, 0.0, None)
    total = np.sum(durations)
    if total <= 0:
        return np.array([], dtype=float)
    probs = durations / total
    choice = rng.choice(len(durations), size=n, p=probs)
    u = rng.random(n)
    return tstart[choice] + u * durations[choice]


def _fold_events_with_exposure(times, epochs, period, nphase):
    phase = (times / period) % 1.0
    bins = np.linspace(0, 1, nphase + 1)
    counts, _ = np.histogram(phase, bins=bins)

    lc = make_even_lightcurve_from_events(times, epochs, dt=max(period / (nphase * 8), 1e-6))
    valid = lc["exposure_per_bin"] > 0
    ph_lc = (lc["t_mid"][valid] / period) % 1.0
    exp_counts, _ = np.histogram(ph_lc, bins=bins, weights=lc["exposure_per_bin"][valid])

    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.divide(counts, exp_counts, out=np.full(nphase, np.nan), where=exp_counts > 0)
        rate_err = np.divide(np.sqrt(counts), exp_counts, out=np.full(nphase, np.nan), where=exp_counts > 0)

    phase_center = 0.5 * (bins[:-1] + bins[1:])
    return phase_center, counts.astype(float), exp_counts.astype(float), rate, rate_err


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True)
    parser.add_argument("--srcid", required=True)
    parser.add_argument("--suffix", default="_p90")
    parser.add_argument("--emin", type=float, default=500)
    parser.add_argument("--emax", type=float, default=8000)
    parser.add_argument("--pmin", type=float, default=50)
    parser.add_argument("--pmax", type=float, default=50000)
    parser.add_argument("--oversample", type=float, default=5)
    parser.add_argument("--mmax", type=int, default=12)
    parser.add_argument("--nphase", type=int, default=16)
    parser.add_argument("--nsim", type=int, default=0)
    args = parser.parse_args()

    t_run0 = time.perf_counter()
    paths = get_default_paths(args.base, args.srcid, suffix=args.suffix)
    outdir = paths["outdir"]
    outdir.mkdir(parents=True, exist_ok=True)

    epochs = load_epoch_txt(paths["epoch"])
    src = filter_time_by_epochs(filter_energy(load_photon_txt(paths["src"]), args.emin, args.emax), epochs)
    bkg = filter_time_by_epochs(filter_energy(load_photon_txt(paths["bkg"]), args.emin, args.emax), epochs)

    summary = summarize_timing_data(src, bkg, epochs)
    summary.update({
        "method": "GL-like",
        "is_exact_gregory_loredo": False,
        "srcid": str(args.srcid),
        "emin": args.emin,
        "emax": args.emax,
        "pmin": args.pmin,
        "pmax": args.pmax,
        "oversample": args.oversample,
        "mmax": args.mmax,
        "nphase": args.nphase,
        "nsim": args.nsim,
        "warnings": ["This is a GL-like period search, not a full Gregory-Loredo odds-ratio implementation."],
    })

    if summary["n_source_events"] == 0:
        raise RuntimeError("No source events remain after filtering; cannot run GL-like search.")

    t = np.asarray(src["time"], dtype=float)
    baseline = summary["baseline_s"]

    if args.pmax > baseline / 2.0:
        msg = "pmax > baseline/2"
        warnings.warn(msg)
        summary["warnings"].append(msg)

    freq = _get_frequency_grid(baseline, args.pmin, args.pmax, args.oversample)
    tref = np.min(epochs["tstart"])

    gl_stat, best_m = _search_gl_like(t, freq, args.mmax, tref)

    i = int(np.nanargmax(gl_stat))
    best_freq = float(freq[i])
    best_period = 1.0 / best_freq
    best_stat = float(gl_stat[i])

    if summary["total_exposure_s"] / best_period < 3:
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

    empirical_fap = None
    if args.nsim > 0:
        rng = np.random.default_rng(12345)
        max_stats = []
        for _ in range(args.nsim):
            tsim = _simulate_uniform_in_epochs(len(t), epochs["tstart"], epochs["tstop"], rng)
            sim_stat, _ = _search_gl_like(tsim, freq, args.mmax, tref)
            max_stats.append(np.nanmax(sim_stat))
        max_stats = np.asarray(max_stats)
        empirical_fap = float(np.mean(max_stats >= best_stat))

    period = 1.0 / freq
    np.savetxt(
        outdir / f"src_{args.srcid}_gl_periodogram.txt",
        np.column_stack([freq, period, gl_stat, best_m]),
        header="frequency period gl_stat best_m",
    )

    ph, counts, exposure, rate, rate_err = _fold_events_with_exposure(t, epochs, best_period, args.nphase)
    ph2 = np.concatenate([ph, ph + 1.0])
    c2 = np.concatenate([counts, counts])
    e2 = np.concatenate([exposure, exposure])
    r2 = np.concatenate([rate, rate])
    re2 = np.concatenate([rate_err, rate_err])

    np.savetxt(
        outdir / f"src_{args.srcid}_gl_folded.txt",
        np.column_stack([ph2, c2, e2, r2, re2]),
        header="phase counts exposure rate rate_err",
    )

    runtime_s = float(time.perf_counter() - t_run0)
    w = 2.0 * np.pi * freq
    wpeak = float(w[i])
    period_from_wpeak = float(2.0 * np.pi / wpeak)

    # Weight-based mean/confidence interval from non-negative GL-like stats
    sshift = gl_stat - np.nanmin(gl_stat)
    sshift = np.clip(sshift, 0.0, None)
    if np.sum(sshift) > 0:
        wmean = float(np.sum(w * sshift) / np.sum(sshift))
        cdf = np.cumsum(sshift / np.sum(sshift))
        i_lo = int(np.searchsorted(cdf, 0.025))
        i_hi = int(np.searchsorted(cdf, 0.975))
        i_lo = min(max(i_lo, 0), len(w) - 1)
        i_hi = min(max(i_hi, 0), len(w) - 1)
        wconf_lo = float(w[i_lo])
        wconf_hi = float(w[i_hi])
    else:
        wmean = None
        wconf_lo = None
        wconf_hi = None

    if empirical_fap is not None:
        prob = float(max(0.0, 1.0 - empirical_fap))
    else:
        prob = None

    result_row = [
        int(args.srcid),
        runtime_s,
        prob,
        period_from_wpeak,
        wpeak,
        wmean,
        int(best_m[i]),
        wconf_lo,
        wconf_hi,
        int(summary["n_source_events"]),
    ]

    summary.update(
        {
            "best_frequency_hz": best_freq,
            "best_period_s": best_period,
            "best_gl_stat": best_stat,
            "best_m": int(best_m[i]),
            "n_freq": int(len(freq)),
            "empirical_fap": empirical_fap,
            "runtime_s": runtime_s,
            "result_vector": result_row,
            "result_vector_definition": [
                "srcid", "runtime_s", "Prob", "2pi_over_wpeak_s", "wpeak_rad_s",
                "wmean_rad_s", "mopt", "wconf_lo_rad_s", "wconf_hi_rad_s", "counts"
            ],
        }
    )
    save_summary_json(summary, outdir / f"src_{args.srcid}_gl_summary.json")

    np.savetxt(
        outdir / f"src_{args.srcid}_gl_result.txt",
        np.array([result_row], dtype=object),
        fmt="%s",
        header="srcid runtime_s Prob 2pi_over_wpeak_s wpeak_rad_s wmean_rad_s mopt wconf_lo_rad_s wconf_hi_rad_s counts",
    )

    fig = plt.figure(figsize=(8, 5))
    plt.plot(period, gl_stat, lw=1)
    plt.xlabel("Period (s)")
    plt.ylabel("GL-like statistic")
    if args.pmax / args.pmin > 20:
        plt.xscale("log")
    plt.tight_layout()
    fig.savefig(outdir / f"src_{args.srcid}_gl_periodogram.png", dpi=150)
    plt.close(fig)

    fig = plt.figure(figsize=(8, 5))
    plt.errorbar(ph2, r2, yerr=re2, fmt="o", ms=4, capsize=2)
    plt.xlim(0, 2)
    plt.xlabel("Phase")
    plt.ylabel("Rate (counts/s)")
    plt.tight_layout()
    fig.savefig(outdir / f"src_{args.srcid}_gl_folded.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
