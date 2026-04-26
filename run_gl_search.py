#!/usr/bin/env python3
"""Run independent unbinned exact Gregory-Loredo odds period search from Chandra TXT timing products."""

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


def _precompute_log_factorial(nmax):
    f = np.zeros(nmax + 1, dtype=float)
    for i in range(2, nmax + 1):
        f[i] = f[i - 1] + np.log(i)
    return f


def _compute_bin_counts(times, m, w, phi):
    phase = np.mod(w * times + phi, 2.0 * np.pi) / (2.0 * np.pi)
    idx = np.floor(m * phase).astype(int)
    idx = np.clip(idx, 0, m - 1)
    return np.bincount(idx, minlength=m).astype(int)


def _compute_factor(N, m, v):
    # log( m^N * (m-1)! / ((N+m-1)! * 2*pi*v) )
    f1 = N * np.log(m)
    f2 = np.sum(np.log(np.arange(1, N + m)))
    f3 = np.sum(np.log(np.arange(1, m))) if m > 1 else 0.0
    return f1 + f3 - f2 - np.log(2.0 * np.pi * v)


def _get_T_in_mbins(epoch_info, w, m, phi):
    tstart = np.asarray(epoch_info["tstart"], dtype=float)
    tstop = np.asarray(epoch_info["tstop"], dtype=float)
    T = 2.0 * np.pi / w
    tbin = T / m
    T_in = np.zeros(m, dtype=float)

    nstart = tstart / tbin + m * phi / (2.0 * np.pi)
    nstop = tstop / tbin + m * phi / (2.0 * np.pi)
    istart = np.floor(nstart).astype(int) + 1
    istop = np.floor(nstop).astype(int)

    for i in range(len(nstart)):
        if istop[i] >= istart[i]:
            T_in += int((istop[i] - istart[i]) / m) * tbin
            T_in[(istart[i] % m) - 1] += (istart[i] - nstart[i]) * tbin
            T_in[istop[i] % m] += (nstop[i] - istop[i]) * tbin
            rest = (istop[i] - istart[i]) % m
            for k in range(rest):
                T_in[(istart[i] + k) % m] += tbin
        else:
            T_in[(istart[i] % m) - 1] += (nstop[i] - nstart[i]) * tbin
    return T_in


def _compute_ln_S(epoch_info, times, w, m, phi):
    n = _compute_bin_counts(times, m, w, phi).astype(float)
    tau = _get_T_in_mbins(epoch_info, w, m, phi)
    tau_norm = tau / (np.sum(tau) / m)
    tau_norm = np.clip(tau_norm, 1e-300, None)
    return -np.sum(n * np.log(tau_norm))


def _compute_om1_single_w(times, epoch_info, m, w, log_factor, log_fact, ni):
    pgrid = np.arange(ni, dtype=float) / float(ni) * 2.0 * np.pi / m
    vals = np.zeros_like(pgrid)
    for i, phi in enumerate(pgrid):
        n = _compute_bin_counts(times, m, w, phi)
        log_mult = np.sum(log_fact[n])
        ln_s = _compute_ln_S(epoch_info, times, w, m, phi)
        vals[i] = np.exp(log_mult + log_factor + ln_s)
    return np.trapz(vals, pgrid) * m


def compute_gl_exact(times, epoch_info, w_range, m_max=12, ni=10):
    N = len(times)
    if N <= 0:
        raise ValueError("No events for GL search")
    if m_max < 2:
        raise ValueError("m_max must be >=2")

    v = m_max - 1
    log_fact = _precompute_log_factorial(N)
    factors = np.array([_compute_factor(N, m, v) for m in range(1, m_max + 1)], dtype=float)

    Om1w = np.zeros((m_max, len(w_range)), dtype=float)
    for mi, m in enumerate(range(1, m_max + 1)):
        for wi, w in enumerate(w_range):
            Om1w[mi, wi] = _compute_om1_single_w(times, epoch_info, m, w, factors[mi], log_fact, ni)

    w_lo = np.min(w_range)
    w_hi = np.max(w_range)
    pw = 1.0 / (w_range * np.log(w_hi / w_lo))

    O1m = np.array([np.trapz(pw * Om1w[mi], w_range) for mi in range(m_max)], dtype=float)
    # periodic vs constant model (exclude m=1 term)
    O_period = float(np.sum(O1m[1:]))
    p_period = float(O_period / (1.0 + O_period))

    m_opt_global = int(np.argmax(O1m)) + 1

    # Spectral probability for optimal m
    S_raw = Om1w[m_opt_global - 1] / w_range
    C = np.trapz(S_raw, w_range)
    S = S_raw / C if C > 0 else np.full_like(S_raw, np.nan)

    cdf = np.zeros_like(S)
    for i in range(len(S)):
        cdf[i] = np.trapz(S[: i + 1], w_range[: i + 1])

    w_peak = float(w_range[int(np.nanargmax(S))])
    w_mean = float(np.trapz(S * w_range, w_range))

    wr = w_range[(cdf > 0.025) & (cdf < 0.975)]
    if len(wr) > 0:
        w_conf = [float(np.min(wr)), float(np.max(wr))]
    else:
        dw = w_range[1] - w_range[0] if len(w_range) > 1 else 0.0
        w_conf = [float(w_peak - dw), float(w_peak + dw)]

    best_m_per_w = np.argmax(Om1w, axis=0) + 1

    return {
        "O_period": O_period,
        "p_period": p_period,
        "m_opt_global": m_opt_global,
        "S": S,
        "w": w_range,
        "w_peak": w_peak,
        "w_mean": w_mean,
        "w_conf": w_conf,
        "cdf": cdf,
        "Om1w": Om1w,
        "best_m_per_w": best_m_per_w,
    }


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
    parser.add_argument("--ni", type=int, default=10, help="phase-offset integration steps for exact GL")
    args = parser.parse_args()

    t_run0 = time.perf_counter()
    paths = get_default_paths(args.base, args.srcid, suffix=args.suffix)
    outdir = paths["outdir"]
    outdir.mkdir(parents=True, exist_ok=True)

    epochs = load_epoch_txt(paths["epoch"])
    src = filter_time_by_epochs(filter_energy(load_photon_txt(paths["src"]), args.emin, args.emax), epochs)
    bkg = filter_time_by_epochs(filter_energy(load_photon_txt(paths["bkg"]), args.emin, args.emax), epochs)

    summary = summarize_timing_data(src, bkg, epochs)
    summary.update(
        {
            "method": "Gregory-Loredo exact odds",
            "is_exact_gregory_loredo": True,
            "srcid": str(args.srcid),
            "emin": args.emin,
            "emax": args.emax,
            "pmin": args.pmin,
            "pmax": args.pmax,
            "oversample": args.oversample,
            "mmax": args.mmax,
            "nphase": args.nphase,
            "nsim": args.nsim,
            "ni": args.ni,
            "warnings": [],
        }
    )

    if summary["n_source_events"] == 0:
        raise RuntimeError("No source events remain after filtering; cannot run GL search.")

    t = np.asarray(src["time"], dtype=float)
    baseline = summary["baseline_s"]

    if args.pmax > baseline / 2.0:
        msg = "pmax > baseline/2"
        warnings.warn(msg)
        summary["warnings"].append(msg)

    freq = _get_frequency_grid(baseline, args.pmin, args.pmax, args.oversample)
    w = 2.0 * np.pi * freq

    gl = compute_gl_exact(t, epochs, w, m_max=args.mmax, ni=args.ni)
    S = gl["S"]
    best_i = int(np.nanargmax(S))
    best_w = float(gl["w_peak"])
    best_period = float(2.0 * np.pi / best_w)
    best_freq = float(best_w / (2.0 * np.pi))

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
        max_s = []
        for _ in range(args.nsim):
            tsim = _simulate_uniform_in_epochs(len(t), epochs["tstart"], epochs["tstop"], rng)
            gl_sim = compute_gl_exact(tsim, epochs, w, m_max=args.mmax, ni=args.ni)
            max_s.append(np.nanmax(gl_sim["S"]))
        empirical_fap = float(np.mean(np.asarray(max_s) >= np.nanmax(S)))

    period = 1.0 / freq
    np.savetxt(
        outdir / f"src_{args.srcid}_gl_periodogram.txt",
        np.column_stack([freq, period, S, gl["best_m_per_w"]]),
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
    prob = float(gl["p_period"])
    result_row = [
        int(args.srcid),
        runtime_s,
        prob,
        float(2.0 * np.pi / gl["w_peak"]),
        float(gl["w_peak"]),
        float(gl["w_mean"]),
        int(gl["m_opt_global"]),
        float(gl["w_conf"][0]),
        float(gl["w_conf"][1]),
        int(summary["n_source_events"]),
    ]

    summary.update(
        {
            "best_frequency_hz": best_freq,
            "best_period_s": best_period,
            "best_gl_stat": float(np.nanmax(S)),
            "best_m": int(gl["best_m_per_w"][best_i]),
            "m_opt_global": int(gl["m_opt_global"]),
            "n_freq": int(len(freq)),
            "empirical_fap": empirical_fap,
            "runtime_s": runtime_s,
            "O_period": float(gl["O_period"]),
            "p_period": float(gl["p_period"]),
            "w_peak": float(gl["w_peak"]),
            "w_mean": float(gl["w_mean"]),
            "w_conf": [float(gl["w_conf"][0]), float(gl["w_conf"][1])],
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
    plt.plot(period, S, lw=1)
    plt.xlabel("Period (s)")
    plt.ylabel("GL odds spectral probability")
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
