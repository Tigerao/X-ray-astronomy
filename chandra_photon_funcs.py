#!/usr/bin/env python3
"""Helper functions for ACIS timing photon extraction."""

from __future__ import annotations

from pathlib import Path
import csv
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def read_sources_csv(csvfile):
    """Read source list CSV with columns: srcid,ra,dec."""
    csv_path = Path(csvfile)
    if not csv_path.exists():
        raise FileNotFoundError(f"Sources CSV not found: {csv_path}")

    srcid, ra, dec = [], [], []
    with csv_path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        required = {"srcid", "ra", "dec"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"CSV must contain columns {sorted(required)}; got {reader.fieldnames}")
        for row in reader:
            srcid.append(int(row["srcid"]))
            ra.append(float(row["ra"]))
            dec.append(float(row["dec"]))

    return np.asarray(srcid, dtype=int), np.asarray(ra, dtype=float), np.asarray(dec, dtype=float)


def read_circle_region(regfile):
    """Read DS9 physical region line like: circle(x,y,r)."""
    reg_path = Path(regfile)
    if not reg_path.exists():
        raise FileNotFoundError(f"Region file not found: {reg_path}")

    lines = [ln.strip() for ln in reg_path.read_text().splitlines() if ln.strip()]
    circle_line = None
    for line in lines:
        if line.startswith("circle("):
            circle_line = line
            break
    if circle_line is None:
        raise ValueError(f"No circle() found in region file: {reg_path}")

    content = circle_line[circle_line.find("(") + 1 : circle_line.rfind(")")]
    parts = [p.strip() for p in content.split(",")]
    if len(parts) != 3:
        raise ValueError(f"Invalid circle region format in {reg_path}: {circle_line}")
    x, y, r = map(float, parts)
    return x, y, r


def _read_annulus_with_exclusions(regfile):
    reg_path = Path(regfile)
    if not reg_path.exists():
        raise FileNotFoundError(f"Background region file not found: {reg_path}")

    lines = [ln.strip() for ln in reg_path.read_text().splitlines() if ln.strip() and not ln.strip().startswith("#")]

    annulus = None
    exclusions = []
    for line in lines:
        if line.startswith("annulus("):
            content = line[line.find("(") + 1 : line.rfind(")")]
            x, y, rin, rout = map(float, [p.strip() for p in content.split(",")])
            annulus = (x, y, rin, rout)
        elif line.startswith("-circle("):
            content = line[line.find("(") + 1 : line.rfind(")")]
            x, y, r = map(float, [p.strip() for p in content.split(",")])
            exclusions.append((x, y, r))

    if annulus is None:
        raise ValueError(f"No annulus() line found in {reg_path}")
    return annulus, exclusions


def _sample_psf_radius(psf_data, xpix, ypix):
    ny, nx = psf_data.shape
    xi = int(np.clip(np.rint(xpix) - 1, 0, nx - 1))
    yi = int(np.clip(np.rint(ypix) - 1, 0, ny - 1))
    r = float(psf_data[yi, xi])
    if not np.isfinite(r) or r <= 0:
        raise ValueError(f"Invalid PSF radius sampled at pixel ({xpix:.2f}, {ypix:.2f}): {r}")
    return r


def _psf_radius_to_pixels(psf_header, radius_value, pixel_scale_arcsec=0.492):
    bunit = str(psf_header.get("BUNIT", "")).strip().lower()
    if "arcsec" in bunit:
        return radius_value / pixel_scale_arcsec
    return radius_value


def write_src_bkg_regions_for_obs(
    xdata_dir,
    timing_reg_dir,
    obsid,
    srcid,
    ra,
    dec,
    ecf=90,
    src_radius_factor=1.0,
    bkg_inner_factor=2.0,
    bkg_outer_factor=4.0,
    all_sources_xy_r=None,
    exclude_radius_factor=1.5,
    pixel_scale_arcsec=0.492,
):
    """Create source/background physical region files for one source in one observation."""
    xdata_dir = Path(xdata_dir)
    timing_reg_dir = Path(timing_reg_dir)

    img_file = xdata_dir / f"img_{obsid}_300_8000.fits"
    psf_file = xdata_dir / f"reproj_psf90_{obsid}_300_8000.fits"

    if not img_file.exists():
        raise FileNotFoundError(f"Missing image file for WCS: {img_file}")
    if not psf_file.exists():
        raise FileNotFoundError(f"Missing PSF map file: {psf_file}")

    with fits.open(img_file) as hdul_img:
        w = WCS(hdul_img[0].header)
        xpix, ypix = w.all_world2pix(ra, dec, 1)

    with fits.open(psf_file) as hdul_psf:
        psf_data = np.asarray(hdul_psf[0].data, dtype=float)
        if psf_data.ndim != 2:
            raise ValueError(f"PSF map is not 2D: {psf_file}")
        r_native = _sample_psf_radius(psf_data, xpix, ypix)
        r_pix = _psf_radius_to_pixels(hdul_psf[0].header, r_native, pixel_scale_arcsec=pixel_scale_arcsec)

    if r_pix <= 0:
        raise ValueError(f"Computed source radius <= 0 for srcid={srcid}, obsid={obsid}")

    src_r = r_pix * src_radius_factor
    bkg_rin = r_pix * bkg_inner_factor
    bkg_rout = r_pix * bkg_outer_factor

    out_dir = timing_reg_dir / f"region_{obsid}" / f"region_{ecf}"
    out_dir.mkdir(parents=True, exist_ok=True)

    src_reg = out_dir / f"{srcid}.reg"
    bkg_reg = out_dir / f"{srcid}_bkg.reg"

    src_reg.write_text(f"circle({xpix:.6f},{ypix:.6f},{src_r:.6f})\n")

    bkg_lines = [f"annulus({xpix:.6f},{ypix:.6f},{bkg_rin:.6f},{bkg_rout:.6f})"]
    if all_sources_xy_r is not None:
        for other_srcid, ox, oy, orad in all_sources_xy_r:
            if int(other_srcid) == int(srcid):
                continue
            dx = ox - xpix
            dy = oy - ypix
            if np.hypot(dx, dy) < (bkg_rout + orad * exclude_radius_factor):
                bkg_lines.append(
                    f"-circle({ox:.6f},{oy:.6f},{(orad * exclude_radius_factor):.6f})"
                )

    bkg_reg.write_text("\n".join(bkg_lines) + "\n")

    return {
        "src_reg": src_reg,
        "bkg_reg": bkg_reg,
        "x": xpix,
        "y": ypix,
        "psf_radius_pix": r_pix,
    }


def _load_event_columns(evtfile):
    evtfile = Path(evtfile)
    if not evtfile.exists():
        raise FileNotFoundError(f"Event file not found: {evtfile}")

    with fits.open(evtfile) as hdul:
        data = hdul[1].data
        time = np.asarray(data["TIME"], dtype=float)
        x = np.asarray(data["X"], dtype=float)
        y = np.asarray(data["Y"], dtype=float)
        energy = np.asarray(data["ENERGY"], dtype=float)
    return time, x, y, energy


def _save_photons(outfile, time, energy, obsid):
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    if time.size == 0:
        outfile.write_text("")
        return

    obs = np.full(time.size, int(obsid), dtype=int)
    arr = np.column_stack((time, energy, obs))
    arr = arr[np.argsort(arr[:, 0])]
    np.savetxt(outfile, arr, fmt="%.7f  %.3f  %d")


def extract_source_photons(evtfile, regfile, outfile, obsid, emin=300, emax=8000):
    """Extract source photons inside circular region and energy range."""
    x0, y0, r = read_circle_region(regfile)
    time, x, y, energy = _load_event_columns(evtfile)

    in_circle = (x - x0) ** 2 + (y - y0) ** 2 <= r**2
    in_energy = (energy >= emin) & (energy <= emax)
    sel = in_circle & in_energy

    _save_photons(outfile, time[sel], energy[sel], obsid)


def extract_background_photons(evtfile, bkg_regfile, outfile, obsid, emin=300, emax=8000):
    """Extract background photons in annulus and outside exclusion circles."""
    (x0, y0, rin, rout), exclusions = _read_annulus_with_exclusions(bkg_regfile)
    time, x, y, energy = _load_event_columns(evtfile)

    r2 = (x - x0) ** 2 + (y - y0) ** 2
    in_annulus = (r2 >= rin**2) & (r2 <= rout**2)

    outside_excl = np.ones_like(in_annulus, dtype=bool)
    for ex, ey, er in exclusions:
        outside_excl &= ((x - ex) ** 2 + (y - ey) ** 2) > er**2

    in_energy = (energy >= emin) & (energy <= emax)
    sel = in_annulus & outside_excl & in_energy

    _save_photons(outfile, time[sel], energy[sel], obsid)


def make_epoch_file(obsids, xdata_dir, outfile):
    """Create epoch file: tstart tstop obsid exposure from barycentered event files."""
    xdata_dir = Path(xdata_dir)
    rows = []

    for obsid in obsids:
        evt = xdata_dir / f"all_bcc_{obsid}_reproj_evt.fits"
        if not evt.exists():
            continue
        with fits.open(evt) as hdul:
            t = np.asarray(hdul[1].data["TIME"], dtype=float)
        if t.size == 0:
            continue
        tstart = float(np.min(t))
        tstop = float(np.max(t))
        rows.append((tstart, tstop, int(obsid), tstop - tstart))

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        outfile.write_text("")
        return

    arr = np.asarray(rows, dtype=float)
    arr = arr[np.argsort(arr[:, 0])]
    np.savetxt(outfile, arr, fmt="%.7f  %.7f  %d  %.7f")


def _read_txt_file(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return np.empty((0, 3), dtype=float)
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr


def merge_txt_for_source(
    srcid,
    obsids,
    timing_txt_dir,
    outdir,
    suffix="_p90",
    include_bkg=True,
):
    """Merge per-observation source/background TXT files into all-observation outputs."""
    timing_txt_dir = Path(timing_txt_dir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    src_rows = []
    bkg_rows = []
    epoch_rows = []

    for obsid in obsids:
        src_file = timing_txt_dir / f"txt_{obsid}{suffix}" / f"{srcid}.txt"
        bkg_file = timing_txt_dir / f"txt_{obsid}{suffix}" / f"{srcid}_bkg.txt"

        src_arr = _read_txt_file(src_file)
        if src_arr.size > 0:
            src_rows.append(src_arr)
            tstart = float(np.min(src_arr[:, 0]))
            tstop = float(np.max(src_arr[:, 0]))
            epoch_rows.append((tstart, tstop, int(obsid), tstop - tstart))

        if include_bkg:
            bkg_arr = _read_txt_file(bkg_file)
            if bkg_arr.size > 0:
                bkg_rows.append(bkg_arr)

    src_out = outdir / f"{srcid}.txt"
    if src_rows:
        src_all = np.vstack(src_rows)
        src_all = src_all[np.argsort(src_all[:, 0])]
        np.savetxt(src_out, src_all, fmt="%.7f  %.3f  %d")
    else:
        src_out.write_text("")

    if include_bkg:
        bkg_out = outdir / f"{srcid}_bkg.txt"
        if bkg_rows:
            bkg_all = np.vstack(bkg_rows)
            bkg_all = bkg_all[np.argsort(bkg_all[:, 0])]
            np.savetxt(bkg_out, bkg_all, fmt="%.7f  %.3f  %d")
        else:
            bkg_out.write_text("")

    epoch_out = outdir / f"epoch_src_{srcid}.txt"
    if epoch_rows:
        earr = np.asarray(epoch_rows, dtype=float)
        earr = earr[np.argsort(earr[:, 0])]
        np.savetxt(epoch_out, earr, fmt="%.7f  %.7f  %d  %.7f")
    else:
        epoch_out.write_text("")
