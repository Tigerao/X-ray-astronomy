#!/usr/bin/env python3
"""Command-line ACIS photon extraction pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from chandra_photon_funcs import (
    read_sources_csv,
    write_src_bkg_regions_for_obs,
    extract_source_photons,
    extract_background_photons,
    make_epoch_file,
    merge_txt_for_source,
)


def _default_paths(base: Path):
    xdata_dir = base / "merge_data" / "xdata"
    timing_dir = base / "merge_data" / "timing"
    reg_dir = timing_dir / "reg"
    txt_dir = timing_dir / "txt"
    return xdata_dir, timing_dir, reg_dir, txt_dir


def _prepare_all_source_geometry(xdata_dir: Path, obsid: int, srcids, ras, decs, pixel_scale_arcsec=0.492):
    """Optional helper: derive per-source (x,y,r_pix) list for exclusion circles."""
    img_file = xdata_dir / f"img_{obsid}_300_8000.fits"
    psf_file = xdata_dir / f"reproj_psf90_{obsid}_300_8000.fits"
    if not (img_file.exists() and psf_file.exists()):
        return None

    with fits.open(img_file) as himg, fits.open(psf_file) as hpsf:
        w = WCS(himg[0].header)
        psf_data = np.asarray(hpsf[0].data, dtype=float)
        bunit = str(hpsf[0].header.get("BUNIT", "")).lower()
        ny, nx = psf_data.shape

        records = []
        for srcid, ra, dec in zip(srcids, ras, decs):
            x, y = w.all_world2pix(ra, dec, 1)
            xi = int(np.clip(np.rint(x) - 1, 0, nx - 1))
            yi = int(np.clip(np.rint(y) - 1, 0, ny - 1))
            r = float(psf_data[yi, xi])
            if not np.isfinite(r) or r <= 0:
                continue
            if "arcsec" in bunit:
                r = r / pixel_scale_arcsec
            records.append((int(srcid), float(x), float(y), float(r)))
        return records


def cmd_make_regions(args):
    base = Path(args.base)
    xdata_dir, _, reg_dir, _ = _default_paths(base)

    srcids, ras, decs = read_sources_csv(args.sources)

    n_ok = 0
    n_fail = 0
    for obsid in args.obsids:
        all_geom = _prepare_all_source_geometry(xdata_dir, obsid, srcids, ras, decs)
        for srcid, ra, dec in zip(srcids, ras, decs):
            try:
                write_src_bkg_regions_for_obs(
                    xdata_dir=xdata_dir,
                    timing_reg_dir=reg_dir,
                    obsid=obsid,
                    srcid=srcid,
                    ra=ra,
                    dec=dec,
                    ecf=args.ecf,
                    src_radius_factor=args.src_radius_factor,
                    bkg_inner_factor=args.bkg_inner_factor,
                    bkg_outer_factor=args.bkg_outer_factor,
                    all_sources_xy_r=all_geom if args.exclude_neighbors else None,
                    exclude_radius_factor=args.exclude_radius_factor,
                    pixel_scale_arcsec=args.pixel_scale_arcsec,
                )
                n_ok += 1
            except Exception as exc:
                n_fail += 1
                print(f"[WARN] make-regions failed for obsid={obsid} srcid={srcid}: {exc}")

    print(f"make-regions summary: wrote {n_ok} source/background region pairs; failures={n_fail}")


def cmd_extract(args):
    base = Path(args.base)
    xdata_dir, _, reg_dir, txt_dir = _default_paths(base)

    srcids, _, _ = read_sources_csv(args.sources)
    suffix = f"_p{args.ecf}"

    n_src = 0
    n_bkg = 0
    n_warn = 0

    for obsid in args.obsids:
        evtfile = xdata_dir / f"all_bcc_{obsid}_reproj_evt.fits"
        per_obs_txt = txt_dir / f"txt_{obsid}{suffix}"
        per_obs_txt.mkdir(parents=True, exist_ok=True)

        for srcid in srcids:
            regfile = reg_dir / f"region_{obsid}" / f"region_{args.ecf}" / f"{srcid}.reg"
            bkg_regfile = reg_dir / f"region_{obsid}" / f"region_{args.ecf}" / f"{srcid}_bkg.reg"
            out_src = per_obs_txt / f"{srcid}.txt"
            out_bkg = per_obs_txt / f"{srcid}_bkg.txt"

            try:
                extract_source_photons(evtfile, regfile, out_src, obsid, emin=args.emin, emax=args.emax)
                n_src += 1
            except Exception as exc:
                n_warn += 1
                out_src.write_text("")
                print(f"[WARN] source extraction failed for obsid={obsid} srcid={srcid}: {exc}")

            try:
                extract_background_photons(evtfile, bkg_regfile, out_bkg, obsid, emin=args.emin, emax=args.emax)
                n_bkg += 1
            except Exception as exc:
                n_warn += 1
                out_bkg.write_text("")
                print(f"[WARN] background extraction failed for obsid={obsid} srcid={srcid}: {exc}")

    print(f"extract summary: source files={n_src}, background files={n_bkg}, warnings={n_warn}")


def cmd_epoch(args):
    base = Path(args.base)
    xdata_dir, timing_dir, _, _ = _default_paths(base)
    epoch_out = timing_dir / "epoch_all_obs.txt"
    make_epoch_file(args.obsids, xdata_dir, epoch_out)
    print(f"epoch summary: wrote {epoch_out}")


def cmd_merge(args):
    base = Path(args.base)
    _, _, _, txt_dir = _default_paths(base)

    srcids, _, _ = read_sources_csv(args.sources)
    outdir = txt_dir / f"txt_all_obs{args.suffix}"
    outdir.mkdir(parents=True, exist_ok=True)

    n = 0
    for srcid in srcids:
        merge_txt_for_source(
            srcid=srcid,
            obsids=args.obsids,
            timing_txt_dir=txt_dir,
            outdir=outdir,
            suffix=args.suffix,
            include_bkg=not args.no_bkg,
        )
        n += 1

    print(f"merge summary: merged {n} sources into {outdir}")


def cmd_all(args):
    cmd_make_regions(args)
    cmd_extract(args)
    cmd_epoch(args)

    merge_args = argparse.Namespace(**vars(args))
    merge_args.suffix = f"_p{args.ecf}"
    merge_args.no_bkg = False
    cmd_merge(merge_args)

    base = Path(args.base)
    _, timing_dir, reg_dir, txt_dir = _default_paths(base)
    print("\nALL summary")
    print(f"  Regions root : {reg_dir}")
    print(f"  Text root    : {txt_dir}")
    print(f"  Epoch file   : {timing_dir / 'epoch_all_obs.txt'}")
    print(f"  Merged txt   : {txt_dir / f'txt_all_obs_p{args.ecf}'}")


def build_parser():
    parser = argparse.ArgumentParser(description="Extract ACIS timing photons from barycentered event files.")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_common(p, need_sources=True):
        p.add_argument("--base", default="/data/home/tiger/chandra/course", help="Pipeline base directory")
        p.add_argument("--obsids", nargs="+", type=int, required=True, help="ObsID list")
        if need_sources:
            p.add_argument("--sources", required=True, help="CSV with columns srcid,ra,dec")

    p_reg = sub.add_parser("make-regions", help="Create source/background region files")
    add_common(p_reg, need_sources=True)
    p_reg.add_argument("--ecf", type=int, default=90)
    p_reg.add_argument("--src-radius-factor", type=float, default=1.0)
    p_reg.add_argument("--bkg-inner-factor", type=float, default=2.0)
    p_reg.add_argument("--bkg-outer-factor", type=float, default=4.0)
    p_reg.add_argument("--pixel-scale-arcsec", type=float, default=0.492)
    p_reg.add_argument("--exclude-neighbors", action="store_true", help="Exclude nearby source circles from background")
    p_reg.add_argument("--exclude-radius-factor", type=float, default=1.5)
    p_reg.set_defaults(func=cmd_make_regions)

    p_ext = sub.add_parser("extract", help="Extract source/background photons into txt files")
    add_common(p_ext, need_sources=True)
    p_ext.add_argument("--ecf", type=int, default=90)
    p_ext.add_argument("--emin", type=float, default=300)
    p_ext.add_argument("--emax", type=float, default=8000)
    p_ext.set_defaults(func=cmd_extract)

    p_epoch = sub.add_parser("epoch", help="Create epoch file for all obsids")
    add_common(p_epoch, need_sources=False)
    p_epoch.set_defaults(func=cmd_epoch)

    p_merge = sub.add_parser("merge", help="Merge per-obs txt files into all-observation products")
    add_common(p_merge, need_sources=True)
    p_merge.add_argument("--suffix", default="_p90")
    p_merge.add_argument("--no-bkg", action="store_true")
    p_merge.set_defaults(func=cmd_merge)

    p_all = sub.add_parser("all", help="Run make-regions -> extract -> epoch -> merge")
    add_common(p_all, need_sources=True)
    p_all.add_argument("--ecf", type=int, default=90)
    p_all.add_argument("--emin", type=float, default=300)
    p_all.add_argument("--emax", type=float, default=8000)
    p_all.add_argument("--src-radius-factor", type=float, default=1.0)
    p_all.add_argument("--bkg-inner-factor", type=float, default=2.0)
    p_all.add_argument("--bkg-outer-factor", type=float, default=4.0)
    p_all.add_argument("--pixel-scale-arcsec", type=float, default=0.492)
    p_all.add_argument("--exclude-neighbors", action="store_true")
    p_all.add_argument("--exclude-radius-factor", type=float, default=1.5)
    p_all.set_defaults(func=cmd_all)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
