# ACIS timing extraction quickstart


## 1) Run CIAO preprocessing

```bash
bash acis_course_modular.sh all 23603
```

This generates required old-style files in:

- `/data/home/tiger/chandra/course/merge_data/xdata/all_bcc_23603_evt.fits`
- `/data/home/tiger/chandra/course/merge_data/xdata/psf90_23603_500_8000.fits`
- `/data/home/tiger/chandra/course/merge_data/xdata/expmap_23603_500_8000.fits`
- `/data/home/tiger/chandra/course/merge_data/xdata/img_23603_500_8000.fits`

## 2) Extract photons

```bash
python extract_chandra_photons.py all \
    --base /data/home/tiger/chandra/course \
    --obsids 23603 \
    --sources sources.csv \
    --ecf 90 \
    --emin 500 \
    --emax 8000
```

## 3) Expected outputs

Region files:

- `timing/reg/region_23603/region_90/1.reg`
- `timing/reg/region_23603/region_90/1_bkg.reg`

Per-ObsID photon text files:

- `timing/txt/txt_23603_p90/1.txt`
- `timing/txt/txt_23603_p90/1_bkg.txt`

Merged all-observation products:

- `timing/txt/txt_all_obs_p90/1.txt`
- `timing/txt/txt_all_obs_p90/1_bkg.txt`
- `timing/txt/txt_all_obs_p90/epoch_src_1.txt`


# Chandra period-search scripts (independent LS and GL-like workflows)

## Required inputs
For each source ID (`srcid`), place these files under:

`{BASE}/merge_data/timing/txt/txt_all_obs_p90/`

- `{srcid}.txt` (source photons: `time  energy  obsid`)
- `{srcid}_bkg.txt` (background photons: `time  energy  obsid`)
- `epoch_src_{srcid}.txt` (epochs: `tstart  tstop  obsid  exposure`)

Default suffix is `_p90`.

## Directory structure
- Inputs:
  - `${BASE}/merge_data/timing/txt/txt_all_obs_p90/{srcid}.txt`
  - `${BASE}/merge_data/timing/txt/txt_all_obs_p90/{srcid}_bkg.txt`
  - `${BASE}/merge_data/timing/txt/txt_all_obs_p90/epoch_src_{srcid}.txt`
- Outputs (both scripts):
  - `${BASE}/merge_data/timing/period_search/src_{srcid}/`

## Example commands (run separately)

### Lomb–Scargle (binned light curve)
```bash
python run_ls_search.py \
    --base /data/home/tiger/chandra/course \
    --srcid 1 \
    --pmin 100 \
    --pmax 20000 \
    --dt 50
```

### Gregory–Loredo exact odds (unbinned events)
```bash
python run_gl_search.py \
    --base /data/home/tiger/chandra/course \
    --srcid 1 \
    --pmin 100 \
    --pmax 20000 \
    --mmax 12 \
    --ni 10
```

## Notes and interpretation
- `run_ls_search.py` uses **binned** light curves (`counts/exposure_per_bin`) and runs `astropy.timeseries.LombScargle`.
- `run_gl_search.py` now uses the **exact Gregory–Loredo odds** style workflow (unbinned phases, m-bin odds integration with phase-offset integration step `--ni`).
- The two scripts are independent and should be run independently (no wrapper needed).
- Compare best periods with literature values for your source before claiming a detection.
- Always inspect possible instrumental aliases near Chandra dither periods **706.96 s** and **999.96 s** (plus harmonics/subharmonics).
- Do not claim a robust detection if only one or two cycles are covered by exposure or baseline.
- Multi-ObsID windows and large inter-observation gaps can produce strong aliases and period ambiguities.
