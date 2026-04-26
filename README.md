# X-ray-astronomy

## ACIS timing extraction quickstart

### 1) Run CIAO preprocessing

```bash
bash acis_course_modular.sh all 23603
```

This generates required files in:

- `/data/home/tiger/chandra/course/merge_data/xdata/all_bcc_23603_reproj_evt.fits`
- `/data/home/tiger/chandra/course/merge_data/xdata/reproj_psf90_23603_500_8000.fits`
- `/data/home/tiger/chandra/course/merge_data/xdata/reproj_expmap_23603_500_8000.fits`
- `/data/home/tiger/chandra/course/merge_data/xdata/img_23603_500_8000.fits`

### 2) Extract photons

```bash
python extract_chandra_photons.py all \
    --base /data/home/tiger/chandra/course \
    --obsids 23603 \
    --sources sources.csv \
    --ecf 90 \
    --emin 500 \
    --emax 8000
```

### 3) Expected outputs

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
