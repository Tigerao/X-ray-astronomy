#!/usr/bin/env bash
# Minimal CIAO ACIS imaging pipeline:
# chandra_repro -> fluximage(expmap) -> mkpsfmap(from expmap) -> axbary(for timing)
# Run after: source /path/to/ciao/bin/ciao.sh

set -euo pipefail

# ===== user settings =====
BASEDIR="/data/home/tiger/chandra/course"
WORKDIR="${BASEDIR}/merge_data/xdata"
OBSIDS="23603"
RA=23.462042
DEC=30.660222
BANDS="broad"       # use broad for minimal result; can change to "broad,soft,medium,hard"
BINSIZE=1
PSF_ENERGY=2.3      # keV
PSF_ECF=0.9         # 90% encircled counts fraction
NPROC=1             # use 1 for stability; change to -1 later if desired
mkdir -p "$WORKDIR"

log="${WORKDIR}/minimal_acis_pipeline.log"
: > "$log"

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found; initialize CIAO first" >&2; exit 1; }; }
for t in chandra_repro fluximage mkpsfmap axbary dmkeypar; do need "$t"; done

for OBSID in $OBSIDS; do
  INDIR="${BASEDIR}/${OBSID}"
  REPRO="${INDIR}/repro"
  OUT="${WORKDIR}/${OBSID}"
  mkdir -p "$OUT"

  echo "===== ${OBSID}: chandra_repro =====" | tee -a "$log"
  punlearn chandra_repro
  # Important: pix_adj value is lowercase in CIAO: default|edser|none|randomize|centroid
  # set_ardlib=no avoids leaving an obs-specific ardlib setting in your session.
  chandra_repro indir="$INDIR" outdir="$REPRO" \
      pix_adj=edser check_vf_pha=yes set_ardlib=no cleanup=yes clobber=yes verbose=1 \
      2>&1 | tee -a "$log"

  EVT=$(ls -1 "$REPRO"/*_repro_evt2.fits | head -n 1)
  ASOL=$(ls -1 "$REPRO"/pcadf*_asol1.fits 2>/dev/null | paste -sd, -)
  BPIX=$(ls -1 "$REPRO"/*_repro_bpix1.fits 2>/dev/null | head -n 1)
  MASK=$(ls -1 "$REPRO"/*_msk1.fits 2>/dev/null | head -n 1)
  EPH=$(find "$INDIR" -name 'orbit*eph1.fits*' | sort | head -n 1)

  [[ -n "$EVT"  ]] || { echo "ERROR: no repro evt2 found" >&2; exit 2; }
  [[ -n "$ASOL" ]] || { echo "ERROR: no asol found in $REPRO" >&2; exit 3; }
  [[ -n "$BPIX" ]] || { echo "ERROR: no repro bpix found in $REPRO" >&2; exit 4; }
  [[ -n "$MASK" ]] || { echo "ERROR: no msk found in $REPRO" >&2; exit 5; }
  [[ -n "$EPH"  ]] || { echo "ERROR: no eph1 found under $INDIR; download eph1" >&2; exit 6; }

  GRATING=$(dmkeypar "$EVT" GRATING echo+ 2>/dev/null || echo "UNKNOWN")
  if [[ "$GRATING" != "NONE" ]]; then
    echo "WARNING: ObsID $OBSID has GRATING=$GRATING; expected non-grating ACIS imaging." | tee -a "$log"
  fi

  echo "EVT=$EVT"   | tee -a "$log"
  echo "ASOL=$ASOL" | tee -a "$log"
  echo "BPIX=$BPIX" | tee -a "$log"
  echo "MASK=$MASK" | tee -a "$log"
  echo "EPH=$EPH"   | tee -a "$log"

  echo "===== ${OBSID}: fluximage on un-barycentered repro EVT =====" | tee -a "$log"
  # Do NOT run fluximage on the axbary output copied to WORKDIR: ASOLFILE keyword may be relative/basename only.
  # Instead use the repro evt2 and explicitly pass asolfile/badpixfile/maskfile.
  fluximage infile="$EVT" outroot="${OUT}/flux" \
      bands="$BANDS" binsize="$BINSIZE" \
      asolfile="$ASOL" badpixfile="$BPIX" maskfile="$MASK" \
      units=default background=none psfecf=none \
      parallel=no nproc="$NPROC" clobber=yes verbose=1 \
      2>&1 | tee -a "$log"

  # Locate broad expmap. CIAO names are normally flux_broad_thresh.expmap or flux_broad.expmap.
  EXPMAP=$(ls -1 "${OUT}"/flux_broad*expmap 2>/dev/null | head -n 1)
  [[ -n "$EXPMAP" ]] || { echo "ERROR: no broad expmap produced in $OUT" >&2; ls -l "$OUT"; exit 7; }

  echo "===== ${OBSID}: mkpsfmap from expmap =====" | tee -a "$log"
  PSFMAP="${OUT}/psf${PSF_ECF}_${OBSID}_from_expmap.fits"
  punlearn mkpsfmap
  mkpsfmap infile="$EXPMAP" outfile="$PSFMAP" energy="$PSF_ENERGY" ecf="$PSF_ECF" \
      units=arcsec clobber=yes verbose=1 \
      2>&1 | tee -a "$log"

  echo "===== ${OBSID}: axbary for timing only =====" | tee -a "$log"
  BARY="${OUT}/evt2_${OBSID}_bcc.fits"
  punlearn axbary
  axbary infile="$EVT" outfile="$BARY" orbitfile="$EPH" ra="$RA" dec="$DEC" \
      refframe=ICRS clobber=yes \
      2>&1 | tee -a "$log"

  echo "DONE ${OBSID}" | tee -a "$log"
  echo "  expmap: $EXPMAP" | tee -a "$log"
  echo "  psfmap: $PSFMAP" | tee -a "$log"
  echo "  bcc evt: $BARY" | tee -a "$log"
done

echo "All done: $(date)" | tee -a "$log"
