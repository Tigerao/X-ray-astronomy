#!/usr/bin/env bash
# ACIS imaging course pipeline, modular version.
# CIAO tested style: chandra_repro -> fluximage(expmap/psfmap) -> optional mkpsfmap-from-expmap -> axbary.
# Run after initializing CIAO, e.g. source /path/to/ciao/bin/ciao.sh
#
# Usage examples:
#   bash acis_course_modular.sh all 23603
#   bash acis_course_modular.sh repro 23603
#   bash acis_course_modular.sh expmap 23603
#   bash acis_course_modular.sh psfmap 23603
#   bash acis_course_modular.sh bary 23603
#   bash acis_course_modular.sh status 23603

set -euo pipefail

# =====================
# USER SETTINGS
# =====================
BASEDIR="/data/home/tiger/chandra/course"      # contains directories like $BASEDIR/23603/primary
WORKDIR="${BASEDIR}/merge_data/xdata"         # outputs go to $WORKDIR/$OBSID
RA=23.462042                                  # deg, target/source RA for axbary
DEC=30.660222                                 # deg, target/source Dec for axbary
BANDS="broad"                                 # broad or broad,soft,medium,hard
BINSIZE=1
PSF_ENERGY=2.3                                # keV
PSF_ECF=0.9                                   # 90% encircled counts fraction
NPROC=1                                       # keep 1 for stability
PIX_ADJ="edser"                               # legal CIAO value: default|edser|none|randomize|centroid

mkdir -p "$WORKDIR"
LOG="${WORKDIR}/acis_course_modular.log"

# =====================
# HELPERS
# =====================
need() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: $1 not found. Initialize CIAO first." >&2
    exit 1
  }
}

check_tools() {
  for t in chandra_repro fluximage mkpsfmap axbary dmkeypar dmlist; do
    need "$t"
  done
}

msg() {
  echo "===== $* =====" | tee -a "$LOG"
}

obs_paths() {
  OBSID="$1"
  INDIR="${BASEDIR}/${OBSID}"
  REPRO="${INDIR}/repro"
  OUT="${WORKDIR}/${OBSID}"
  mkdir -p "$OUT"
}

find_repro_products() {
  EVT=$(ls -1 "${REPRO}"/*_repro_evt2.fits 2>/dev/null | head -n 1 || true)
  ASOL=$(ls -1 "${REPRO}"/pcadf*_asol1.fits 2>/dev/null | paste -sd, - || true)
  BPIX=$(ls -1 "${REPRO}"/*_repro_bpix1.fits 2>/dev/null | head -n 1 || true)
  MASK=$(ls -1 "${REPRO}"/*_msk1.fits 2>/dev/null | head -n 1 || true)
  EPH=$(find "${INDIR}" -name 'orbit*eph1.fits*' 2>/dev/null | sort | head -n 1 || true)

  [[ -n "${EVT:-}"  ]] || { echo "ERROR: no repro evt2 found in ${REPRO}. Run: bash $0 repro ${OBSID}" >&2; exit 2; }
  [[ -n "${ASOL:-}" ]] || { echo "ERROR: no asol found in ${REPRO}." >&2; exit 3; }
  [[ -n "${BPIX:-}" ]] || { echo "ERROR: no repro badpix found in ${REPRO}." >&2; exit 4; }
  [[ -n "${MASK:-}" ]] || { echo "ERROR: no mask file found in ${REPRO}." >&2; exit 5; }
}

find_expmap() {
  EXPMAP=$(ls -1 "${OUT}"/flux_broad_thresh.expmap "${OUT}"/flux_broad.expmap "${OUT}"/flux_broad*expmap 2>/dev/null | head -n 1 || true)
  [[ -n "${EXPMAP:-}" ]] || { echo "ERROR: no broad expmap found in ${OUT}. Run: bash $0 expmap ${OBSID}" >&2; exit 7; }
}

# =====================
# TASKS
# =====================
do_repro() {
  obs_paths "$1"
  [[ -d "$INDIR" ]] || { echo "ERROR: missing ObsID directory: $INDIR" >&2; exit 10; }

  msg "${OBSID}: chandra_repro"
  punlearn chandra_repro
  chandra_repro indir="$INDIR" outdir="$REPRO" \
    pix_adj="$PIX_ADJ" check_vf_pha=yes set_ardlib=no cleanup=yes clobber=yes verbose=1 \
    2>&1 | tee -a "$LOG"

  find_repro_products
  GRATING=$(dmkeypar "$EVT" GRATING echo+ 2>/dev/null || echo "UNKNOWN")
  if [[ "$GRATING" != "NONE" ]]; then
    echo "WARNING: ObsID $OBSID has GRATING=$GRATING; expected ACIS imaging/non-grating." | tee -a "$LOG"
  fi
  echo "EVT=$EVT" | tee -a "$LOG"
  echo "ASOL=$ASOL" | tee -a "$LOG"
  echo "BPIX=$BPIX" | tee -a "$LOG"
  echo "MASK=$MASK" | tee -a "$LOG"
}

do_expmap() {
  obs_paths "$1"
  find_repro_products

  msg "${OBSID}: fluximage from un-barycentered repro_evt2"
  # Keep psfecf=$PSF_ECF so fluximage itself creates flux_broad_thresh.psfmap.
  # This is usually enough; the separate psfmap task below is optional and uses the expmap as input.
  punlearn fluximage
  fluximage infile="$EVT" outroot="${OUT}/flux" \
    bands="$BANDS" binsize="$BINSIZE" \
    asolfile="$ASOL" badpixfile="$BPIX" maskfile="$MASK" \
    units=default background=none psfecf="$PSF_ECF" \
    parallel=no nproc="$NPROC" clobber=yes verbose=1 \
    2>&1 | tee -a "$LOG"

  find_expmap
  echo "EXPMAP=$EXPMAP" | tee -a "$LOG"
  if [[ -f "${OUT}/flux_broad_thresh.psfmap" ]]; then
    echo "FLUXIMAGE_PSFMAP=${OUT}/flux_broad_thresh.psfmap" | tee -a "$LOG"
  fi
}

do_psfmap() {
  obs_paths "$1"
  find_expmap

  msg "${OBSID}: mkpsfmap from expmap"
  PSFMAP="${OUT}/psf_ecf${PSF_ECF}_e${PSF_ENERGY}keV_from_expmap.fits"

  # IMPORTANT for CIAO 4.15/4.18:
  # - mkpsfmap has no 'verbose' parameter, so do not pass verbose.
  # - set every required/query parameter explicitly and force mode=h to avoid prompts for ecf.
  # - spectrum is only used if energy=INDEF; with numeric energy it can be blank.
  punlearn mkpsfmap
  pset mkpsfmap infile="$EXPMAP"
  pset mkpsfmap outfile="$PSFMAP"
  pset mkpsfmap energy="$PSF_ENERGY"
  pset mkpsfmap spectrum=""
  pset mkpsfmap ecf="$PSF_ECF"
  pset mkpsfmap units="arcsec"
  pset mkpsfmap clobber=yes
  pset mkpsfmap mode=h
  mkpsfmap 2>&1 | tee -a "$LOG"

  echo "PSFMAP=$PSFMAP" | tee -a "$LOG"
}

do_bary() {
  obs_paths "$1"
  find_repro_products
  [[ -n "${EPH:-}" ]] || { echo "ERROR: no eph1 found under ${INDIR}; download eph1." >&2; exit 6; }

  msg "${OBSID}: axbary for timing only"
  BARY="${OUT}/evt2_${OBSID}_bcc.fits"
  punlearn axbary
  axbary infile="$EVT" outfile="$BARY" orbitfile="$EPH" \
    ra="$RA" dec="$DEC" refframe=ICRS clobber=yes \
    2>&1 | tee -a "$LOG"
  echo "BARY_EVT=$BARY" | tee -a "$LOG"
}

do_status() {
  obs_paths "$1"
  echo "INDIR=$INDIR"
  echo "REPRO=$REPRO"
  echo "OUT=$OUT"
  echo "--- repro products ---"
  ls -1 "$REPRO"/*_repro_evt2.fits "$REPRO"/pcadf*_asol1.fits "$REPRO"/*_repro_bpix1.fits "$REPRO"/*_msk1.fits 2>/dev/null || true
  echo "--- output products ---"
  ls -1 "$OUT"/* 2>/dev/null || true
}

usage() {
  cat <<USAGE
Usage:
  bash $0 <task> <obsid>

Tasks:
  all      = repro + expmap + psfmap + bary
  repro    = run chandra_repro only
  expmap   = run fluximage only, using repro_evt2 and explicit ASOL/BPIX/MASK
  psfmap   = run mkpsfmap only, using flux_broad_thresh.expmap
  bary     = run axbary only, using repro_evt2 and eph1
  status   = print detected files

Edit USER SETTINGS at the top before running.
USAGE
}

# =====================
# MAIN
# =====================
check_tools
TASK="${1:-}"
OBSID_ARG="${2:-}"
[[ -n "$TASK" && -n "$OBSID_ARG" ]] || { usage; exit 1; }
: >> "$LOG"
msg "START task=${TASK} obsid=${OBSID_ARG} date=$(date)"

case "$TASK" in
  all)
    do_repro "$OBSID_ARG"
    do_expmap "$OBSID_ARG"
    do_psfmap "$OBSID_ARG"
    do_bary "$OBSID_ARG"
    ;;
  repro)  do_repro "$OBSID_ARG" ;;
  expmap) do_expmap "$OBSID_ARG" ;;
  psfmap) do_psfmap "$OBSID_ARG" ;;
  bary)   do_bary "$OBSID_ARG" ;;
  status) do_status "$OBSID_ARG" ;;
  *) usage; exit 1 ;;
esac

msg "DONE task=${TASK} obsid=${OBSID_ARG} date=$(date)"
