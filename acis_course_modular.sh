#!/usr/bin/env bash
# Modular ACIS preprocessing pipeline (direct filename outputs; no symlinks).

set -euo pipefail

BASE="/data/home/tiger/chandra/course"
XDATA_DIR="${BASE}/merge_data/xdata"
LOGFILE="${XDATA_DIR}/acis_course_modular.log"

# barycenter coordinates (deg)
: "${BARY_RA:=85.83383}"
: "${BARY_DEC:=-41.03175}"

PIX_ADJ="edser"
PSF_ENERGY_KEV="2.3"
PSF_ECF="0.9"
BAND_MIN=500
BAND_MAX=8000

err(){ echo "ERROR: $*" >&2; exit 1; }
msg(){ echo "[$(date -u +"%Y-%m-%dT%H:%M:%SZ")] $*" | tee -a "$LOGFILE"; }

need_tool(){ command -v "$1" >/dev/null 2>&1 || err "Missing CIAO tool: $1"; }
check_tools(){ for t in chandra_repro fluximage mkpsfmap axbary dmcopy punlearn pset; do need_tool "$t"; done; }

setup_paths(){
  OBSID="$1"
  OBS_DIR="${BASE}/${OBSID}"
  REPRO_DIR="${OBS_DIR}/repro"
  OBS_WORK_DIR="${XDATA_DIR}/${OBSID}"
  mkdir -p "$REPRO_DIR" "$OBS_WORK_DIR" "$XDATA_DIR"
  : >> "$LOGFILE"
  [[ -d "${OBS_DIR}/primary" ]] || err "Missing ${OBS_DIR}/primary"
  [[ -d "${OBS_DIR}/secondary" ]] || err "Missing ${OBS_DIR}/secondary"
}

find_repro_inputs(){
  EVT_REPRO=$(find "$REPRO_DIR" -maxdepth 1 -name '*_repro_evt2.fits' | head -n 1 || true)
  [[ -n "$EVT_REPRO" ]] || err "No *_repro_evt2.fits in $REPRO_DIR"
  ASOL=$(find "$REPRO_DIR" -maxdepth 1 -name 'pcadf*_asol1.fits' | sort | paste -sd, - || true)
  [[ -n "$ASOL" ]] || err "No pcadf*_asol1.fits in $REPRO_DIR"
  BADPIX=$(find "$REPRO_DIR" -maxdepth 1 -name '*_repro_bpix1.fits' | head -n 1 || true)
  [[ -n "$BADPIX" ]] || err "No *_repro_bpix1.fits in $REPRO_DIR"
  MASK=$(find "$REPRO_DIR" -maxdepth 1 -name '*_msk1.fits' | head -n 1 || true)
  [[ -n "$MASK" ]] || err "No *_msk1.fits in $REPRO_DIR"
  ORBIT=$(find "$OBS_DIR" -type f -name 'orbit*eph1.fits*' | head -n 1 || true)
}

do_repro(){
  msg "ObsID $OBSID: chandra_repro pix_adj=${PIX_ADJ}"
  punlearn chandra_repro
  chandra_repro indir="$OBS_DIR" outdir="$REPRO_DIR" pix_adj="$PIX_ADJ" check_vf_pha=yes cleanup=yes clobber=yes
}

do_expmap(){
  find_repro_inputs
  msg "ObsID $OBSID: fluximage on un-barycentered repro evt2 (500-8000 eV band)"

  local outroot="${OBS_WORK_DIR}/flux_${OBSID}"
  punlearn fluximage
  fluximage infile="$EVT_REPRO" outroot="$outroot" \
    bands="0.5:8.0:2.3" binsize=1 \
    asolfile="$ASOL" badpixfile="$BADPIX" maskfile="$MASK" \
    units=default psfecf="$PSF_ECF" background=none clobber=yes

  local expmap_tmp
  expmap_tmp=$(find "$OBS_WORK_DIR" -maxdepth 1 -name "flux_${OBSID}_*.expmap" | head -n 1 || true)
  [[ -n "$expmap_tmp" ]] || err "fluximage did not produce expmap in $OBS_WORK_DIR"

  IMG_OUT="${XDATA_DIR}/img_${OBSID}_500_8000.fits"
  EXPMAP_OUT="${XDATA_DIR}/expmap_${OBSID}_500_8000.fits"

  dmcopy "${EVT_REPRO}[energy=${BAND_MIN}:${BAND_MAX}][bin x=::1,y=::1]" "$IMG_OUT" clobber=yes
  dmcopy "$expmap_tmp" "$EXPMAP_OUT" clobber=yes

  msg "image: $IMG_OUT"
  msg "expmap: $EXPMAP_OUT"
}

do_psfmap(){
  EXPMAP_OUT="${XDATA_DIR}/expmap_${OBSID}_500_8000.fits"
  [[ -f "$EXPMAP_OUT" ]] || err "Missing $EXPMAP_OUT. Run expmap first."
  PSFMAP_OUT="${XDATA_DIR}/psf90_${OBSID}_500_8000.fits"

  msg "ObsID $OBSID: mkpsfmap (non-interactive), infile=$EXPMAP_OUT"
  punlearn mkpsfmap
  pset mkpsfmap infile="$EXPMAP_OUT"
  pset mkpsfmap outfile="$PSFMAP_OUT"
  pset mkpsfmap energy="$PSF_ENERGY_KEV"
  pset mkpsfmap ecf="$PSF_ECF"
  pset mkpsfmap units=arcsec
  pset mkpsfmap clobber=yes
  pset mkpsfmap mode=h
  mkpsfmap

  [[ -f "$PSFMAP_OUT" ]] || err "mkpsfmap failed: $PSFMAP_OUT"
  msg "psfmap: $PSFMAP_OUT"
}

do_bary(){
  find_repro_inputs
  [[ -n "$ORBIT" ]] || err "No orbit*eph1.fits found under $OBS_DIR"
  BARY_OUT="${XDATA_DIR}/all_bcc_${OBSID}_evt.fits"

  msg "ObsID $OBSID: axbary for timing (output uses required final filename)"
  axbary infile="$EVT_REPRO" outfile="$BARY_OUT" orbitfile="$ORBIT" \
    ra="$BARY_RA" dec="$BARY_DEC" refframe=ICRS clobber=yes

  [[ -f "$BARY_OUT" ]] || err "axbary failed: $BARY_OUT"
  msg "bary evt: $BARY_OUT"
}

status(){
  echo "repro evt2: $(find "$REPRO_DIR" -maxdepth 1 -name '*_repro_evt2.fits' | head -n 1 || true)"
  echo "xdata final products:"
  echo "  ${XDATA_DIR}/all_bcc_${OBSID}_evt.fits"
  echo "  ${XDATA_DIR}/psf90_${OBSID}_500_8000.fits"
  echo "  ${XDATA_DIR}/expmap_${OBSID}_500_8000.fits"
  echo "  ${XDATA_DIR}/img_${OBSID}_500_8000.fits"
}

usage(){
  cat <<USAGE
Usage: bash $0 <task> <obsid>

Tasks:
  all     : repro + expmap + psfmap + bary
  repro   : run chandra_repro only
  expmap  : run fluximage (un-bary evt2) + create img/expmap final files (500-8000 eV)
  psfmap  : run mkpsfmap using final expmap filename
  bary    : run axbary to final bary event filename
  status  : show expected final file paths
USAGE
}

main(){
  local task="${1:-}" obsid="${2:-}"
  [[ -n "$task" && -n "$obsid" ]] || { usage; exit 1; }
  check_tools
  setup_paths "$obsid"
  msg "START task=$task obsid=$obsid"

  case "$task" in
    all) do_repro; do_expmap; do_psfmap; do_bary ;;
    repro) do_repro ;;
    expmap) do_expmap ;;
    psfmap) do_psfmap ;;
    bary) do_bary ;;
    status) status ;;
    *) usage; exit 1 ;;
  esac

  msg "DONE task=$task obsid=$obsid"
}

main "$@"
