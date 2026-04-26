#!/usr/bin/env bash
# Modular ACIS preprocessing pipeline.
#
# Produces (under BASE/merge_data/xdata):
#   all_bcc_{obsid}_reproj_evt.fits
#   reproj_psf90_{obsid}_500_8000.fits
#   reproj_expmap_{obsid}_500_8000.fits
#   img_{obsid}_500_8000.fits
#
# Usage:
#   bash acis_course_modular.sh all 23603
#   bash acis_course_modular.sh repro 23603
#   bash acis_course_modular.sh expmap 23603
#   bash acis_course_modular.sh psfmap 23603
#   bash acis_course_modular.sh bary 23603
#   bash acis_course_modular.sh links 23603
#   bash acis_course_modular.sh status 23603

set -euo pipefail

BASE="/data/home/tiger/chandra/course"
XDATA_DIR="${BASE}/merge_data/xdata"
LOGFILE="${XDATA_DIR}/acis_course_modular.log"

# Default source position used by axbary; override with environment if desired.
BARY_RA="${BARY_RA:-23.462042}"
BARY_DEC="${BARY_DEC:-30.660222}"

PSF_ENERGY_KEV="2.3"
PSF_ECF="0.9"
PIX_ADJ="edser"
BINSIZE="1"

err() {
  echo "ERROR: $*" >&2
  exit 1
}

msg() {
  echo "[$(date -u +"%Y-%m-%dT%H:%M:%SZ")] $*" | tee -a "$LOGFILE"
}

need_tool() {
  command -v "$1" >/dev/null 2>&1 || err "Required CIAO tool '$1' not found in PATH. Did you source ciao.sh?"
}

check_tools() {
  local tools=(chandra_repro fluximage mkpsfmap axbary dmcopy dmcoords pset punlearn)
  for t in "${tools[@]}"; do
    need_tool "$t"
  done
}

setup_paths() {
  OBSID="$1"
  OBS_DIR="${BASE}/${OBSID}"
  PRIMARY_DIR="${OBS_DIR}/primary"
  SECONDARY_DIR="${OBS_DIR}/secondary"
  REPRO_DIR="${OBS_DIR}/repro"
  OUT_DIR="${XDATA_DIR}/${OBSID}"

  [[ -d "$OBS_DIR" ]] || err "ObsID directory not found: $OBS_DIR"
  [[ -d "$PRIMARY_DIR" ]] || err "Missing primary directory: $PRIMARY_DIR"
  [[ -d "$SECONDARY_DIR" ]] || err "Missing secondary directory: $SECONDARY_DIR"

  mkdir -p "$REPRO_DIR" "$OUT_DIR" "$XDATA_DIR"
  : >> "$LOGFILE"
}

find_repro_inputs() {
  EVT_REPRO=$(find "$REPRO_DIR" -maxdepth 1 -name '*_repro_evt2.fits' | head -n 1 || true)
  [[ -n "$EVT_REPRO" ]] || err "No *_repro_evt2.fits found in $REPRO_DIR. Run repro first."

  ASOL=$(find "$REPRO_DIR" -maxdepth 1 -name 'pcadf*_asol1.fits' | sort | paste -sd, - || true)
  [[ -n "$ASOL" ]] || err "No pcadf*_asol1.fits found in $REPRO_DIR"

  BADPIX=$(find "$REPRO_DIR" -maxdepth 1 -name '*_repro_bpix1.fits' | head -n 1 || true)
  [[ -n "$BADPIX" ]] || err "No *_repro_bpix1.fits found in $REPRO_DIR"

  MASK=$(find "$REPRO_DIR" -maxdepth 1 -name '*_msk1.fits' | head -n 1 || true)
  [[ -n "$MASK" ]] || err "No *_msk1.fits found in $REPRO_DIR"

  ORBIT=$(find "$OBS_DIR" -type f -name 'orbit*eph1.fits*' | head -n 1 || true)
}

do_repro() {
  msg "ObsID $OBSID: running chandra_repro (pix_adj=${PIX_ADJ})"
  punlearn chandra_repro
  chandra_repro \
    indir="$OBS_DIR" \
    outdir="$REPRO_DIR" \
    pix_adj="$PIX_ADJ" \
    check_vf_pha=yes \
    cleanup=yes \
    clobber=yes

  find_repro_inputs
  msg "repro evt2: $EVT_REPRO"
}

do_expmap() {
  find_repro_inputs
  msg "ObsID $OBSID: running fluximage on UN-BARYCENTERED repro evt2"

  local outroot="${OUT_DIR}/flux"
  punlearn fluximage
  fluximage \
    infile="$EVT_REPRO" \
    outroot="$outroot" \
    bands="broad" \
    binsize="$BINSIZE" \
    asolfile="$ASOL" \
    badpixfile="$BADPIX" \
    maskfile="$MASK" \
    units=default \
    psfecf="$PSF_ECF" \
    background=none \
    clobber=yes

  EXPMAP_BROAD=$(find "$OUT_DIR" -maxdepth 1 -name 'flux_broad*.expmap' | head -n 1 || true)
  IMG_BROAD=$(find "$OUT_DIR" -maxdepth 1 -name 'flux_broad*.img' | head -n 1 || true)

  [[ -n "$EXPMAP_BROAD" ]] || err "fluximage did not produce broad-band exposure map in $OUT_DIR"

  IMG_500_8000="${OUT_DIR}/img_${OBSID}_500_8000.fits"
  msg "ObsID $OBSID: creating 0.3-8.0 keV image: $IMG_500_8000"
  dmcopy "${EVT_REPRO}[energy=500:8000][bin x=::1,y=::1]" "$IMG_500_8000" clobber=yes

  msg "broad expmap: $EXPMAP_BROAD"
  [[ -n "$IMG_BROAD" ]] && msg "broad img from fluximage: $IMG_BROAD"
}

do_psfmap() {
  local expmap_in
  expmap_in=$(find "$OUT_DIR" -maxdepth 1 -name 'flux_broad*.expmap' | head -n 1 || true)
  [[ -n "$expmap_in" ]] || err "No broad expmap found in $OUT_DIR. Run expmap first."

  PSFMAP_OUT="${OUT_DIR}/psfmap_ecf90_${OBSID}.fits"

  msg "ObsID $OBSID: running mkpsfmap (non-interactive), infile=$expmap_in"
  punlearn mkpsfmap
  pset mkpsfmap infile="$expmap_in"
  pset mkpsfmap outfile="$PSFMAP_OUT"
  pset mkpsfmap energy="$PSF_ENERGY_KEV"
  pset mkpsfmap ecf="$PSF_ECF"
  pset mkpsfmap units="arcsec"
  pset mkpsfmap clobber=yes
  pset mkpsfmap mode=h
  mkpsfmap

  [[ -f "$PSFMAP_OUT" ]] || err "mkpsfmap failed to create $PSFMAP_OUT"
  msg "psfmap: $PSFMAP_OUT"
}

do_bary() {
  find_repro_inputs
  [[ -n "$ORBIT" ]] || err "No orbit*eph1.fits found below $OBS_DIR; required for axbary"

  BARY_EVT="${OUT_DIR}/evt2_${OBSID}_bcc.fits"
  msg "ObsID $OBSID: running axbary for timing products only"
  axbary \
    infile="$EVT_REPRO" \
    outfile="$BARY_EVT" \
    orbitfile="$ORBIT" \
    ra="$BARY_RA" \
    dec="$BARY_DEC" \
    refframe=ICRS \
    clobber=yes

  [[ -f "$BARY_EVT" ]] || err "axbary failed to produce $BARY_EVT"
  msg "bary evt2: $BARY_EVT"
}

make_oldstyle_links() {
  local bary_evt="${OUT_DIR}/evt2_${OBSID}_bcc.fits"
  local psfmap="${OUT_DIR}/psfmap_ecf90_${OBSID}.fits"
  local expmap
  expmap=$(find "$OUT_DIR" -maxdepth 1 -name 'flux_broad*.expmap' | head -n 1 || true)
  local img500="${OUT_DIR}/img_${OBSID}_500_8000.fits"

  [[ -f "$bary_evt" ]] || err "Missing barycentered evt2: $bary_evt (run bary first)"
  [[ -f "$psfmap" ]] || err "Missing PSF map: $psfmap (run psfmap first)"
  [[ -n "$expmap" && -f "$expmap" ]] || err "Missing broad exposure map in $OUT_DIR (run expmap first)"
  [[ -f "$img500" ]] || err "Missing 0.3-8.0 keV image: $img500 (run expmap first)"

  local old_evt="${XDATA_DIR}/all_bcc_${OBSID}_reproj_evt.fits"
  local old_psf="${XDATA_DIR}/reproj_psf90_${OBSID}_500_8000.fits"
  local old_exp="${XDATA_DIR}/reproj_expmap_${OBSID}_500_8000.fits"
  local old_img="${XDATA_DIR}/img_${OBSID}_500_8000.fits"

  ln -sfn "$bary_evt" "$old_evt"
  ln -sfn "$psfmap" "$old_psf"
  ln -sfn "$expmap" "$old_exp"
  ln -sfn "$img500" "$old_img"

  msg "Created/updated old-style links:"
  msg "  $old_evt -> $bary_evt"
  msg "  $old_psf -> $psfmap"
  msg "  $old_exp -> $expmap"
  msg "  $old_img -> $img500"
}

status() {
  msg "Status for ObsID $OBSID"
  echo "OBS_DIR=$OBS_DIR"
  echo "REPRO_DIR=$REPRO_DIR"
  echo "OUT_DIR=$OUT_DIR"
  echo "XDATA_DIR=$XDATA_DIR"
  echo
  echo "Repro files:"
  find "$REPRO_DIR" -maxdepth 1 -type f | sort || true
  echo
  echo "Obs output files:"
  find "$OUT_DIR" -maxdepth 1 -type f | sort || true
  echo
  echo "Old-style products in xdata/:"
  find "$XDATA_DIR" -maxdepth 1 -type l -o -maxdepth 1 -type f | rg "${OBSID}|all_bcc_|reproj_psf90_|reproj_expmap_|img_" || true
}

usage() {
  cat <<USAGE
Usage: bash $0 <task> <obsid>

Tasks:
  all     : repro + expmap + psfmap + bary + links
  repro   : run chandra_repro
  expmap  : run fluximage (broad 0.5-7.0 keV) and create img_{obsid}_500_8000.fits
  psfmap  : run mkpsfmap using fluximage broad exposure map
  bary    : run axbary on repro evt2 for timing products
  links   : create/update required old-style filenames in xdata/
  status  : print file status for this obsid
USAGE
}

main() {
  local task="${1:-}"
  local obsid="${2:-}"
  [[ -n "$task" && -n "$obsid" ]] || { usage; exit 1; }

  check_tools
  setup_paths "$obsid"

  msg "START task=$task obsid=$obsid"
  case "$task" in
    all)
      do_repro
      do_expmap
      do_psfmap
      do_bary
      make_oldstyle_links
      ;;
    repro)
      do_repro
      ;;
    expmap)
      do_expmap
      ;;
    psfmap)
      do_psfmap
      ;;
    bary)
      do_bary
      ;;
    links)
      make_oldstyle_links
      ;;
    status)
      status
      ;;
    *)
      usage
      exit 1
      ;;
  esac

  msg "Final product paths for ObsID $obsid:"
  echo "  ${XDATA_DIR}/all_bcc_${obsid}_reproj_evt.fits"
  echo "  ${XDATA_DIR}/reproj_psf90_${obsid}_500_8000.fits"
  echo "  ${XDATA_DIR}/reproj_expmap_${obsid}_500_8000.fits"
  echo "  ${XDATA_DIR}/img_${obsid}_500_8000.fits"
  msg "DONE task=$task obsid=$obsid"
}

main "$@"
