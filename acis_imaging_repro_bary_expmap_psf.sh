#!/usr/bin/env bash
# Modern CIAO ACIS imaging pipeline replacing most of mulred410.e + bcc.e
# CIAO >= 4.18 recommended. Run after initializing CIAO:  source /path/to/ciao/bin/ciao.sh
# Directory layout expected: $BASEDIR/$OBSID/primary and $BASEDIR/$OBSID/secondary

set -euo pipefail

# =========================
# USER SETTINGS
# =========================
BASEDIR="/data/home/tiger/chandra/course"             # parent directory containing ObsID directories
WORKDIR="${BASEDIR}/merge_data/xdata"  # output directory
RA=10.685                              # barycenter/reference coordinate, deg
DEC=41.268972                          # barycenter/reference coordinate, deg

# ACIS non-grating only; edit this list
OBSIDS="23603"

# Analysis choices
CCD_FILTER="ccd_id=0,1,2,3,6,7"
BANDS="broad,soft,medium,hard"         # CSC bands: broad=0.5-7, soft=0.5-1.2, medium=1.2-2, hard=2-7 keV
BINSIZE=1                              # ACIS pixels; use 1 for source detection, larger for quick-look
DO_DEFLARE=0                           # 1 = make background light curve and sigma-clip GTI
DEFLARE_BAND="9500:12000"              # PI/energy filter in eV for high-energy background flare screening
DEFLARE_BIN=324.014                    # seconds, close to old script
DO_WAVDETECT=0
DO_MERGE=0
DO_EXPMAP=1
DO_PSF=1
DO_BARY=1
NPROC=-1                               # -1 means all but one core for CIAO scripts that support it

mkdir -p "${WORKDIR}"
LOG="${WORKDIR}/acis_pipeline.log"
: > "${LOG}"

echo "Pipeline started: $(date)" | tee -a "${LOG}"
echo "CIAO version:" | tee -a "${LOG}"
ciaover 2>&1 | tee -a "${LOG}" || true

require_tool() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found. Initialize CIAO first." >&2; exit 1; }
}
for t in chandra_repro dmcopy dmextract deflare axbary fluximage merge_obs mkpsfmap wavdetect dmkeypar dmlist; do
  require_tool "$t"
done

# Build comma-separated stacks for merge_obs/reproject_obs later
REPRO_STACK=""
BARY_STACK=""

# =========================
# 1. chandra_repro: replaces manual acis_clear_status_bits + badpix + acis_process_events + grade/status/GTI filtering
# =========================
for OBSID in ${OBSIDS}; do
  INDIR="${BASEDIR}/${OBSID}"
  REPRO="${INDIR}/repro"
  echo "===== ObsID ${OBSID}: chandra_repro =====" | tee -a "${LOG}"

  if [[ ! -d "${INDIR}" ]]; then
    echo "ERROR: missing directory ${INDIR}" | tee -a "${LOG}"
    exit 2
  fi

  punlearn chandra_repro
  # pix_adj=default currently uses CIAO-recommended processing; use pix_adj=EDSER explicitly for ACIS imaging.
  # check_vf_pha=yes cleans VFAINT background where applicable.
  chandra_repro indir="${INDIR}" outdir="${REPRO}" \
      pix_adj=EDSER check_vf_pha=yes cleanup=yes clobber=yes verbose=1 \
      2>&1 | tee -a "${LOG}"

  EVT=$(ls -1 "${REPRO}"/*_repro_evt2.fits | head -n 1)
  echo "Reprocessed EVT2: ${EVT}" | tee -a "${LOG}"

  GRATING=$(dmkeypar "${EVT}" GRATING echo+ || echo "UNKNOWN")
  if [[ "${GRATING}" != "NONE" ]]; then
    echo "WARNING: ObsID ${OBSID} has GRATING=${GRATING}; this pipeline is intended for ACIS imaging/non-grating." | tee -a "${LOG}"
  fi

  # Keep only useful ACIS event columns and selected chips, analogous to old evt2file_pps_xyet.fits.
  EVT_MIN="${WORKDIR}/evt2_${OBSID}_repro_xyet.fits"
  dmcopy "${EVT}[${CCD_FILTER}][cols time,x,y,chipx,chipy,tdetx,tdety,detx,dety,ccd_id,energy,pi,grade,status]" \
      "${EVT_MIN}" clobber=yes

  # =========================
  # 2. Optional flare filtering, replacing old IDL flarefilter block
  # =========================
  EVT_CLEAN="${EVT_MIN}"
  if [[ "${DO_DEFLARE}" == "1" ]]; then
    echo "===== ObsID ${OBSID}: deflare =====" | tee -a "${LOG}"
    LC="${WORKDIR}/bkg_lc_${OBSID}.fits"
    GTI="${WORKDIR}/deflare_${OBSID}.gti"
    EVT_CLEAN="${WORKDIR}/evt2_${OBSID}_clean.fits"

    dmextract "${EVT_MIN}[energy=${DEFLARE_BAND}][bin time=::${DEFLARE_BIN}]" "${LC}" opt=ltc1 clobber=yes
    deflare "${LC}" "${GTI}" method=sigma nsigma=3 plot=no save="${WORKDIR}/deflare_${OBSID}.png" verbose=1
    dmcopy "${EVT_MIN}[@${GTI}]" "${EVT_CLEAN}" clobber=yes
  fi

  # =========================
  # 3. Barycentric correction: replaces bcc.e axbary block
  # =========================
  if [[ "${DO_BARY}" == "1" ]]; then
    echo "===== ObsID ${OBSID}: axbary =====" | tee -a "${LOG}"
    EPH=$(find "${INDIR}" -name 'orbit*eph1.fits*' | sort | head -n 1)
    if [[ -z "${EPH}" ]]; then
      echo "ERROR: no orbit ephemeris file found for ObsID ${OBSID}. Download eph1." | tee -a "${LOG}"
      exit 3
    fi
    # axbary accepts .gz in modern CIAO/CFITSIO in most cases, but use funpack/gunzip if your system fails.
    EVT_BARY="${WORKDIR}/evt2_${OBSID}_bcc.fits"
    punlearn axbary
    axbary infile="${EVT_CLEAN}" orbitfile="${EPH}" outfile="${EVT_BARY}" ra="${RA}" dec="${DEC}" refframe=ICRS clobber=yes
    echo "Barycentered EVT2: ${EVT_BARY}" | tee -a "${LOG}"
    USE_EVT="${EVT_BARY}"
  else
    USE_EVT="${EVT_CLEAN}"
  fi

  # link/copy-like products matching old naming convention
  ln -sf "${USE_EVT}" "${WORKDIR}/evt2_final_${OBSID}.fits"

  if [[ -z "${REPRO_STACK}" ]]; then
    REPRO_STACK="${USE_EVT}"
  else
    REPRO_STACK="${REPRO_STACK},${USE_EVT}"
  fi

  # =========================
  # 4. Per-ObsID exposure maps/images: replaces manual images + adds exposure maps
  # =========================
  if [[ "${DO_EXPMAP}" == "1" ]]; then
    echo "===== ObsID ${OBSID}: fluximage =====" | tee -a "${LOG}"
    fluximage "${USE_EVT}" "${WORKDIR}/flux_${OBSID}/" bands="${BANDS}" binsize="${BINSIZE}" \
        units=default parallel=yes nproc="${NPROC}" clobber=yes verbose=1 \
        2>&1 | tee -a "${LOG}"
  fi

  # =========================
  # 5. Per-ObsID PSF maps, similar to mkpsfmap block in mulred/bcc
  # =========================
  if [[ "${DO_PSF}" == "1" ]]; then
    echo "===== ObsID ${OBSID}: mkpsfmap =====" | tee -a "${LOG}"
    IMG="${WORKDIR}/flux_${OBSID}/broad_thresh.img"
    if [[ ! -f "${IMG}" ]]; then
      # fallback image for old CIAO naming or if thresholded image absent
      IMG="${WORKDIR}/flux_${OBSID}/broad.img"
    fi
    if [[ ! -f "${IMG}" ]]; then
      IMG="${WORKDIR}/img_${OBSID}_0p5_7.fits"
      dmcopy "${USE_EVT}[energy=500:7000][bin x=::${BINSIZE},y=::${BINSIZE}]" "${IMG}" clobber=yes
    fi
    mkpsfmap "${IMG}" "${WORKDIR}/psf90_${OBSID}.fits" energy=2.3 ecf=0.9 units=arcsec clobber=yes
    mkpsfmap "${IMG}" "${WORKDIR}/psf50_${OBSID}.fits" energy=2.3 ecf=0.5 units=arcsec clobber=yes
  fi

  # =========================
  # 6. Per-ObsID wavdetect, now with exposure/PSF support when available
  # =========================
  if [[ "${DO_WAVDETECT}" == "1" ]]; then
    echo "===== ObsID ${OBSID}: wavdetect =====" | tee -a "${LOG}"
    IMG="${WORKDIR}/img_${OBSID}_2_8.fits"
    dmcopy "${USE_EVT}[energy=2000:8000][bin x=::1,y=::1]" "${IMG}" clobber=yes
    PSF="${WORKDIR}/psf90_${OBSID}.fits"
    punlearn wavdetect
    wavdetect infile="${IMG}" outfile="${WORKDIR}/src_${OBSID}.fits" \
        scellfile="${WORKDIR}/src_${OBSID}_scell.fits" \
        imagefile="${WORKDIR}/src_${OBSID}_imgfile.fits" \
        defnbkgfile="${WORKDIR}/src_${OBSID}_nbgd.fits" \
        regfile="${WORKDIR}/src_${OBSID}.reg" \
        psffile="${PSF}" scales="1 2 4 8" clobber=yes \
        2>&1 | tee -a "${LOG}"
  fi

done

# =========================
# 7. Merge products: replaces old reproject_obs/merge_all quick-look, but with modern merge_obs exposure-map handling
# =========================
if [[ "${DO_MERGE}" == "1" ]]; then
  echo "===== merge_obs all observations =====" | tee -a "${LOG}"
  merge_obs "${REPRO_STACK}" "${WORKDIR}/merged/" bands="${BANDS}" binsize="${BINSIZE}" \
      refcoord="${RA},${DEC}" units=default psfecf=0.9 psfmerge=exptime \
      parallel=yes nproc="${NPROC}" clobber=yes verbose=1 \
      2>&1 | tee -a "${LOG}"

  if [[ "${DO_WAVDETECT}" == "1" ]]; then
    echo "===== wavdetect merged broad image =====" | tee -a "${LOG}"
    MERGE_IMG="${WORKDIR}/merged/broad_thresh.img"
    MERGE_EXP="${WORKDIR}/merged/broad_thresh.expmap"
    MERGE_PSF="${WORKDIR}/merged/broad_thresh.psfmap"
    # If merge_obs did not make a psfmap with this exact name, make one from the merged image.
    if [[ ! -f "${MERGE_PSF}" ]]; then
      mkpsfmap "${MERGE_IMG}" "${WORKDIR}/merged/psf90_merged.fits" energy=2.3 ecf=0.9 units=arcsec clobber=yes
      MERGE_PSF="${WORKDIR}/merged/psf90_merged.fits"
    fi
    punlearn wavdetect
    wavdetect infile="${MERGE_IMG}" outfile="${WORKDIR}/merged/src_merged.fits" \
      scellfile="${WORKDIR}/merged/src_merged_scell.fits" \
      imagefile="${WORKDIR}/merged/src_merged_imgfile.fits" \
      defnbkgfile="${WORKDIR}/merged/src_merged_nbgd.fits" \
      regfile="${WORKDIR}/merged/src_merged.reg" \
      expfile="${MERGE_EXP}" psffile="${MERGE_PSF}" scales="1 2 4 8" clobber=yes \
      2>&1 | tee -a "${LOG}"
  fi
fi

echo "Pipeline finished: $(date)" | tee -a "${LOG}"
