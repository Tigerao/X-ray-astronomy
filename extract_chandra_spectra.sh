#!/usr/bin/env bash
#
# extract_chandra_spectra_auto.sh
#
# Extract Chandra/ACIS source and background spectra with CIAO specextract.
#
# This version is designed for the directory layout used by the course/GitHub
# timing-extraction products:
#
#   ${BASE}/merge_data/xdata/all_bcc_${OBSID}_evt.fits
#   ${BASE}/merge_data/timing/reg/region_${OBSID}/region_${ECF}/${SRCID}.reg
#   ${BASE}/merge_data/timing/reg/region_${OBSID}/region_${ECF}/${SRCID}_bkg.reg
#
# It does NOT rely on relative ASOLFILE/BPIXFILE/MSKFILE keywords stored in the
# event header.  Instead, it searches under --ancillary-root, writes an asol list
# when multiple aspect-solution files exist, and explicitly passes:
#
#   asp=...
#   mskfile=...
#   badpixfile=...
#
# to specextract.
#
# Example:
#   bash extract_chandra_spectra_auto.sh all \
#     --base /data/home/tiger/chandra/course \
#     --obsids 914 \
#     --srcids 1 \
#     --ecf 90
#
# Dry run:
#   bash extract_chandra_spectra_auto.sh all --base /data/home/tiger/chandra/course \
#     --obsids 914 --srcids 1 --dry-run

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  bash extract_chandra_spectra_auto.sh MODE --base BASE --obsids OBSIDS --srcids SRCIDS [options]

Required:
  MODE                  A label used in the output directory, e.g. all, flare, nonflare.
  --base BASE           Project base directory, e.g. /data/home/tiger/chandra/course.
  --obsids LIST         ObsID list, comma or space separated: "914,23603" or "914 23603".
  --srcids LIST         Source ID list, comma or space separated: "1,2,3" or "1 2 3".

Common options:
  --ecf N               Region ECF label. Default: 90.
  --emin EV             Minimum event energy in eV. Default: 500.
  --emax EV             Maximum event energy in eV. Default: 8000.
  --outdir DIR          Output root. Default: ${BASE}/merge_data/spec/spec_${MODE}_p${ECF}.
  --ancillary-root DIR  Where to search for asol/msk/bpix. Default: ${BASE}.
  --dry-run             Print commands and paths but do not run specextract.

Advanced path templates:
  --evt-template T      Event path relative to BASE. Use {obsid} and/or {obsid5}.
                        Default: merge_data/xdata/all_bcc_{obsid}_evt.fits
  --src-reg-template T  Source region path relative to BASE. Use {obsid}, {obsid5}, {ecf}, {srcid}.
                        Default: merge_data/timing/reg/region_{obsid}/region_{ecf}/{srcid}.reg
  --bkg-reg-template T  Background region path relative to BASE. Use {obsid}, {obsid5}, {ecf}, {srcid}.
                        Default: merge_data/timing/reg/region_{obsid}/region_{ecf}/{srcid}_bkg.reg

specextract options:
  --weight yes|no       specextract weight parameter. Default: no.
  --correctpsf yes|no   specextract correctpsf parameter. Default: no.
  --bkgresp yes|no      Generate background response. Default: yes.
  --grouptype VALUE     Source grouping. Default: NONE.
  --binspec VALUE       Source grouping value. Default: NONE.
  --verbose N           specextract verbose level. Default: 1.

Outputs:
  ${OUTDIR}/src_${SRCID}/obs_${OBSID}/src_${SRCID}_obs_${OBSID}.pi
  ${OUTDIR}/src_${SRCID}/obs_${OBSID}/src_${SRCID}_obs_${OBSID}.arf
  ${OUTDIR}/src_${SRCID}/obs_${OBSID}/src_${SRCID}_obs_${OBSID}.rmf
  plus background and log files created by specextract.
USAGE
}

# ------------------------- default parameters -------------------------
MODE="${1:-}"
if [[ -z "$MODE" || "$MODE" == "-h" || "$MODE" == "--help" ]]; then
    usage
    exit 0
fi
shift

BASE=""
OBSIDS_RAW=""
SRCIDS_RAW=""
ECF="90"
EMIN="500"
EMAX="8000"
OUTDIR=""
ANC_ROOT=""
DRY_RUN=0

EVT_TEMPLATE='merge_data/xdata/all_bcc_{obsid}_evt.fits'
SRC_REG_TEMPLATE='merge_data/timing/reg/region_{obsid}/region_{ecf}/{srcid}.reg'
BKG_REG_TEMPLATE='merge_data/timing/reg/region_{obsid}/region_{ecf}/{srcid}_bkg.reg'

WEIGHT="no"
CORRECTPSF="no"
BKGRESP="yes"
GROUPTYPE="NONE"
BINSPEC="NONE"
VERBOSE="1"

# ------------------------- parse command line -------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --base) BASE="$2"; shift 2 ;;
        --obsids) OBSIDS_RAW="$2"; shift 2 ;;
        --srcids) SRCIDS_RAW="$2"; shift 2 ;;
        --ecf) ECF="$2"; shift 2 ;;
        --emin) EMIN="$2"; shift 2 ;;
        --emax) EMAX="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --ancillary-root) ANC_ROOT="$2"; shift 2 ;;
        --evt-template) EVT_TEMPLATE="$2"; shift 2 ;;
        --src-reg-template) SRC_REG_TEMPLATE="$2"; shift 2 ;;
        --bkg-reg-template) BKG_REG_TEMPLATE="$2"; shift 2 ;;
        --weight) WEIGHT="$2"; shift 2 ;;
        --correctpsf) CORRECTPSF="$2"; shift 2 ;;
        --bkgresp) BKGRESP="$2"; shift 2 ;;
        --grouptype) GROUPTYPE="$2"; shift 2 ;;
        --binspec) BINSPEC="$2"; shift 2 ;;
        --verbose) VERBOSE="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done

if [[ -z "$BASE" || -z "$OBSIDS_RAW" || -z "$SRCIDS_RAW" ]]; then
    echo "ERROR: --base, --obsids, and --srcids are required." >&2
    usage
    exit 2
fi

BASE="${BASE%/}"
ANC_ROOT="${ANC_ROOT:-$BASE}"
ANC_ROOT="${ANC_ROOT%/}"
OUTDIR="${OUTDIR:-${BASE}/merge_data/spec/spec_${MODE}_p${ECF}}"

# Convert comma-separated or space-separated strings into arrays.
parse_list() {
    local raw="$1"
    raw="${raw//,/ }"
    # shellcheck disable=SC2206
    local arr=( $raw )
    printf '%s\n' "${arr[@]}"
}

mapfile -t OBSIDS < <(parse_list "$OBSIDS_RAW")
mapfile -t SRCIDS < <(parse_list "$SRCIDS_RAW")

# Replace placeholders in a relative template and prepend BASE.
expand_template() {
    local template="$1"
    local obsid="$2"
    local obsid5="$3"
    local srcid="$4"
    local out="$template"
    out="${out//\{obsid\}/$obsid}"
    out="${out//\{obsid5\}/$obsid5}"
    out="${out//\{ecf\}/$ECF}"
    out="${out//\{srcid\}/$srcid}"
    if [[ "$out" = /* ]]; then
        printf '%s\n' "$out"
    else
        printf '%s/%s\n' "$BASE" "$out"
    fi
}

require_file() {
    local label="$1"
    local path="$2"
    if [[ ! -f "$path" ]]; then
        echo "ERROR: missing ${label}: $path" >&2
        return 1
    fi
}

# Return sorted matches for a find expression.  This intentionally searches
# beneath ANC_ROOT rather than trusting possibly-stale relative header keywords.
find_matches() {
    local pattern1="$1"
    local pattern2="${2:-}"
    if [[ -n "$pattern2" ]]; then
        find "$ANC_ROOT" -type f \( -name "$pattern1" -o -name "$pattern2" \) 2>/dev/null | sort -u
    else
        find "$ANC_ROOT" -type f -name "$pattern1" 2>/dev/null | sort -u
    fi
}

# Prefer files in repro directories, but fall back to any sorted match.
choose_first_prefer_repro() {
    local file
    local first=""
    while IFS= read -r file; do
        [[ -z "$file" ]] && continue
        [[ -z "$first" ]] && first="$file"
        if [[ "$file" == *repro* ]]; then
            printf '%s\n' "$file"
            return 0
        fi
    done
    [[ -n "$first" ]] && printf '%s\n' "$first"
}

find_asol_stack() {
    local obsid="$1"
    local obsid5="$2"
    local listfile="$3"
    local -a hits=()

    # pcadf00914_000N001_asol1.fits is the common pattern.
    mapfile -t hits < <(find_matches "pcadf${obsid5}*asol*.fits" "pcad*${obsid5}*asol*.fits")

    # Fallback for less-standard names that include unpadded ObsID.
    if [[ ${#hits[@]} -eq 0 ]]; then
        mapfile -t hits < <(find_matches "*${obsid}*asol*.fits")
    fi

    if [[ ${#hits[@]} -eq 0 ]]; then
        echo "ERROR: cannot find aspect solution file for ObsID ${obsid}." >&2
        echo "Tried under: ${ANC_ROOT}" >&2
        echo "Useful manual check: find ${ANC_ROOT} -name '*${obsid5}*asol*.fits'" >&2
        return 1
    fi

    mkdir -p "$(dirname "$listfile")"
    printf '%s\n' "${hits[@]}" > "$listfile"

    if [[ ${#hits[@]} -eq 1 ]]; then
        printf '%s\n' "${hits[0]}"
    else
        printf '@%s\n' "$listfile"
    fi
}

find_msk_file() {
    local obsid="$1"
    local obsid5="$2"
    local hits chosen
    hits=$(find_matches "acisf${obsid5}*msk*.fits" "*${obsid5}*msk*.fits") || true
    if [[ -z "$hits" ]]; then
        hits=$(find_matches "*${obsid}*msk*.fits") || true
    fi
    chosen=$(printf '%s\n' "$hits" | choose_first_prefer_repro || true)
    if [[ -z "$chosen" ]]; then
        echo "ERROR: cannot find mask file for ObsID ${obsid}." >&2
        echo "Useful manual check: find ${ANC_ROOT} -name '*${obsid5}*msk*.fits'" >&2
        return 1
    fi
    printf '%s\n' "$chosen"
}

find_bpix_file() {
    local obsid="$1"
    local obsid5="$2"
    local hits chosen
    hits=$(find_matches "acisf${obsid5}*bpix*.fits" "*${obsid5}*bpix*.fits") || true
    if [[ -z "$hits" ]]; then
        hits=$(find_matches "*${obsid}*bpix*.fits") || true
    fi
    chosen=$(printf '%s\n' "$hits" | choose_first_prefer_repro || true)
    if [[ -z "$chosen" ]]; then
        echo "ERROR: cannot find bad-pixel file for ObsID ${obsid}." >&2
        echo "Useful manual check: find ${ANC_ROOT} -name '*${obsid5}*bpix*.fits'" >&2
        return 1
    fi
    printf '%s\n' "$chosen"
}

run_or_print() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf 'DRY-RUN:'
        printf ' %q' "$@"
        printf '\n'
    else
        "$@"
    fi
}

# Confirm CIAO command exists early.  In dry-run mode, do not require it.
if [[ "$DRY_RUN" -eq 0 ]] && ! command -v specextract >/dev/null 2>&1; then
    echo "ERROR: specextract not found. Did you activate CIAO?" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# ----------------------------- main loop ------------------------------
for OBSID in "${OBSIDS[@]}"; do
    OBSID5=$(printf '%05d' "$OBSID")

    echo "============================================================"
    echo "ObsID: ${OBSID}  padded: ${OBSID5}"
    echo "Searching ancillary files under: ${ANC_ROOT}"

    # Per-ObsID ancillary files are independent of source region.
    ANC_DIR="${OUTDIR}/ancillary/obs_${OBSID}"
    ASOL_LIST="${ANC_DIR}/asol_${OBSID}.lis"
    ASP=$(find_asol_stack "$OBSID" "$OBSID5" "$ASOL_LIST")
    MSK=$(find_msk_file "$OBSID" "$OBSID5")
    BPIX=$(find_bpix_file "$OBSID" "$OBSID5")

    echo "ASP        = ${ASP}"
    echo "MSKFILE    = ${MSK}"
    echo "BADPIXFILE = ${BPIX}"

    for SRCID in "${SRCIDS[@]}"; do
        EVT=$(expand_template "$EVT_TEMPLATE" "$OBSID" "$OBSID5" "$SRCID")
        REG=$(expand_template "$SRC_REG_TEMPLATE" "$OBSID" "$OBSID5" "$SRCID")
        BKGREG=$(expand_template "$BKG_REG_TEMPLATE" "$OBSID" "$OBSID5" "$SRCID")

        require_file "event file" "$EVT"
        require_file "source region" "$REG"
        require_file "background region" "$BKGREG"
        require_file "mask file" "$MSK"
        require_file "bad-pixel file" "$BPIX"

        SRC_OUTDIR="${OUTDIR}/src_${SRCID}/obs_${OBSID}"
        mkdir -p "$SRC_OUTDIR"
        OUTROOT="${SRC_OUTDIR}/src_${SRCID}_obs_${OBSID}"
        LOG="${OUTROOT}_specextract.log"

        INFILE="${EVT}[sky=region(${REG})][energy=${EMIN}:${EMAX}]"
        BKGFILE="${EVT}[sky=region(${BKGREG})][energy=${EMIN}:${EMAX}]"

        echo "------------------------------------------------------------"
        echo "Source ID  = ${SRCID}"
        echo "EVT        = ${EVT}"
        echo "SRC REG    = ${REG}"
        echo "BKG REG    = ${BKGREG}"
        echo "OUTROOT    = ${OUTROOT}"
        echo "LOG        = ${LOG}"

        # punlearn prevents hidden parameter values from previous interactive runs.
        run_or_print punlearn specextract

        if [[ "$DRY_RUN" -eq 1 ]]; then
            run_or_print specextract \
                infile="$INFILE" \
                bkgfile="$BKGFILE" \
                outroot="$OUTROOT" \
                asp="$ASP" \
                mskfile="$MSK" \
                badpixfile="$BPIX" \
                bkgresp="$BKGRESP" \
                weight="$WEIGHT" \
                correctpsf="$CORRECTPSF" \
                grouptype="$GROUPTYPE" \
                binspec="$BINSPEC" \
                bkg_grouptype="$GROUPTYPE" \
                bkg_binspec="$BINSPEC" \
                clobber=yes \
                verbose="$VERBOSE"
        else
            specextract \
                infile="$INFILE" \
                bkgfile="$BKGFILE" \
                outroot="$OUTROOT" \
                asp="$ASP" \
                mskfile="$MSK" \
                badpixfile="$BPIX" \
                bkgresp="$BKGRESP" \
                weight="$WEIGHT" \
                correctpsf="$CORRECTPSF" \
                grouptype="$GROUPTYPE" \
                binspec="$BINSPEC" \
                bkg_grouptype="$GROUPTYPE" \
                bkg_binspec="$BINSPEC" \
                clobber=yes \
                verbose="$VERBOSE" 2>&1 | tee "$LOG"
        fi
    done
done

echo "============================================================"
echo "Done. Output root: ${OUTDIR}"
