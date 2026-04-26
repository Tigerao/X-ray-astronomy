#!/usr/bin/env bash
# Minimal CIAO pipeline for ACIS: chandra_repro -> fluximage -> mkpsfmap -> axbary
# Usage:
#   ./min_acis_repro_flux_psf_bary.sh <obsid_dir> <evt2_file_or_auto> <orbit_file> <ra_deg> <dec_deg> [outdir]
#
# Example:
#   ./min_acis_repro_flux_psf_bary.sh ./12345 auto ./orbitf12345.fits 83.6331 22.0145 ./products

set -euo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  sed -n '2,9p' "$0"
  exit 0
fi

if [[ $# -lt 5 || $# -gt 6 ]]; then
  echo "[ERROR] 参数数量不正确。请使用: $0 <obsid_dir> <evt2_file_or_auto> <orbit_file> <ra_deg> <dec_deg> [outdir]" >&2
  exit 2
fi

obsid_dir="$1"
evt2_input="$2"
orbit_file="$3"
ra_deg="$4"
dec_deg="$5"
outdir="${6:-./ciao_products}"

logdir="${outdir}/logs"
mkdir -p "$outdir" "$logdir"
logfile="${logdir}/min_acis_repro_flux_psf_bary.$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$logfile") 2>&1

echo "[INFO] $(date -u '+%F %T UTC') pipeline start"

# ===== 环境检查 =====
if ! command -v punlearn >/dev/null 2>&1; then
  echo "[ERROR] 未检测到 CIAO 命令 punlearn。请先执行 ciao 环境初始化（如 source /path/to/ciao/bin/ciao.sh）。" >&2
  exit 3
fi

if ! command -v chandra_repro >/dev/null 2>&1 || ! command -v fluximage >/dev/null 2>&1 || ! command -v mkpsfmap >/dev/null 2>&1 || ! command -v axbary >/dev/null 2>&1; then
  echo "[ERROR] 缺少 CIAO 核心工具（chandra_repro/fluximage/mkpsfmap/axbary）。请确认 CIAO 已完整安装。" >&2
  exit 3
fi

if [[ ! -d "$obsid_dir" ]]; then
  echo "[ERROR] obsid_dir 不存在: $obsid_dir" >&2
  exit 4
fi

if [[ ! -f "$orbit_file" ]]; then
  echo "[ERROR] orbit file 不存在: $orbit_file" >&2
  exit 4
fi

if ! [[ "$ra_deg" =~ ^-?[0-9]+([.][0-9]+)?$ ]] || ! [[ "$dec_deg" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
  echo "[ERROR] RA/DEC 必须为十进制度浮点数。当前: RA=$ra_deg DEC=$dec_deg" >&2
  exit 5
fi

# ===== 目录约定 =====
repro_dir="${outdir}/repro"
flux_dir="${outdir}/flux"
psf_dir="${outdir}/psf"
bary_dir="${outdir}/bary"
mkdir -p "$repro_dir" "$flux_dir" "$psf_dir" "$bary_dir"

# ===== 1) 重处理 =====
echo "[INFO] Step 1: chandra_repro"
punlearn chandra_repro
chandra_repro "$obsid_dir" "$repro_dir" verbose=1 clobber=yes check_vf_pha=yes

# ===== 2) 事件文件定位 =====
if [[ "$evt2_input" == "auto" ]]; then
  evt2_file="$(find "$repro_dir" -maxdepth 1 -type f -name '*_evt2.fits' | head -n 1)"
  if [[ -z "$evt2_file" ]]; then
    echo "[ERROR] auto 模式未找到 *_evt2.fits。请检查 chandra_repro 输出。" >&2
    exit 6
  fi
else
  evt2_file="$evt2_input"
fi

if [[ ! -f "$evt2_file" ]]; then
  echo "[ERROR] evt2 file 不存在: $evt2_file" >&2
  exit 6
fi

# ===== 3) 能段筛选后生成 flux image =====
# CIAO syntax: energy in eV, bands with name:min:max
bands="broad:500:7000"
echo "[INFO] Step 2: fluximage (bands=${bands})"
punlearn fluximage
fluximage "$evt2_file" "$flux_dir" bands="$bands" binsize=1 verbose=1 clobber=yes

# 自动找 broad thresh.img 作为 PSF map 输入
img_file="$(find "$flux_dir" -maxdepth 1 -type f -name '*broad*thresh.img' | head -n 1)"
if [[ -z "$img_file" ]]; then
  img_file="$(find "$flux_dir" -maxdepth 1 -type f -name '*broad*.img' | head -n 1)"
fi
if [[ -z "$img_file" ]]; then
  echo "[ERROR] 未找到 fluximage 生成的 broad 图像。" >&2
  exit 7
fi

# ===== 4) 生成 PSF map =====
# fluximage 通常会直接生成 *psfmap；若存在则直接复用，避免重复调用 mkpsfmap。
existing_psfmap="$(find "$flux_dir" -maxdepth 1 -type f -name '*broad*psfmap*' | head -n 1)"
if [[ -n "$existing_psfmap" ]]; then
  psfmap_file="$existing_psfmap"
  echo "[INFO] Step 3: 复用 fluximage 生成的 PSF map: $psfmap_file"
else
  # CIAO 4.15 的 mkpsfmap 不接受 verbose 参数；传入会触发 parammatch 错误。
  psfmap_file="${psf_dir}/broad_psfmap.fits"
  echo "[INFO] Step 3: mkpsfmap"
  punlearn mkpsfmap
  mkpsfmap infile="$img_file" outfile="$psfmap_file" energy=1.496 ecf=0.9 clobber=yes
fi

# ===== 5) 重心改正（barycentric correction） =====
bary_evt_file="${bary_dir}/$(basename "${evt2_file%.fits}")_bary.fits"
echo "[INFO] Step 4: axbary"
punlearn axbary
axbary "$evt2_file" "$orbit_file" "$bary_evt_file" ra="$ra_deg" dec="$dec_deg" clobber=yes

echo "[INFO] pipeline done"
echo "[INFO] Repro EVT2 : $evt2_file"
echo "[INFO] Flux dir   : $flux_dir"
echo "[INFO] PSF map    : $psfmap_file"
echo "[INFO] Bary evt   : $bary_evt_file"
echo "[INFO] Log file   : $logfile"
