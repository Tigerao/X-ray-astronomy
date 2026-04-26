# Chandra ACIS 时序 + 光谱分析流程（按任务顺序）

下面给出**可直接落地**、并与本仓库脚本一致的工作流：

1. 预处理数据（CIAO）
2. 抽取光子（本仓库脚本）
3. 周期搜寻（LS / GL）
4. 抽取光谱并用 `ftgrouppha` 分组（bin）

> 下面示例统一使用：
>
> - `BASE=/data/home/tiger/chandra/course`
> - `OBSID=23603`
> - `SRCID=1`

---

## 0. 环境准备（建议）

```bash
# CIAO 环境
ciao

# 如果要用 ftgrouppha，还需 HEASoft 环境（按你的机器实际安装方式加载）
# 例如：source /path/to/heasoft-init.sh
```

确保以下工具可用：`chandra_repro`、`fluximage`、`mkpsfmap`、`axbary`、`specextract`、`ftgrouppha`。

---

## 1) 预处理数据（CIAO）

本仓库已提供模块化脚本：`acis_course_modular.sh`。

### 1.1 一步跑完（推荐）

```bash
bash acis_course_modular.sh all 23603
```

### 1.2 产物检查

```bash
bash acis_course_modular.sh status 23603
```

关键产物（供后续步骤使用）：

- `${BASE}/merge_data/xdata/all_bcc_23603_evt.fits`
- `${BASE}/merge_data/xdata/psf90_23603_500_8000.fits`
- `${BASE}/merge_data/xdata/expmap_23603_500_8000.fits`
- `${BASE}/merge_data/xdata/img_23603_500_8000.fits`

---

## 2) 抽取光子（按脚本真实接口）

`extract_chandra_photons.py` 是**子命令式**接口，建议直接用 `all`（内部顺序：`make-regions -> extract -> epoch -> merge`）。

```bash
python extract_chandra_photons.py all \
  --base /data/home/tiger/chandra/course \
  --obsids 23603 \
  --sources sources.csv \
  --ecf 90 \
  --emin 500 \
  --emax 8000
```

### 2.1 你会得到的核心文件

- 区域文件：
  - `timing/reg/region_23603/region_90/1.reg`
  - `timing/reg/region_23603/region_90/1_bkg.reg`
- 单个 ObsID 光子：
  - `timing/txt/txt_23603_p90/1.txt`
  - `timing/txt/txt_23603_p90/1_bkg.txt`
- 合并后（供周期搜索）：
  - `timing/txt/txt_all_obs_p90/1.txt`
  - `timing/txt/txt_all_obs_p90/1_bkg.txt`
  - `timing/txt/txt_all_obs_p90/epoch_src_1.txt`

---

## 3) 周期搜寻

本仓库提供两条**独立**流程：

- `run_ls_search.py`：Lomb–Scargle（基于分箱光变）
- `run_gl_search.py`：Gregory–Loredo 风格（基于非分箱事件）

### 3.1 LS 搜寻

```bash
python run_ls_search.py \
  --base /data/home/tiger/chandra/course \
  --srcid 1 \
  --pmin 100 \
  --pmax 20000 \
  --dt 50
```

### 3.2 GL 搜寻

```bash
python run_gl_search.py \
  --base /data/home/tiger/chandra/course \
  --srcid 1 \
  --pmin 100 \
  --pmax 20000 \
  --mmax 12 \
  --ni 10
```

输出目录（两者一致）：

- `${BASE}/merge_data/timing/period_search/src_1/`

实践建议：

- 对比文献周期；
- 重点检查 Chandra 抖动相关别名：706.96 s、999.96 s 及其谐波；
- 若有效覆盖周期数很少，不宜做强结论。

---

## 4) 抽光谱 + `ftgrouppha` bin 光谱

> 这里使用第 2 步生成的源/背景 region，确保与时序分析口径一致。

### 4.1 用 `specextract` 抽谱（每个 ObsID、每个源各跑一次）

```bash
BASE=/data/home/tiger/chandra/course
OBSID=23603
SRCID=1

EVT="${BASE}/merge_data/xdata/all_bcc_${OBSID}_evt.fits"
REG="${BASE}/merge_data/timing/reg/region_${OBSID}/region_90/${SRCID}.reg"
BKGREG="${BASE}/merge_data/timing/reg/region_${OBSID}/region_90/${SRCID}_bkg.reg"
OUTDIR="${BASE}/merge_data/spectra/obs_${OBSID}"
OUTROOT="${OUTDIR}/src_${SRCID}_${OBSID}"

mkdir -p "${OUTDIR}"

specextract \
  infile="${EVT}[sky=region(${REG})]" \
  bkgfile="${EVT}[sky=region(${BKGREG})]" \
  outroot="${OUTROOT}" \
  bkgresp=yes \
  weight=no \
  correctpsf=yes \
  clobber=yes
```

常见输出：

- `${OUTROOT}.pi`（源谱）
- `${OUTROOT}_bkg.pi`（背景谱）
- `${OUTROOT}.arf`
- `${OUTROOT}.rmf`

### 4.2 用 `ftgrouppha` 分组（bin）

```bash
ftgrouppha \
  infile="${OUTROOT}.pi" \
  outfile="${OUTROOT}_grp20.pi" \
  grouptype=optmin \
  groupscale=20 \
  clobber=yes
```

含义：每个 bin 至少约 20 counts（常用于后续 χ² 拟合；若用 C-stat 可考虑更细分组或不分组）。

### 4.3 批处理示例（多个源）

```bash
BASE=/data/home/tiger/chandra/course
OBSID=23603

for SRCID in $(awk -F, 'NR>1{print $1}' sources.csv); do
  EVT="${BASE}/merge_data/xdata/all_bcc_${OBSID}_evt.fits"
  REG="${BASE}/merge_data/timing/reg/region_${OBSID}/region_90/${SRCID}.reg"
  BKGREG="${BASE}/merge_data/timing/reg/region_${OBSID}/region_90/${SRCID}_bkg.reg"
  OUTDIR="${BASE}/merge_data/spectra/obs_${OBSID}"
  OUTROOT="${OUTDIR}/src_${SRCID}_${OBSID}"

  mkdir -p "${OUTDIR}"

  specextract \
    infile="${EVT}[sky=region(${REG})]" \
    bkgfile="${EVT}[sky=region(${BKGREG})]" \
    outroot="${OUTROOT}" \
    bkgresp=yes weight=no correctpsf=yes clobber=yes

  ftgrouppha \
    infile="${OUTROOT}.pi" \
    outfile="${OUTROOT}_grp20.pi" \
    grouptype=optmin groupscale=20 clobber=yes
done
```

---

## 5) 最小执行清单（从头到尾）

```bash
# 1) 预处理
bash acis_course_modular.sh all 23603

# 2) 抽光子
python extract_chandra_photons.py all \
  --base /data/home/tiger/chandra/course \
  --obsids 23603 \
  --sources sources.csv \
  --ecf 90 --emin 500 --emax 8000

# 3) 周期搜寻（示例：srcid=1）
python run_ls_search.py --base /data/home/tiger/chandra/course --srcid 1 --pmin 100 --pmax 20000 --dt 50
python run_gl_search.py --base /data/home/tiger/chandra/course --srcid 1 --pmin 100 --pmax 20000 --mmax 12 --ni 10

# 4) 抽光谱并分组（示例：srcid=1）
#    先 specextract，再 ftgrouppha（见上文 4.1 / 4.2）
```
