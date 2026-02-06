## 桑吉号（Sanchi）溢油事件：可复现数据集构建流程（Python）

这个仓库提供一个**可复现**的 Python 流程，用于围绕桑吉号事故时段自动构建“用于数值模拟验证”的观测数据集：

- 自动检索 Sentinel‑1 GRD 场景清单（ASF Search，输出 CSV/GeoJSON）
- （可选）按清单批量下载原始 GRD `.zip`（通常需要 NASA Earthdata 账号）
- 对**已预处理的 GeoTIFF**做油膜候选区域提取（输出 GeoJSON 多边形）
- 生成时间序列指标（面积、质心、外接框、像素面积等，输出 CSV）

> 说明：原始 Sentinel‑1 GRD（SAFE/ZIP）到可分析 GeoTIFF（标定/去斑/地形校正/裁剪）的预处理，工程上最稳定的是 ESA SNAP `gpt`。本项目提供 Python 封装与 SNAP graph 模板；如果你已经有预处理后的 GeoTIFF，也可以直接跳过该步。

---

## 1) 安装

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

---

## 2) 自动拉取场景清单（推荐从这里开始）

默认使用一个覆盖事故附近海域的 bbox 和 2018-01 的时间窗（可改）：

```bash
python3 -m sanchi_pipeline search \
  --start 2018-01-06T00:00:00Z \
  --end   2018-01-20T23:59:59Z \
  --bbox  123.5 26.5 126.5 29.5 \
  --outdir data/scenes
```

输出：

- `data/scenes/scenes.geojson`：每个场景 footprint + 元数据（含下载 URL）
- `data/scenes/scenes.csv`：便于筛选/排序的表格清单

---

## 3) （可选）批量下载 GRD ZIP

ASF Datapool 的下载通常需要 Earthdata 登录。推荐用环境变量传入：

```bash
export EARTHDATA_USERNAME='你的用户名'
export EARTHDATA_PASSWORD='你的密码'

python3 -m sanchi_pipeline download \
  --scenes data/scenes/scenes.geojson \
  --outdir data/raw_zips
```

---

## 4) 预处理（两种方式）

### A. 使用 SNAP `gpt`（推荐、最可复现）

安装 SNAP 后，执行：

```bash
python3 -m sanchi_pipeline snap-preprocess \
  --input-zip data/raw_zips/S1A_*.zip \
  --outdir data/geotiff \
  --bbox 123.5 26.5 126.5 29.5
```

会输出裁剪后的 GeoTIFF（示例为 VV Sigma0 dB）。你也可以自行修改 `snap_graphs/` 下的 graph。

### B. 已有 GeoTIFF

如果你已经从其他流程获得了 GeoTIFF（例如 RTC 产品或自制 SNAP 输出），直接进入下一步即可。

---

## 5) 油膜多边形提取（GeoTIFF → GeoJSON）

```bash
python3 -m sanchi_pipeline extract \
  --input data/geotiff \
  --outdir data/oil_polygons \
  --min-area-m2 50000
```

输出：每景一个 `*.geojson`（FeatureCollection），并生成汇总表（可选见下一步）。

---

## 6) 时间序列指标（GeoJSON → CSV）

```bash
python3 -m sanchi_pipeline timeseries \
  --polygons data/oil_polygons \
  --out data/oil_timeseries.csv
```

指标包含：`sceneName, startTime, stopTime, area_m2, centroid_lon, centroid_lat, bbox_*` 等。

---

## 重要提醒（科学与工程层面）

- SAR “暗斑”并不总是油膜：低风区、内波、雨、海面膜等都会产生假阳性。本文提供的是**可复现的基线提取**，真实科研验证建议叠加风场/海况、人工审阅或更强模型。
- 如果你的目标是验证数值模拟的漂移扩散，最常用的量化指标是：
  - 质心漂移误差（km）
  - 油膜面积误差（m²）
  - 形状重叠（IoU / Dice）

