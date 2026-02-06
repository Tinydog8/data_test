# 中国近海溢油公开数据资源——用于数值模拟验证

## 概述

以下整理了可公开获取的中国近海溢油事故观测数据和相关环境数据资源，可用于验证溢油数值模拟（油膜漂移扩散、风化过程等）的准确性。

---

## 1. 重大溢油事故案例数据

### 1.1 蓬莱 19-3 油田溢油事故（2011年）

- **事故概况**：2011年6月，渤海蓬莱19-3油田（康菲石油中国有限公司运营）发生溢油事故，是中国近海最严重的溢油事件之一。溢油面积累计超过 5500 km²，影响了渤海大面积海域。
- **位置**：渤海中部，约 38.3°N, 120.1°E
- **公开数据来源**：
  - **国家海洋局（现自然资源部）通报**：事故期间国家海洋局发布了多期监测通报，包含油膜分布范围、面积和漂移方向等信息。
  - **遥感影像**：
    - MODIS（Terra/Aqua）卫星影像可从 NASA Worldview / LAADS DAAC 免费获取。
    - SAR（合成孔径雷达）影像：ESA Envisat ASAR 数据可通过 ESA Earth Online 申请获取，覆盖了事故期间的渤海区域。
    - 中国环境卫星 HJ-1 数据可通过中国资源卫星应用中心（CRESDA）获取。
  - **文献参考**：
    - Guo, W., et al. (2013). "Monitoring and assessment of the Penglai 19-3 oil spill in the Bohai Sea using satellite SAR." *Marine Pollution Bulletin*, 71(1-2), 150-162.
    - Liu, X., et al. (2015). "Trajectory and fate of oil spill in the Bohai Sea: a numerical simulation study." *Journal of Ocean University of China*, 14(5), 751-758.
    - Pan, G., et al. (2015). "Numerical simulation of the oil spill transport in Bohai Sea." *Acta Oceanologica Sinica*, 34(5), 46-55.

### 1.2 大连新港输油管道爆炸溢油事故（2010年）

- **事故概况**：2010年7月16日，大连新港一条输油管道在作业中发生爆炸起火，导致约 1500 吨原油泄漏入海，污染面积约 430 km²。
- **位置**：大连新港附近，约 38.97°N, 121.63°E
- **公开数据来源**：
  - **遥感影像**：
    - MODIS 卫星可见光影像（NASA LAADS DAAC 免费获取）。
    - Landsat 影像（USGS EarthExplorer 免费获取）。
    - SAR 影像：ESA Envisat ASAR、RADARSAT-2（需申请）。
  - **文献参考**：
    - Guo, J., et al. (2012). "Monitoring of oil spill in Dalian using SAR imagery." *Proceedings of SPIE*, 8532.
    - Li, Y., et al. (2012). "Numerical simulation of the Dalian oil spill accident." *Marine Pollution Bulletin*, 64(12), 2757-2762.
    - Liu, S.Y., et al. (2011). "Remote sensing monitoring and numerical simulation of the Dalian oil spill." *Acta Oceanologica Sinica*.

### 1.3 青岛黄岛输油管道爆炸事故（2013年）

- **事故概况**：2013年11月22日，中石化东黄输油管道在青岛经济技术开发区段发生泄漏爆炸，部分原油流入胶州湾。
- **位置**：青岛市黄岛区，胶州湾附近，约 36.05°N, 120.19°E
- **公开数据来源**：
  - 相关遥感影像同样可从 MODIS、Landsat 获取。
  - **文献参考**：
    - 相关中文期刊论文可通过中国知网（CNKI）检索关键词"青岛溢油 数值模拟"获取。

---

## 2. 环境驱动数据（用于模拟输入）

验证数值模拟不仅需要溢油观测数据，还需要准确的环境场数据作为模拟输入：

### 2.1 海流数据

| 数据集 | 来源 | 分辨率 | 获取方式 |
|--------|------|--------|----------|
| HYCOM (HYbrid Coordinate Ocean Model) | HYCOM.org | 1/12° (~8km), 逐日/逐3小时 | 免费下载 |
| CMEMS (Copernicus Marine) 全球海流再分析 | Copernicus Marine Service | 1/12°, 逐日 | 免费注册下载 |
| FVCOM 渤海模型输出 | 相关文献 | 非结构网格，高分辨率 | 联系作者 |
| ROMS 中国近海区域模型 | 相关文献 | 可变分辨率 | 联系作者 |

### 2.2 风场数据

| 数据集 | 来源 | 分辨率 | 获取方式 |
|--------|------|--------|----------|
| ERA5 再分析风场 | ECMWF / CDS | 0.25° (~25km), 逐小时 | 免费注册下载 |
| NCEP/NCAR 再分析 | NOAA | 2.5°, 逐6小时 | 免费下载 |
| CCMP 卫星融合风场 | NASA PO.DAAC | 0.25°, 逐6小时 | 免费下载 |
| 中国气象局（CMA）再分析 | CMA | 多种分辨率 | 部分公开 |

### 2.3 海浪数据

| 数据集 | 来源 | 分辨率 | 获取方式 |
|--------|------|--------|----------|
| ERA5 海浪再分析 | ECMWF / CDS | 0.5°, 逐小时 | 免费注册下载 |
| WaveWatch III 全球模型 | NOAA/NCEP | 0.5°, 逐3小时 | 免费下载 |

### 2.4 海表温度（SST）

| 数据集 | 来源 | 分辨率 | 获取方式 |
|--------|------|--------|----------|
| OSTIA SST | CMEMS / Met Office | 0.05° (~5km), 逐日 | 免费下载 |
| MODIS SST | NASA OceanColor | 1km/4km, 逐日 | 免费下载 |

---

## 3. 遥感油膜监测数据来源

| 数据/平台 | 类型 | 用途 | 获取方式 |
|-----------|------|------|----------|
| Sentinel-1 SAR | C波段 SAR | 油膜检测（全天候） | ESA Copernicus Open Access Hub，免费 |
| RADARSAT-2 | C波段 SAR | 油膜检测 | CSA/MDA，需商业申请或科研合作 |
| Envisat ASAR（历史数据，至2012年） | C波段 SAR | 历史溢油事件监测 | ESA Earth Online，免费申请 |
| MODIS (Terra/Aqua) | 光学多光谱 | 大面积油膜可见光检测 | NASA LAADS DAAC，免费 |
| Landsat 8/9 OLI | 光学多光谱 | 高分辨率油膜检测 | USGS EarthExplorer，免费 |
| Sentinel-2 MSI | 光学多光谱 | 高分辨率油膜检测 | ESA Copernicus Open Access Hub，免费 |
| 高分系列卫星 (GF-1/2/3) | 光学/SAR | 油膜检测 | 中国资源卫星应用中心（CRESDA），部分公开 |

---

## 4. 中国海洋观测与数据共享平台

| 平台名称 | 网址 | 数据内容 |
|----------|------|----------|
| 国家海洋科学数据中心 | https://mds.nmdis.org.cn/ | 海洋环境观测数据、历史海洋资料 |
| 中国海洋信息网 | https://www.nmdis.org.cn/ | 海洋环境公报、海洋灾害公报 |
| 自然资源部海洋预警监测司 | http://www.mnr.gov.cn/ | 海洋环境状况公报、溢油事故通报 |
| 中国资源卫星应用中心（CRESDA） | https://www.cresda.com/ | 高分卫星、环境卫星影像 |
| 中国知网（CNKI） | https://www.cnki.net/ | 中文学术文献、溢油观测与模拟研究论文 |

---

## 5. 推荐验证策略

针对数值模拟验证，建议采用以下策略：

### 5.1 油膜轨迹验证
- 使用 SAR 或光学遥感影像提取油膜位置和范围，与模拟的油膜漂移轨迹进行对比。
- 推荐使用 **蓬莱19-3** 或 **大连新港** 事故数据，因为这两起事故持续时间较长，遥感覆盖较好，且已有大量文献可参考。

### 5.2 油膜面积验证
- 对比遥感提取的油膜面积与模拟的扩散面积随时间的变化。

### 5.3 到岸时间验证
- 对比模拟预测的油膜到达海岸线的时间与实际观测报告。

### 5.4 常用验证指标
- **油粒子命中率 (Hit Rate)**：模拟油粒子落入遥感观测油膜范围的比例。
- **重叠面积比 (Overlap Index)**：模拟油膜与观测油膜的重叠面积与并集面积之比。
- **质心偏差 (Centroid Deviation)**：模拟油膜质心与观测油膜质心的距离偏差。
- **RMSE**：轨迹点之间的均方根误差。

---

## 6. 关键参考文献

1. Guo, W., & Wang, Y. (2009). "A numerical oil spill model based on a hybrid method." *Marine Pollution Bulletin*, 58(5), 726-734.
2. Li, P., et al. (2016). "An overview of the oil spill models for the oil spill emergency response in China's seas." *Acta Oceanologica Sinica*, 35(9), 1-14.
3. Chen, J., et al. (2019). "A review of oil spill remote sensing." *Remote Sensing*, 11(11), 1305.
4. Guo, W., et al. (2013). "Monitoring the Penglai 19-3 oil spill using SAR." *Marine Pollution Bulletin*, 71(1-2), 150-162.
5. 国家海洋局. (2012). "蓬莱19-3油田溢油事故联合调查报告."
6. Berry, A., et al. (2012). "The oil spill model OILTRANS and its application to the Celtic Sea." *Marine Pollution Bulletin*, 64(11), 2489-2501.

---

## 总结

中国近海溢油数值模拟验证最为推荐的公开数据来源是 **2011年蓬莱19-3溢油事故** 和 **2010年大连新港溢油事故**。这两起事故：
- 影响范围大、持续时间长，遥感数据覆盖较完整；
- 国家海洋局发布了详细的监测通报；
- 已有大量中英文学术文献进行了分析和数值模拟研究，可作为对比基准；
- 环境驱动数据（海流、风场、海浪等）可从多个国际公开数据平台免费获取。

对于海流等环境场输入数据，推荐使用 **HYCOM** 或 **CMEMS** 海流再分析产品，配合 **ERA5** 风场再分析数据。遥感验证数据可优先选用 **Sentinel-1 SAR** 和 **MODIS** 光学影像。
