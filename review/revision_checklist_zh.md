# 论文修改执行清单

> 使用方式：建议作者按“高优先级 -> 中优先级 -> 低优先级”顺序逐项核对。  
> 状态栏可自行填写：`未开始 / 进行中 / 已完成 / 不采用`

---

## 一、高优先级：决定论文是否成立的修改项

| 优先级 | 修改项 | 具体动作 | 状态 |
|---|---|---|---|
| 高 | 明确研究创新点 | 在引言末尾用 1 段话明确：本文不是单纯研究 overlap，而是研究**日变化分层**与**潮流-风浪相对位相**如何共同决定 overlap 的时机、持续时间和重复次数 | 未开始 |
| 高 | 给出 overlap/merge 的操作性定义 | 在方法或结果开头明确：如何定义 SBL/BBL 厚度；何时判定 overlap；merge 与 overlap 是否同义；判据是否对所有 case 一致 | 未开始 |
| 高 | 收紧对 Langmuir supercells 的表述 | 决定采用以下其一：1）降调为 “Langmuir-supercell-like structures”；2）补充直接流场图证明为 Langmuir supercells | 未开始 |
| 高 | 处理“无法解释”的关键现象 | 将文中“no reasonable explanation...”改写为“稳健现象 + 可能机制 + 未来工作”，避免主文中直接承认无法解释 | 未开始 |
| 高 | 增强 Ri / shear / Pk 的定量证据 | 至少增加 1 组补充诊断：代表性相位 shear 剖面、Pk 时间序列、Ri-Pk-k 联合对比、aligned vs misaligned 对比图 | 未开始 |
| 高 | 重写结论 | 将结论由“结果重复”改为“机制总结 + 参数化启示 + 适用范围 + 局限性” | 未开始 |

---

## 二、中优先级：显著影响论文质量与说服力

| 优先级 | 修改项 | 具体动作 | 状态 |
|---|---|---|---|
| 中 | 统一术语体系 | 统一 SBL/OSBL、merge/overlap、full-column/full-depth、positive/negative polarity、alignment 等术语 | 未开始 |
| 中 | 规范三类 regime 命名 | 全文统一为 `fully separated regime`、`once-daily overlap regime`、`twice-daily overlap regime` | 未开始 |
| 中 | 强化理想化设置的合理性说明 | 在方法或讨论中说明：为什么采用恒定风、恒定波、圆形潮流；这些简化支持哪些结论、限制哪些外推 | 未开始 |
| 中 | 增加总结图或 regime 示意图 | 建议增加一个 summary figure / regime diagram，概括热分层与 alignment 如何决定三类型态 | 未开始 |
| 中 | 精简结果段重复表述 | 第 3 节讲现象、第 4 节讲能量学、第 5 节讲谱、第 6 节讲各向异性；避免同一结论在多节重复 | 未开始 |
| 中 | 限定 NP/PP 结论的适用条件 | 在相关段落明确“该结论是在本文参数组合与北半球设定下成立” | 未开始 |

---

## 三、语言与写作层面：建议全文统一修正

### 1. 标题、Key Points、摘要

| 修改项 | 具体动作 | 状态 |
|---|---|---|
| 调整标题表达 | 可考虑突出“boundary-layer overlap”或“joint control” | 未开始 |
| 重写 Key Points | 避免 `once-overlap`, `twice-overlap`, `Boundary layers merge is...` 等不地道表达 | 未开始 |
| 压缩摘要句长 | 将过长句拆分，减少套话，提高机制表达密度 | 未开始 |
| 降低结论性措辞强度 | 对证据链尚不够直接的结论使用 `consistent with`, `suggest`, `indicate` 等表述 | 未开始 |

### 2. 共性语言问题

| 常见问题 | 修改要求 | 状态 |
|---|---|---|
| 连字符 | 统一 `wave-driven`, `current-driven`, `large-eddy`, `12-hour` 等写法 | 未开始 |
| 主谓一致 | 全文检查 `boundary layers merge is`、`results shows` 等问题 | 未开始 |
| 冠词与复数 | 检查 `the large eddy simulation`、`the turbulence` 等不自然表达 | 未开始 |
| 时态 | 方法通常用一般现在时或一般过去时之一，全文统一 | 未开始 |
| 介词搭配 | 如 `overlap once per day`、`transition from ... to ...`、`consistent with` | 未开始 |
| 中式英语 | 如 `have chances to`, `is the primary role`, `is being set` 等应全部改为自然英文 | 未开始 |

### 3. 代表性语句重点修正

| 原问题类型 | 建议动作 | 状态 |
|---|---|---|
| `Our simulations using Oceananigans ...` 句子不完整 | 重写整句，使主谓结构完整 | 未开始 |
| `the vertical temperature gradient is always satisfying...` | 改为简洁客观表达，如 `remains` / `stays` | 未开始 |
| `Tk is the primary role` | 改为 `Tk plays the primary role` | 未开始 |
| `These powerful turbulence events` | 改为 `These intense turbulent events` | 未开始 |
| `no satisfactory explanation` | 改成“可能机制 + scope limitation”结构 | 未开始 |

---

## 四、分节修改核对表

### 1. Introduction

- [ ] 引言末尾单独成段写清楚本文知识缺口  
- [ ] 说明与 Yan et al. (2022)、Valcarcel et al. (2025)、Shrestha et al. (2019) 的区别  
- [ ] 将“poorly understood”改为更具体的问题陈述  
- [ ] 明确本文要回答的是“什么时候 overlap、为什么 overlap、什么时候形成全水柱湍流”  

### 2. Method

- [ ] 检查所有方程变量首次出现是否定义  
- [ ] 检查单位格式：数字与单位之间留空格  
- [ ] 检查 `g`、盐度、温度、热通量等物理量写法是否规范  
- [ ] 明确 spin-up 长度与统计稳态判据  
- [ ] 明确 case 命名规则及控制实验含义  
- [ ] 说明圆形潮流和固定风浪方向的物理目的  

### 3. Results: Boundary-layer evolution

- [ ] 给出三类 regime 的明确判据  
- [ ] 区分 `overlap`、`merge`、`interaction` 的含义  
- [ ] 对 NP/PP 差异加入条件限定  
- [ ] 如果可能，增加 regime 总结图  

### 4. Results: Energy budget

- [ ] 每个 case 按统一逻辑写：主导项 -> 中层变化 -> 结论  
- [ ] 删除过细但对主结论帮助不大的枝节叙述  
- [ ] 对 pressure transport 与 internal waves 的关系使用更谨慎表述  
- [ ] 检查 budget 各项符号、正负号、归一化说明是否一致  

### 5. Results: Spectral analysis

- [ ] 将章节标题改为 `Spectral Analysis of Turbulence` 或类似更地道表达  
- [ ] 明确区分“观察结果”和“机制推断”  
- [ ] 对 Th/8 周期现象提出可检验的假说，而不是简单说尚无解释  
- [ ] 检查图 7-9 图注是否能独立读懂  

### 6. Results: Anisotropy analysis

- [ ] 统一 `anisotropy barycentric map` 等术语  
- [ ] 统一 `cigar-like` / `pancake-like` 写法  
- [ ] 在节首说明该分析如何服务主科学问题  
- [ ] 若使用该部分作为 LSC 证据，请说明证据等级和局限性  

### 7. Discussion and Conclusion

- [ ] 不再逐节复述结果  
- [ ] 明确区分机制性结论与配置依赖结论  
- [ ] 提出更具体的参数化启示  
- [ ] 增加局限性说明：真实风浪谱、复杂地形、时变风、非圆潮流等  

---

## 五、图表与排版检查

| 修改项 | 具体动作 | 状态 |
|---|---|---|
| 图 1 图注 | 重写为完整句，避免碎片化表达 | 未开始 |
| 表 1 | 优化标题、注释和大小写格式，使其符合期刊风格 | 未开始 |
| 图 2-3 | 在图注中补充 regime 含义或关键读图信息 | 未开始 |
| 图 4-6 | 统一 budget term 命名与顺序 | 未开始 |
| 图 7-9 | 解释参考频率线的物理含义 | 未开始 |
| 图 10-12 | 统一 marker 命名，避免 `xcross` 等非标准表述 | 未开始 |
| 全文图表 | 统一大小写、缩写首次定义、归一化说明 | 未开始 |

---

## 六、投稿规范与期刊要求

| 优先级 | 修改项 | 具体动作 | 状态 |
|---|---|---|---|
| 高 | 替换掉模板式 Open Research / Inclusion 文本 | 删除模板说明，改为作者真实内容；若不适用，按期刊要求放到正确位置 | 未开始 |
| 高 | 补全 Acknowledgments | 增加资助、合作、计算资源等说明 | 未开始 |
| 高 | 补全 Data Availability Statement | 给出数据、代码、脚本的获取链接或说明 | 未开始 |
| 中 | 核查参考文献格式 | 统一大小写、DOI、期刊名、断行、斜体等 | 未开始 |
| 中 | 核查图表引用顺序 | 确保正文中按出现顺序引用图表 | 未开始 |

---

## 七、如果只做“最小可投稿版本”，建议至少完成以下 10 项

- [ ] 引言末尾重写，明确知识缺口与创新点  
- [ ] 明确定义 overlap/merge 判据  
- [ ] 收紧 Langmuir supercells 的表述强度  
- [ ] 补充一个定量图支撑 Ri/shear/Pk 机制  
- [ ] 处理文中“无法解释”的现象  
- [ ] 重写摘要与 Key Points  
- [ ] 重写结论，增加参数化启示  
- [ ] 替换 Open Research 模板文本  
- [ ] 补全 Acknowledgments 与 Data Availability  
- [ ] 做一轮系统性英文润色与术语统一  

---

## 八、建议的改稿流程

1. **先改逻辑，不先改语言**  
2. **先补判据和证据，再补图**  
3. **结论收紧后，再统一摘要和 Key Points**  
4. **最后做全文语言润色、图注和格式清理**

