# 论文详细修改意见（审稿式）

## 基本信息

- 论文题目：**Diurnal Modulation of Coastal Dual-Boundary Layer Turbulence Under Periodic Forcing**
- 目标期刊：**JGR: Oceans（推定）**
- 评阅定位：**海洋湍流/边界层过程方向的审稿式修改建议**
- 结论性建议：**建议“大修（major revision）”后再投**

---

## 一、总体评价

这篇论文围绕**浅海表层边界层（SBL/OSBL）与底边界层（BBL）在日变化强迫下的耦合、重叠与全水柱湍流爆发**展开，主题具有明显的学术价值，也契合当前沿海湍流、Langmuir turbulence、boundary-layer overlap 及参数化研究的热点。文章的主要亮点包括：

1. **科学问题重要**：讨论了热通量、风浪、潮流共同作用下双边界层日变化调制问题；
2. **结果具有结构性**：提出了三类边界层演变型态（fully separated / once overlap / twice overlap）；
3. **诊断维度较完整**：结合了 TKE、Ri、TKE budget、wavelet spectrum、anisotropy barycentric map；
4. **潜在应用清晰**：对近岸垂向混合、悬沙再悬浮、物质输运和参数化方案有启发意义。

但从达到 **JGR: Oceans** 论文发表标准的角度看，当前稿件仍存在几类关键问题：

- **创新点与已有工作的边界尚不够尖锐**；
- **“边界层重叠/合并”的判据仍偏定性**；
- **对 Langmuir supercells 的识别证据尚不充分，部分表述偏超前**；
- **若干核心现象作者自己也承认“尚无合理解释”，削弱了说服力**；
- **语言层面存在系统性的语法、搭配、冠词、时态、主谓一致和术语表达问题**；
- **期刊规范性要素仍不完整**，尤其是开放研究声明、致谢、数据/代码可用性等。

因此，本稿**具备发表潜力，但需要系统性重写和实质性补强**，而不只是润色英语。

---

## 二、总体修改建议：按优先级排序

### A. 必须优先解决的重大问题（决定稿件能否成立）

#### 1. 创新点与知识缺口（knowledge gap）界定还不够锋利

**涉及位置**：引言第 1 节，尤其是 L92-L154 一带。

**问题**：

当前引言虽然回顾了 Langmuir supercells、overlapping boundary layers、stratification 和 misalignment 的已有研究，但“本研究究竟比 Yan et al. (2022)、Valcarcel et al. (2025)、Shrestha et al. (2019) 新在哪里”并没有被非常清晰地钉住。现有写法容易给审稿人一种印象：  
你们是在把“已有因素”放到同一个 LES 框架里再做一次组合实验，但**核心新意**尚未被一句话概括清楚。

**建议修改方向**：

将知识缺口明确压缩为以下三点中的一到两点主线，不要平均用力：

1. **不是简单研究 overlap 是否发生，而是研究 overlap 的日周期位相锁定机制（phase locking）**；
2. **不是简单研究 stratification 或 misalignment 的单因素效应，而是研究二者在同一日周期中的交替主控关系**；
3. **不是简单分类三种型态，而是说明“什么时候 thermal control 主导、什么时候 directional control 主导”。**

**建议强化句式（可直接改写引言结尾）**：

> Existing studies have examined boundary-layer overlap, Langmuir turbulence, and directional misalignment largely in isolation. What remains unclear is how diurnally varying stratification and rotating tidal-current alignment jointly determine the timing, duration, and recurrence of overlap events within a repeating daily cycle.

**修改理由**：  
JGR 审稿人首先会问“新意是否足够集中”。如果创新点不够聚焦，后面的结果再丰富，也容易被评价为“interesting but descriptive”。

---

#### 2. “boundary-layer overlap / merge” 缺少可复现、可量化的判据

**涉及位置**：结果第 3 节，特别是 L375-L381、L422-L430 等。

**问题**：

目前三类边界层型态的分类主要基于图 2 的 TKE 分布进行视觉判断，但文中没有给出一个明确的、可复现的定义，例如：

- 以 TKE 阈值定义上下边界层是否连通？
- 以混合层深度 + 底边界层厚度之和是否超过水深定义？
- 以中层 Ri、shear production 或湍流通量阈值定义？
- “merge” 与 “overlap” 是否完全同义？

如果判据不明，审稿人很可能会质疑三类 regime 是**定性描述**而非严格分类。

**建议修改方向**：

建议在方法或结果开头补充一个**操作性定义（operational definition）**。至少应说明以下内容：

1. 你们如何定义 surface boundary layer depth；
2. 你们如何定义 bottom boundary layer depth；
3. 何种条件下判定 overlap/merge event 已发生；
4. 该判据是否对所有图中的结论稳健。

**建议可采用的判据方案（给出备选）**：

- **方案 A：TKE 连通阈值法**  
  若在同一相位内，从表层到近底存在一条连续深度带满足  
  `k / u_*^2 > k_crit`，则定义为 full-column overlap。

- **方案 B：边界层厚度和法**  
  若 `h_sbl + h_bbl >= H`，则定义为 overlap。

- **方案 C：动力学连通法**  
  若中层 `Pk + Bk + Tk` 中至少一个净正项持续超过阈值，并伴随 TKE 显著增长，则定义为 dynamical merging event。

**建议**：正文保留一个主判据，附录中给出敏感性检查。

**修改理由**：  
没有判据，分类就不够严谨；有判据，文章的“机制归纳”才能站住。

---

#### 3. 对 Langmuir supercells 的识别证据仍偏间接，建议降调或补证据

**涉及位置**：L490-L506、L819-L839、结论部分 L851-L858。

**问题**：

文中多处将 full-column turbulence 与 Langmuir supercells 紧密关联，甚至近乎将两者等同。但目前主要证据来自：

- TKE 全水柱增强；
- 低频窄带谱峰；
- 各向异性应力状态变化。

这些证据**支持存在大尺度、有组织的全水柱结构**，但若要更强地宣称其为 **Langmuir supercells**，通常还需要更直接的流场证据，例如：

- 横截面涡对结构；
- 沿风向条带/汇聚-辐散图样；
- streamwise roll coherence；
- spanwise spacing 与已知 LSC 标度的一致性。

**建议修改方向**：

两种方案二选一：

- **方案 A（更稳妥）**：将正文中的强判断改为  
  “full-column turbulent structures consistent with Langmuir supercells”  
  或  
  “Langmuir-supercell-like structures”。

- **方案 B（更强）**：补充一幅或两幅代表性流场快照/剖面图，直接展示 overlap 事件中的成对旋涡、汇聚带或卷轴结构。

**修改理由**：  
这是典型的“证据链强度与结论强度不匹配”问题。保守表述比过度命名更安全。

---

#### 4. 存在“作者自己也解释不了”的关键现象，必须处理

**涉及位置**：

- L579-L580：*no reasonable explanation can be proposed for the oscillations of TKE observed in the bottom boundary layer*  
- L663-L664：关于 Th/8 周期“there is no satisfactory explanation for it yet”

**问题**：

这类句子如果直接保留在主文结果段落中，会显著削弱稿件说服力。审稿人往往会据此认为：

- 关键机制尚未阐明；
- 某些解释可能仍是经验性描述；
- 结果是否稳健尚不清楚。

**建议修改方向**：

不要简单写“无法解释”。建议改为以下三层表达：

1. **先界定现象**：该振荡/周期确实存在且在多个 case 中重复出现；
2. **再给出最可能机制假说**：例如与潮流旋转、底边界层内 shear-production modulation、惯性/半日潮谐波耦合有关；
3. **最后将其作为 future work / limitation**，而不是直接宣告无法解释。

**建议表述模板**：

> The recurring Th/8 signal suggests an internally generated modulation within the bottom boundary layer, possibly associated with nonlinear interactions among tidal rotation, near-bottom shear adjustment, and the phase-dependent stress response. A dedicated diagnosis is beyond the scope of the present study, but the robustness of this signal across cases indicates that it is unlikely to be numerical noise.

**修改理由**：  
可以承认尚未完全解决，但不能把核心现象以“无解释”结束。

---

#### 5. 理想化强迫设置需要更充分地讨论适用范围与局限性

**涉及位置**：方法第 2 节，特别是 L204-L259、L317-L323、L324-L340。

**问题**：

本研究采用了高度理想化的配置：

- 恒定风应力；
- 恒定单色波；
- 24 h 热通量循环；
- 12 h 圆形潮流；
- 热通量与潮流周期被人为设置为简单可比较关系；
- 无复杂地形、无水平非均匀性、无盐度变化、无波轨道速度与底边界相互作用。

这些简化对于“机制实验”是合理的，但当前稿件对**为什么必须这样简化、这些简化不会影响哪些结论、又会限制哪些外推**，讨论还不够。

**建议修改方向**：

在 2.3 或 Discussion 中增加一个专门小段，明确区分：

- **机制层面的普适结论**：例如 stratification 与 alignment 的双重控制；
- **依赖理想化配置的定量结果**：例如 once/twice overlap 的具体分界、NP/PP 差异幅度；
- **仍需真实海况验证的结论**：例如在 broadband waves、非圆潮流、时变风下是否仍保留三类 regime。

**修改理由**：  
越理想化，越需要主动说明“哪些能外推，哪些不能外推”。

---

#### 6. Ri 与 TKE 之间的机制解释需要更多定量支撑，而不仅是叙述性说明

**涉及位置**：L437-L503。

**问题**：

你们已经提出一个很有价值的观点：  
**中层是否出现 overlap，不仅由 stratification 决定，也受 shear minimum 控制。**

但当前写法主要是文字解释，缺少更直接的定量支撑，例如：

- 代表性相位的 mean shear profile；
- 中层 `Pk` 的相位变化曲线；
- `Rig` 与 `Pk` / `k` 的联合相图；
- 关键深度处的时间序列比较。

尤其在 Case 4(PP) 中，你们强调 `Ri_g < 0.25` 仍不足以触发强湍流，这实际上是全文的关键科学判断之一，值得配更强的证据。

**建议修改方向**：

至少补充一种额外定量展示：

1. 一个代表深度处的 `Ri_g`, `Pk`, `k` 三联时间序列；
2. 两个相位的剖面比较（aligned vs oppositely directed）；
3. 一个散点图显示中层 turbulence onset 与 shear-production threshold 的关系。

**修改理由**：  
这会让“alignment regulates local energy supply” 从定性论断升级为定量结论。

---

#### 7. 讨论与结论部分仍偏“结果重述”，需要提升机制整合与外推价值

**涉及位置**：第 7 节，L840-L896。

**问题**：

结论部分目前更像结果摘要，而不是高水平期刊常见的“机制归纳 + 适用范围 + 对参数化的启发”。  
尤其以下三点可进一步增强：

- 与已有研究相比，你们推进了什么；
- 结果最可能如何进入 parameterization 思路；
- 哪些结论可能依赖北半球、圆形潮流、固定 La_t、固定深度 H=45 m。

**建议修改方向**：

建议把结论压缩成 4 个层次：

1. **最核心发现**：三种 regime + 触发条件；
2. **机制归纳**：strong stratification control -> mixed control -> directional control；
3. **参数化启示**：仅靠 bulk stratification 不能判断 overlap，需要引入 phase-dependent alignment metric；
4. **适用范围与局限性**：哪些结论是机制性的，哪些需要更真实海况验证。

---

### B. 需要重点修改的重要问题（影响论文质量和可读性）

#### 8. 文章结构较长，结果段落之间存在重复解释

**问题表现**：

- 第 3 节、第 4 节、第 5 节、第 6 节之间存在一定重复：
  - “三类 regime”被多次重新解释；
  - “stratification + alignment”的二元控制逻辑反复出现；
  - 个别句子更像结论，提前出现在结果段中。

**建议**：

- 第 3 节只回答“现象是什么”；
- 第 4 节只回答“能量学机制是什么”；
- 第 5 节只回答“结构尺度如何变化”；
- 第 6 节只回答“各向异性如何证明结构性质变化”；
- 将总括性判断尽量收束到第 7 节。

---

#### 9. 个别图注和正文的衔接不够严谨

**例子**：

- L330-L335 一带图 1 图注语法不完整；
- 表 1 标题和注释排版略显仓促；
- 图 7-9 的“white line / red line / yellow line”说明需要更明确说明物理意义；
- 图 10-12 中“xcross”这种表述不规范，建议统一为 “cross markers” 或 “x-shaped markers”。

**建议**：

- 图注不要只描述颜色和线型，更要说明“读图结论是什么”；
- 图题首句尽量写清：数据类型、归一化方式、坐标含义、代表的 case。

---

#### 10. 术语需要统一

**建议统一项**：

- `surface boundary layer (SBL)` 与 `ocean surface boundary layer (OSBL)` 二选一，并全文统一；
- `merge` / `overlap` / `interaction` / `coupling` 含义应区分；
- `full-column turbulence` 与 `full-depth turbulence` 选择其一为主；
- `positive polarity tide / negative polarity tide` 建议固定缩写首次定义后统一使用；
- `wind-wave-current alignment` 与 `alignment of wind, waves, and currents` 二者择一。

**修改理由**：  
术语不统一会让稿件显得尚未定稿，尤其在多诊断量并存时更明显。

---

## 三、分节详细修改意见

### 1. Title / Key Points / Abstract / Plain Language Summary

#### 1.1 标题

当前标题总体可用，但略显平直。建议增强“机制性”和“结果性”。

**可选标题 1（较保守）**  
**Diurnal Modulation of Turbulence in Coastal Dual Boundary Layers under Periodic Forcing**

**可选标题 2（更突出机制）**  
**Diurnal Controls on Boundary-Layer Overlap and Full-Column Turbulence in Shallow Coastal Oceans**

**可选标题 3（更突出核心发现）**  
**Stratification and Flow Alignment Jointly Control Diurnal Overlap of Coastal Dual Boundary Layers**

---

#### 1.2 Key Points

当前 Key Points 的语言存在明显语法问题，例如：

- “once-overlap” / “twice-overlap” 不自然；
- “Boundary layers merge is controlled by ...” 主谓不一致；
- “wind-waves” 是否连字符写法需统一。

**建议改写版本**：

1. **Diurnal forcing organizes coastal dual-boundary layers into three regimes: fully separated, overlapping once per day, and overlapping twice per day.**
2. **Full-column turbulence develops when the surface and bottom boundary layers merge and decays after they separate.**
3. **Boundary-layer overlap is jointly controlled by water-column stratification and the alignment among tides, wind, and waves.**

---

#### 1.3 Abstract

**优点**：问题、方法、结果、结论基本齐全。  
**主要问题**：

- 句子偏长；
- 局部搭配不地道；
- “reveal three distinct regimes”之后，机制解释和结论强度略超前；
- “Negative polarity tides ... are more likely to induce boundary layer overlap” 需要说明是在当前参数设置下。

**建议**：

1. 首句可再压缩，避免“critical role... with significant implications...”这类套话叠加；
2. 明确研究对象是 idealized LES；
3. 结尾建议区分“our simulations show”与“we infer”；
4. 对“full-column turbulence fundamentally altering the stress state of coastal oceans”这类较强表达略微降调。

**建议性改写框架**：

> Turbulent mixing in shallow coastal oceans is often confined to the surface and bottom boundary layers, but under favorable conditions these layers can interact and even merge. The mechanisms governing such interactions under diurnal forcing remain insufficiently understood. Here we use idealized large-eddy simulations to examine how diurnally varying surface heat flux, steady wind-wave forcing, and rotating semidiurnal tides jointly modulate turbulence in coastal dual boundary layers. The simulations identify three recurring regimes: fully separated boundary layers, boundary layers overlapping once per day, and boundary layers overlapping twice per day. The regime transitions are controlled by the combined effects of water-column stratification and the phase-dependent alignment among wind, waves, and tidal currents. Even after stratification weakens, strong current-wave misalignment can suppress overlap. Under otherwise identical forcing, negative-polarity tides in the Northern Hemisphere favor thicker bottom boundary layers and more frequent overlap. Full-column turbulence emerges during overlap events and is accompanied by substantial changes in turbulence structure and stress anisotropy.

---

#### 1.4 Plain Language Summary

**问题**：

- “have chances to interact” 过于口语；
- “These powerful turbulence events” 语法不对，应为 “these intense turbulent events” 或 “such energetic turbulence events”；
- 需要更强调“为什么公众或跨学科读者应该关心”。

**建议**：

Plain Language Summary 应少用术语堆砌，多强调：

- 近岸海洋为什么会突然全水柱混合；
- 这对泥沙、营养盐、缺氧、水质和生态有何影响；
- 你们用数值模拟发现了“什么时候会发生、什么时候不会发生”。

---

### 2. Introduction

#### 2.1 文献综述较全面，但“本文的独特贡献”需要单独成段

建议在引言末尾增加一个非常明确的段落，例如：

> In contrast to previous studies that focused on either stratification, directional misalignment, or overlap events separately, this study investigates how these controls interact within a repeating diurnal cycle. We further distinguish whether overlap is suppressed by residual stratification or by phase-dependent reductions in local shear production.

#### 2.2 若要强调“poorly understood”，需更精确

建议避免笼统说“remain poorly understood”，而改为：

- the recurrence of overlap events within a diurnal cycle remains unclear
- the relative roles of stratification and phase-dependent alignment remain unresolved
- it remains unclear when overlap leads to full-column turbulence rather than weak interior mixing

这样更符合审稿人对“精确缺口陈述”的预期。

---

### 3. Method

#### 3.1 第 2.1 节存在若干英语句法问题

例如 L158-L167 一段中：

- “Our simulations using Oceananigans ...” 不是完整句；
- “a modern computational fluid dynamics model distinguished by...” 结构有断裂；
- 方程说明中的标点、空格、变量定义需要统一。

**建议改写**：

> We perform the simulations using Oceananigans (Ramadhan et al., 2020), a modern computational fluid dynamics framework with a user-friendly interface and efficient GPU acceleration. Oceananigans has been widely applied to studies of turbulence and mixing in oceanic boundary layers.

#### 3.2 方程与符号说明需要更规范

建议逐项检查：

- 变量首次出现是否定义完整；
- 是否使用一致的斜体/正体；
- 单位前是否留空格；
- `g = -9.81 m/s` 建议检查符号与量纲表达，通常应写为加速度大小 `g = 9.81 m s^-2`，浮力方向由方程符号体现；
- `S = 35 psu` 若期刊风格允许，最好考虑使用更规范的盐度单位说明。

#### 3.3 热通量和潮流设置的解释应更服务于“机制问题”

建议强调：

- 为什么采用 12 h heating / 12 h cooling 的理想化循环；
- 为什么选择 circular tide；
- 为什么 Umag = Umin 有利于隔离方向效应；
- 为什么把风浪固定为 x 方向有助于研究 alignment 而不是 forcing magnitude variability。

#### 3.4 自旋时间与统计稳态定义需要更量化

L337-L340 当前写得偏口语化。建议说明：

- spin-up 持续多少个周期；
- 统计稳态的判据是什么；
- 最后六天数据是否足以收敛；
- 是否对不同 case 一致。

---

### 4. Results: Boundary-layer evolution

#### 4.1 三类 regime 的命名建议规范化

当前写法：

- fully separated
- overlap once a day
- overlap twice a day

建议在文中首次定义后，统一写成：

- **fully separated regime**
- **once-daily overlap regime**
- **twice-daily overlap regime**

这会比正文中反复写完整描述更利落。

#### 4.2 关于 NP/PP 差异的讨论需要避免过度泛化

文中指出 negative polarity tides 在北半球更易导致 overlap，这一结果很有意义，但建议限定条件：

- 在本文参数设定下；
- 在相同风浪、热通量、深度与旋转参数下；
- 该结论本质上源于底边界层发展时间尺度差异。

否则容易被理解为“普适海洋结论”。

#### 4.3 建议增加一张归纳图或 regime diagram

如果篇幅允许，建议增加一个总结图，横轴表示净热损失或 stratification strength，纵轴表示 alignment favorability / polarity，标出三类 regime 所在区域。  
即便只是概念图，也会显著增强文章结构感。

---

### 5. Results: Energy budget

#### 5.1 TKE budget 分析总体扎实，但表达可更简洁

当前问题：

- 段落偏长；
- 每个 budget term 都解释得过细，容易淹没主结论；
- “pressure transport feeds internal waves” 等表述可以更谨慎。

**建议结构**：

每个 case 用统一模板：

1. dominant balance near boundaries；
2. what changes in the mid-layer；
3. which term distinguishes merged from separated states；
4. one takeaway sentence。

#### 5.2 “Tk is the primary role” 等表达需修改

这类中式英语较多，建议统一替换为：

- `Tk plays the primary role in...`
- `Tk provides the dominant positive contribution to...`
- `Tk bridges the gap between local production and dissipation in the mid-layer`

---

### 6. Results: Spectral analysis

#### 6.1 章节标题建议修改

当前标题：**Spectra analysis of turbulence**  
建议改为：

- **Spectral Analysis of Turbulence**
- 或 **Wavelet-Based Spectral Analysis of Turbulence**

#### 6.2 对频谱结果的解释需注意“观察”与“推断”分开

例如：

- “The reason for this frequency shift may be attributed to...” 可以，但后面要说明这是 interpretation；
- 对 internal waves 的角色应避免下结论过强；
- 对 convective supercells 的比较应点到为止，避免引入新的未充分展开的话题。

---

### 7. Results: Anisotropy barycentric map

#### 7.1 术语与数学表述应再规范

例如：

- `anisotropic barycentric map` 更常见表达是 **anisotropy barycentric map** 或 **barycentric map of turbulence anisotropy**；
- `principal invariants of the anisotropy stress tensor` 建议改为 `principal invariants of the anisotropy tensor`；
- `cigarlike` / `pancakelike` 建议统一写为 `cigar-like` / `pancake-like`。

#### 7.2 这一节很有特色，但需要更明确地服务主结论

建议在本节开头增加一句：

> This analysis is used to diagnose whether overlap events merely intensify turbulence or fundamentally reorganize its structure across the full water column.

这样读者更容易理解为什么需要 barycentric map，而不是把它看作附加展示。

---

### 8. Discussion and Conclusion

#### 8.1 建议重写为“机制综合”而不是“逐节回顾”

推荐结论结构：

1. 三类 regime；
2. 决定性机制：stratification vs alignment；
3. NP/PP 差异的动力学解释；
4. full-column turbulence 的条件性发生；
5. 对参数化和沿海过程研究的意义；
6. 限制与未来工作。

#### 8.2 参数化启示应再具体

当前结论说“parameterizations must account for diurnal variability”，这还偏泛。  
可进一步指出：

- 仅用 bulk mixed-layer depth 可能不足以判断 overlap；
- 需要考虑 phase-dependent alignment metric；
- BBL/SBL 连通性可能需要事件式参数化而非连续参数化。

---

## 四、语言与行文问题：代表性例子

以下不是全文逐句校改，而是**最典型的一批共性问题**。建议作者据此对全文统一修订。

| 位置 | 原表达 | 建议表达 | 修改理由 |
|---|---|---|---|
| Key Points | `once-overlap`, `twice-overlap` | `overlapping once per day`, `overlapping twice per day` | 英语搭配不自然 |
| L17-L18 | `Boundary layers merge is controlled by...` | `Boundary-layer merging is controlled by...` | 主谓一致 + 名词化更自然 |
| L29-L31 | `wave driven` / `current driven` | `wave-driven` / `current-driven` | 复合形容词需加连字符 |
| L33-L34 | `the large eddy simulation is employed` | `large-eddy simulations are employed` | LES 通常复数表述更自然 |
| L36-L37 | `once overlap a day` | `overlap once per day` | 介词搭配错误 |
| L41-L43 | `are more likely to induce` | `more readily induce` / `are more likely to favor` | 学术表达更自然 |
| L46-L47 | `have chances to interact` | `can interact` / `may interact` | 避免口语化 |
| L56-L58 | `These powerful turbulence events` | `These intense turbulent events` | 形容词与名词搭配错误 |
| L92-L93 | `lead to significant changes in the turbulence` | `substantially modify turbulence intensity and structure` | 更具体 |
| L103-L104 | `depends not only on water depth but also on a subtle matching` | `depends not only on water depth but also on a favorable combination` | `subtle matching` 不自然 |
| L116 | `vertical transport(e.g.,` | `vertical transport (e.g.,` | 括号前缺空格 |
| L164-L166 | `Our simulations using...` | `We perform the simulations using...` | 原句不完整 |
| L191-L193 | `here, g = -9.81m/s` | `Here, g = 9.81 m s^-2` | 量纲和格式规范 |
| L205-L206 | `As we mentioned in section 1` | `As discussed in Section 1` | 更正式 |
| L240-L241 | `the driving force is therefore being set` | `the driving force is therefore prescribed` | 被动式更简洁 |
| L256 | `using the polarity` | `using tidal polarity` | 术语更清晰 |
| L264 | `assumption of horizontally homogeneous` | `assumption of horizontal homogeneity` | 名词形式正确 |
| L309-L315 | `12hour` | `12-hour` | 连字符与格式 |
| L338-L339 | `which is the statistical characteristics of turbulence...` | 重写整句 | 原句语法不成立 |
| L377-L380 | `the pattern can be classified into three types` | `the simulated diurnal evolution can be classified into three regimes` | 更准确 |
| L401-L402 | `transit from` | `transition from` | 动词搭配 |
| L428-L429 | `persistent cooling rules out thermal forcing as the only determinant` | 可保留，但建议重写 | 逻辑较绕 |
| L458-L460 | `Fig. 2d) The difference` | `Fig. 2d). This mismatch between ...` | 标点和衔接 |
| L471-L472 | `depends on the magnitude of the shear is also important` | `also depends on the magnitude of the shear` | 语法错误 |
| L497-L499 | `the vertical temperature gradient is always satisfying` | `the vertical temperature gradient remains` | 避免中式进行体 |
| L509-L510 | `was used here` | `is used here` 或 `is used to diagnose` | 时态统一 |
| L588-L589 | `Tk is the primary role` | `Tk plays the primary role` | 固定搭配 |
| L640-L641 | `non-stationary processes We still focus` | 句号缺失，需断句 | 标点错误 |
| L693-L694 | `The cooling is much intense` | `Cooling is much stronger` | 比较级错误 |
| L707 | `frequency` lines in figure caption | 建议补“reference frequencies for interpretation” | 图注信息需更完整 |
| L744-L748 | 引号和标签格式 | 统一使用英文半角引号或不加引号 | 排版统一 |
| L759-L762 | `shows the anisotropy by scalar` | `represents anisotropy using scalar measures` | 更准确 |
| L887-L889 | `better assess the robustness` | 可保留，但建议补“of the identified regimes” | 指向更明确 |
| 全文 | `Langmuir turbulence` 有时小写/大写不统一 | 全文统一 | 术语规范 |

---

## 五、学术规范与期刊格式问题

### 1. Open Research / Inclusion Statement 目前显然还是模板文本

**涉及位置**：Appendix A, L897-L926。

**问题**：

这部分仍是期刊模板提示语，尚未替换成作者自己的内容。若直接投稿，极可能被编辑部退回。

**建议**：

- 若需要该声明，请按期刊要求填写真实内容；
- 若本稿不适用 Appendix A 形式，请按期刊模板改放至正确位置；
- 不要保留编辑说明和模板指导语。

---

### 2. Acknowledgments 与 Data/Code Availability 需补全

提取文本中只看到 **Acknowledgments** 和 **References** 标题，未见完整内容。  
JGR/AGU 通常非常重视：

- 资助信息；
- 数据可获取性；
- 代码或模拟配置可获取性；
- 图表复现说明。

**建议**：

至少准备以下内容：

1. funding sources；
2. model version（Oceananigans version）；
3. simulation scripts / post-processing scripts 是否公开；
4. output data 是否上传 Zenodo、figshare、institutional repository 或 GitHub release；
5. 数据可用性声明与链接。

---

### 3. 参考文献格式建议统一全面复核

当前参考文献整体较完整，但从提取文本看，存在以下潜在问题：

- 专有名词大小写不统一，如 `langmuir turbulence`；
- 个别期刊名/标题大小写可能不完全符合 AGU 风格；
- DOI 前后空格、大小写与断行需要统一；
- 个别文献看起来来自非正式抓取源，需核对最终 BibTeX 输出质量。

---

## 六、我建议作者优先新增或补强的内容

如果作者准备进行实质性大修，最值得新增的内容有：

1. **一个明确的 overlap 判据**；
2. **一张 regime summary / conceptual schematic**；
3. **至少一组 aligned vs misaligned 的 mean shear / Pk / Ri 对比图**；
4. **若坚持使用“Langmuir supercells”术语，则补一张直接流场结构图**；
5. **一个专门讨论理想化设置适用范围的小节**；
6. **将无法解释的现象改写成“稳健现象 + 假说 + 局限”三段式表述**。

---

## 七、推荐修改顺序

为提高效率，建议按以下顺序改稿：

1. **先重写引言最后 3 段和结论**，把创新点与主结论钉牢；
2. **补定义**：overlap / merge / regime；
3. **补机制证据**：Ri、Pk、shear 的量化支撑；
4. **决定是否保留“Langmuir supercells”强表述**；
5. **补规范性内容**：Open Research、Acknowledgments、Data/Code Availability；
6. **最后再做全文英语统一润色**。

---

## 八、给作者的最终判断

### 优势

- 选题好；
- 结果有层次；
- 模拟设计具有机制研究价值；
- 多种诊断方法互相支撑，具备发展成高质量论文的基础。

### 主要风险

- 分类定义不够定量；
- 核心机制部分仍偏叙述化；
- 术语和结论强度存在不匹配；
- 语言问题较系统，不是简单润色可完全解决；
- 期刊模板信息尚未完成。

### 最终建议

**建议作为一篇有潜力的稿件进行“大修”，重点不是“修英语”，而是“强化机制判据 + 收紧结论边界 + 完善期刊规范”。**

