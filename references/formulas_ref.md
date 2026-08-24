# 新版 SPC 精确公式速查（技能自带 grounding）

> 本文件为技能发布时自包含的精简版公式速查，不依赖任何外部路径。

## 1. 计量值控制图——精确分布控制限

### 1.1 X̄-R 图（w 分布，α=0.0027）

基础量：
- $\bar{X}_i = \frac{1}{n}\sum x$，$\bar{\bar{X}} = \frac{1}{k}\sum \bar{X}_i$
- $R_i = x_{\max}-x_{\min}$，$\bar{R} = \frac{1}{k}\sum R_i$
- $\hat{\sigma} = \bar{R}/d_2$（$d_2$ 为常数，见末表）

精确控制限：
- $UCL_R = w_{n;1-\alpha/2}\cdot\hat{\sigma}$，$LCL_R = w_{n;\alpha/2}\cdot\hat{\sigma}$（$LCL_R\le0$ 取 0）
- $w_{n;p}$ = 标准化极差 $R/\sigma$ 的 $p$ 分位数 = `scipy.stats.studentized_range.ppf(p, n, np.inf)`
- $UCL_{\bar{X}} = \bar{\bar{X}} + A_2\bar{R}$，$LCL_{\bar{X}} = \bar{\bar{X}} - A_2\bar{R}$，$A_2 = 3/(d_2\sqrt{n})$
- （手册精确式含 $c_n$ 修正：$UCL_{\bar{X}} = \bar{\bar{X}} + u_{1-\alpha/2}\frac{1}{\sqrt{n}}\frac{c_n}{d_2}\bar{R}$；引擎当前用 A2 等价，留作精化）

### 1.2 X̄-s 图（χ² 分布）

- $\hat{\sigma} = \bar{s}/c_4$
- $UCL_s = \hat{\sigma}\sqrt{\chi^2_{\nu;1-\alpha/2}/\nu}$，$LCL_s = \hat{\sigma}\sqrt{\chi^2_{\nu;\alpha/2}/\nu}$

### 1.3 I-MR 图（n=1）

- $MR_i=|x_i-x_{i-1}|$，$\overline{MR}=\frac{1}{k-1}\sum_{i=2}^k MR_i$
- $\hat{\sigma}=\overline{MR}/d_2(n=2)$
- $UCL_X=\bar{x}+E_2\overline{MR}$，$UCL_{MR}=D_4\overline{MR}$（移动极差点相关，仅越限为有效信号）

## 2. 计数值控制图——精确分位数限

### 2.1 p / np 图（二项分布）
- $\hat{p}=\sum x_i/\sum n_i$
- 精确限：$G(x_{lo};n_i,\hat{p})\le\alpha/2$，$G(x_{up};n_i,\hat{p})\ge1-\alpha/2$ → $LCL_i=x_{lo}/n_i$，$UCL_i=x_{up}/n_i$
- 实现：`scipy.stats.binom.ppf` / 线性搜索满足 CDF 边界

### 2.2 u / c 图（泊松分布）
- $\bar{u}=\sum x_i/\sum n_i$（u）；$\hat{\mu}=\sum x_i/k$（c）
- 精确限：$F(x_{lo};\bar{u}\cdot n_i)\le\alpha/2$，$F(x_{up};\bar{u}\cdot n_i)\ge1-\alpha/2$
- 实现：`scipy.stats.poisson.ppf`

## 3. 能力指数四族

### 3.1 正态近似（用 $s$ 总体标准差）
- $Cp=(U-L)/(6s)$，$Cpk=\min[(U-\bar{x})/(3s),(\bar{x}-L)/(3s)]$
- $Pp=(U-L)/(6s)$，$Ppk=\min[(U-\bar{x})/(3s),(\bar{x}-L)/(3s)]$（overall）
- 组内：$\hat{\sigma}_w=\bar{R}/d_2$ 或 $\bar{s}/c_4$；$C_w=(U-L)/(6\hat{\sigma}_w)$，$C_{wk}=\min[\cdot]$

### 3.2 通用几何法（分位数，不假设正态，ISO 22514）
- $C_{pG}=(U-L)/(X_{99.865\%}-X_{0.135\%})$
- $C_{pU.G}=(U-X_{50\%})/(X_{99.865\%}-X_{50\%})$，$C_{pL.G}=(X_{50\%}-L)/(X_{50\%}-X_{0.135\%})$
- $C_{pk.G}=\min(C_{pU.G},C_{pL.G})$

### 3.3 z-分数法 / 超标比例法（非正态）
- $p_L=\int_{-\infty}^L f(x)dx$，$p_U=\int_U^\infty f(x)dx$
- $z_L=\Phi^{-1}(p_L)$，$z_U=\Phi^{-1}(1-p_U)$
- $P_{pk.z}=\min(-z_L,z_U)/3$，$P_{p.z}=(-z_L+z_U)/6$

### 3.4 机器性能 Pm/Pmk（第 7 章）
- $P_m=(U-L)/(6s_m)$，$P_{mk}=\min[(U-\bar{x})/(3s_m),(\bar{x}-L)/(3s_m)]$
- 目标值（N≥50）：关键 Pm=2.33/Pmk=2.00；主要 2.00/1.67；次要 1.67/1.33；其他 1.00/1.00
- 减样本量修正（χ² 分位数，见知识库 §2.4）

## 4. 时间相关分布模型（ISO 22514-2）

| 模型 | 位置 | 变差 | 含义 |
|------|------|------|------|
| A | 恒定 | 恒定 | 经典稳定（A1 正态 / A2 单峰非正态） |
| B | 恒定 | 变化 | 散布随时间变 |
| C | 变化 | 恒定 | 位置变（C1 随机/C2 趋势/C3 批间/C4 系统+随机） |
| D | 变化 | 变化 | 最复杂/多峰 |

## 5. 稳态硬规则

Nelson 1-8 任意触发 → 过程不稳态 → Cp/Cpk/Cpm/Cpmk 输出 NA（仅 Pp/Ppk 可报）。

## 6. d₂ 常数表（n=2..25）

| n | d₂ | n | d₂ | n | d₂ |
|---|----|---|----|---|----|
| 2 | 1.128 | 10 | 3.078 | 18 | 3.640 |
| 3 | 1.693 | 11 | 3.173 | 19 | 3.689 |
| 4 | 2.059 | 12 | 3.258 | 20 | 3.735 |
| 5 | 2.326 | 13 | 3.336 | 21 | 3.778 |
| 6 | 2.534 | 14 | 3.407 | 22 | 3.819 |
| 7 | 2.704 | 15 | 3.472 | 23 | 3.858 |
| 8 | 2.847 | 16 | 3.532 | 24 | 3.895 |
| 9 | 2.970 | 17 | 3.588 | 25 | 3.931 |

## 7. 页码索引（知识库）

- 能力指数四族 / 通用几何法：知识库 §1（p.108-116）
- 机器性能 Pm/Pmk：知识库 §2（p.125-145）
- 时间模型 A/B/C/D：知识库 §3（p.对应 ISO 22514-2）
- 控制图精确限：知识库 §4（p.238 w 分布；X̄-s χ²；p/c 二项/泊松）
- 旧版基线对照：知识库 §5
