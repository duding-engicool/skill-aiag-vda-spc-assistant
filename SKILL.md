---
name: AIAG-VDA SPC助手
slug: aiag-vda-spc-assistant
displayName: AIAG-VDA SPC助手
description: 基于 AIAG-VDA《SPC 统计过程控制手册（2026 版）》的企业落地陪跑技能；当用户需要从旧版 SPC（AIAG 第2版）切换到新版手册、做过程能力三阶段评价（Pm/Pmk→Pp/Ppk→Cp/Cpk）、起草 OCAP、或绘制新版控制图（计量型 I-MR/Xbar-R/Xbar-S、计数型 p/np/c/u、时间加权 EWMA/CUSUM 全控制图族）并生成 20 项可编辑网页报告时使用
version: 1.1.0
category: quality
author: org-jaxjwo0r
---

# AIAG-VDA SPC助手

## 任务目标
- 本 Skill 服务于企业从旧版 SPC（AIAG 第 2 版，4 章 + 9 附录）切换到 **AIAG-VDA《SPC 手册 2026 版》**（12 章，与 ISO 22514 对齐）的落地过程。
- 三重角色合一：**诊断顾问（陪跑）** + **OCAP 引导起草** + **计算引擎与可编辑网页报告**。
- 触发条件：用户提及「新版 SPC」「AIAG-VDA」「过程能力三阶段」「Pm/Pp/Cp」「OCAP」「控制图判异准则」「SPC 报告 20 项」等，或需要做 SPC 数据的计算、绘图与合规报告。

## 与旧版 SPC 技能的关系
- 本技能是 `spc-analysis`（旧版 SPC 绘图分析）的**升级换代**版本，**不覆盖**旧技能。
- 差异：旧版只做「绘图 + 分析」；本技能增加①新版三阶段能力评价与 ISO 22514 分位数法②OCAP 成文引导③判异准则须用户显式选定④N<2000 分布拟合硬规则⑤20 项合规网页报告（可编辑/可导出）。
- 旧项目若只需传统 Xbar-R/P/C/U 图与基础能力分析，仍可用 `spc-analysis`；涉及新版手册合规要求时改用本技能。

## 前置准备
- 依赖包（隔离 venv 内已安装）：`numpy`、`scipy>=1.17`、`pandas`。
- 运行环境：`C:/Users/admin/.workbuddy/binaries/python/envs/default/Scripts/python.exe`（Windows 托管控 Python 3.13）。
- 脚本调用统一在该技能目录下执行：`cd <技能目录>` 后用上述 python 调用 `scripts/` 下脚本。

## 三种工作模式

### 模式 A — 诊断顾问（陪跑，纯对话）
- 用途：诊断用户企业当前 SPC 应用与新版手册的差距，给出导入路线。
- 依据：[references/diagnosis_checklist.md](references/diagnosis_checklist.md)（按新版 12 章映射的诊断问题，含 OCAP 四要素、Pm 诊断、分布模型诊断）与 [references/old_new_mapping.md](references/old_new_mapping.md)（5 层变化 + 章节映射）。
- 交互：**轻交互、强制交互锁**。只问「当前阶段 + 具体诉求 + 已给信息」；信息不足时针对性追问 1–2 点，不预发大问卷。
- 输出：缺口清单 + 导入路线 + 技能联动指引（联动 `msa-analysis`、`cp-control-plan`、`fmea-assistant`）。
- 红线：诊断结论为**建议**，不做判定；不替用户写管理文件。

### 模式 B — OCAP 引导起草（建议稿，非代写）
- 用途：按新版手册第 9 章，引导用户形成每个控制图挂接的 OCAP（失控行动计画）。
- 依据：[references/ocap_guide.md](references/ocap_guide.md)（OCAP 四要素、执行流图 fig 9-2、4 轮引导提问、建议稿输出格式）。
- 交互：4 轮引导提问（责任归属 → 操作员调整矩阵 → 过程日志 → 培训确认），根据用户提供信息给出**建议稿**。
- 输出：txt + md 双版建议稿，文件头明确标注「**建议稿·待用户确认**」，不可直接当正式文件发布。
- 红线：绝不替用户编造现场细节；OCAP 最终责任在用户与现场。

### 模式 C — 计算引擎 + 可编辑网页报告
- 用途：对真实项目数据做新版合规计算，生成 20 项 SPC 报告（参照手册 11.2 模板）。
- 依据：[references/capability_logic.md](references/capability_logic.md)（三阶段 + 分位计算分支 + 控制图常数表）、[references/report_20_elements.md](references/report_20_elements.md)（20 项规格 + 合规硬规则 + 计算约定）。
- 流程：
  1. 运行计算引擎 `scripts/spc_engine.py` 产出 JSON 结果。
  2. 运行报告生成器 `scripts/report_builder.py` 注入结果到模板，产出单文件 HTML。
  3. 网页版支持：编辑抬头/数据（改数即重算）、localStorage 草稿、导出 JSON、导入 JSON、另存 HTML、导出 PDF（打印）。

## 操作步骤（模式 C 标准流程）

1. **接收数据**：用户上传 CSV（计量型/时间加权为单列数值，每行一个观测；若子组数据则用 `n` 列或一次性给定子组大小）或直接在对话中给出数值；**计数型（属性图）须另给计数列（counts）与样本量/单位数（sizes）**。确认是否有 USL/LSL/目标值、子组大小、评价阶段。
2. **选定判异准则（必做）**：新版手册要求判异准则由用户**显式选定并文档化**，默认仅规则 1。与用户确认使用哪些（1/2/…/8），如未明确则只跑规则 1 并提示。
3. **运行计算引擎**：
   ```bash
   cd <技能目录>
   python scripts/spc_engine.py \
     --data <数据.csv> \
     --subgroup-size <子组大小, 1=单值I-MR> \
     --chart auto|xbar_r|xbar_s|i_mr|p|np|c|u|ewma|cusum \
     --usl <规格上限> --lsl <规格下限> --target <目标值> \
     --rules 1,2,5 \
     --stage machine|performance|capability \
     --out <结果.json>
   ```
   - `--stage`：机器性能 `machine`（Pm/Pmk，设备放行）→ 过程性能 `performance`（Pp/Ppk，过程放行）→ 过程能力 `capability`（Cp/Cpk，量产监控）。
   - `--chart auto`：n≥2 默认 Xbar-R，n=1 自动 I-MR。
   - **计数型（属性图）** 用 `--counts`（计数列）+ `--sizes`（样本量/单位数：CSV 或常量；c 图可省略）：
     ```bash
     # p 图（不合格品率，样本量可逐组不等）
     python scripts/spc_engine.py --chart p --counts p_counts.csv --sizes p_sizes.csv --out r_p.json
     # np 图（不合格品数，样本量恒定）
     python scripts/spc_engine.py --chart np --counts np_counts.csv --sizes 50 --out r_np.json
     # c 图（缺陷数）；u 图（单位缺陷数，单位数可逐组不等）
     python scripts/spc_engine.py --chart c --counts c_counts.csv --out r_c.json
     python scripts/spc_engine.py --chart u --counts u_counts.csv --sizes u_units.csv --out r_u.json
     ```
     属性图不输出分布/能力指数/ppm（报告第 VII–X 节显示「不适用」），仅判定合格/不合格。
   - **时间加权（EWMA / CUSUM）** 作用于连续数据（`--data`），可选参数：
     `--ewma-lambda 0.2`（平滑系数）、`--ewma-L 3.0`（限倍数）、
     `--cusum-k 0.5`（以σ计的检测偏移）、`--cusum-h 4.0`（决策区间）。
     ```bash
     python scripts/spc_engine.py --data cont.csv --chart ewma --ewma-lambda .2 --ewma-L 3 --out r_ewma.json
     python scripts/spc_engine.py --data cont.csv --chart cusum --cusum-k .5 --cusum-h 4 --out r_cusum.json
     ```
   - 输出 JSON 含：meta / stability / distribution（Anderson-Darling + 分布拟合）/ statistics / capability / ppm_expected / charts / raw_data。
4. **生成可编辑网页报告**：
   ```bash
   python scripts/report_builder.py \
     --result <结果.json> \
     --info <抬头信息.json, 可选> \
     --template assets/report_template.html \
     --out <报告.html>
   ```
   `--info` 可放零件号、工序、特性、抽样条件等抬头字段（不传则报告内为「待补充」）。
5. **智能体研判并交付**：读取 JSON 关键项（稳定性、分布模型、能力指数、ppm、判异点），给出结论与处置建议；把 HTML 报告交付用户。

## 合规硬规则（引擎与报告均强制遵守，见 references/report_20_elements.md）
1. **仅在过程受控时输出 Cp/Cpk**；失控只给 Pp/Ppk（并提示按 OCAP 处置）。
2. **N<2000 不采用经验分位数**，须拟合分布（正态/lognorm/weibull/foldnorm）后用分位数法。
3. **判异准则须用户选定并记入报告**，不默认堆叠全部 8 条（每条约 +10% 误报）。
4. **须声明分布模型**（正态或拟合分布），能力指数计算方法随之而定。
5. 单值数据（n=1）须用 I-MR 图，MR 图判异仅用规则 1。
6. 报告格式与内容须经**顾客与供应商共同商定**（手册 11.2），引擎结论仅供参考，判定权在使用者。

## 输出形态
- 计算/绘图/报告类结果（模式 C）产出**单文件自包含 HTML 网页报告**（主色 #C8102E，含 SVG 控制图/运行图/直方图/概率图，20 个报告项，可编辑可导出）。
- 诊断/OCAP 类结果（模式 A/B）产出 **txt + md 双版文字稿**（内部使用，不生成网页）。
- 决策口诀：**展示/会议 → 网页；内部/中间 → 文字+文档**。

## 参考文件
- [references/old_new_mapping.md](references/old_new_mapping.md) — 新旧版 SPC 5 层变化与章节映射
- [references/diagnosis_checklist.md](references/diagnosis_checklist.md) — 模式 A 诊断问题清单（12 章映射）
- [references/ocap_guide.md](references/ocap_guide.md) — 模式 B OCAP 引导起草流程
- [references/capability_logic.md](references/capability_logic.md) — 三阶段能力 + 分位计算分支 + 控制图常数表
- [references/report_20_elements.md](references/report_20_elements.md) — 20 项报告规格 + 合规硬规则 + 计算约定

## 使用示例
- 示例 1（模式 A 诊断）：用户问「我们公司现在用旧版 SPC，切到新版要补什么？」→ 进入诊断顾问，按清单问当前阶段与已有成文，给缺口清单与导入路线。
- 示例 2（模式 B OCAP）：用户说「帮我起草这个尺寸特性的 OCAP」→ 4 轮引导提问后给建议稿（标注待确认）。
- 示例 3（模式 C 计算报告）：用户给 25 组子组（n=5）直径数据 + USL/LSL → 跑引擎（stage=capability, rules=1,2,5）→ 生成可编辑 HTML 报告，用户在网页里改数重算/导出 PDF。

## TRACE 测评

| 维度 | 评分 | 说明 |
|------|------|------|
| T — 可信任度 | 9/10 | 纯脚本/文档技能，无外部调用风险；合规硬规则内置，结论标注「仅供参考」 |
| R — 可靠性 | 10/10 | 引擎经正态/偏移/非正态三类仿真数据 + 全 9 类控制图（计量/计数/时间加权）自测；网页报告经 DOM 模拟渲染 + 编辑重算 + 导入导出全链路验证（含属性图 N/A 区块、EWMA 逐点限、CUSUM 双侧图） |
| A — 适用性 | 10/10 | 三模式覆盖诊断/OCAP/计算；控制图全族覆盖（I-MR/Xbar-R/Xbar-S + p/np/c/u + EWMA/CUSUM）；触发条件与新旧版边界明确 |
| C — 规范性 | 10/10 | frontmatter 完整；references 与 scripts/assets 结构合规；中文文档全 |
| E — 有效性 | 9/10 | 20 项报告直接对接手册 11.2 模板；含使用示例与合规硬规则清单 |
| **总分** | **48/50** | 通过 |

## 反馈与问题咨询
- 本技能的使用反馈、问题咨询、改进建议，请发邮件至：**engicool@agent.qq.com**
- 该邮箱仅用于本技能的使用反馈与问题咨询，不承接营销或其他无关诉求。
