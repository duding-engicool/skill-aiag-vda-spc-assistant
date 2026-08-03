# -*- coding: utf-8 -*-
"""
spc_engine.py — 新版 AIAG-VDA SPC 计算引擎
按新版手册合规主线计算：稳定性判定 → 分布判定 → 三阶段能力指数 → ppm。
合规硬规则（见 references/report_20_elements.md）在本引擎强制执行：
  1) 过程受控才输出 Cp/Cpk；不受控只输出 Pp/Ppk
  2) N<2000 禁用经验分位数，非正态必须先拟合分布
  3) 判异准则按用户选定子集执行（默认仅准则1），报告注明
  4) 指数注明计算方法（正态经典 / ISO 22514 分位数法）

控制图类型（--chart）：
  计量型 Shewhart：auto | xbar_r | xbar_s | i_mr
  计数型（属性图）：p | np | c | u
  时间加权：        ewma | cusum

用法:
  # 计量型
  python spc_engine.py --data <csv> --subgroup-size <n> --chart auto|xbar_r|xbar_s|i_mr ...
  # 计数型 p/u 图（counts=不合格数/缺陷数, sizes=子组样本量/单位数）
  python spc_engine.py --chart p --counts <csv> --sizes <csv或数字> ...
  # 计数型 c 图（counts=缺陷数, 单位面积恒定）
  python spc_engine.py --chart c --counts <csv> ...
  # 时间加权（作用于计量数据）
  python spc_engine.py --data <csv> --chart ewma [--ewma-lambda .2 --ewma-L 3]
  python spc_engine.py --data <csv> --chart cusum [--cusum-k .5 --cusum-h 4]

数据格式: CSV 单列数值（计量型）；属性图 counts/sizes 各为单列 CSV（表头可有可无）。
输出: JSON（供 report_builder.py 消费），含控制图数据、判异结果、分布检验、统计量、
      能力指数、ppm。
"""
import argparse
import csv
import json
import math
import os
import sys

import warnings

import numpy as np
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------- 控制图常数表 (n=2..10，计量型 Shewhart 用) ----------
CONST = {
    2:  dict(A2=1.880, A3=2.659, d2=1.128, D3=0.0,   D4=3.267, B3=0.0,   B4=3.267, c4=0.7979),
    3:  dict(A2=1.023, A3=1.954, d2=1.693, D3=0.0,   D4=2.574, B3=0.0,   B4=2.568, c4=0.8862),
    4:  dict(A2=0.729, A3=1.628, d2=2.059, D3=0.0,   D4=2.282, B3=0.0,   B4=2.266, c4=0.9213),
    5:  dict(A2=0.577, A3=1.427, d2=2.326, D3=0.0,   D4=2.114, B3=0.0,   B4=2.089, c4=0.9400),
    6:  dict(A2=0.483, A3=1.287, d2=2.534, D3=0.0,   D4=2.004, B3=0.030, B4=1.970, c4=0.9515),
    7:  dict(A2=0.419, A3=1.182, d2=2.704, D3=0.076, D4=1.924, B3=0.118, B4=1.882, c4=0.9594),
    8:  dict(A2=0.373, A3=1.099, d2=2.847, D3=0.136, D4=1.864, B3=0.185, B4=1.815, c4=0.9650),
    9:  dict(A2=0.337, A3=1.032, d2=2.970, D3=0.184, D4=1.816, B3=0.239, B4=1.761, c4=0.9693),
    10: dict(A2=0.308, A3=0.975, d2=3.078, D3=0.223, D4=1.777, B3=0.284, B4=1.716, c4=0.9727),
}

RULE_DESC = {
    1: "准则1: 1点落在±3σ控制限外",
    2: "准则2: 连续9点在中心线同一侧",
    3: "准则3: 连续6点持续上升或下降",
    4: "准则4: 连续14点上下交替",
    5: "准则5: 3点中有2点落在同侧±2σ~±3σ区",
    6: "准则6: 5点中有4点落在同侧超±1σ区",
    7: "准则7: 连续15点落在±1σ区内(中心化)",
    8: "准则8: 连续8点落在±1σ区外(两侧均可)",
}

ATTR_CHARTS = {"p", "np", "c", "u"}
TW_CHARTS = {"ewma", "cusum"}
CONT_CHARTS = {"xbar_r", "xbar_s", "i_mr"}


# ---------- 数据读取 ----------
def read_data(path):
    """读取 CSV 单列数值。"""
    vals = []
    with open(path, "r", encoding="utf-8-sig") as f:
        for row in csv.reader(f):
            if not row:
                continue
            cell = row[0].strip()
            if not cell:
                continue
            try:
                vals.append(float(cell))
            except ValueError:
                continue  # 跳过表头/非数值
    if len(vals) < 10:
        raise SystemExit(f"数据量过少({len(vals)})，至少需要10个观测值")
    return np.array(vals, dtype=float)


def read_pairs(counts_path, sizes_arg, chart):
    """属性图读取：counts 单列；sizes 为 CSV 文件路径或常量数字（c 图可省略）。"""
    counts = read_data(counts_path)
    if chart == "c":
        sizes = None
    else:
        if sizes_arg is None:
            raise SystemExit(f"{chart} 图必须为每个子组提供样本量（--sizes <文件或数字>）")
        if os.path.exists(str(sizes_arg)):
            sizes = read_data(sizes_arg)
            if len(sizes) != len(counts):
                raise SystemExit(f"counts({len(counts)}) 与 sizes({len(sizes)}) 行数不一致")
        else:
            try:
                sval = float(sizes_arg)
            except ValueError:
                raise SystemExit("--sizes 须为 CSV 文件路径或数字")
            sizes = np.full(len(counts), sval)
    return counts, sizes


def make_subgroups(data, n):
    """按行序切子组，末尾不足一组的丢弃并记录。"""
    k = len(data) // n
    dropped = len(data) - k * n
    sub = data[: k * n].reshape(k, n)
    return sub, dropped


# ---------- 判异准则 (作用于点序列；sigma 可为标量或逐点数组) ----------
def apply_rules(points, cl, sigma, rules):
    """返回 {rule_no: [违规点索引...]}，索引0基。sigma 标量或逐点数组均可。"""
    pts = np.asarray(points, dtype=float)
    sa = np.asarray(sigma, dtype=float)
    if sa.ndim == 0:
        sa = np.full_like(pts, float(sa))
    safe = np.where(sa > 0, sa, 1.0)
    z = np.where(sa > 0, (pts - cl) / safe, 0.0)
    viol = {}

    def mark(rno, idxs):
        if idxs:
            viol.setdefault(rno, set()).update(idxs)

    if 1 in rules:
        mark(1, [i for i, v in enumerate(z) if abs(v) > 3])
    if 2 in rules:
        for i in range(len(z) - 8):
            w = z[i:i + 9]
            if np.all(w > 0) or np.all(w < 0):
                mark(2, list(range(i, i + 9)))
    if 3 in rules:
        d = np.diff(pts)
        for i in range(len(d) - 5):
            w = d[i:i + 6]
            if np.all(w > 0) or np.all(w < 0):
                mark(3, list(range(i, i + 7)))
    if 4 in rules:
        d = np.diff(pts)
        for i in range(len(d) - 13):
            w = d[i:i + 14]
            if np.all(w[:-1] * w[1:] < 0):
                mark(4, list(range(i, i + 15)))
    if 5 in rules:
        for i in range(len(z) - 2):
            w = z[i:i + 3]
            if np.sum(w > 2) >= 2 or np.sum(w < -2) >= 2:
                mark(5, [i + j for j in range(3) if abs(w[j]) > 2])
    if 6 in rules:
        for i in range(len(z) - 4):
            w = z[i:i + 5]
            if np.sum(w > 1) >= 4 or np.sum(w < -1) >= 4:
                mark(6, [i + j for j in range(5) if abs(w[j]) > 1])
    if 7 in rules:
        for i in range(len(z) - 14):
            w = z[i:i + 15]
            if np.all(np.abs(w) < 1):
                mark(7, list(range(i, i + 15)))
    if 8 in rules:
        for i in range(len(z) - 7):
            w = z[i:i + 8]
            if np.all(np.abs(w) > 1):
                mark(8, list(range(i, i + 8)))
    return {k: sorted(v) for k, v in viol.items()}


# ---------- 计量型 Shewhart 控制图 ----------
def build_charts(data, sub, chart, rules):
    charts = {}
    n = sub.shape[1] if sub is not None else 1
    if chart == "i_mr":
        mr = np.abs(np.diff(data))
        mrbar = float(np.mean(mr))
        xbar = float(np.mean(data))
        sig = mrbar / 1.128
        cl = dict(center=xbar, ucl=xbar + 2.66 * mrbar, lcl=xbar - 2.66 * mrbar)
        viol = apply_rules(data, xbar, sig, rules)
        charts["primary"] = dict(type="I", points=[round(float(v), 5) for v in data],
                                 **{k: round(v, 5) for k, v in cl.items()}, violations=viol)
        charts["secondary"] = dict(type="MR", points=[round(float(v), 5) for v in mr],
                                   center=round(mrbar, 5), ucl=round(3.267 * mrbar, 5), lcl=0.0,
                                   violations=apply_rules(mr, mrbar, mrbar * 0.8525, {1} & set(rules)))
    else:
        means = np.mean(sub, axis=1)
        xbar = float(np.mean(means))
        c = CONST[n]
        if chart == "xbar_r":
            rng = np.ptp(sub, axis=1)
            rbar = float(np.mean(rng))
            sig_x = (c["A2"] * rbar) / 3.0
            charts["primary"] = dict(type="Xbar", points=[round(float(v), 5) for v in means],
                                     center=round(xbar, 5),
                                     ucl=round(xbar + c["A2"] * rbar, 5),
                                     lcl=round(xbar - c["A2"] * rbar, 5),
                                     violations=apply_rules(means, xbar, sig_x, rules))
            charts["secondary"] = dict(type="R", points=[round(float(v), 5) for v in rng],
                                       center=round(rbar, 5), ucl=round(c["D4"] * rbar, 5),
                                       lcl=round(c["D3"] * rbar, 5),
                                       violations=apply_rules(rng, rbar, (c["D4"] - 1) * rbar / 3, {1} & set(rules)))
        else:  # xbar_s
            sds = np.std(sub, axis=1, ddof=1)
            sbar = float(np.mean(sds))
            sig_x = (c["A3"] * sbar) / 3.0
            charts["primary"] = dict(type="Xbar", points=[round(float(v), 5) for v in means],
                                     center=round(xbar, 5),
                                     ucl=round(xbar + c["A3"] * sbar, 5),
                                     lcl=round(xbar - c["A3"] * sbar, 5),
                                     violations=apply_rules(means, xbar, sig_x, rules))
            charts["secondary"] = dict(type="S", points=[round(float(v), 5) for v in sds],
                                       center=round(sbar, 5), ucl=round(c["B4"] * sbar, 5),
                                       lcl=round(c["B3"] * sbar, 5),
                                       violations=apply_rules(sds, sbar, (c["B4"] - 1) * sbar / 3, {1} & set(rules)))
    return charts


# ---------- 计数型（属性）控制图 ----------
def build_attribute_charts(counts, sizes, chart, rules):
    """p/np/c/u 图。p/u 为比率（逐点 UCL/LCL）；np/c 为计数值（恒定限）。"""
    k = len(counts)
    if chart in ("p", "u"):
        if sizes is None:
            sizes = np.ones(k)
        n_sum = float(np.sum(sizes))
        pbar = float(np.sum(counts)) / n_sum
        center = pbar
        pts, ucls, lcls, sigmas = [], [], [], []
        for i in range(k):
            ni = float(sizes[i])
            if chart == "p":
                sigma_i = math.sqrt(pbar * (1 - pbar) / ni)
                pts.append(float(counts[i] / ni))
                ucls.append(min(1.0, pbar + 3 * sigma_i))
            else:  # u
                sigma_i = math.sqrt(pbar / ni)
                pts.append(float(counts[i] / ni))
                ucls.append(pbar + 3 * sigma_i)
            lcls.append(max(0.0, pbar - 3 * sigma_i))
            sigmas.append(sigma_i)
        primary = dict(type=chart.upper(), points=[round(v, 6) for v in pts],
                       center=round(center, 6),
                       ucl=[round(v, 6) for v in ucls],
                       lcl=[round(v, 6) for v in lcls],
                       violations=apply_rules(pts, center, sigmas, rules))
    elif chart == "np":
        n = float(sizes[0])
        pbar = float(np.sum(counts)) / (n * k)
        npbar = n * pbar
        sigma = math.sqrt(npbar * (1 - pbar))
        ucl = npbar + 3 * sigma
        lcl = max(0.0, npbar - 3 * sigma)
        primary = dict(type="np", points=[round(float(v), 5) for v in counts],
                       center=round(npbar, 5), ucl=round(ucl, 5), lcl=round(lcl, 5),
                       violations=apply_rules(counts, npbar, sigma, rules))
    else:  # c
        cbar = float(np.mean(counts))
        sigma = math.sqrt(cbar)
        ucl = cbar + 3 * sigma
        lcl = max(0.0, cbar - 3 * sigma)
        primary = dict(type="c", points=[round(float(v), 5) for v in counts],
                       center=round(cbar, 5), ucl=round(ucl, 5), lcl=round(lcl, 5),
                       violations=apply_rules(counts, cbar, sigma, rules))
    return {"primary": primary}


# ---------- 时间加权控制图（EWMA / CUSUM，作用于连续数据） ----------
def build_time_weighted(data, sub, chart, rules, target, params):
    """EWMA / CUSUM。n=1 应用于单值，n>1 应用于子组均值；sigma 取相应尺度。"""
    n = sub.shape[1] if sub is not None else 1
    if n == 1:
        series = data
        mr = np.abs(np.diff(data))
        sigma = float(np.mean(mr)) / 1.128
        sigma_within = sigma
    else:
        means = np.mean(sub, axis=1)
        series = means
        sigma = float(np.std(means, ddof=1))          # 子组均值变异
        sigma_within = sigma * math.sqrt(n)
    tgt = float(target) if target is not None else float(np.mean(series))
    charts = {}
    if chart == "ewma":
        lam = params.get("lambda", 0.2)
        L = params.get("L", 3.0)
        z = tgt
        Z = []
        for x in series:
            z = lam * x + (1 - lam) * z
            Z.append(z)
        ucls, lcls = [], []
        for i in range(1, len(series) + 1):
            f = math.sqrt(lam / (2 - lam) * (1 - (1 - lam) ** (2 * i)))
            ucls.append(tgt + L * sigma * f)
            lcls.append(tgt - L * sigma * f)
        sig_e = sigma * math.sqrt(lam / (2 - lam))
        charts["primary"] = dict(type="EWMA", points=[round(float(v), 5) for v in Z],
                                 center=round(tgt, 5),
                                 ucl=[round(v, 5) for v in ucls],
                                 lcl=[round(v, 5) for v in lcls],
                                 violations=apply_rules(Z, tgt, sig_e, rules))
        charts["secondary"] = None
    else:  # cusum（tabular 双侧）
        kk = params.get("k", 0.5)
        h = params.get("h", 4.0)
        shift = kk * sigma
        cplus = cminus = 0.0
        CP, CM = [], []
        for x in series:
            cplus = max(0.0, x - (tgt + shift) + cplus)
            cminus = max(0.0, (tgt - shift) - x + cminus)
            CP.append(cplus)
            CM.append(cminus)
        hsigma = h * sigma
        charts["primary"] = dict(type="CUSUM+", points=[round(float(v), 5) for v in CP],
                                 center=0.0, ucl=round(hsigma, 5), lcl=0.0,
                                 violations=apply_rules(CP, 0.0, hsigma, {1} & set(rules)))
        charts["secondary"] = dict(type="CUSUM-", points=[round(float(v), 5) for v in CM],
                                   center=0.0, ucl=round(hsigma, 5), lcl=0.0,
                                   violations=apply_rules(CM, 0.0, hsigma, {1} & set(rules)))
    return charts, sigma_within


# ---------- 分布拟合 ----------
def fit_distributions(data):
    """AD 正态检验 + 备选分布拟合，返回判定结果。"""
    res = {"normal": {}, "candidates": [], "selected": None}
    ad = stats.anderson(data, dist="norm")
    crit_5 = None
    for sl, cv in zip(ad.significance_level, ad.critical_values):
        if abs(sl - 5.0) < 1e-9:
            crit_5 = cv
    is_normal = bool(ad.statistic < crit_5)
    res["normal"] = {
        "method": "Anderson-Darling",
        "statistic": round(float(ad.statistic), 4),
        "critical_5pct": round(float(crit_5), 4),
        "is_normal": is_normal,
    }
    if is_normal:
        mu, s = float(np.mean(data)), float(np.std(data, ddof=1))
        res["selected"] = {
            "name": "normal", "label": "正态分布",
            "params": {"mu": mu, "sigma": s},
            "quantiles": {
                "x_00135": mu - 3 * s, "x_50": mu, "x_99865": mu + 3 * s,
            },
        }
        return res

    candidates = []
    shift = 0.0
    d = data.copy()
    if np.min(d) <= 0:
        shift = float(np.min(d)) - 1e-6
        d = d - shift
    try:
        p = stats.lognorm.fit(d, floc=0)
        a2 = _ad_stat(d, stats.lognorm, p)
        candidates.append(("lognorm", "对数正态分布", p, a2))
    except Exception:
        pass
    try:
        p = stats.weibull_min.fit(d, floc=0)
        a2 = _ad_stat(d, stats.weibull_min, p)
        candidates.append(("weibull_min", "威布尔分布", p, a2))
    except Exception:
        pass
    try:
        p = stats.foldnorm.fit(d, floc=0)
        a2 = _ad_stat(d, stats.foldnorm, p)
        candidates.append(("foldnorm", "折叠正态分布", p, a2))
    except Exception:
        pass

    if not candidates:
        raise SystemExit("非正态且所有备选分布拟合失败，请人工检查数据")
    candidates.sort(key=lambda c: c[3])
    for name, label, p, a2 in candidates:
        res["candidates"].append({"name": name, "label": label, "ad_stat": round(a2, 4)})
    name, label, p, a2 = candidates[0]
    dist = getattr(stats, name)
    q = lambda pr: float(dist.ppf(pr, *p)) + shift
    res["selected"] = {
        "name": name, "label": label,
        "params": [round(float(x), 6) for x in p], "shift": shift,
        "ad_stat": round(a2, 4),
        "quantiles": {"x_00135": q(0.00135), "x_50": q(0.5), "x_99865": q(0.99865)},
    }
    return res


def _ad_stat(data, dist, params):
    """通用 AD 统计量（对拟合分布）。"""
    x = np.sort(data)
    n = len(x)
    cdf = np.clip(dist.cdf(x, *params), 1e-10, 1 - 1e-10)
    i = np.arange(1, n + 1)
    return float(-n - np.mean((2 * i - 1) * (np.log(cdf) + np.log(1 - cdf[::-1]))))


# ---------- ppm ----------
def calc_ppm(sel, usl, lsl):
    """按所选分布模型计算规格外比例 ppm。"""
    name = sel["name"]
    if name == "normal":
        mu, s = sel["params"]["mu"], sel["params"]["sigma"]
        p_out = 0.0
        if usl is not None:
            p_out += 1 - stats.norm.cdf(usl, mu, s)
        if lsl is not None:
            p_out += stats.norm.cdf(lsl, mu, s)
    else:
        dist = getattr(stats, name)
        params, shift = sel["params"], sel.get("shift", 0.0)
        p_out = 0.0
        if usl is not None:
            p_out += 1 - float(dist.cdf(usl - shift, *params))
        if lsl is not None:
            p_out += float(dist.cdf(lsl - shift, *params))
    return round(p_out * 1e6, 2)


# ---------- 能力指数 ----------
def capability_indices(data, sub, sel, in_control, usl, lsl, stage, chart):
    """按合规分叉输出指数。chart 用于选组内σ估计方式。"""
    out = {"stage": stage, "method": None, "indices": {}, "notes": []}
    n = sub.shape[1] if sub is not None else 1
    s_total = float(np.std(data, ddof=1))
    mu = float(np.mean(data))

    # 组内σ
    if chart in ("xbar_r", "i_mr"):
        if chart == "i_mr":
            mr = np.abs(np.diff(data))
            sigma_w = float(np.mean(mr) / 1.128)
        else:
            rbar = float(np.mean(np.ptp(sub, axis=1)))
            sigma_w = rbar / CONST[n]["d2"]
    else:  # xbar_s
        sbar = float(np.mean(np.std(sub, axis=1, ddof=1)))
        sigma_w = sbar / CONST[n]["c4"]

    is_norm = sel["name"] == "normal"
    q = sel["quantiles"]
    spread = q["x_99865"] - q["x_00135"]

    def two_sided(lo_s, hi_s, span):
        v = {}
        if usl is not None and lsl is not None:
            v["p"] = (usl - lsl) / span
        if usl is not None:
            v["pu"] = (usl - q["x_50"]) / hi_s
        if lsl is not None:
            v["pl"] = (q["x_50"] - lsl) / lo_s
        ks = [v[k] for k in ("pu", "pl") if k in v]
        if ks:
            v["pk"] = min(ks)
        return v

    if usl is None and lsl is None:
        out["notes"].append("未提供规格限，仅做稳定性与分布分析，不输出能力指数")
        return out, sigma_w, s_total

    prefix_perf = {"machine": "Pm", "performance": "Pp", "capability": "Pp"}[stage]

    if is_norm:
        out["method"] = "正态经典公式（ISO 22514 分位数法在正态下的特例）"
        vp = two_sided(3 * s_total, 3 * s_total, 6 * s_total)
        for k, val in vp.items():
            out["indices"][prefix_perf + ("" if k == "p" else "k" if k == "pk" else k.upper()[-1])] = round(val, 3)
        if stage == "capability":
            if in_control:
                vc = two_sided(3 * sigma_w, 3 * sigma_w, 6 * sigma_w)
                for k, val in vc.items():
                    out["indices"]["C" + ("p" if k == "p" else "pk" if k == "pk" else "p" + k[-1].upper())] = round(val, 3)
            else:
                out["notes"].append("过程不受控：按新版硬规则不输出 Cp/Cpk，仅输出 Pp/Ppk；请先按 OCAP 处置并恢复受控")
    else:
        out["method"] = f"ISO 22514 分位数法（{sel['label']}，X0.135%/X50%/X99.865%）"
        lo_s = q["x_50"] - q["x_00135"]
        hi_s = q["x_99865"] - q["x_50"]
        vp = two_sided(lo_s, hi_s, spread)
        for k, val in vp.items():
            out["indices"][prefix_perf + ("" if k == "p" else "k" if k == "pk" else k.upper()[-1])] = round(val, 3)
        if stage == "capability":
            if in_control:
                out["notes"].append("非正态受控过程：Cp/Cpk 同样按分位数法（组内变差的分位数估计较复杂，本引擎以总体分位数近似并注明）")
                for k, val in vp.items():
                    key = "C" + ("p" if k == "p" else "pk" if k == "pk" else "p" + k[-1].upper())
                    out["indices"][key] = round(val, 3)
            else:
                out["notes"].append("过程不受控：不输出 Cp/Cpk")
        if len(data) < 2000:
            out["notes"].append(f"N={len(data)}<2000：按新版硬规则未使用经验分位数，分位数来自拟合分布（{sel['label']}）")
    return out, sigma_w, s_total


# ---------- 属性图 NA 占位（分布/能力不适用） ----------
def attr_na_blocks(stage):
    return (
        {"normal": {"method": "—", "statistic": None, "critical_5pct": None, "is_normal": None},
         "selected": {"name": "attribute", "label": "计数型控制图（不适用连续分布分析）"},
         "candidates": []},
        {"mean": None, "std_total": None, "std_within": None, "min": None, "max": None,
         "x_00135": None, "x_50": None, "x_99865": None, "mean_ci95": [None, None]},
        {"stage": stage, "method": "计数型控制图不适用过程能力指数（Cp/Cpk/Pp/Ppk）", "indices": {},
         "notes": ["计数型控制图判定合格/不合格，不计算连续型过程能力指数"]},
        None,
    )


def main():
    ap = argparse.ArgumentParser(description="AIAG-VDA 新版 SPC 计算引擎")
    ap.add_argument("--data", default=None, help="CSV 数据文件（计量型/时间加权，单列数值）")
    ap.add_argument("--counts", default=None, help="属性图计数列 CSV（不合格数/缺陷数）")
    ap.add_argument("--sizes", default=None, help="属性图样本量：CSV 文件路径或常量数字（c 图可省略）")
    ap.add_argument("--subgroup-size", type=int, default=1, help="子组大小(1=单值)")
    ap.add_argument("--chart", default="auto",
                    choices=["auto", "xbar_r", "xbar_s", "i_mr", "p", "np", "c", "u", "ewma", "cusum"])
    ap.add_argument("--usl", type=float, default=None)
    ap.add_argument("--lsl", type=float, default=None)
    ap.add_argument("--target", type=float, default=None)
    ap.add_argument("--rules", default="1", help="判异准则子集，如 1,2,5,6（默认仅1，须经用户选定）")
    ap.add_argument("--stage", default="capability", choices=["machine", "performance", "capability"])
    ap.add_argument("--ewma-lambda", type=float, default=0.2, help="EWMA 平滑系数 λ(默认0.2)")
    ap.add_argument("--ewma-L", type=float, default=3.0, help="EWMA 控制限倍数 L(默认3)")
    ap.add_argument("--cusum-k", type=float, default=0.5, help="CUSUM  slack k(以σ计,默认0.5)")
    ap.add_argument("--cusum-h", type=float, default=4.0, help="CUSUM 决策区间 h(以σ计,默认4)")
    ap.add_argument("--out", default=None, help="输出 JSON 路径（默认当前工作目录 spc_result.json）")
    args = ap.parse_args()

    rules = set(int(r) for r in args.rules.split(",") if r.strip())
    bad = rules - set(RULE_DESC)
    if bad:
        raise SystemExit(f"无效判异准则: {bad}")

    chart = args.chart
    is_attr = chart in ATTR_CHARTS
    is_tw = chart in TW_CHARTS

    # ---------- 属性图 ----------
    if is_attr:
        counts, sizes = read_pairs(args.counts, args.sizes, chart)
        k = len(counts)
        charts = build_attribute_charts(counts, sizes, chart, rules)
        in_control = not any(charts[c]["violations"] for c in charts if charts.get(c))
        dist, stats_blk, cap, ppm = attr_na_blocks(args.stage)
        mu = float(np.mean(counts))
        result = {
            "meta": {
                "engine": "spc_engine v1.1 (AIAG-VDA 新版合规主线)",
                "data_file": os.path.basename(args.counts or "?"),
                "chart_type": chart, "is_attribute": True,
                "n_total": int(k), "subgroup_size": (float(sizes[0]) if sizes is not None else 1.0),
                "n_subgroups": int(k),
                "dropped_tail": 0, "stage": args.stage,
                "rules_used": sorted(rules), "rules_desc": [RULE_DESC[r] for r in sorted(rules)],
                "rule_note": "判异准则须按新版9.2.2由用户选定成文；每增加一条准则第一类错误概率约增10%",
                "usl": args.usl, "lsl": args.lsl, "target": args.target,
                "attr_chart": chart,
            },
            "stability": {
                "in_control": in_control,
                "statement": "过程受控（所选判异准则下无违规）" if in_control
                             else "过程不受控（存在判异违规，详见控制图），应按 OCAP 处置",
                "violations_summary": {c: {str(r): idx for r, idx in charts[c]["violations"].items()}
                                        for c in charts if charts.get(c)},
            },
            "distribution": dist,
            "statistics": stats_blk,
            "capability": cap,
            "ppm_expected": ppm,
            "hints": ["计数型控制图判定合格/不合格，不涉及连续型过程能力指数与 ppm"],
            "charts": charts,
            "attr": {
                "chart": chart,
                "counts": [round(float(v), 5) for v in counts],
                "sizes": None if sizes is None else [round(float(v), 5) for v in sizes],
            },
            "raw_data": [round(float(v), 5) for v in counts],
        }
        out = args.out or os.path.join(os.getcwd(), "spc_result.json")
        with open(out, "w", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False, indent=1)
        print(f"[spc_engine] 完成(属性图): {out}")
        print(f"  图类型: {chart.upper()} | 受控性: {result['stability']['statement']}")
        print(f"  子组数: {k} | 计数合计: {int(np.sum(counts))}")
        return

    # ---------- 计量型 / 时间加权（连续数据） ----------
    if args.data is None:
        raise SystemExit("计量型/时间加权图必须提供 --data（CSV 单列数值）")
    data = read_data(args.data)
    n = max(1, args.subgroup_size)
    if chart == "auto":
        chart = "i_mr" if n == 1 else ("xbar_r" if n <= 8 else "xbar_s")

    if chart == "i_mr":
        sub, dropped = data.reshape(-1, 1), 0
    elif is_tw:
        if n == 1:
            sub, dropped = None, 0
        else:
            if n < 2 or n > 10:
                raise SystemExit("子组大小须在2-10之间(超出请用 xbar_s 并扩充常数表)")
            sub, dropped = make_subgroups(data, n)
    else:  # xbar_r / xbar_s
        if n < 2 or n > 10:
            raise SystemExit("子组大小须在2-10之间(超出请用 xbar_s 并扩充常数表)")
        sub, dropped = make_subgroups(data, n)

    if is_tw:
        params = {"lambda": args.ewma_lambda, "L": args.ewma_L,
                  "k": args.cusum_k, "h": args.cusum_h} if chart == "ewma" else \
                 {"k": args.cusum_k, "h": args.cusum_h}
        charts, sigma_w = build_time_weighted(data, sub if chart != "i_mr" else None,
                                              chart, rules, args.target, params)
        in_control = not any(charts[c]["violations"] for c in charts if charts.get(c))
        # 能力分析仍基于底层连续数据（组内σ）
        dist = fit_distributions(data)
        cap, _, s_total = capability_indices(
            data, sub, dist["selected"], in_control, args.usl, args.lsl, args.stage,
            "i_mr" if sub is None else "xbar_s")
        ppm = calc_ppm(dist["selected"], args.usl, args.lsl) if (args.usl is not None or args.lsl is not None) else None
    else:
        charts = build_charts(data, sub if chart != "i_mr" else None, chart, rules)
        in_control = not any(charts[c]["violations"] for c in charts if charts.get(c))
        dist = fit_distributions(data)
        cap, sigma_w, s_total = capability_indices(
            data, sub if chart != "i_mr" else None, dist["selected"],
            in_control, args.usl, args.lsl, args.stage, chart)
        ppm = calc_ppm(dist["selected"], args.usl, args.lsl) if (args.usl is not None or args.lsl is not None) else None

    # 时间相关模型提示（σb/σw 检测）
    hints = []
    if chart in ("xbar_r", "xbar_s"):
        means = np.mean(sub, axis=1)
        sigma_b = float(np.std(means, ddof=1))
        if sigma_w > 0 and sigma_b / (sigma_w / math.sqrt(sub.shape[1])) > 2.0:
            hints.append("子组间变差显著大于组内变差：疑似时间相关分布模型(ISO 22514-2 C/D类)，建议核查漂移因素(刀具磨损/批间差异)")

    mu = float(np.mean(data))
    q = dist["selected"]["quantiles"]
    se = s_total / math.sqrt(len(data))
    result = {
        "meta": {
            "engine": "spc_engine v1.1 (AIAG-VDA 新版合规主线)",
            "data_file": os.path.basename(args.data),
            "n_total": len(data), "subgroup_size": n,
            "n_subgroups": int(sub.shape[0]) if (sub is not None and chart != "i_mr") else len(data),
            "dropped_tail": dropped, "chart_type": chart, "is_attribute": False,
            "stage": args.stage,
            "rules_used": sorted(rules), "rules_desc": [RULE_DESC[r] for r in sorted(rules)],
            "rule_note": "判异准则须按新版9.2.2由用户选定成文；每增加一条准则第一类错误概率约增10%",
            "usl": args.usl, "lsl": args.lsl, "target": args.target,
            "ewma_lambda": args.ewma_lambda, "ewma_L": args.ewma_L,
            "cusum_k": args.cusum_k, "cusum_h": args.cusum_h,
        },
        "stability": {
            "in_control": in_control,
            "statement": "过程受控（所选判异准则下无违规）" if in_control
                         else "过程不受控（存在判异违规，详见控制图），应按 OCAP 处置",
            "violations_summary": {c: {str(r): idx for r, idx in charts[c]["violations"].items()}
                                    for c in charts if charts.get(c)},
        },
        "distribution": dist,
        "statistics": {
            "mean": round(mu, 5), "std_total": round(s_total, 5), "std_within": round(sigma_w, 5),
            "min": round(float(np.min(data)), 5), "max": round(float(np.max(data)), 5),
            "x_00135": round(q["x_00135"], 5), "x_50": round(q["x_50"], 5), "x_99865": round(q["x_99865"], 5),
            "mean_ci95": [round(mu - 1.96 * se, 5), round(mu + 1.96 * se, 5)],
        },
        "capability": cap,
        "ppm_expected": ppm,
        "hints": hints,
        "charts": charts,
        "raw_data": [round(float(v), 5) for v in data],
    }

    out = args.out or os.path.join(os.getcwd(), "spc_result.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=1)
    print(f"[spc_engine] 完成: {out}")
    print(f"  受控性: {result['stability']['statement']}")
    dist_label = dist["selected"]["label"] if dist["selected"] else "—"
    print(f"  分布: {dist_label} | 指数: {cap['indices']} | ppm: {ppm}")


if __name__ == "__main__":
    main()
