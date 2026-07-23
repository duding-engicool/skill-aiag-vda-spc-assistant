# -*- coding: utf-8 -*-
"""
spc_engine.py — 新版 AIAG-VDA SPC 计算引擎
按新版手册合规主线计算：稳定性判定 → 分布判定 → 三阶段能力指数 → ppm。
合规硬规则（见 references/report_20_elements.md）在本引擎强制执行：
  1) 过程受控才输出 Cp/Cpk；不受控只输出 Pp/Ppk
  2) N<2000 禁用经验分位数，非正态必须先拟合分布
  3) 判异准则按用户选定子集执行（默认仅准则1），报告注明
  4) 指数注明计算方法（正态经典 / ISO 22514 分位数法）

用法:
  python spc_engine.py --data <csv文件> --subgroup-size <n> \
      [--chart auto|xbar_r|xbar_s|i_mr] [--usl U] [--lsl L] [--target T] \
      [--rules 1,2,5] [--stage capability|performance|machine] \
      [--out result.json]

数据格式: CSV，单列数值（表头可有可无）；按行序依次分组为子组。
输出: JSON（供 report_builder.py 消费），含控制图数据、判异结果、
      分布检验、统计量/分位数、能力指数、ppm。
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

# ---------- 控制图常数表 (n=2..10) ----------
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


def make_subgroups(data, n):
    """按行序切子组，末尾不足一组的丢弃并记录。"""
    k = len(data) // n
    dropped = len(data) - k * n
    sub = data[: k * n].reshape(k, n)
    return sub, dropped


# ---------- 判异准则 (作用于均值图/单值图点序列) ----------
def apply_rules(points, cl, sigma, rules):
    """返回 {rule_no: [违规点索引...]}，索引0基。"""
    pts = np.asarray(points)
    z = (pts - cl) / sigma if sigma > 0 else np.zeros_like(pts)
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


# ---------- 分布拟合 ----------
def fit_distributions(data):
    """AD 正态检验 + 备选分布拟合，返回判定结果。"""
    res = {"normal": {}, "candidates": [], "selected": None}
    ad = stats.anderson(data, dist="norm")
    # α=0.05 临界值位于 significance_level==5.0
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

    # 非正态：拟合候选分布，按 AD 统计量最小选取
    candidates = []
    shift = 0.0
    d = data.copy()
    if np.min(d) <= 0:  # 对数正态/威布尔要求正值，做平移记录
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
    """按合规分叉输出指数。"""
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
        # 性能指数（总体 s）
        vp = two_sided(3 * s_total, 3 * s_total, 6 * s_total)
        for k, val in vp.items():
            out["indices"][prefix_perf + ("" if k == "p" else "k" if k == "pk" else k.upper()[-1])] = round(val, 3)
        # 受控才允许 Cp/Cpk（组内σ）——仅 capability 阶段
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


# ---------- 控制图 ----------
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


def main():
    ap = argparse.ArgumentParser(description="AIAG-VDA 新版 SPC 计算引擎")
    ap.add_argument("--data", required=True, help="CSV 数据文件（单列数值）")
    ap.add_argument("--subgroup-size", type=int, default=1, help="子组大小(1=单值)")
    ap.add_argument("--chart", default="auto", choices=["auto", "xbar_r", "xbar_s", "i_mr"])
    ap.add_argument("--usl", type=float, default=None)
    ap.add_argument("--lsl", type=float, default=None)
    ap.add_argument("--target", type=float, default=None)
    ap.add_argument("--rules", default="1", help="判异准则子集，如 1,2,5,6（默认仅1，须经用户选定）")
    ap.add_argument("--stage", default="capability", choices=["machine", "performance", "capability"])
    ap.add_argument("--out", default=None, help="输出 JSON 路径（默认当前工作目录 spc_result.json）")
    args = ap.parse_args()

    data = read_data(args.data)
    n = max(1, args.subgroup_size)
    rules = set(int(r) for r in args.rules.split(",") if r.strip())
    bad = rules - set(RULE_DESC)
    if bad:
        raise SystemExit(f"无效判异准则: {bad}")

    chart = args.chart
    if chart == "auto":
        chart = "i_mr" if n == 1 else ("xbar_r" if n <= 8 else "xbar_s")
    if chart == "i_mr":
        sub, dropped = data.reshape(-1, 1), 0
    else:
        if n < 2 or n > 10:
            raise SystemExit("子组大小须在2-10之间(超出请用 xbar_s 并扩充常数表)")
        sub, dropped = make_subgroups(data, n)

    charts = build_charts(data, sub if chart != "i_mr" else None, chart, rules)
    in_control = not any(charts[c]["violations"] for c in charts)

    dist = fit_distributions(data)
    cap, sigma_w, s_total = capability_indices(
        data, sub if chart != "i_mr" else None, dist["selected"],
        in_control, args.usl, args.lsl, args.stage, chart)
    ppm = calc_ppm(dist["selected"], args.usl, args.lsl) if (args.usl is not None or args.lsl is not None) else None

    # 时间相关模型提示（σb/σw 检测）
    hints = []
    if chart != "i_mr":
        means = np.mean(sub, axis=1)
        sigma_b = float(np.std(means, ddof=1))
        if sigma_w > 0 and sigma_b / (sigma_w / math.sqrt(sub.shape[1])) > 2.0:
            hints.append("子组间变差显著大于组内变差：疑似时间相关分布模型(ISO 22514-2 C/D类)，建议核查漂移因素(刀具磨损/批间差异)")

    mu = float(np.mean(data))
    q = dist["selected"]["quantiles"]
    se = s_total / math.sqrt(len(data))
    result = {
        "meta": {
            "engine": "spc_engine v1.0 (AIAG-VDA 新版合规主线)",
            "data_file": os.path.basename(args.data),
            "n_total": len(data), "subgroup_size": n, "n_subgroups": int(sub.shape[0]) if chart != "i_mr" else len(data),
            "dropped_tail": dropped, "chart_type": chart, "stage": args.stage,
            "rules_used": sorted(rules), "rules_desc": [RULE_DESC[r] for r in sorted(rules)],
            "rule_note": "判异准则须按新版9.2.2由用户选定成文；每增加一条准则第一类错误概率约增10%",
            "usl": args.usl, "lsl": args.lsl, "target": args.target,
        },
        "stability": {
            "in_control": in_control,
            "statement": "过程受控（所选判异准则下无违规）" if in_control
                         else "过程不受控（存在判异违规，详见控制图），应按 OCAP 处置",
            "violations_summary": {c: {str(r): idx for r, idx in charts[c]["violations"].items()} for c in charts},
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
    print(f"  分布: {dist['selected']['label']} | 指数: {cap['indices']} | ppm: {ppm}")
    for note in cap["notes"] + hints:
        print(f"  注: {note}")


if __name__ == "__main__":
    main()
