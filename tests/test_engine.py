#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
新版 SPC 引擎自测（对齐开发期验收铁律：核心公式 + 端到端）
- 合成数据校验 Cp/Cpk/Cw 公式、w 分布精确限、Nelson、Anderson-Darling
- 端到端用技能内置样本数据（tests 内联，不依赖外部文件），核对关键结论（稳态/Form/能力/图表类型）
"""
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

SKILL_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(SKILL_ROOT, "scripts"))
import spc_v2_engine as eng  # noqa: E402

PY = r"C:\Users\admin\.workbuddy\binaries\python\envs\default\Scripts\python.exe"
ENGINE = os.path.join(SKILL_ROOT, "scripts", "spc_v2_engine.py")

PASS = 0
FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  PASS  {name}  {detail}")
    else:
        FAIL += 1
        print(f"  FAIL  {name}  {detail}")


def synth_matrix(mean, sigma, k, n, seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(mean, sigma, size=(k, n)).round(4).tolist()


def test_capability_formula():
    print("[1] 能力指数公式（合成正态数据）")
    # μ=25, σ=1, tol=6 -> Cp 理论 = 1.0
    mat = synth_matrix(25, 1.0, 20, 5, seed=1)
    r = eng.xbar_r_analysis(mat, usl=28, lsl=22, target=25)
    cap = r["capability"]
    check("Cp ≈ 1.0", abs(cap["Cp"] - 1.0) < 0.25, f"Cp={cap['Cp']}")
    check("Cp≈Cpk（中心重合）", abs(cap["Cp"] - cap["Cpk"]) < 0.2, f"Cpk={cap['Cpk']}")
    check("Pp≈Cp（稳定正态）", abs(cap["Pp"] - cap["Cp"]) < 0.2, f"Pp={cap['Pp']}")
    check("Cw 存在且为正", isinstance(cap["Cw"], (int, float)) and cap["Cw"] > 0, f"Cw={cap['Cw']}")
    check("CpG 通用几何法存在", cap["CpG"] is not None and cap["CpG"] > 0, f"CpG={cap['CpG']}")


def test_w_distribution_limit():
    print("[2] w 分布精确 R 限（n=5, σ=1）")
    mat = synth_matrix(25, 1.0, 20, 5, seed=2)
    r = eng.xbar_r_analysis(mat, usl=28, lsl=22, target=25)
    ucl_r = r["control_limits"]["R"]["UCL"]
    # 标准极差 99.865% 分位数（n=5，蒙特卡洛精确值）对照内嵌表
    w5 = eng._W_HI[5]
    check("UCL_R ≈ w_hi*sigma", abs(ucl_r / r["sigma_est"] - w5) < 0.15,
          f"UCL_R/sig={ucl_r / r['sigma_est']:.3f} w5={w5:.3f}")
    lcl_r = r["control_limits"]["R"]["LCL"]
    check("LCL_R ≥ 0", lcl_r >= 0, f"LCL_R={lcl_r}")


def test_nelson():
    print("[3] Nelson 规则触发")
    # 构造 6 点单调上升 -> 规则3
    rising = [[i + 0.01 * j for j in range(5)] for i in range(20)]
    r = eng.xbar_r_analysis(rising, usl=None, lsl=None, target=None)
    check("规则3（连续6升）触发", len(r["nelson_violations"]["3"]) > 0,
          f"v3={r['nelson_violations']['3']}")
    # 纯随机稳定数据不应触发规则1
    mat = synth_matrix(25, 1.0, 20, 5, seed=3)
    r2 = eng.xbar_r_analysis(mat, usl=None, lsl=None, target=None)
    check("随机数据规则1不误报", len(r2["nelson_violations"]["1"]) == 0)


def test_anderson():
    print("[4] Anderson-Darling 正态性")
    normal = synth_matrix(25, 1.0, 30, 5, seed=4)
    flat = np.array(normal).flatten()
    _, rej_n, concl_n = eng.anderson_test(flat)
    check("正态数据不拒绝", not rej_n, concl_n)
    # 指数（右偏）应拒绝
    rng = np.random.default_rng(5)
    expo = rng.exponential(1.0, size=150)
    _, rej_e, concl_e = eng.anderson_test(expo)
    check("右偏数据拒绝正态", rej_e, concl_e)


def test_stable_hard_rule():
    print("[5] 稳态硬规则（不稳态 Cp/Cpk -> NA）")
    rising = [[i + 0.01 * j for j in range(5)] for i in range(20)]
    r = eng.xbar_r_analysis(rising, usl=28, lsl=22, target=25)
    check("不稳态判定", r["stable"] is False)
    check("Cp -> NA", r["capability"]["Cp"] == "NA")
    check("Cpk -> NA", r["capability"]["Cpk"] == "NA")
    check("Pp 仍可报（整体）", isinstance(r["capability"]["Pp"], (int, float)))


def test_attribute():
    print("[6] 计数值精确分位数限")
    # c 图：20 点缺陷数，均值≈4
    rng = np.random.default_rng(7)
    c = rng.poisson(4.0, size=20).tolist()
    res = eng.attribute_analysis(c, [1] * 20, chart='c')
    check("c 图类型", res["chart_type"] == "c")
    check("c 图 UCL>LCL", res["control_limits"][0]["UCL"] > res["control_limits"][0]["LCL"])
    # p 图：20 批 n=50，p≈0.1
    x = rng.binomial(50, 0.1, size=20).tolist()
    res2 = eng.attribute_analysis(x, [50] * 20, chart='p')
    check("p 图类型", res2["chart_type"] == "p")
    check("p 图 UCL> phat", res2["control_limits"][0]["UCL"] > res2["rate"])


def _embedded_samples():
    """技能内置样本（不依赖外部文件）。固定种子保证可复现。"""
    rng = np.random.default_rng(20260822)
    # 1) 正态稳态、能力充足：25 子组 ×5，中心 25，sd≈0.8，USL=28 LSL=22
    stable = np.round(25 + rng.normal(0, 0.8, size=(25, 5)), 3).tolist()
    # 2) 不稳定（单调趋势）：20 子组 ×5，均值 20→30
    trend = [np.round((20 + i * 0.5) + rng.normal(0, 0.5, size=5), 3).tolist()
             for i in range(20)]
    # 3) c 图：20 点缺陷数 ~ Poisson(4)
    c = rng.poisson(4.0, size=20).tolist()
    # 4) p 图：20 批 n=50，缺陷 ~ Binomial(50,0.10)
    p_x = rng.binomial(50, 0.10, size=20).tolist()
    p_n = [50] * 20
    return stable, trend, c, p_x, p_n


def _write_csv(rows, path):
    with open(path, "w", encoding="utf-8") as f:
        for r in rows:
            if isinstance(r, (list, tuple)):
                f.write(",".join(str(v) for v in r) + "\n")
            else:
                f.write(str(r) + "\n")


def _run_engine(path, extra=None):
    cmd = [PY, ENGINE, "--data-file", path, "--format", "json"]
    if extra:
        cmd += extra
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        return None, out.stderr
    return json.loads(out.stdout), None


def test_e2e_internal():
    print("[7] 端到端：技能内置样本数据（不依赖外部文件）")
    stable, trend, c, p_x, p_n = _embedded_samples()

    # 7.1 计量值 Xbar-R：正态稳态
    p1 = os.path.join(tempfile.gettempdir(), "spc_e2e_stable.csv")
    _write_csv(stable, p1)
    res, err = _run_engine(p1, ["--usl", "28", "--lsl", "22"])
    if err:
        check("e2e 稳态数据运行", False, err[:200])
    else:
        check("e2e 稳态判定", res.get("stable") is True)
        check("e2e Form=A", res.get("form") == "A", f"form={res.get('form')}")
        check("e2e 图表=Xbar-R", res.get("chart_type") == "Xbar-R")

    # 7.2 计量值 Xbar-R：趋势不稳态
    p2 = os.path.join(tempfile.gettempdir(), "spc_e2e_trend.csv")
    _write_csv(trend, p2)
    res2, err2 = _run_engine(p2)
    if err2:
        check("e2e 趋势数据运行", False, err2[:200])
    else:
        check("e2e 趋势判异触发", res2.get("any_violation") is True)
        check("e2e 趋势 Cp->NA", res2["capability"]["Cp"] == "NA")

    # 7.3 计数值 c 图（单列缺陷数，自动判定）
    p3 = os.path.join(tempfile.gettempdir(), "spc_e2e_c.csv")
    _write_csv(c, p3)
    res3, err3 = _run_engine(p3)
    if err3:
        check("e2e c 数据运行", False, err3[:200])
    else:
        check("e2e c 图类型", res3.get("chart_type") == "c")

    # 7.4 计数值 p 图（两列 n,x，自动判定）
    p4 = os.path.join(tempfile.gettempdir(), "spc_e2e_p.csv")
    _write_csv(list(zip(p_n, p_x)), p4)
    res4, err4 = _run_engine(p4)
    if err4:
        check("e2e p 数据运行", False, err4[:200])
    else:
        check("e2e p 图类型", res4.get("chart_type") == "p")


if __name__ == "__main__":
    test_capability_formula()
    test_w_distribution_limit()
    test_nelson()
    test_anderson()
    test_stable_hard_rule()
    test_attribute()
    test_e2e_internal()
    print(f"\n=== 结果：{PASS} PASS / {FAIL} FAIL ===")
    sys.exit(1 if FAIL > 0 else 0)
