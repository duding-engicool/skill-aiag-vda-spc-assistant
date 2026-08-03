# -*- coding: utf-8 -*-
"""
report_builder.py — 将 spc_engine.py 的计算结果注入可编辑 HTML 报告模板。

用法:
  python report_builder.py --result spc_result.json [--info info.json] \
      [--template <模板路径>] [--out SPC报告.html]

info.json（可选）为报告抬头/条件信息（要素1-10、20），键名见模板 IDENT/COND/MEAS：
  part_no, part_name, process, characteristic, unit, special_char,
  location, period, sampling, operator, gauge, environment,
  resolution, accuracy, grr, conclusion
未提供的字段在报告中显示【待补充】，用户可在网页编辑态直接填写。

输出: 单文件自包含 HTML（可编辑/保存草稿/导出JSON/另存HTML/打印PDF），
      写入用户当前工作目录（或 --out 指定路径）。
"""
import argparse
import json
import os
import sys
import time
import warnings

import numpy as np
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

DEFAULT_TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "assets", "report_template.html")


def pdf_curve(result, n_pts=80):
    """按引擎所选分布模型生成拟合曲线点 [(x, pdf)]，供直方图初始态使用。
    属性图（计数型）不适用连续分布，返回空列表。"""
    if result.get("meta", {}).get("is_attribute"):
        return []
    sel = result["distribution"]["selected"]
    data = np.array(result["raw_data"])
    lo, hi = float(np.min(data)), float(np.max(data))
    xs = np.linspace(lo, hi, n_pts)
    if sel["name"] == "normal":
        mu, s = sel["params"]["mu"], sel["params"]["sigma"]
        ys = stats.norm.pdf(xs, mu, s)
    else:
        dist = getattr(stats, sel["name"])
        params, shift = sel["params"], sel.get("shift", 0.0)
        ys = dist.pdf(xs - shift, *params)
    return [[round(float(x), 6), round(float(y), 8)] for x, y in zip(xs, ys)]


def prob_theo(result):
    """按所选分布生成概率图理论分位数（与排序数据一一对应）。
    属性图（计数型）不适用概率图，返回空列表。"""
    if result.get("meta", {}).get("is_attribute"):
        return []
    sel = result["distribution"]["selected"]
    n = len(result["raw_data"])
    p = (np.arange(1, n + 1) - 0.5) / n
    if sel["name"] == "normal":
        q = stats.norm.ppf(p)  # 标准正态分位数
    else:
        dist = getattr(stats, sel["name"])
        params, shift = sel["params"], sel.get("shift", 0.0)
        q = dist.ppf(p, *params) + shift
    return [round(float(v), 6) for v in q]


def main():
    ap = argparse.ArgumentParser(description="SPC 可编辑报告生成器")
    ap.add_argument("--result", required=True, help="spc_engine 输出的 JSON")
    ap.add_argument("--info", default=None, help="抬头/条件信息 JSON（可选）")
    ap.add_argument("--template", default=DEFAULT_TEMPLATE)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    with open(args.result, "r", encoding="utf-8") as f:
        result = json.load(f)
    info = {}
    if args.info and os.path.exists(args.info):
        with open(args.info, "r", encoding="utf-8") as f:
            info = json.load(f)

    payload = {
        "report_id": time.strftime("%Y%m%d_%H%M%S"),
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "result": result,
        "info": info,
        "pdf_curve": pdf_curve(result),
        "prob_theo": prob_theo(result),
    }

    with open(args.template, "r", encoding="utf-8") as f:
        tpl = f.read()
    if "__PAYLOAD__" not in tpl:
        sys.exit("模板缺少 __PAYLOAD__ 占位符")
    # </script> 防注入转义
    html = tpl.replace("__PAYLOAD__",
                       json.dumps(payload, ensure_ascii=False).replace("</", "<\\/"))

    out = args.out or os.path.join(os.getcwd(), f"SPC报告_{payload['report_id']}.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[report_builder] 报告已生成: {out}")
    print("  功能: 编辑态修改抬头/数据(改数即重算) | localStorage草稿 | 导出JSON/另存HTML | 导出PDF(打印)")


if __name__ == "__main__":
    main()
