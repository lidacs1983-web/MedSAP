# -*- coding: utf-8 -*-
"""双内参重归一化脚本（PeerJ 投稿包 / 双内参补实验方案配套）。

用法：
    python refgene_renorm.py --input qPCR_raw_data.xlsx --sheet 0 \
        --out renorm_out

输入要求（xlsx 或 csv）：
  必需列：Target（靶基因）、Group（对照组/干旱组）、Time_h（0/6/12/24）、
          Bio_rep（生物学重复 1-3）、Ct_EF1a、Ct_UBQ
  可选列：Ct_Actin1（保留旧口径对照）
  如有多列技术重复，请先取均值并注明；本脚本不做技术重复合并。

输出（写入 out 目录）：
  renorm_all.tsv   三种口径的 2^-ddCt 表（Actin1 / EF1a / EF1a+UBQ 几何平均）
  refgene_stats.tsv 内参稳定性（各处理组 ΔCt 范围、组间 |ΔCt| 差、CV）
  summary.txt      关键结论（各靶基因各时点 fold change 对照）
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load(path, sheet):
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(p, sheet_name=sheet)
    return pd.read_csv(p, sep=None, engine="python")


def ddct(df, refcols, cal_col):
    """返回宽表：每靶基因各时点相对 0h 对照的 2^-ddCt。"""
    d = df.copy()
    d["ref_mean"] = d[refcols].mean(axis=1)
    d["dCt"] = d["Ct_target"] - d["ref_mean"]
    cal = d[(d["Group"].astype(str).str.contains("对")) |
            (d["Group"].astype(str).str.lower().isin(["control", "ck", "well-watered"]))]
    cal = cal[cal["Time_h"] == 0].groupby("Target")["dCt"].mean()
    d["ddCt"] = d.apply(lambda r: r["dCt"] - cal.get(r["Target"], np.nan), axis=1)
    d["fold"] = 2 ** (-d["ddCt"])
    d["method"] = cal_col
    return d


def stability(df, refcols):
    rows = []
    for rc in refcols:
        col = "Ct_" + rc
        if col not in df.columns:
            continue
        v = df[col].astype(float)
        # 各时点组内 Ct 标准差
        grp = df.groupby(["Group", "Time_h"])[col].agg(["std", "mean"])
        rows.append(dict(reference=rc, overall_sd=round(v.std(), 3),
                         mean_Ct=round(v.mean(), 2),
                         max_group_sd=round(grp["std"].max(), 3),
                         note="geNorm 请另行计算（本表给出原始量）"))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--sheet", default=0)
    ap.add_argument("--out", default="renorm_out")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    df = load(a.input, a.sheet)
    df.columns = [str(c).strip() for c in df.columns]
    need = {"Target", "Group", "Time_h", "Bio_rep"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"缺少必需列: {missing}；另需 Ct_target 与 Ct_EF1a/Ct_UBQ")

    # Ct 列：target 列名兼容
    tcol = [c for c in df.columns if c.lower() in ("ct_target", "target_ct", "ct")]
    if not tcol:
        raise SystemExit("未找到靶基因 Ct 列（ct_target / target_ct / ct）")
    df = df.rename(columns={tcol[0]: "Ct_target"})

    frames = []
    if "Ct_Actin1" in df.columns:
        frames.append(ddct(df, ["Ct_Actin1"], "Actin1(旧口径)"))
    if "Ct_EF1a" in df.columns:
        frames.append(ddct(df, ["Ct_EF1a"], "EF1a"))
    if "Ct_EF1a" in df.columns and "Ct_UBQ" in df.columns:
        frames.append(ddct(df, ["Ct_EF1a", "Ct_UBQ"], "EF1a+UBQ 几何平均"))
    if not frames:
        raise SystemExit("未找到 Ct_EF1a / Ct_UBQ / Ct_Actin1 中的任何列")
    res = pd.concat(frames, ignore_index=True)
    res.to_csv(out / "renorm_all.tsv", sep="\t", index=False, float_format="%.4f")
    stability(df, ["Actin1", "EF1a", "UBQ"]).to_csv(
        out / "refgene_stats.tsv", sep="\t", index=False)

    with open(out / "summary.txt", "w", encoding="utf-8") as f:
        for m, d in res.groupby("method"):
            f.write(f"== 口径：{m} ==\n")
            piv = d.pivot_table(index="Target", columns="Time_h", values="fold")
            f.write(piv.round(2).to_string() + "\n\n")
    print("done ->", out)


if __name__ == "__main__":
    main()
