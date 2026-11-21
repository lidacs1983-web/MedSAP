# -*- coding: utf-8 -*-
"""地菍 SAP 跨物种 Ka/Ks —— 出图 + 结果表（中文版）"""
import json, statistics, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = "C:/Users/lidac/WorkBuddy/2026-09-29-20-40-35"
R = json.load(open(BASE + "/tmp_kaks/sap_kaks_cross.json"))

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

SP = {
    "Mcan": dict(cn="野牡丹", c="#C0392B", m="o", z=6, mew=1.0, s=78),
    "Egr":  dict(cn="巨桉",   c="#2E86AB", m="^", z=4, mew=0.8, s=52),
    "Vvi":  dict(cn="葡萄",   c="#6A4C93", m="s", z=4, mew=0.8, s=46),
    "Vco":  dict(cn="蓝莓",   c="#3B7A57", m="D", z=4, mew=0.8, s=46),
    "Csi":  dict(cn="茶树",   c="#B8860B", m="P", z=4, mew=0.8, s=50),
    "Ngr":  dict(cn="显齿蛇葡萄", c="#1B998B", m="X", z=4, mew=0.8, s=52),
    "Nta":  dict(cn="烟草",   c="#E07A5F", m="v", z=4, mew=0.8, s=48),
    "At":   dict(cn="拟南芥", c="#7F8C8D", m="*", z=4, mew=0.8, s=70),
}
ORDER = ["Mcan", "Ngr", "Egr", "Vvi", "Vco", "Csi", "Nta", "At"]

fig = plt.figure(figsize=(13.4, 6.6), dpi=200)
fig.patch.set_facecolor("white")

# ---------------- Panel A: Ka–Ks 散点（log-log） ----------------
axA = fig.add_axes([0.055, 0.170, 0.40, 0.660])
LO, HI = 1.2e-3, 40.0
axA.set_xscale("log"); axA.set_yscale("log")
axA.set_xlim(LO, HI); axA.set_ylim(1.5e-3, 3.0)
axA.axvspan(2.0, HI, color="#F2F2F2", zorder=0)
axA.text(9.0, 2.5, "Ks > 2：同义位点饱和\n（Ka/Ks 不可估）", fontsize=8.2,
         color="#888888", ha="center", va="top", linespacing=1.5)
axA.plot([LO, HI], [LO, HI], ls="--", lw=1.0, color="#B0B0B0", zorder=1)
axA.plot([LO, HI], [0.5 * LO, 0.5 * HI], ls=":", lw=1.0, color="#D0D0D0", zorder=1)
axA.text(0.11, 0.165, "Ka = Ks（中性）", fontsize=7.8, color="#999999",
         rotation=42, ha="center", va="center")
axA.text(0.30, 0.165, "Ka = 0.5Ks", fontsize=7.6, color="#BBBBBB",
         rotation=42, ha="center", va="center")

for sp in ORDER:
    rs = [r for r in R if r["sp"] == sp and r["ds"] is not None
          and math.isfinite(r["ds"]) and r["ds"] > 0 and r["dn"] > 0]
    if not rs:
        continue
    cfg = SP[sp]
    gen = (sp == "Mcan")
    axA.scatter([r["ds"] for r in rs], [r["dn"] for r in rs],
                s=cfg["s"], marker=cfg["m"], color=cfg["c"],
                edgecolors=("#7B241C" if gen else "white"),
                linewidths=(cfg["mew"] if gen else 0.7),
                alpha=(1.0 if gen else 0.85), zorder=cfg["z"],
                label=f"{cfg['cn']}（n={len(rs)}）")

# 标注几个属级代表点
ann = {"MedSAP8": (-38, 22), "MedSAP12": (34, 8), "MedSAP1": (28, 20),
       "MedSAP21": (30, 10)}
for r in R:
    if r["sp"] == "Mcan" and r["med_name"] in ann:
        dx, dy = ann[r["med_name"]]
        axA.annotate(r["med_name"], xy=(r["ds"], r["dn"]), xytext=(dx, dy),
                     textcoords="offset points", fontsize=7.0, color="#7B241C",
                     ha="center",
                     arrowprops=dict(arrowstyle="-", color="#C0392B", lw=0.5,
                                     shrinkA=0, shrinkB=2))

axA.set_xlabel("Ks（同义替换率, dS）", fontsize=11)
axA.set_ylabel("Ka（非同义替换率, dN）", fontsize=11)
axA.set_title("A  MedSAP 与近缘物种同源基因对的 Ka–Ks 分布", fontsize=12.5,
              fontweight="bold", pad=9, color="#1F3A5F")
h, l = axA.get_legend_handles_labels()
order2 = sorted(range(len(l)), key=lambda i: ORDER.index([k for k in ORDER if SP[k]['cn'] in l[i]][0]))
axA.legend([h[i] for i in order2], [l[i] for i in order2], frameon=False,
           loc="lower right", prop={"family": "Microsoft YaHei", "size": 8.2},
           handletextpad=0.4, labelspacing=0.35)
axA.tick_params(labelsize=9)
for s in ("top", "right"):
    axA.spines[s].set_visible(False)

# ---------------- Panel B: 各物种 Ka/Ks（可靠配对） ----------------
axB = fig.add_axes([0.585, 0.170, 0.385, 0.660])
rows = []
for sp in ORDER:
    rs = [r for r in R if r["sp"] == sp and r.get("reliable")]
    if rs:
        rows.append((sp, sorted(rs, key=lambda r: r["ratio"]), False))
    else:
        allr = [r for r in R if r["sp"] == sp]
        if allr and all(r.get("saturated") for r in allr):
            rows.append((sp, allr, True))   # 全部饱和：占位行
rows = rows[::-1]     # Mcan 在顶部，全部饱和的物种置于底部

ymap = {}
y = 0.0
for sp, rs, _ in rows:
    ymap[sp] = y
    y += 1.0
ytop = y

for sp, rs, is_sat in rows:
    cfg = SP[sp]; yy = ymap[sp]
    gen = (sp == "Mcan")
    if is_sat:
        axB.plot([0, 1.45], [yy, yy], color="#DCDCDC", lw=1.3,
                 ls=(0, (4, 3)), zorder=1)
        axB.text(0.60, yy, "全部配对饱和（Ks > 2），Ka/Ks 不可估",
                 fontsize=7.8, color="#909090", ha="center", va="center",
                 zorder=2)
        continue
    for r in rs:
        axB.plot([0, r["ratio"]], [yy, yy], color=cfg["c"], lw=1.0,
                 alpha=0.18 if gen else 0.10, zorder=1)
    jit = np.linspace(-0.22, 0.22, len(rs)) if len(rs) > 1 else [0]
    axB.scatter([r["ratio"] for r in rs], [yy + j for j in jit],
                s=(66 if gen else 42), color=cfg["c"],
                edgecolors=("#7B241C" if gen else "white"),
                linewidths=(0.9 if gen else 0.6), zorder=3)
    med = statistics.median([r["ratio"] for r in rs])
    axB.plot([med, med], [yy - 0.30, yy + 0.30], color="#111111", lw=1.8, zorder=4)
    axB.annotate(f"{med:.2f}", xy=(med, yy - 0.34), fontsize=8.0, color="#111111",
                 ha="center", va="top", fontweight="bold")

axB.axvline(1.0, color="#C0392B", ls="--", lw=1.1, zorder=1)
axB.text(1.03, ytop - 0.45, "Ka/Ks = 1（中性）", fontsize=8.0, color="#C0392B",
         va="top")
axB.set_yticks([ymap[sp] for sp, _, _ in rows])
axB.set_yticklabels(
    [(f"{SP[sp]['cn']}  (n={len(rs)}，全部饱和)" if sat
      else f"{SP[sp]['cn']}  (n={len(rs)})") for sp, rs, sat in rows],
    fontsize=9.2)
for lab, (sp, rs, sat) in zip(axB.get_yticklabels(), rows):
    lab.set_color("#A3A3A3" if sat else "#1A1A1A")
axB.set_ylim(-0.62, ytop - 0.28)
axB.set_xlim(0, 1.45)
axB.set_xlabel("Ka/Ks（dN/dS）", fontsize=11)
axB.set_title("B  各物种同源配对的 Ka/Ks（黑线 = 中位数）", fontsize=12.5,
              fontweight="bold", pad=9, color="#1F3A5F")
axB.tick_params(labelsize=9)
for s in ("top", "right"):
    axB.spines[s].set_visible(False)

fig.text(0.055, 0.013,
         "注：同源配对 = 地菍 22 个 MedSAP 与各物种 SAP 家族成员间的互惠最佳命中（RBH）；Ka/Ks 主算法为 ML（ML 失败时回退 NG86）。\n"
         "属级配对（地菍 vs 野牡丹，同属 Melastoma）Ks 仅 0.014–0.056、未饱和，Ka/Ks 可靠；烟草（8 对）与拟南芥（6 对）的全部配对同义位点已饱和（Ks > 2），\n"
         "Ka/Ks 无法估计，故在图 B 中保留物种行并以灰虚线占位（A、B 面板物种一致），其配对仍全部显示于图 A 的饱和区（灰底）。",
         fontsize=7.6, color="#555555", ha="left", va="bottom", linespacing=1.5)

out = BASE + "/outputs/地菍SAP_跨物种KaKs分析图.png"
fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
print("saved", out)

# ---------------- 结果表 ----------------
def fmt(v, n=4):
    return "" if v is None or (isinstance(v, float) and not math.isfinite(v)) else f"{v:.{n}f}"

hdr = ["物种", "物种中文", "配对类型", "MedSAP", "对方同源基因", "蛋白一致性(%)",
       "比对长度(aa)", "主算法", "Ka", "Ks", "Ka/Ks", "Ka(ML)", "Ka(NG86)",
       "Ks(ML)", "Ks(NG86)", "可靠性"]
lines = ["\t".join(hdr)]
for r in sorted(R, key=lambda r: (ORDER.index(r["sp"]),
                                  int(r["med_name"][6:]) if r["med_name"][6:].isdigit() else 999)):
    typ = ("属级(同属)" if r["genus"] else ("参照(拟南芥)" if r["ref"] else "跨属/科"))
    rel = "饱和" if r["saturated"] else ("可靠" if r.get("reliable") else "Ks过小")
    lines.append("\t".join([
        r["sp"], r["cn"], typ, r["med_name"], r["sp_gene"],
        f"{r['ident']*100:.1f}", str(r["aln_len"]), r["method"],
        fmt(r["dn"]), fmt(r["ds"]),
        fmt(r["ratio"], 3), fmt(r["dn_ml"]), fmt(r["dn_ng"]),
        fmt(r["ds_ml"]), fmt(r["ds_ng"]), rel]))
open(BASE + "/outputs/地菍SAP_跨物种KaKs分析结果.tsv", "w", encoding="utf-8-sig").write("\n".join(lines) + "\n")

# 属级专表
gl = ["\t".join(["MedSAP", "野牡丹同源基因", "蛋白一致性(%)", "比对长度(aa)", "主算法",
                 "Ka", "Ks", "Ka/Ks", "可靠性"])]
for r in sorted([r for r in R if r["sp"] == "Mcan"],
                key=lambda r: int(r["med_name"][6:])):
    rel = "饱和" if r["saturated"] else ("可靠" if r.get("reliable") else "Ks过小(参考)")
    gl.append("\t".join([r["med_name"], r["sp_gene"], f"{r['ident']*100:.1f}",
                         str(r["aln_len"]), r["method"], fmt(r["dn"]), fmt(r["ds"]),
                         fmt(r["ratio"], 3), rel]))
open(BASE + "/outputs/地菍SAP_属级KaKs_地菍vs野牡丹.tsv", "w", encoding="utf-8-sig").write("\n".join(gl) + "\n")
print("tables written")

# ---------------- 摘要 ----------------
print("\n--- 属级（Med–Mcan）可靠配对 ---")
ok = [r for r in R if r["sp"] == "Mcan" and r.get("reliable")]
rr = [r["ratio"] for r in ok]
print(f"n={len(ok)}  Ka/Ks {min(rr):.3f}–{max(rr):.3f}  中位 {statistics.median(rr):.3f}  "
      f"均值 {statistics.mean(rr):.3f}")
print(f"Ks {min(r['ds'] for r in ok):.4f}–{max(r['ds'] for r in ok):.4f}  "
      f"Ka {min(r['dn'] for r in ok):.4f}–{max(r['dn'] for r in ok):.4f}")
print("> 1 的配对:", [f"{r['med_name']}({r['ratio']:.2f})" for r in ok if r["ratio"] > 1])
