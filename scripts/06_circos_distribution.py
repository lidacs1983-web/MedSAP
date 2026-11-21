# -*- coding: utf-8 -*-
"""地菍基因组 MedSAP 家族分布与共线性圈图（circos 风格，极坐标原生绘制）
外圈：12 条假染色体 + MedSAP 成员标记（红=家族内共线/WGD/片段重复，黑=无共线/分散重复）
中间两环：基因密度（外=1 Mb 窗口热图，内=1 Mb 窗口柱状图）
连线：灰色=假染色体间共线性背景，红色=MedSAP 成员间共线关系（5 对，标 Ka/Ks）
"""
import json, sys, math, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, PathPatch, Patch
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib import colors as mcolors

sys.stdout.reconfigure(encoding="utf-8")
BASE = "C:/Users/lidac/WorkBuddy/2026-09-29-20-40-35"

plt.rcParams.update({
    "font.family": ["Arial", "Microsoft YaHei", "SimHei", "DejaVu Sans"],
    "axes.unicode_minus": False,
})

# ---------------- 数据 ----------------
med = json.load(open(f"{BASE}/tmp_synteny/med_genes.json"))
CHRS = [f"GWHBCLA000000{i:02d}" for i in range(1, 13)]
chr_len = {c: max(e for g, (cc, s, e, st) in med.items() if cc == c) for c in CHRS}

dr2sap = {}
for line in open(f"{BASE}/tmp_prop/med_sap_proteins.fa", encoding="utf-8"):
    if line.startswith(">"):
        nm, dr = line[1:].strip().split("|")
        dr2sap[dr] = nm
sap_pos = {dr2sap[dr]: med["Med|" + dr] for dr in dr2sap}

dup = json.load(open(f"{BASE}/tmp_kaks/dup_classification.json"))
wgd_pairs = [(p["g1"], p["g2"], p["KaKs"]) for p in dup["pairs"] if p["类型"] == "WGD/片段重复"]
red_genes = {g for pr in wgd_pairs for g in pr[:2]}

blocks = json.load(open(f"{BASE}/tmp_synteny/blocks_intra_Med.json"))

# 基因密度（1 Mb 窗口）
WIN = 1_000_000
dens, maxdens = {}, 0
for c in CHRS:
    n = int(chr_len[c] // WIN) + 1
    arr = np.zeros(n)
    for g, (cc, s, e, st) in med.items():
        if cc == c:
            for w in range(int(s // WIN), int(e // WIN) + 1):
                if w < n:
                    arr[w] += 1
    dens[c] = arr
    maxdens = max(maxdens, arr.max())

# ---------------- 角度体系（度） ----------------
GAP = 2.4
tot = sum(chr_len.values())
avail = 360 - GAP * 12
span = {c: chr_len[c] / tot * avail for c in CHRS}
cum = 90.0 - (GAP * 6 + span[CHRS[0]] / 2)   # Chr1 中点在正上方
theta = {}
for c in CHRS:
    theta[c] = (cum + GAP / 2, cum + GAP / 2 + span[c])
    cum += GAP + span[c]

def ang(c, pos):
    """基因组位置 -> 弧度"""
    t0, t1 = theta[c]
    return math.radians(t0 + (pos / chr_len[c]) * (t1 - t0))

def arc_poly(t0d, t1d, r0, r1, step=1.5):
    """(theta_deg, r) 数据坐标下的环形扇面多边形"""
    ts = np.arange(t0d, t1d + step / 2, step)
    ts = np.append(ts, t1d) if ts[-1] < t1d else ts
    outer = [(math.radians(t), r1) for t in ts]
    inner = [(math.radians(t), r0) for t in ts[::-1]]
    return outer + inner

# ---------------- 半径体系 ----------------
R_CHR0, R_CHR1 = 1.00, 1.065
R_HM0, R_HM1 = 0.905, 0.970
R_BAR1, R_BAR0 = 0.890, 0.775
R_LINK = 0.752
R_LBL = 1.155

MED_RED, MED_BLK = "#D7263D", "#1A1A1A"
CHR_FC, CHR_EC = "#7FB3D5", "#2E5E8C"
GRAY_LNK = "#B4B8BC"

cmap = matplotlib.colormaps["Blues"]
norm = mcolors.Normalize(vmin=0, vmax=maxdens)

L = {
 "cn": dict(ttl="地菍基因组 MedSAP 基因家族分布与共线性圈图",
            sub="Genome-wide distribution and intra-family colinearity of the MedSAP family in Melastoma dodecandrum",
            ring1="基因密度热图（1 Mb）", ring2="基因密度柱状图（1 Mb）",
            red_lab="红色标记：存在家族内共线性（WGD/片段重复，n=10）",
            blk_lab="黑色标记：无家族内共线性（分散重复，n=12）",
            gray_lab="灰色连线：全基因组共线性背景",
            redlnk="红色连线：MedSAP 成员间共线关系（5 对）",
            out=f"{BASE}/outputs/地菍SAP_基因组分布与共线性圈图.png"),
 "en": dict(ttl="Genomic distribution and colinearity of the MedSAP family",
            sub="Circos view of M. dodecandrum pseudochromosomes; red / black MedSAP = with / without intra-family colinearity",
            ring1="Gene density heatmap (1 Mb)", ring2="Gene density histogram (1 Mb)",
            red_lab="Red MedSAP: intra-family colinear (WGD/segmental, n=10)",
            blk_lab="Black MedSAP: no intra-family colinearity (dispersed, n=12)",
            gray_lab="Gray links: genome-wide colinear background",
            redlnk="Red links: MedSAP-MedSAP colinear pairs (5)",
            out=f"{BASE}/outputs/MedSAP_genome_distribution_circos_EN.png"),
}

n_gray = 0
for lang, T in L.items():
    fig = plt.figure(figsize=(12.4, 12.9))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1.46)
    ax.axis("off")

    # ---- 灰色共线背景连线 ----
    for blk in blocks:
        spanpos = collections.defaultdict(list)
        for anc in blk["anchors"]:
            for gid in anc:
                c, s, e, st = med[gid]
                spanpos[c].append((s + e) / 2)
        cs = [c for c in spanpos if c in theta]
        if len(cs) < 2:
            continue
        c1, c2 = cs[0], cs[1]
        t1, t2 = ang(c1, np.mean(spanpos[c1])), ang(c2, np.mean(spanpos[c2]))
        dt = abs(t2 - t1)
        rc = R_LINK * math.cos(min(dt, math.pi) / 2)
        pth = Path([(t1, R_LINK), (t1, rc), (t2, rc), (t2, R_LINK)],
                   [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        ax.add_patch(PathPatch(pth, fill=False, edgecolor=GRAY_LNK,
                               lw=0.5, alpha=0.30, zorder=2))
        n_gray += 1

    # ---- 红色 MedSAP 成员间连线 ----
    for g1, g2, kk in wgd_pairs:
        c1, s1, e1, _ = sap_pos[g1]
        c2, s2, e2, _ = sap_pos[g2]
        t1, t2 = ang(c1, (s1 + e1) / 2), ang(c2, (s2 + e2) / 2)
        dt = abs(t2 - t1)
        rc = R_LINK * math.cos(min(dt, math.pi) / 2)
        pth = Path([(t1, R_LINK), (t1, rc), (t2, rc), (t2, R_LINK)],
                   [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        ax.add_patch(PathPatch(pth, fill=False, edgecolor=MED_RED,
                               lw=2.3, alpha=0.92, zorder=6))
        if kk:
            tm = (t1 + t2) / 2 + (math.pi if dt > math.pi else 0)
            ax.text(tm, rc * 0.60, f"Ka/Ks={float(kk):.3f}",
                    ha="center", va="center", fontsize=8.4, style="italic",
                    color=MED_RED, zorder=7,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white",
                              ec="none", alpha=0.78))

    # ---- 染色体带 ----
    for i, c in enumerate(CHRS, 1):
        t0, t1 = theta[c]
        ax.add_patch(Polygon(arc_poly(t0, t1, R_CHR0, R_CHR1), closed=True,
                             facecolor=CHR_FC, edgecolor=CHR_EC, lw=0.8, zorder=4))
        # Chr 标签：内嵌带内，切向白字
        mid = (t0 + t1) / 2
        mid = mid % 360
        rot = -mid if mid <= 90 or mid >= 270 else -mid + 180
        ax.text(math.radians(mid), (R_CHR0 + R_CHR1) / 2, f"Chr{i}",
                ha="center", va="center", fontsize=9.2, fontweight="bold",
                color="white", rotation=rot, rotation_mode="anchor", zorder=5)

    # ---- 基因密度热图环 ----
    for c in CHRS:
        arr, t0, t1 = dens[c], theta[c][0], theta[c][1]
        n = len(arr)
        for w in range(n):
            a0 = t0 + w * WIN / chr_len[c] * (t1 - t0)
            a1 = t0 + min((w + 1) * WIN, chr_len[c]) / chr_len[c] * (t1 - t0)
            ax.add_patch(Polygon(arc_poly(a0, a1, R_HM0, R_HM1, step=3),
                                 closed=True, facecolor=cmap(norm(arr[w])),
                                 edgecolor="none", zorder=3))

    # ---- 基因密度柱状环 ----
    for c in CHRS:
        arr, t0, t1 = dens[c], theta[c][0], theta[c][1]
        n = len(arr)
        wds = WIN / chr_len[c] * (t1 - t0)          # 度 -> bar 需弧度
        angs = [math.radians(t0 + (w + 0.5) * WIN / chr_len[c] * (t1 - t0)) for w in range(n)]
        hts = [(R_BAR1 - R_BAR0) * arr[w] / maxdens for w in range(n)]
        ax.bar(angs, hts, width=[math.radians(wds)] * n, bottom=R_BAR0,
               color="#2171B5", edgecolor="none", alpha=0.78, zorder=3)

    # 环分隔细线
    for rr in (R_CHR0 - 0.006, R_HM0 - 0.006, R_BAR1 + 0.004, R_BAR0 - 0.004):
        ax.plot(np.linspace(0, 2 * math.pi, 721), [rr] * 721,
                color="#DDDDDD", lw=0.7, zorder=1)

    # ---- MedSAP 标记与标签（贪心分层防碰撞） ----
    marks = sorted((ang(c, (s + e) / 2), g) for g, (c, s, e, st) in sap_pos.items())
    TEXTW = 0.19           # 标签文字近似角宽（rad）
    STEP, MAXLV = 0.075, 4
    placed = []            # (theta, level)
    levels = []
    for t, g in marks:
        lv = 0
        while lv < MAXLV and any(abs(t - tp) < TEXTW and lp == lv for tp, lp in placed):
            lv += 1
        placed.append((t, lv))
        levels.append(lv)
    for (t, g), lv in zip(marks, levels):
        col = MED_RED if g in red_genes else MED_BLK
        rt = R_LBL + lv * STEP
        ax.plot([t, t], [R_CHR1, R_CHR1 + 0.030], color=col, lw=1.6, zorder=6)
        ax.plot(t, R_CHR1 + 0.040, "o", ms=5.4, color=col, zorder=7,
                markeredgecolor="white", markeredgewidth=0.7)
        ax.plot([t, t], [R_CHR1 + 0.056, rt - 0.014], color=col, lw=0.6,
                alpha=0.7, zorder=5)
        # 按所在象限选择文字对齐，避免圆点/引线压住自己的标签
        st_, ct_ = math.sin(t), math.cos(t)
        if st_ > 0.06:            # 右半：文字向右外延伸
            ha, va = "left", "center"
        elif st_ < -0.06:         # 左半：文字向左外延伸
            ha, va = "right", "center"
        elif ct_ > 0:             # 正上：文字在锚点上方
            ha, va = "center", "bottom"
        else:                     # 正下
            ha, va = "center", "top"
        ax.text(t, rt, g, ha=ha, va=va, fontsize=8.4,
                fontweight="bold" if g in red_genes else "normal",
                color=col, zorder=8)

    # ---- 环注释 ----
    ax.text(math.radians(197), (R_HM0 + R_HM1) / 2, T["ring1"],
            ha="center", va="center", fontsize=8.6, color="#1C4E80",
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#BBBBBB", lw=0.6, alpha=0.92))
    ax.text(math.radians(197), (R_BAR0 + R_BAR1) / 2, T["ring2"],
            ha="center", va="center", fontsize=8.6, color="#1C4E80",
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#BBBBBB", lw=0.6, alpha=0.92))

    # ---- 标题与图例 ----
    fig.text(0.5, 0.978, T["ttl"], ha="center", fontsize=17,
             fontweight="bold", color="#1F3A5F", fontfamily="Microsoft YaHei")
    fig.text(0.5, 0.9565, T["sub"], ha="center", fontsize=10,
             style="italic", color="#666666")
    handles = [
        Line2D([], [], marker="o", ls="none", ms=7, mfc=MED_RED, mec="white", label=T["red_lab"]),
        Line2D([], [], marker="o", ls="none", ms=7, mfc=MED_BLK, mec="white", label=T["blk_lab"]),
        Line2D([], [], color=GRAY_LNK, lw=1.4, alpha=0.8, label=T["gray_lab"]),
        Line2D([], [], color=MED_RED, lw=2.4, label=T["redlnk"]),
        Patch(facecolor=cmap(0.75), label=T["ring1"]),
        Patch(facecolor="#2171B5", alpha=0.78, label=T["ring2"]),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.058),
              ncol=3, fontsize=9.6, frameon=True, framealpha=0.95,
              borderpad=0.7, columnspacing=1.6, handlelength=1.6)

    fig.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.085)
    fig.savefig(T["out"], dpi=300, facecolor="white")
    plt.close(fig)
    print("saved:", T["out"], " gray links drawn:", n_gray)
