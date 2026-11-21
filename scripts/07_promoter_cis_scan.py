# -*- coding: utf-8 -*-
"""地菍 MedSAP 家族启动子顺式作用元件分布分析（专项）

口径：启动子 = TSS 上游 2000 bp（负链取下游反向互补，坐标以 TSS 为 0，仅正值）
输出：
  1. outputs/地菍SAP_启动子顺式作用元件图.png        分布图 + 类别组成 + 元件排名
  2. outputs/地菍SAP_启动子顺式作用元件_计数矩阵.tsv  22 x 27
  3. outputs/地菍SAP_启动子顺式作用元件_位置明细.tsv  长表（含 TSS 上游 bp）
  4. outputs/地菍SAP_启动子顺式作用元件_类别统计.tsv  类别/元件汇总
"""
import os, re, json
from collections import defaultdict
from Bio import Phylo, SeqIO
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

BASE = "C:/Users/lidac/WorkBuddy/2026-09-29-20-40-35"
GST  = os.path.join(BASE, "tmp_gstruct")
PHY  = os.path.join(BASE, "tmp_phylo")
OUTD = os.path.join(BASE, "outputs")
GENOME = r"D:/博士资料/13.M. dodecandrum,基因组文件/GWHBCLA00000000.genome.fasta"
UP = 2000

# ---------- PlantCARE 一致性序列 ----------
IUPAC = {"R": "AG", "Y": "CT", "W": "AT", "S": "CG", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def iupac_re(p):
    return "".join("[%s]" % IUPAC[c] if c in IUPAC else c for c in p)


ELEMS = [
    ("Box 4", "ATTAAT", "light"), ("G-Box", "CACGTG", "light"),
    ("GT1-motif", "GGTTAA", "light"), ("TCT-motif", "TCTTAC", "light"),
    ("Sp1", "GGGCGG", "light"), ("ACE", "ACGTG", "light"),
    ("AE-box", "AGAAACAA", "light"), ("TCCC-motif", "TCTCCCT", "light"),
    ("ARE", "AAACAAA", "stress"), ("GC-motif", "CCCCCG", "stress"),
    ("LTR", "CCGAAA", "stress"), ("MBS", "CAACTG", "stress"),
    ("TC-rich repeats", "GTTTTCTTAC", "stress"), ("W box", "TTGAC[CT]", "stress"),
    ("ABRE", "ACGTG", "hormone"), ("CGTCA-motif", "CGTCA", "hormone"),
    ("TGACG-motif", "TGACG", "hormone"), ("TCA-element", "CCATCTTTTT", "hormone"),
    ("TGA-element", "AACGAC", "hormone"), ("GARE-motif", "AAACAGA", "hormone"),
    ("P-box", "CCTTTTG", "hormone"), ("TATC-box", "TATCCCA", "hormone"),
    ("ERE", "ATTTCAAA", "hormone"), ("CAT-box", "GCCACT", "growth"),
    ("O2-site", "GATGATATGG", "growth"), ("Skn-1-like", "GTCAT", "growth"),
    ("circadian", "CAANNNNATC", "growth"),
]
ELEM_NAME = [e[0] for e in ELEMS]
ELEM_CAT  = {e[0]: e[2] for e in ELEMS}
PATTERNS  = {e[0]: re.compile(iupac_re(e[1])) for e in ELEMS}
CAT_ZH = {"light": "光响应", "stress": "胁迫相关", "hormone": "激素响应", "growth": "生长发育"}
CAT_ORDER = ["light", "stress", "hormone", "growth"]
CAT_COLOR = {"light": "#E8A33D", "stress": "#C94F4F",
             "hormone": "#5B8FD6", "growth": "#5FA86A"}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
revcomp = lambda s: s.translate(COMP)[::-1]

# ---------- 取 SAP 基因与树序 ----------
models = json.load(open(os.path.join(GST, "med_sap_models.json")))
meta   = json.load(open(os.path.join(PHY, "sap_all9_meta.json")))
n2o    = {m["new"]: m["orig"] for m in meta if m["species"] == "Med"}
t = Phylo.read(os.path.join(PHY, "sap_all9.nwk"), "newick"); t.ladderize()
ORDER = [lf.name.split("|")[0] for lf in t.get_terminals()
         if lf.name.split("|")[0].startswith("MedSAP")]
print("树序(%d):" % len(ORDER), ORDER)

gidx = SeqIO.index(GENOME, "fasta")


def upstream(gid):
    m = models[gid]
    gs, ge, strand = m["gene"]
    rec = gidx[m["chr"]]
    L = len(rec.seq)
    if strand == "+":
        a, b = max(0, gs - 1 - UP), gs - 1
        s = str(rec.seq[a:b]).upper()
    else:
        a, b = ge, min(L, ge + UP)
        s = revcomp(str(rec.seq[a:b]).upper())
    if len(s) < UP:
        s = "N" * (UP - len(s)) + s
    return s


results, gids = {}, {}
for new in ORDER:
    gid = n2o[new]
    gids[new] = gid
    s = upstream(gid)
    hits = {}
    for name, pat in PATTERNS.items():
        hs = [(mm.start(), mm.end()) for mm in pat.finditer(s)]
        if hs:
            hits[name] = hs
    results[new] = hits
    print("%-9s %-9s elements=%d" % (new, gid, sum(len(v) for v in hits.values())))

# ---------- 统计 ----------
tot_elem = sum(sum(len(v) for v in results[n].values()) for n in ORDER)
elem_cnt = defaultdict(int)
for n in ORDER:
    for name, hs in results[n].items():
        elem_cnt[name] += len(hs)
print("总元件数 = %d，均值 %.1f/基因" % (tot_elem, tot_elem / len(ORDER)))

# ---------- 输出 1：计数矩阵 ----------
with open(os.path.join(OUTD, u"地菍SAP_启动子顺式作用元件_计数矩阵.tsv"), "w",
          encoding="utf-8-sig", newline="") as f:
    f.write("蛋白\tGeneID\t" + "\t".join(ELEM_NAME) + "\t总数\n")
    for new in ORDER:
        cnt = [len(results[new].get(nm, [])) for nm in ELEM_NAME]
        f.write("%s\t%s\t%s\t%d\n" % (new, gids[new], "\t".join(map(str, cnt)), sum(cnt)))

# ---------- 输出 2：位置明细 ----------
with open(os.path.join(OUTD, u"地菍SAP_启动子顺式作用元件_位置明细.tsv"), "w",
          encoding="utf-8-sig", newline="") as f:
    f.write("蛋白\tGeneID\t元件\t类别\t元件数\tTSS上游位置(bp,0=转录起始位点)\t距TSS最近(bp)\n")
    for new in ORDER:
        for name in ELEM_NAME:
            hs = results[new].get(name)
            if not hs:
                continue
            pos = ";".join("%d-%d" % h for h in hs)
            nearest = min(UP - b for _, b in hs)
            f.write("%s\t%s\t%s\t%s\t%d\t%s\t%d\n" %
                    (new, gids[new], name, ELEM_CAT[name], len(hs), pos, nearest))

# ---------- 输出 3：类别统计 ----------
with open(os.path.join(OUTD, u"地菍SAP_启动子顺式作用元件_类别统计.tsv"), "w",
          encoding="utf-8-sig", newline="") as f:
    f.write("类别\t元件\t总数量\t占全部(%)\t基因数(含该元件)\n")
    for c in CAT_ORDER:
        for nm in [e[0] for e in ELEMS if e[2] == c]:
            v = elem_cnt.get(nm, 0)
            if v == 0:
                continue
            ng = sum(1 for n in ORDER if nm in results[n])
            f.write("%s\t%s\t%d\t%.1f\t%d\n" % (CAT_ZH[c], nm, v, 100.0 * v / tot_elem, ng))
    for c in CAT_ORDER:
        names = [e[0] for e in ELEMS if e[2] == c]
        sub = sum(elem_cnt.get(nm, 0) for nm in names)
        f.write("%s（合计）\t-\t%d\t%.1f\t-\n" % (CAT_ZH[c], sub, 100.0 * sub / tot_elem))

# ================= 绘图 =================
N = len(ORDER)
fig = plt.figure(figsize=(12.8, 10.4))

fig.text(0.5, 0.972, u"地菍 MedSAP 基因家族启动子顺式作用元件分布",
         fontsize=15, fontweight="bold", color="#111111", ha="center", va="center")
fig.text(0.5, 0.944,
         u"启动子 = TSS 上游 2000 bp；PlantCARE 一致性序列扫描，27 种元件 / 4 大类；"
         u"基因按系统发育树顺序排列；★ = RT-PCR 验证基因",
         fontsize=8.8, color="#555555", ha="center", va="center")

# ---- 主图：x = TSS 上游距离（0 在右 = TSS，2000 在左）
LAX, RAX = 0.150, 0.760
ax = fig.add_axes([LAX, 0.255, RAX - LAX, 0.615])
ax.set_xlim(2060, -80)                       # 反转：上游 2000 在左，TSS 在右
ax.set_ylim(N - 0.5, -1.15)
ax.set_facecolor("white")

for i in range(N):
    if i % 2 == 1:
        ax.axhspan(i - 0.5, i + 0.5, color="#F3F5F7", zorder=0)

# 每行 4 条类别基线（自上而下：光响应/胁迫/激素/生长发育），重叠错层
LANE_H, LANE_PITCH = 0.108, 0.135
for i, new in enumerate(ORDER):
    for ci, c in enumerate(CAT_ORDER):
        ybase = i - 0.28 + ci * LANE_PITCH
        lane = []
        for nm in [e[0] for e in ELEMS if e[2] == c]:
            for a, b in results[new].get(nm, []):
                lane.append((a, b))
        lane.sort()
        rows = []                             # 贪心车道分配
        for a, b in lane:
            placed = False
            for r in range(len(rows)):
                if a > rows[r][-1][1] + 4:
                    rows[r].append((a, b)); placed = True; break
            if not placed:
                rows.append([(a, b)])
        for ri, seg in enumerate(rows):
            y = ybase + ri * 0.075
            for a, b in seg:
                ax.add_patch(Rectangle((UP - b, y - 0.048), b - a, 0.096,
                                       facecolor=CAT_COLOR[c], edgecolor="none",
                                       alpha=0.92, zorder=3))

# 基因名（左缘外，两行：名字 + GeneID）+ RT-PCR 红星
for i, new in enumerate(ORDER):
    ax.text(2120, i - 0.08, new, ha="right", va="center", fontsize=8.6,
            color="#1A1A1A", fontweight="bold", transform=ax.transData)
    ax.text(2120, i + 0.24, gids[new], ha="right", va="center",
            fontsize=6.6, color="#8A8A8A", transform=ax.transData)
    if new in ("MedSAP8", "MedSAP16"):
        ax.text(2113, i - 0.08, u"★", ha="left", va="center", fontsize=10,
                color="#D0342C", transform=ax.transData)

# TSS 线与坐标轴
ax.axvline(0, color="#333333", lw=1.4, zorder=4)
ax.set_xticks([0, 500, 1000, 1500, 2000])
ax.set_xticklabels(["0\n(TSS)", "500", "1000", "1500", "2000"], fontsize=9.0)
ax.set_xlabel(u"TSS 上游距离 (bp)", fontsize=10.6, labelpad=5)
ax.xaxis.set_label_coords(0.28, -0.085)          # 左移，避开右下面板标题
ax.set_yticks([])
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_linewidth(0.9)
ax.tick_params(axis="x", length=3, width=0.8)

# ---- 左下：类别组成
axc = fig.add_axes([LAX, 0.095, 0.27, 0.105])
left = 0.0
for c in CAT_ORDER:
    names = [e[0] for e in ELEMS if e[2] == c]
    sub = sum(elem_cnt.get(nm, 0) for nm in names)
    axc.barh(0, sub, left=left, height=0.6, color=CAT_COLOR[c],
             edgecolor="white", lw=0.8)
    if sub / tot_elem > 0.04:
        axc.text(left + sub / 2, 0, "%.0f%%" % (100.0 * sub / tot_elem),
                 ha="center", va="center", fontsize=8.4, color="white",
                 fontweight="bold")
    left += sub
axc.set_xlim(0, tot_elem); axc.set_ylim(-0.5, 0.5)
axc.set_yticks([]); axc.set_xticks([0, 100, 200, 300, 400])
axc.tick_params(labelsize=7.6)
axc.set_xlabel(u"元件总数（%d 个）" % tot_elem, fontsize=8.6, labelpad=1)
for s in ("top", "right", "left"):
    axc.spines[s].set_visible(False)
axc.set_title(u"四大类别组成", fontsize=9.8, fontweight="bold", loc="left", pad=4)

# ---- 右下：Top12 元件
axe = fig.add_axes([0.485, 0.062, 0.42, 0.158])
top = sorted(elem_cnt.items(), key=lambda kv: -kv[1])[:12]
names = [k for k, _ in top][::-1]
vals = [v for _, v in top][::-1]
cols = [CAT_COLOR[ELEM_CAT[k]] for k in names]
axe.barh(range(len(names)), vals, color=cols, edgecolor="none", height=0.72)
axe.set_yticks(range(len(names)))
axe.set_yticklabels(names, fontsize=7.4)
axe.tick_params(axis="x", labelsize=7.6)
for i, v in enumerate(vals):
    axe.text(v + 0.8, i, str(v), va="center", fontsize=7.2, color="#444444")
axe.set_xlim(0, max(vals) * 1.13)
for s in ("top", "right", "left"):
    axe.spines[s].set_visible(False)
axe.set_title(u"各元件出现次数（前 12）", fontsize=9.6, fontweight="bold",
              loc="left", pad=5)

# ---- 图例（右上，独立区域）
handles = [Patch(facecolor=CAT_COLOR[c], label=CAT_ZH[c]) for c in CAT_ORDER]
handles.append(plt.Line2D([0], [0], marker="*", color="#D0342C", linestyle="",
                          markersize=11, label=u"RT-PCR 验证"))
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.775, 0.885),
           frameon=False, fontsize=9.2, title=u"元件类别", title_fontsize=9.6)

foot = (u"注：启动子定义为各基因转录起始位点 (TSS) 上游 2000 bp（负链基因取下游 2000 bp 并反向互补）；"
        u"元件按 PlantCARE 一致性序列扫描（ACE 与 ABRE 共用 ACGTG 核心，同一位点各自计一次）。\n"
        u"横坐标为 TSS 上游距离：0 = 转录起始位点，2000 = 上游 2000 bp；"
        u"每行 4 条色带自上而下依次为光响应 / 胁迫相关 / 激素响应 / 生长发育，同一类别内相邻元件自动错层显示。")
fig.text(LAX, 0.014, foot, fontsize=7.2, color="#777777", ha="left",
         va="bottom", linespacing=1.7)

out_png = os.path.join(OUTD, u"地菍SAP_启动子顺式作用元件图.png")
fig.savefig(out_png, dpi=300, facecolor="white")
plt.close(fig)
print("saved:", out_png)
