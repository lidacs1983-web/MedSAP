# -*- coding: utf-8 -*-
"""地菍 MedSAP 基因重复事件类型分析（五分类：WGD/片段、串联、近距离、分散、单拷贝）
判据（方法学说明将写入论文 Methods）：
  旁系同源对 = 家族内蛋白局部比对 一致性>=40% 且覆盖短序列>=60%（BLOSUM62, gap -11/-1）
  WGD/片段重复：同源对位于基因组内共线性块（blocks_intra_Med.json）中
  串联重复    ：同染色体、两基因间间隔基因数 <= 5
  近距离重复  ：同染色体、间隔基因数 6-25
  分散重复    ：其余同源对（异染色体无共线性 / 同染色体间隔>25 基因）
  单拷贝      ：家族内无同源对
每基因按优先级 WGD/片段 > 串联 > 近距离 > 分散 归类一次。
"""
import csv
import json
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Align

sys.stdout.reconfigure(encoding="utf-8")
BASE = "C:/Users/lidac/WorkBuddy/2026-09-29-20-40-35"

# ---------- 1. 基因位置 ----------
pos = {}
for row in csv.DictReader(open(f"{BASE}/outputs/地菍SAP染色体定位信息.tsv", encoding="utf-8-sig"), delimiter="\t"):
    pos[row["基因名"]] = (row["染色体ID"], int(row["起始(bp)"]), int(row["终止(bp)"]), row["链"])
assert len(pos) == 22

# 全基因组基因序（用于计算间隔基因数）
med_genes = json.load(open(f"{BASE}/tmp_synteny/med_genes.json"))
chr_order = defaultdict(list)          # chr -> [(start,end,geneID)]
for gid, (c, s, e, st) in med_genes.items():
    chr_order[c].append((s, e, gid.split("|")[1]))
for c in chr_order:
    chr_order[c].sort()

def n_intervening(a, b):
    c1, s1, e1, _ = pos[a]
    c2, s2, e2, _ = pos[b]
    if c1 != c2:
        return None
    lo, hi = sorted((min(s1, e1), min(s2, e2)))[0], max(max(s1, e1), max(s2, e2))
    lo, hi = min(s1, e1), max(s1, e1)
    lo2, hi2 = min(s2, e2), max(s2, e2)
    left, right = sorted(((lo, hi), (lo2, hi2)))
    return sum(1 for s, e, g in chr_order[c1] if s > left[1] and e < right[0])

# ---------- 2. 共线性块中的 SAP 旁系对 ----------
dr2sap = {}
for line in open(f"{BASE}/tmp_prop/med_sap_proteins.fa", encoding="utf-8"):
    if line.startswith(">"):
        nm, dr = line[1:].strip().split("|")
        dr2sap[dr] = nm
blocks = json.load(open(f"{BASE}/tmp_synteny/blocks_intra_Med.json"))
block_pairs, block_meta = set(), {}
for blk in blocks:
    ids = {x.split("|")[1] for anc in blk["anchors"] for x in anc}
    hit = [dr2sap[i] for i in ids if i in dr2sap]
    for i in range(len(hit)):
        for j in range(i + 1, len(hit)):
            pr = frozenset((hit[i], hit[j]))
            block_pairs.add(pr)
            block_meta[pr] = (blk["c1"], blk["c2"], len(blk["anchors"]))
print("共线性块内 SAP 旁系对:", [tuple(sorted(p)) for p in sorted(block_pairs, key=sorted)])

# ---------- 3. 蛋白两两比对 ----------
pep, nm, buf = {}, None, []
for line in open(f"{BASE}/tmp_prop/med_sap_proteins.fa", encoding="utf-8"):
    if line.startswith(">"):
        if nm: pep[nm] = "".join(buf)
        nm = line[1:].strip().split("|")[0]; buf = []
    else:
        buf.append(line.strip())
pep[nm] = "".join(buf)

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.mode = "local"

HOMO_ID, HOMO_COV = 0.40, 0.60
genes = sorted(pep)
homolog_pairs = {}
for i in range(len(genes)):
    for j in range(i + 1, len(genes)):
        a, b = genes[i], genes[j]
        aln = aligner.align(pep[a], pep[b])[0]
        if not aln.length:
            continue
        ident = aln.counts().identities / aln.length
        cov = aln.length / min(len(pep[a]), len(pep[b]))
        if ident >= HOMO_ID and cov >= HOMO_COV:
            homolog_pairs[frozenset((a, b))] = (ident, cov, aln.length)
print(f"同源判据 id>={HOMO_ID}, cov>={HOMO_COV} → {len(homolog_pairs)} 对")

# ---------- 4. Ka/Ks ----------
kaks = json.load(open(f"{BASE}/tmp_kaks/sap_kaks.json"))
pair_kaks = {frozenset((x["g1"], x["g2"])): x for x in kaks
             if x["type"].startswith("物种内")}

# ---------- 5. 分类 ----------
ORDER = ["WGD/片段重复", "串联重复", "近距离重复", "分散重复"]
gene_type, gene_partners = defaultdict(set), defaultdict(set)
pair_rows = []
for p in homolog_pairs:
    a, b = sorted(p)
    ident, cov, L = homolog_pairs[p]
    n_iv = n_intervening(a, b)
    if p in block_pairs:
        typ = "WGD/片段重复"
    elif n_iv is not None and n_iv <= 5:
        typ = "串联重复"
    elif n_iv is not None and n_iv <= 25:
        typ = "近距离重复"
    else:
        typ = "分散重复"
    kk = pair_kaks.get(p)
    row = dict(g1=a, g2=b, 类型=typ, 间隔基因数=(n_iv if n_iv is not None else "—"),
               蛋白一致性=round(ident * 100, 1),
               KaKs=(round(kk["ratio"], 3) if kk and kk.get("ratio") else ""),
               Ks=(round(kk["ds"], 4) if kk else ""),
               证据=("基因组内共线性块" if p in block_pairs else "家族内同源比对"))
    pair_rows.append(row)
    gene_type[a].add(typ); gene_type[b].add(typ)
    gene_partners[a].add(b); gene_partners[b].add(a)

gene_final = {g: next((t for t in ORDER if t in gene_type.get(g, ())), "单拷贝")
              for g in sorted(pos)}
cnt = defaultdict(int)
for t in gene_final.values():
    cnt[t] += 1

print("\n每基因分类:")
for g in sorted(pos, key=lambda x: (pos[x][0], pos[x][1])):
    print(f"  {g}  {gene_final[g]:<8}  伴侣={sorted(gene_partners.get(g, []))}")
print("\n类别计数:", dict(cnt))

# ---------- 6. 输出 ----------
def chr_lab(c):
    n = int(c.replace("GWHBCLA000000", ""))
    return (f"Chr{n}" if n <= 12 else f"Scaf{n}")

with open(f"{BASE}/outputs/地菍SAP_基因重复事件类型.tsv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["序列ID", "基因ID", "染色体", "起始(bp)", "终止(bp)", "重复事件类型", "家族内重复伴侣"])
    for g in sorted(pos, key=lambda x: (pos[x][0], pos[x][1])):
        c, s, e, st = pos[g]
        dr = [k for k, v in dr2sap.items() if v == g][0]
        w.writerow([g, dr, chr_lab(c), s, e, gene_final[g],
                    "、".join(sorted(gene_partners.get(g, []))) or "—"])
    w.writerow([])
    w.writerow(["配对明细：g1", "g2", "类型", "染色体1", "染色体2", "间隔基因数",
                "蛋白一致性(%)", "Ka/Ks", "Ks", "证据"])
    for r in sorted(pair_rows, key=lambda x: (ORDER.index(x["类型"]), x["g1"])):
        w.writerow([r["g1"], r["g2"], r["类型"], chr_lab(pos[r["g1"]][0]), chr_lab(pos[r["g2"]][0]),
                    r["间隔基因数"], r["蛋白一致性"], r["KaKs"], r["Ks"], r["证据"]])
print("TSV 已输出")

json.dump(dict(gene_final=gene_final,
               pairs=pair_rows,
               partners={g: sorted(v) for g, v in gene_partners.items()}),
          open(f"{BASE}/tmp_kaks/dup_classification.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1, default=str)
