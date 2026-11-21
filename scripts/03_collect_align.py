# -*- coding: utf-8 -*-
"""构建 9 物种 SAP 家族数据集（Stress-Associated Protein, A20/AN1 锌指）

流程:
 1. 读取 sap9_hits.json（E≤1e-10 的 PF01754/PF01428 命中）
 2. 保留含完整结构域的成员（A20 或 AN1 覆盖 ≥0.85）
 3. 物种内剔除完全相同的蛋白序列（同一基因的多转录本/重复注释，保留一条代表）
 4. 按物种重编号写出 fasta + meta（记录域类型：A20-AN1 / A20 only / AN1 only）
"""
import os, json, collections
from Bio import SeqIO

BASE = "C:/Users/lidac/WorkBuddy/2026-09-29-20-40-35"
PHY = os.path.join(BASE, "tmp_phylo")
TEA = os.path.join(BASE, "tmp_tea")

SPNAME = {"Mcan": "Mca", "Egr": "Egr", "Vvi": "Vvi",
          "Nta": "Nta", "Vco": "Vco", "Csi": "Csi", "Ngr": "Ngr"}
ORDER = ["Med", "Mcan", "Egr", "Vvi", "Vco", "Nta", "At", "Csi", "Ngr"]
PROT = {
    "Med":  os.path.join(PHY, "med.pep.faa"),
    "Mcan": os.path.join(PHY, "mcan.pep.faa"),
    "Egr":  os.path.join(PHY, "egrandis.pep.faa"),
    "Vvi":  os.path.join(PHY, "vvinifera.pep.faa"),
    "Nta":  os.path.join(PHY, "ntabacum.pep.faa"),
    "Vco":  os.path.join(PHY, "vcorymbosum.pep.faa"),
    "At":   os.path.join(PHY, "ath.pep.faa"),
    "Csi":  os.path.join(TEA, "ncbi_dataset/data/GCA_004153795.2/protein.faa"),
    "Ngr":  os.path.join(TEA, "Ngr.protein.faa"),
}
COV_MIN = 0.85

hits = json.load(open(os.path.join(PHY, "sap9_hits.json")))
bysp = collections.defaultdict(list)
dropped_cov, dropped_dup = [], []

for sp in ORDER:
    seqs = {r.id.split()[0]: str(r.seq)
            for r in SeqIO.parse(PROT[sp], "fasta")}
    seen = {}
    for g, x in sorted(hits[sp]["hits"].items()):
        a = x["dom"].get("PF01754", {}).get("cov", 0)
        b = x["dom"].get("PF01428", {}).get("cov", 0)
        if max(a, b) < COV_MIN:
            dropped_cov.append((sp, g)); continue
        s = seqs.get(g)
        if s is None:
            continue
        if s in seen:                       # 完全相同序列 -> 只保留一条代表
            dropped_dup.append((sp, g, seen[s])); continue
        seen[s] = g
        dtype = ("A20-AN1" if (a >= COV_MIN and b >= COV_MIN)
                 else ("A20 only" if a >= COV_MIN else "AN1 only"))
        bysp[sp].append({"orig": g, "seq": s, "len": len(s),
                         "a20": round(a, 3), "an1": round(b, 3),
                         "type": dtype, "evalue": x["evalue"],
                         "both": bool(a >= COV_MIN and b >= COV_MIN)})

out = open(os.path.join(PHY, "sap_all9.fasta"), "w", encoding="utf-8")
meta = []
for sp in ORDER:
    items = sorted(bysp[sp], key=lambda z: z["orig"])
    for i, it in enumerate(items, 1):
        nm = (f"MedSAP{i}" if sp == "Med"
              else (f"AtSAP{i}" if sp == "At" else f"{SPNAME[sp]}_SAP{i:02d}"))
        out.write(f">{nm}|{sp}|{it['orig']}\n{it['seq']}\n")
        m = dict(it); m.pop("seq")
        m.update({"new": nm, "species": sp, "data": "sap_all9"})
        meta.append(m)
out.close()
json.dump(meta, open(os.path.join(PHY, "sap_all9_meta.json"), "w"), indent=1)

cnt = collections.Counter(m["species"] for m in meta)
typ = collections.Counter(m["type"] for m in meta)
print(f"sap_all9.fasta 写出 {len(meta)} 条")
for sp in ORDER:
    print(f"  {sp:<5}{cnt[sp]:>3}  " + ", ".join(
        m["new"] for m in meta if m["species"] == sp))
print("\n域类型:", dict(typ))
print(f"剔除覆盖不足 {len(dropped_cov)} 条；剔除完全相同冗余 {len(dropped_dup)} 条")
for sp, g, keep in dropped_dup:
    print(f"    {sp:<5}{g:<22} == {keep}")
