# -*- coding: utf-8 -*-
"""9 物种 SAP 家族统一鉴定（Stress-Associated Protein, A20/AN1 锌指）

模型: PF01754 (zf-A20, M=29) + PF01428 (zf-AN1, M=40)
对每个物种完整蛋白组做 hmmsearch，逐基因记录：
  两条模型的 domain 覆盖比例、E 值、序列长度
输出 sap9_hits.json
"""
import os, json
import pyhmmer
from pyhmmer.easel import SequenceFile

BASE = "C:/Users/lidac/WorkBuddy/2026-09-29-20-40-35"
PHY = os.path.join(BASE, "tmp_phylo")
TEA = os.path.join(BASE, "tmp_tea")
HMM = os.path.join(PHY, "sap_models.hmm")
E_CUT = 1e-10

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


def dec(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


with pyhmmer.plan7.HMMFile(HMM) as hf:
    hmms = list(hf)
for h in hmms:
    print("model", dec(h.accession), "M=", h.M,
          dec(h.description) if h.description else "")

out = {}
for sp, fa in PROT.items():
    if not os.path.exists(fa):
        print(f"{sp}: 缺失 {fa}")
        continue
    with SequenceFile(fa, digital=True) as sf:
        seqs = sf.read_block()
    ln = {}
    for s in seqs:
        ln[dec(s.name).split()[0]] = len(s)
    rec = {}
    for res in pyhmmer.hmmsearch(hmms, seqs, E=E_CUT, cpus=4):
        acc = dec(res.query.accession).split(".")[0]
        M = res.query.M
        for hit in res:
            nm = dec(hit.name).split()[0]
            d = rec.setdefault(nm, {"len": ln.get(nm, 0), "evalue": hit.evalue,
                                    "dom": {}})
            if hit.evalue < d["evalue"]:
                d["evalue"] = hit.evalue
            for dom in hit.domains:
                al = dom.alignment
                c = (al.hmm_to - al.hmm_from + 1) / M
                prev = d["dom"].get(acc)
                if prev is None or c > prev["cov"]:
                    d["dom"][acc] = {"cov": round(c, 4),
                                     "e": float(f"{dom.i_evalue:.2e}"),
                                     "score": round(dom.score, 1),
                                     "hmm_f": al.hmm_from, "hmm_t": al.hmm_to}
    out[sp] = {"file": os.path.relpath(fa, BASE), "n_proteome": len(seqs),
               "n_hits": len(rec), "hits": rec}
    both = sum(1 for v in rec.values()
               if len(v["dom"]) == 2
               and v["dom"]["PF01754"]["cov"] >= 0.5
               and v["dom"]["PF01428"]["cov"] >= 0.5)
    a20only = sum(1 for v in rec.values()
                  if v["dom"].get("PF01754", {}).get("cov", 0) >= 0.5
                  and v["dom"].get("PF01428", {}).get("cov", 0) < 0.5)
    an1only = sum(1 for v in rec.values()
                  if v["dom"].get("PF01428", {}).get("cov", 0) >= 0.5
                  and v["dom"].get("PF01754", {}).get("cov", 0) < 0.5)
    print(f"{sp:5s} 蛋白组 {len(seqs):6d} -> 命中 {len(rec):5d} "
          f"(A20+AN1 {both}, A20-only {a20only}, AN1-only {an1only})")
    if rec:
        L = sorted(v["len"] for v in rec.values())
        print(f"       长度 min={L[0]} med={L[len(L)//2]} max={L[-1]}")

json.dump(out, open(os.path.join(PHY, "sap9_hits.json"), "w"), indent=1)
print("\nsaved tmp_phylo/sap9_hits.json")
