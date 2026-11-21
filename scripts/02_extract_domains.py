# -*- coding: utf-8 -*-
"""从本地 Pfam-A.hmm.gz 抽取 SAP 家族相关 HMM 模型

SAP = Stress-Associated Protein (A20/AN1 锌指)
  PF01754  zf-A20
  PF01428  zf-AN1
"""
import gzip, os

BASE = "C:/Users/lidac/WorkBuddy/2026-09-29-20-40-35"
PHY = os.path.join(BASE, "tmp_phylo")
SRC = os.path.join(PHY, "Pfam-A.hmm.gz")
DST = os.path.join(PHY, "sap_models.hmm")
WANT = {"PF01754", "PF01428"}

kept, buf, hit_acc = [], [], None
with gzip.open(SRC, "rt", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        buf.append(line)
        if line.startswith("ACC "):
            acc = line.split()[1].split(".")[0]
            if acc in WANT:
                hit_acc = acc
                print("found", acc, line.split()[1])
        elif line.strip() == "//":
            if hit_acc:
                kept.append("".join(buf))
                hit_acc = None
            buf = []

with open(DST, "w", encoding="utf-8") as out:
    out.write("".join(kept))
print("wrote", DST, len(kept), "models")

# 校验
import pyhmmer
with pyhmmer.plan7.HMMFile(DST) as hf:
    for h in hf:
        print(" check:", h.accession.decode() if isinstance(h.accession, bytes)
              else h.accession,
              "M=", h.M,
              (h.description.decode() if isinstance(h.description, bytes)
               else h.description))
