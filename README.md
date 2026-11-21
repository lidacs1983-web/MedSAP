# MedSAP — Analysis scripts for the stress-associated protein (SAP) gene family of *Melastoma dodecandrum*

Companion code for the manuscript "Genome-Wide Identification, Evolution, and Screening of
Drought-Response Candidate Genes of the Stress-Associated Protein (SAP) Gene Family in
*Melastoma dodecandrum*" (submitted to PeerJ).

The reference genome used throughout is GWHBCLA00000000 (NGDC Genome Warehouse,
BioProject PRJCA004869; Hao et al. 2022).

## Repository layout

| Script | Purpose |
|---|---|
| `scripts/01_identify_sap.py` | pyhmmer search of the Pfam A20 (PF01754.24) and AN1 (PF01428.24) profiles against each proteome |
| `scripts/02_extract_domains.py` | domain-level classification (A20-AN1 / AN1-only / A20-only) and quality filtering |
| `scripts/03_collect_align.py` | collection of the nine-species SAP set and multiple-sequence alignment |
| `scripts/04_kaks_duplication.py` | duplication-type classification and intra-specific Ka/Ks with the Biopython codonalign module |
| `scripts/05_kaks_cross_species.py` | genus-level orthologous Ka/Ks (M. dodecandrum vs M. candidum) |
| `scripts/06_circos_distribution.py` | chromosomal distribution and intra-genomic collinearity Circos plot |
| `scripts/07_promoter_cis_scan.py` | 2-kb promoter cis-acting element scan and summary |
| `scripts/08_qpcr_renorm.py` | RT-qPCR normalisation with single or geometric-mean reference genes |

## Environment

```bash
python >= 3.10
pip install -r requirements.txt
# external tools: cd-hit (95% redundancy removal), pyhmmer is pip-installable
```

## Reference commands (as stated in the manuscript, Code availability section)

1. Domain search (pyhmmer):

   ```bash
   python scripts/01_identify_sap.py --pfam Pfam-A.hmm --proteome proteome.fa \
       --evalue 1e-10 --domT 1.0
   ```

2. Redundancy removal (CD-HIT, 95% pairwise identity):

   ```bash
   cd-hit -i all_sap.fa -o sap_nr95.fa -c 0.95 -n 5 -M 0
   ```

3. Pairwise codon alignment and Ka/Ks (Biopython codonalign; Ks > 2 treated as saturated):

   ```python
   from Bio import AlignIO
   from Bio.codonalign import build, cal_dn_ds
   from Bio.phylo.treeio import read as read_tree
   aln = AlignIO.read("prot.aln", "fasta")
   tree = read_tree("species.nwk", "newick")
   for n1, n2 in pairs:
       ccaln = build(aln, [cds1, cds2], tree)
       dn, ds = cal_dn_ds(ccaln)
   ```

4. Promoter scan, duplication classification and plotting: `scripts/07_promoter_cis_scan.py`,
   `scripts/04_kaks_duplication.py`, `scripts/06_circos_distribution.py`.

## Data availability

- Reference genome: https://ngdc.cncb.ac.cn/gwh/ (GWHBCLA00000000, PRJCA004869)
- Raw RT-qPCR data: Figshare (DOI: 10.6084/m9.figshare.XXXXXXX — see manuscript Data availability)

## License

MIT — see [LICENSE](LICENSE).
