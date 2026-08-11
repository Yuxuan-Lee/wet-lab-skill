---
name: sanger-stage1-seq-matrix
description: >-
  Two-stage Sanger QC for binder/primer batches: stage-1 .seq does identity +
  indel/SNV screening; stage-2 .ab1 adds peak purity ≥75%. Use for 测序筛选,
  seq/ab1 分析, HER2 clone QC, flank+CDS+flank validation, or gene×测序文件矩阵.
---

# Sanger Stage-1/2 Seq+AB1 Matrix

Two-stage Sanger QC. At **every** stage maintain one **gene × sequencing-file matrix** and its two views (files worth keeping; genes not fully failed).

## Division of labor (important)

| Stage | Input | Primary job |
|-------|-------|-------------|
| **1** | `.seq` only | Full-insert identity + **indel / SNV** screen (simpler; do this first) |
| **2** | `.ab1` of stage-1 PASS only | Same union check + **peak purity ≥ 75%** |
| **Final** | clone CSVs | Gene × **sample** matrix for 返样质粒 (T7/T7T = one sample) |

Indel/SNV detection belongs primarily in **stage-1**. Stage-2 may still report the same insert fields, but its unique value is chromatogram peak purity.

## Canonical package (one directory)

Project package (source of truth to edit):

`primer_design/sanger-stage1-seq-matrix/`

Cursor skill mirror (keep in sync):

`~/.cursor/skills/sanger-stage1-seq-matrix/`

```
sanger-stage1-seq-matrix/
├── SKILL.md
├── reference.md
├── README.md
└── toolkit/
    ├── requirements.txt
    ├── aliases.example.json
    ├── stage1_seq_matrix.py
    ├── stage2_ab1_matrix.py
    └── export_final_sample_matrix.py
```

```bash
pip install -r toolkit/requirements.txt
```

## Normalized pipeline

```
Stage-1 (.seq: identity + indel/SNV)
  → stage1_out/
Stage-2 (.ab1: peak ≥75%, worth files only)
  → stage2_out/
Final (gene × sample for plasmid return)
  → final_out/
```

### Stage-1 — `.seq` only (identity + indel/SNV)

```bash
python toolkit/stage1_seq_matrix.py \
  --seq-dir <测序结果目录> \
  --insert-csv <pcr_insert_flanks.csv> \
  --aliases <optional aliases.json> \
  --out-dir <stage1_out>
```

Do **not** open `.ab1` in stage-1.

Stage-1 catches:

- **SNV / mismatch**: mapped base ≠ reference
- **Deletion / missing coverage**: reference base uncovered by FW∪RV
- **Insertion**: trusted internal query-only bases (`internal_insert_bp > 0`)

### Stage-2 — `.ab1` of stage-1 worth files only (peak purity)

```bash
python toolkit/stage2_ab1_matrix.py \
  --ab1-dir <测序结果目录> \
  --insert-csv <pcr_insert_flanks.csv> \
  --worth-csv <stage1_out/stage1_files_worth_ab1.csv> \
  --aliases <optional aliases.json> \
  --out-dir <stage2_out> \
  --peak-min 0.75
```

Default / standard threshold: **`--peak-min 0.75`**.

## Non-negotiable rules

1. **Reference** = `full_insert_dna` = flank + CDS + flank.
2. Bidirectional T7/T7T; try as-is and RC; keep better orientation.
3. Per reference base, FW∪RV **union**:
   - Stage-1 SUCCESS: any reaction calls the reference base
   - Stage-2 SUCCESS: any reaction calls the reference base **and** correct-base peak fraction ≥ **0.75**
   - FAILURE: uncovered, mismatch, or (stage-2) peak too low
4. **Indel / SNV — primary gate in stage-1** (also recorded in stage-2):
   - Substitution → mapped mismatch → FAIL
   - Deletion / gap in coverage → uncovered ref base → FAIL
   - Insertion → on a **high-confidence** alignment (match ≥90% of ref **and** identity ≥0.95), query bases inside the reference span that do not map to any reference position count as `internal_insert_bp`. Leading/trailing vector sequence and low-quality opposite-primer artifacts are ignored. Trusted `internal_insert_bp > 0` → FAIL
5. Clone PASS ⇔ every reference base SUCCESS **and** no trusted internal inserts (stage-2 also requires peak gate).
6. Matrix: rows=genes (`vendor_id`), cols=sequencing files, cell=`PASS`/`FAIL`/empty (clone-union status of the owning clone).
7. Always refresh both views with the matrix:
   - Stage-1: `stage1_files_worth_ab1.csv`, `stage1_genes_not_failed.csv`
   - Stage-2: `stage2_files_pass.csv`, `stage2_genes_not_failed.csv`

### Peak fraction (stage-2 only)

At each basecall `PLOC`, read `DATA9–12` ordered by `FWO_1` (usually `GATC`):

`peak_frac = height[correct_dye] / sum(G+A+T+C)`

RC orientation: correct dye = complement of the reference base on the raw chromatogram.

### High-confidence (for trusting insert calls)

Both required:

- `match_bp ≥ 0.90 × ref_len`
- `identity ≥ 0.95`

## Agent workflow checklist

```
Progress:
- [ ] Locate seq/ab1 dir + pcr_insert_flanks.csv (+ aliases if needed)
- [ ] Stage-1: run stage1_seq_matrix.py → stage1_out/  (identity + indel/SNV)
- [ ] Report stage-1: clone PASS/FAIL, genes_not_failed, files_worth_ab1
- [ ] Stage-2: run stage2_ab1_matrix.py --peak-min 0.75 on worth files only → stage2_out/  (peak purity)
- [ ] Report stage-2: clone PASS/FAIL, genes_not_failed, files_pass
- [ ] Final export: gene × sample (T7/T7T merged) for plasmid return
- [ ] Keep stage1_out/, stage2_out/, final_out/ as the task record
```

### Final export (plasmid return)

Stage-1/2 matrices are file-level. **Final delivery is sample-level**: one sample = one clone (FW+RV / T7+T7T merged). Mapping gene↔sample comes from the sample token in the filename (plus aliases); you do not need a separate manual 1:1 table beyond that.

```bash
python toolkit/export_final_sample_matrix.py \
  --stage1-clones <stage1_out/stage1_clone_union.csv> \
  --stage2-clones <stage2_out/stage2_clone_union.csv> \
  --genes-csv <pcr_insert_flanks.csv> \
  --out-dir <final_out>
```

Outputs:

- `final_gene_sample_matrix.csv` / `.html` — rows=genes, cols=samples, `CORRECT` = 可收集返样质粒
- `final_gene_sample_long.csv` — tidy form
- `final_plasmid_return_list.csv` — only CORRECT samples

Status labels: `CORRECT` (stage-2 PASS), `SEQ_OK_AB1_FAIL`, `FAIL`.

### HER2 sample mapping

- Solo `h6277-2` → samples `6277-*`
- Pair `h2950-1`/`h2950-2` → `2950-*` / `29502-*`
- Typos via `--aliases` (see `toolkit/aliases.example.json`)

### pairwise2 pitfall

Walk gapped strings from `ref_pos=0`. Do **not** add `alignment.start`.

## Outputs checklist

**Stage-1** (`stage1_out/`):

- [ ] `stage1_matrix.csv`
- [ ] `stage1_matrix_long.csv`
- [ ] `stage1_clone_union.csv`
- [ ] `stage1_files_worth_ab1.csv`
- [ ] `stage1_genes_not_failed.csv`
- [ ] `stage1_per_read.csv`

**Stage-2** (`stage2_out/`):

- [ ] `stage2_matrix.csv`
- [ ] `stage2_matrix_long.csv`
- [ ] `stage2_clone_union.csv`
- [ ] `stage2_files_pass.csv`
- [ ] `stage2_genes_not_failed.csv`
- [ ] `stage2_per_read.csv`

**Final** (`final_out/`):

- [ ] `final_gene_sample_matrix.csv`
- [ ] `final_gene_sample_matrix.html`
- [ ] `final_gene_sample_long.csv`
- [ ] `final_plasmid_return_list.csv`

## Additional resources

- CLI / install: [toolkit/README.md](toolkit/README.md)
- Scoring & schema: [reference.md](reference.md)
