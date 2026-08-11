# Sanger Stage-1/2 toolkit

Gene × sequencing-file matrix pipeline for binder/primer batches.

**Stage-1 (`.seq`)**: identity + indel/SNV screen (primary place for insert/delete/mismatch).  
**Stage-2 (`.ab1`)**: peak purity ≥ 75% on stage-1 worth files.  
**Final**: gene × sample matrix for plasmid return (T7/T7T merged).

## Install

```bash
pip install -r requirements.txt
```

Dependency: `biopython`.

## Stage-1 (`.seq` — identity + indel/SNV)

```bash
python stage1_seq_matrix.py \
  --seq-dir /path/to/F0xxxxx测序结果 \
  --insert-csv /path/to/pcr_insert_flanks.csv \
  --aliases aliases.json \
  --out-dir /path/to/stage1_out
```

### Inputs

| Input | Role |
|-------|------|
| `*.seq` in `--seq-dir` | Stage-1 reads |
| `pcr_insert_flanks.csv` | `vendor_id`, `full_insert_dna` (flank+CDS+flank) |

Filename: `S########_<sample>_<T7|T7T>_<well>.seq`

HER2 mapping: solo `h6277-2` → `6277-*`; pair `h2950-1`/`h2950-2` → `2950-*`/`29502-*`.

### Outputs (maintain together)

| File | Meaning |
|------|---------|
| `stage1_matrix.csv` | rows=genes, cols=seq files, cell=`PASS`/`FAIL`/empty |
| `stage1_matrix_long.csv` | tidy long form |
| `stage1_clone_union.csv` | per-clone FW∪RV stats (+ `internal_insert_bp`) |
| `stage1_files_worth_ab1.csv` | files for stage-2 |
| `stage1_genes_not_failed.csv` | genes with ≥1 PASS clone |
| `stage1_per_read.csv` | per-reaction detail |

### Stage-1 scoring

1. Align each reaction as-is and RC; keep better orientation.
2. SUCCESS if FW **or** RV calls the reference base (union).
3. FAILURE if uncovered or no match (covers deletions / SNVs).
4. Trusted internal inserts (`match≥90% ref` and `identity≥0.95`) → FAIL.
5. Clone PASS iff every reference base succeeds **and** no trusted inserts.

## Stage-2 (`.ab1`, peak ≥ 75%)

Only analyze files listed in `stage1_files_worth_ab1.csv`. Main added check: peak purity.

```bash
python stage2_ab1_matrix.py \
  --ab1-dir /path/to/F0xxxxx测序结果 \
  --insert-csv /path/to/pcr_insert_flanks.csv \
  --worth-csv /path/to/stage1_out/stage1_files_worth_ab1.csv \
  --aliases aliases.json \
  --out-dir /path/to/stage2_out \
  --peak-min 0.75
```

Standard threshold: **`--peak-min 0.75`** (also the script default).

A site counts only if:

- call == reference base, **and**
- `height[correct] / sum(GATC) ≥ 0.75` at `PLOC` (`FWO_1` dye order)

### Stage-2 outputs

| File | Meaning |
|------|---------|
| `stage2_matrix.csv` | rows=genes, cols=AB1 files |
| `stage2_matrix_long.csv` | tidy long form |
| `stage2_clone_union.csv` | per-clone union + peak_fail / insert counts |
| `stage2_files_pass.csv` | AB1 files whose clone PASS |
| `stage2_genes_not_failed.csv` | genes with ≥1 PASS clone |
| `stage2_per_read.csv` | per-reaction detail |

## Final export (gene × sample, plasmid return)

T7/T7T of the same clone = **one sample**. This is the delivery table for collecting return plasmids.

```bash
python export_final_sample_matrix.py \
  --stage1-clones /path/to/stage1_out/stage1_clone_union.csv \
  --stage2-clones /path/to/stage2_out/stage2_clone_union.csv \
  --genes-csv /path/to/pcr_insert_flanks.csv \
  --out-dir /path/to/final_out
```

| File | Meaning |
|------|---------|
| `final_gene_sample_matrix.csv` | rows=genes, cols=sample_1..n (`CORRECT` / `SEQ_OK_AB1_FAIL` / `FAIL`) |
| `final_gene_sample_matrix.html` | same matrix as webpage |
| `final_gene_sample_long.csv` | tidy long form |
| `final_plasmid_return_list.csv` | CORRECT samples only → 收集返样质粒 |
