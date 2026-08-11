---
name: sanger-stage1-seq-matrix
description: >-
  General two-stage Sanger QC: user-provided target sequences + read assignments
  (any number of primers per clone), stage-1 .seq identity/indel screen, stage-2
  .ab1 peak purity ≥75%, final target×sample plasmid-return matrix. Use for
  Sanger review, seq/ab1 QC, multi-primer coverage, or return-plasmid picking.
---

# Sanger Stage-1/2 QC (general)

## What the user must provide

1. **Target sequences** (`--targets` CSV/FASTA)  
   The DNA interval you want verified. It **must be coverable** by the sequencing primers/reads you assign. Do not include regions nobody sequenced.

2. **Read assignments** (`--assignments` CSV)  
   Explicitly map each read file → `target_id` + `clone_id` + `primer`.  
   - Off-target / wrong-locus reads: **leave them out** of that target’s rows  
   - Long inserts needing 3–4 primers: multiple rows, **same `clone_id`**, different `primer`  
   - Agent may draft assignments from filenames; **user confirms** before trusting results

Example assignment:

```csv
file,target_id,clone_id,primer
read1_T7.seq,plasmidA,clone-1,T7
read1_mid.seq,plasmidA,clone-1,MID
read1_T7T.seq,plasmidA,clone-1,T7T
```

## Pipeline

```
targets + assignments + reads
  → Stage-1 (.seq): identity + indel/SNV union over all primers of a clone
  → Stage-2 (.ab1 of stage-1 PASS): + peak purity ≥ 0.75
  → Final: target × sample (all primers merged) for return plasmids
```

## Commands

```bash
pip install -r toolkit/requirements.txt

python toolkit/stage1_seq_matrix.py \
  --targets targets.csv \
  --assignments assignments.csv \
  --reads-dir READS_DIR \
  --out-dir stage1_out

python toolkit/stage2_ab1_matrix.py \
  --targets targets.csv \
  --assignments assignments.csv \
  --reads-dir READS_DIR \
  --worth-csv stage1_out/stage1_files_worth_ab1.csv \
  --out-dir stage2_out \
  --peak-min 0.75

python toolkit/export_final_sample_matrix.py \
  --stage1-clones stage1_out/stage1_clone_union.csv \
  --stage2-clones stage2_out/stage2_clone_union.csv \
  --targets targets.csv \
  --out-dir final_out
```

## Scoring (summary)

- Reference = user target sequence
- Per clone: **union of all assigned primers** (2–N)
- Stage-1 site SUCCESS if any primer calls the reference base; uncovered/mismatch = fail
- Trusted internal inserts (high-confidence alignment only) → fail
- Stage-2 also requires correct-base peak fraction ≥ **0.75**
- Clone PASS ⇔ every target base succeeds and no trusted inserts
- Final `CORRECT` sample = stage-2 PASS clone (collect return plasmid)

High-confidence for insert trust: `match_bp ≥ 0.90 × ref_len` and `identity ≥ 0.95`.

## Agent checklist

```
- [ ] Get targets file (or write one from user sequences)
- [ ] Confirm targets are within primer coverage
- [ ] Draft/confirm assignments.csv (which reads belong where; multi-primer OK)
- [ ] Run stage-1 → report
- [ ] Run stage-2 on worth files → report
- [ ] Export final target×sample return list
```

## Package layout

```
SKILL.md
reference.md
README.md
examples/
toolkit/
  common.py
  stage1_seq_matrix.py
  stage2_ab1_matrix.py
  export_final_sample_matrix.py
  draft_her2_inputs.py   # optional naming helper only
```

See [reference.md](reference.md) and [toolkit/README.md](toolkit/README.md).
