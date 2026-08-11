---
name: sanger-sequencing-qc
description: >-
  General two-stage Sanger QC: user-provided target sequences + read assignments
  (any number of primers per clone), stage-1 .seq identity/indel screen, stage-2
  .ab1 peak purity ≥75%, boundary insertion QC, final target×sample plasmid-return
  matrix. Use for Sanger review, seq/ab1 QC, multi-primer coverage, insertion
  review, or return-plasmid picking.
---

# Sanger Sequencing QC

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
- **Base correctness** = reference-**position** union (any good primer OK)
- **Insertion** = reference-**boundary** event (separate merge; local flank gate)
- Stage-1: discover HC inserts → `INSERT_CANDIDATE` / `STRONG_INSERT` / `INSERT_CONFLICT`
- Stage-2: peak ≥ **0.75** for bases; AB1 validates inserts → `CONFIRMED` / `POSSIBLE` / `CONFLICT`
- Final `CORRECT` = stage-2 PASS (bases OK + `NO_INSERT_EVIDENCE`); `POSSIBLE`/`CONFLICT` → `REVIEW`

Insert HC defaults (`toolkit/insertion.py`): local identity ≥0.95, flanks ≥20 bp each, not within 25 bp of read edge. **No** whole-target 90% coverage gate.

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
  insertion.py           # boundary insertion logic + thresholds
  stage1_seq_matrix.py
  stage2_ab1_matrix.py
  export_final_sample_matrix.py
  draft_her2_inputs.py   # optional naming helper only
  tests/test_insertion.py
```

See [reference.md](reference.md) and [toolkit/README.md](toolkit/README.md).
