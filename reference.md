# Scoring reference (general)

## Inputs

| File | Role |
|------|------|
| `targets` | Sequences to verify (CSV/FASTA). Must be coverable by assigned reads. |
| `assignments` | `file,target_id,clone_id,primer` — which reads belong to which clone. |
| reads dir | `.seq` / `.ab1` files named as in `assignments.file` (stem match allowed). |

Any number of primers per `clone_id` are allowed. Their calls form a **union** over target bases.

## Stage roles

| Stage | Primary job |
|-------|-------------|
| 1 `.seq` | Identity + SNV + deletion (uncovered) + trusted insertion |
| 2 `.ab1` | Same union + peak purity ≥ 0.75 |

## Union rule

For each target position of a clone: SUCCESS if **any** assigned primer covers it with the reference base (stage-2: and peak_frac ≥ 0.75). Uncovered or wrong = FAILURE for that position.

## Insertions

On high-confidence alignments (`match_bp ≥ 0.90×ref_len` and `identity ≥ 0.95`), query-only bases inside the target span count as `internal_insert_bp`. Leading/trailing vector sequence and low-quality primer artifacts are ignored. Trusted inserts > 0 → clone FAIL.

## Peak fraction (stage-2)

`height[correct_dye] / sum(GATC)` at `PLOC`, dye order from `FWO_1`. Threshold default **0.75**.

## Final sample

One sample = one `clone_id` (all its primers). `CORRECT` = stage-2 PASS → collect return plasmid.
