# Sanger stage-1/2 reference

## Division of labor

| Stage | Primary detections |
|-------|--------------------|
| Stage-1 `.seq` | Identity, SNV, deletion (uncovered), insertion (`internal_insert_bp`) |
| Stage-2 `.ab1` | Same union + **peak purity ≥ 0.75** (main added value) |

Indel/SNV screening is intentionally simpler and primary in stage-1.

## Matrix model

Every task keeps one matrix:

- **Rows** = genes to synthesize (`vendor_id`)
- **Columns** = sequencing files (`.seq` in stage-1, `.ab1` in stage-2)
- **Cell** = clone-union status of the clone owning that file (`PASS` / `FAIL` / empty)

Two views are projections of the same matrix:

1. Files worth keeping / pass
2. Genes not fully failed (`gene_not_failed=YES` ⇔ ≥1 PASS clone)

## Union scoring

For clone \(C\) with reference \(R\) and reactions \(\{FW, RV\}\), each position \(i\):

| Stage | SUCCESS | FAILURE |
|-------|---------|---------|
| 1 | ∃ reaction with call \(= R[i]\) | uncovered or no matching call |
| 2 | ∃ reaction with call \(= R[i]\) **and** peak_frac ≥ **0.75** | uncovered, mismatch, or peak too low |

### Variant classes

| Event | How it fails |
|-------|----------------|
| Substitution (SNV) | Mapped base ≠ reference |
| Deletion / missing span | Reference base uncovered by FW∪RV |
| Insertion | Trusted internal query-only bases inside ref span |

**Trusted internal insertions** (primarily stage-1): on alignments with match ≥90% of ref length **and** identity ≥0.95, query-only bases inside the reference span (`ref` gap + `query` base) count as `internal_insert_bp`. Leading/trailing vector sequence and low-quality opposite-primer gapped artifacts are ignored.

Clone PASS ⟺ all positions SUCCESS **and** trusted `internal_insert_bp == 0` (stage-2 also requires the peak gate).

## High-confidence (insert trust)

Both required:

- `match_bp ≥ 0.90 × ref_len`
- `identity ≥ 0.95`

## Stage-2 peak gate

At basecall index `i` with `PLOC[i]`:

- Read `DATA9..12` ordered by `FWO_1` (typically `GATC`)
- `peak_frac = height[want] / sum(heights)`
- `want = ref_base` (as_is) or `complement(ref_base)` (RC_of_read)
- **Standard threshold: `0.75`**

## Alignment parameters

`pairwise2.align.localms(ref, query, 2, -1, -7, -1)`

- Walk gapped strings from `ref_pos = 0` (do **not** add `alignment.start`)
- Stage-1 orientation: maximize `(match_bp, -trusted_inserts, aligned_len, identity)`
- Stage-2 orientation: maximize `(support_bp, -trusted_inserts, aligned_len, identity)`

## Input CSV

`pcr_insert_flanks.csv` minimum columns:

- `vendor_id`
- `full_insert_dna`
