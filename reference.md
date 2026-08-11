# Scoring reference (general)

## Inputs

| File | Role |
|------|------|
| `targets` | Sequences to verify (CSV/FASTA). Must be coverable by assigned reads. |
| `assignments` | `file,target_id,clone_id,primer` — which reads belong to which clone. |
| reads dir | `.seq` / `.ab1` files named as in `assignments.file` (stem match allowed). |

Any number of primers per `clone_id` are allowed.

## Two evidence tracks (do not mix)

| Track | Unit | Merge rule |
|-------|------|------------|
| **Base correctness** | reference **position** | **Union**: any reliable read calling the correct base → SUCCESS |
| **Insertion** | reference **boundary** (between bases) | Separate merge by boundary / interval |

A bad base call on one primer must **not** veto a correct call from another primer.
Insertion evidence is **not** judged by whole-target coverage.

## Stage roles

| Stage | Primary job |
|-------|-------------|
| 1 `.seq` | Base union + discover high-confidence insertion candidates |
| 2 `.ab1` | Base union + peak ≥ 0.75; validate insertions on chromatogram |

## Base union

For each target position of a clone: SUCCESS if **any** assigned primer covers it with the reference base (stage-2: and peak_frac ≥ 0.75). Uncovered or wrong = FAILURE for that position.

## Insertions (boundary events)

From pairwise alignment, extract REF-gap + QUERY-base runs → `insertion after ref position N`.

**High-confidence** uses **local** flank quality only (constants in `toolkit/insertion.py`):

| Constant | Default |
|----------|---------|
| `LOCAL_IDENTITY_MIN` | 0.95 |
| `FLANK_MATCH_MIN_BP` | 20 |
| `EDGE_MARGIN_BP` | 25 |
| `LOCAL_WINDOW_BP` | 40 |

There is **no** `match_bp ≥ 0.90 × ref_len` gate (that caused miss on long targets with partial Sanger coverage).

### Stage-1 merge statuses

| Status | Meaning |
|--------|---------|
| `NO_INSERT_EVIDENCE` | no HC insertion |
| `INSERT_CANDIDATE` | one HC read supports insert; others uncovered / low quality |
| `STRONG_INSERT` | ≥2 HC reads support the same insert |
| `INSERT_CONFLICT` | HC insert vs HC continuous NO_INSERT spanning the boundary |

Homopolymer / repeat: boundaries may widen to an interval (e.g. positions 120–124).

Stage-2 AB1 peak checks are **per insertion event** `(evidence_id, interval, inserted_seq)`, not a single bool per primer. Same-primer replicate files are separate evidence units (`primer::filename`).

Targets: whitespace stripped; **non-ACGT bases (N/IUPAC) raise a loud error** — never silently deleted.

Stage-1 clone status: `PASS` (bases OK, no insert) or `PASS_INSERT_REVIEW` (bases OK, insert needs AB1) — both go to stage-2 worth list.

### Stage-2 validation statuses

| Status | Meaning |
|--------|---------|
| `NO_INSERT_EVIDENCE` | none |
| `CONFIRMED_INSERT` | seq + clean AB1 peaks |
| `POSSIBLE_INSERT` | seq candidate but peaks inconclusive (never silently dropped) |
| `INSERT_CONFLICT` | conflicting HC evidence |

Stage-2: `PASS` only if bases perfect **and** `NO_INSERT_EVIDENCE`. Otherwise `REVIEW` / `FAIL`.

## Peak fraction (stage-2 bases)

`height[correct_dye] / sum(GATC)` at `PLOC`, dye order from `FWO_1`. Threshold default **0.75**.

## Final sample

One sample = one `clone_id` (all its primers).

| Final | Rule |
|-------|------|
| `CORRECT` | stage-2 `PASS` → collect return plasmid |
| `REVIEW` | `POSSIBLE_INSERT` / `INSERT_CONFLICT` (manual check) |
| `FAIL` | base fail or `CONFIRMED_INSERT` |
