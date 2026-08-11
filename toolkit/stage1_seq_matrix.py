#!/usr/bin/env python3
"""Stage-1 Sanger screen (.seq) → target × read matrix.

User provides:
  --targets      sequences to verify (must be coverable by the assigned reads)
  --assignments  which read files belong to which target/clone/primer
  --reads-dir    folder with .seq files

Base correctness: UNION over primers at each reference position.
Insertions: independent reference-boundary events (see insertion.py).
"""

from __future__ import annotations

import argparse
import csv
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

from Bio.Seq import Seq

warnings.filterwarnings("ignore", category=DeprecationWarning)
from Bio import pairwise2  # noqa: E402

from common import (
    group_key,
    index_assignments,
    load_assignments,
    load_targets,
    resolve_read_path,
    validate_assignments,
)
from insertion import (
    STATUS_CANDIDATE,
    STATUS_CONFLICT,
    STATUS_NO_INSERT,
    STATUS_STRONG,
    InsertionCandidate,
    NoInsertSpan,
    candidates_from_alignment,
    merge_clone_insertions,
    no_insert_spans_from_alignment,
)


WORTH_STATUSES = {"PASS", "PASS_INSERT_REVIEW"}


def read_seq_file(p: Path) -> str:
    lines = [
        x.strip().upper()
        for x in p.read_text(encoding="utf-8", errors="ignore").splitlines()
    ]
    return "".join(x for x in lines if x and not x.startswith(">"))


def align_calls(
    read_seq: str, ref_seq: str, primer: str = ""
) -> Tuple[Dict[int, str], float, int, str, List[InsertionCandidate], List[NoInsertSpan]]:
    def one(q: str) -> Tuple[Dict[int, str], float, int, List[InsertionCandidate], List[NoInsertSpan]]:
        aln = pairwise2.align.localms(
            ref_seq, q, 2.0, -1.0, -7.0, -1.0, one_alignment_only=True
        )
        if not aln:
            return {}, 0.0, 0, [], []
        a = aln[0]
        calls: Dict[int, str] = {}
        ref_pos = 0
        match = aligned = 0
        for rc, qc in zip(a.seqA, a.seqB):
            if rc != "-" and qc != "-" and 0 <= ref_pos < len(ref_seq):
                calls[ref_pos] = qc
                aligned += 1
                if rc == qc:
                    match += 1
            if rc != "-":
                ref_pos += 1
        ident = (match / aligned) if aligned else 0.0
        cands = candidates_from_alignment(ref_seq, a.seqA, a.seqB, len(q), primer=primer)
        spans = no_insert_spans_from_alignment(a.seqA, a.seqB, primer=primer)
        return calls, ident, aligned, cands, spans

    c1, i1, a1, n1, s1 = one(read_seq)
    c2, i2, a2, n2, s2 = one(str(Seq(read_seq).reverse_complement()))
    m1 = sum(1 for pos, b in c1.items() if b == ref_seq[pos])
    m2 = sum(1 for pos, b in c2.items() if b == ref_seq[pos])
    # Prefer orientation with more correct bases; then fewer HC inserts (cleaner)
    hc1 = sum(1 for c in n1 if c.high_confidence)
    hc2 = sum(1 for c in n2 if c.high_confidence)
    if (m2, -hc2, a2, i2) > (m1, -hc1, a1, i1):
        return c2, i2, a2, "RC_of_read", n2, s2
    return c1, i1, a1, "as_is", n1, s1


def evaluate_base_union(ref_seq: str, primer_calls: Dict[str, Dict[int, str]]) -> dict:
    """Reference-position union only — insertions handled separately."""
    ref_len = len(ref_seq)
    pos_bases: Dict[int, Set[str]] = defaultdict(set)
    covered: Set[int] = set()
    for calls in primer_calls.values():
        for pos, b in calls.items():
            covered.add(pos)
            pos_bases[pos].add(b)

    success = mismatch_only = 0
    for pos in range(ref_len):
        bases = pos_bases.get(pos, set())
        if ref_seq[pos] in bases:
            success += 1
        elif bases:
            mismatch_only += 1

    return {
        "ref_len": ref_len,
        "success_bp": success,
        "fail_bp": ref_len - success,
        "covered_bp": len(covered),
        "uncovered_bp": ref_len - len(covered),
        "mismatch_only_bp": mismatch_only,
        "success_ratio": round(success / ref_len, 4) if ref_len else 0.0,
        "coverage_ratio": round(len(covered) / ref_len, 4) if ref_len else 0.0,
        "n_reads_used": len(primer_calls),
        "bases_perfect": success == ref_len,
    }


def stage1_status(bases_perfect: bool, insert_status: str) -> str:
    if not bases_perfect:
        return "FAIL"
    if insert_status == STATUS_NO_INSERT:
        return "PASS"
    return "PASS_INSERT_REVIEW"


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage-1 .seq screen with user targets + assignments")
    ap.add_argument("--targets", type=Path, required=True, help="CSV/FASTA of sequences to verify")
    ap.add_argument("--assignments", type=Path, required=True, help="CSV mapping reads → target/clone")
    ap.add_argument("--reads-dir", type=Path, required=True, help="folder containing .seq files")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    targets = load_targets(args.targets)
    assignments = load_assignments(args.assignments)
    validate_assignments(targets, assignments)
    asg_by_file = index_assignments(assignments)

    by_clone: Dict[str, dict] = {}
    per_read: List[dict] = []
    file_to_meta: Dict[str, dict] = {}
    missing: List[str] = []

    seen_paths: Set[Path] = set()
    for asg in assignments:
        path = resolve_read_path(args.reads_dir, asg.file, ".seq")
        if path is None:
            missing.append(asg.file)
            continue
        if path in seen_paths:
            continue
        seen_paths.add(path)

        asg_use = asg_by_file.get(path.name) or asg_by_file.get(Path(asg.file).name) or asg
        if asg_use.target_id not in targets:
            continue
        ref = targets[asg_use.target_id]
        seq = read_seq_file(path)
        primer = asg_use.primer
        calls, ident, aln_len, ori, cands, spans = align_calls(seq, ref.sequence, primer=primer)
        match_n = sum(1 for pos, b in calls.items() if b == ref.sequence[pos])
        mism_n = len(calls) - match_n
        hc_bp = sum(c.insert_len for c in cands if c.high_confidence)
        key = group_key(asg_use.target_id, asg_use.clone_id)

        slot = by_clone.setdefault(
            key,
            {
                "clone_id": asg_use.clone_id,
                "target_id": asg_use.target_id,
                "ref_seq": ref.sequence,
                "primer_calls": {},
                "primer_files": {},
                "primer_candidates": {},
                "primer_no_insert": {},
            },
        )
        prev = slot["primer_calls"].get(primer)
        if prev is None or len(calls) > len(prev):
            slot["primer_calls"][primer] = calls
            slot["primer_files"][primer] = path.name
            slot["primer_candidates"][primer] = cands
            slot["primer_no_insert"][primer] = spans

        meta = {
            "file": path.name,
            "clone_id": asg_use.clone_id,
            "target_id": asg_use.target_id,
            "primer": primer,
            "orientation": ori,
            "aligned_len": aln_len,
            "identity": round(ident, 4),
            "match_bp": match_n,
            "mismatch_bp": mism_n,
            "internal_insert_bp": hc_bp,
            "insert_status": ",".join(
                sorted({STATUS_CANDIDATE for c in cands if c.high_confidence}) or [STATUS_NO_INSERT]
            ),
            "insert_evidence": ";".join(c.summary() for c in cands if c.high_confidence),
            "read_len": len(seq),
        }
        per_read.append(meta)
        file_to_meta[path.name] = meta

    clone_rows: List[dict] = []
    clone_status: Dict[str, str] = {}
    for key, slot in sorted(by_clone.items()):
        base = evaluate_base_union(slot["ref_seq"], slot["primer_calls"])
        ins = merge_clone_insertions(slot["primer_candidates"], slot["primer_no_insert"])
        status = stage1_status(base["bases_perfect"], ins.status)
        clone_status[key] = status
        primers = sorted(slot["primer_files"])
        clone_rows.append(
            {
                "target_id": slot["target_id"],
                "clone_id": slot["clone_id"],
                "primers": "|".join(primers),
                "files": "|".join(slot["primer_files"][p] for p in primers),
                **base,
                "internal_insert_bp": ins.internal_insert_bp,
                "insert_status": ins.status,
                "insert_evidence": ins.evidence,
                "perfect": base["bases_perfect"] and ins.status == STATUS_NO_INSERT,
                "stage1_status": status,
            }
        )

    seq_names = sorted({m["file"] for m in per_read})
    target_ids = sorted(targets)
    matrix_rows: List[dict] = []
    long_rows: List[dict] = []
    for tid in target_ids:
        row: dict = {"target_id": tid, "ref_len": len(targets[tid].sequence)}
        for fname in seq_names:
            meta = file_to_meta.get(fname)
            if meta is None or meta["target_id"] != tid:
                row[fname] = ""
                continue
            key = group_key(meta["target_id"], meta["clone_id"])
            status = clone_status[key]
            row[fname] = status
            long_rows.append(
                {
                    "target_id": tid,
                    "file": fname,
                    "clone_id": meta["clone_id"],
                    "primer": meta["primer"],
                    "clone_union_status": status,
                    "orientation": meta["orientation"],
                    "aligned_len": meta["aligned_len"],
                    "identity": meta["identity"],
                    "match_bp": meta["match_bp"],
                    "mismatch_bp": meta["mismatch_bp"],
                    "internal_insert_bp": meta["internal_insert_bp"],
                    "insert_status": meta["insert_status"],
                    "insert_evidence": meta["insert_evidence"],
                    "read_len": meta["read_len"],
                }
            )
        gene_clones = [r for r in clone_rows if r["target_id"] == tid]
        n_pass = sum(1 for r in gene_clones if r["stage1_status"] == "PASS")
        n_review = sum(1 for r in gene_clones if r["stage1_status"] == "PASS_INSERT_REVIEW")
        row["n_clones"] = len(gene_clones)
        row["n_pass_clones"] = n_pass
        row["target_not_failed"] = "YES" if (n_pass + n_review) > 0 else "NO"
        matrix_rows.append(row)

    worth_keys = {
        group_key(r["target_id"], r["clone_id"])
        for r in clone_rows
        if r["stage1_status"] in WORTH_STATUSES
    }
    files_worth = []
    for meta in per_read:
        key = group_key(meta["target_id"], meta["clone_id"])
        if key not in worth_keys:
            continue
        stem = Path(meta["file"]).stem
        ab1_name = stem + ".ab1"
        ab1_path = resolve_read_path(args.reads_dir, ab1_name, ".ab1")
        files_worth.append(
            {
                "seq_file": meta["file"],
                "ab1_file": ab1_path.name if ab1_path else ab1_name,
                "ab1_exists": "YES" if ab1_path else "NO",
                "clone_id": meta["clone_id"],
                "target_id": meta["target_id"],
                "primer": meta["primer"],
                "clone_union_status": clone_status[key],
                "worth_ab1": "YES",
            }
        )

    targets_view = []
    for tid in target_ids:
        gene_clones = [r for r in clone_rows if r["target_id"] == tid]
        n_pass = sum(1 for r in gene_clones if r["stage1_status"] == "PASS")
        n_review = sum(1 for r in gene_clones if r["stage1_status"] == "PASS_INSERT_REVIEW")
        targets_view.append(
            {
                "target_id": tid,
                "ref_len": len(targets[tid].sequence),
                "n_clones": len(gene_clones),
                "n_pass_clones": n_pass,
                "target_not_failed": "YES" if (n_pass + n_review) > 0 else "NO",
                "best_success_ratio": max((r["success_ratio"] for r in gene_clones), default=0.0),
                "pass_clone_ids": "|".join(
                    r["clone_id"] for r in gene_clones if r["stage1_status"] == "PASS"
                ),
            }
        )

    out = args.out_dir
    write_csv(
        out / "stage1_matrix.csv",
        matrix_rows,
        ["target_id", "ref_len", "n_clones", "n_pass_clones", "target_not_failed", *seq_names],
    )
    write_csv(
        out / "stage1_matrix_long.csv",
        long_rows,
        [
            "target_id",
            "file",
            "clone_id",
            "primer",
            "clone_union_status",
            "orientation",
            "aligned_len",
            "identity",
            "match_bp",
            "mismatch_bp",
            "internal_insert_bp",
            "insert_status",
            "insert_evidence",
            "read_len",
        ],
    )
    write_csv(
        out / "stage1_clone_union.csv",
        clone_rows,
        [
            "target_id",
            "clone_id",
            "primers",
            "files",
            "ref_len",
            "success_bp",
            "fail_bp",
            "covered_bp",
            "uncovered_bp",
            "mismatch_only_bp",
            "internal_insert_bp",
            "insert_status",
            "insert_evidence",
            "success_ratio",
            "coverage_ratio",
            "n_reads_used",
            "bases_perfect",
            "perfect",
            "stage1_status",
        ],
    )
    write_csv(
        out / "stage1_files_worth_ab1.csv",
        files_worth,
        [
            "seq_file",
            "ab1_file",
            "ab1_exists",
            "clone_id",
            "target_id",
            "primer",
            "clone_union_status",
            "worth_ab1",
        ],
    )
    write_csv(
        out / "stage1_targets_not_failed.csv",
        targets_view,
        [
            "target_id",
            "ref_len",
            "n_clones",
            "n_pass_clones",
            "target_not_failed",
            "best_success_ratio",
            "pass_clone_ids",
        ],
    )
    write_csv(
        out / "stage1_genes_not_failed.csv",
        [
            {
                "vendor_id": r["target_id"],
                "ref_len": r["ref_len"],
                "n_clones": r["n_clones"],
                "n_pass_clones": r["n_pass_clones"],
                "gene_not_failed": r["target_not_failed"],
                "best_success_ratio": r["best_success_ratio"],
                "pass_clone_ids": r["pass_clone_ids"],
            }
            for r in targets_view
        ],
        [
            "vendor_id",
            "ref_len",
            "n_clones",
            "n_pass_clones",
            "gene_not_failed",
            "best_success_ratio",
            "pass_clone_ids",
        ],
    )
    write_csv(
        out / "stage1_per_read.csv",
        per_read,
        [
            "file",
            "clone_id",
            "target_id",
            "primer",
            "orientation",
            "aligned_len",
            "identity",
            "match_bp",
            "mismatch_bp",
            "internal_insert_bp",
            "insert_status",
            "insert_evidence",
            "read_len",
        ],
    )

    n_pass = sum(1 for r in clone_rows if r["stage1_status"] == "PASS")
    n_review = sum(1 for r in clone_rows if r["stage1_status"] == "PASS_INSERT_REVIEW")
    n_ok = sum(1 for r in targets_view if r["target_not_failed"] == "YES")
    print(f"reads_dir={args.reads_dir}")
    print(f"targets={len(targets)} clones={len(clone_rows)} seq_reads={len(per_read)}")
    print(f"clone_PASS={n_pass} PASS_INSERT_REVIEW={n_review} FAIL={len(clone_rows) - n_pass - n_review}")
    print(f"targets_not_failed={n_ok}/{len(targets)}")
    print(f"files_worth_ab1={len(files_worth)}")
    print(f"out_dir={out}")
    if missing:
        print(f"missing_reads={len(missing)} e.g. {missing[:5]}")


if __name__ == "__main__":
    main()
