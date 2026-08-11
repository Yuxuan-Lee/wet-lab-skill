#!/usr/bin/env python3
"""Stage-2 Sanger screen (.ab1) → target × read matrix with peak purity.

Uses the same --targets / --assignments model as stage-1.
Only AB1 files listed in stage1_files_worth_ab1.csv are analyzed.
"""

from __future__ import annotations

import argparse
import csv
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

from Bio.Seq import Seq

warnings.filterwarnings("ignore", category=DeprecationWarning)
from Bio import SeqIO, pairwise2  # noqa: E402

from common import (
    group_key,
    index_assignments,
    load_assignments,
    load_targets,
    resolve_read_path,
    validate_assignments,
)

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def load_worth_ab1(path: Path) -> List[str]:
    names: List[str] = []
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            ab1 = (row.get("ab1_file") or "").strip()
            if ab1:
                names.append(Path(ab1).name)
            else:
                seq = (row.get("seq_file") or "").strip()
                if seq:
                    names.append(Path(seq).with_suffix(".ab1").name)
    return sorted(set(names))


def extract_ab1(path: Path) -> Tuple[str, List[Dict[str, int]]]:
    rec = SeqIO.read(str(path), "abi")
    raw = rec.annotations.get("abif_raw", {})
    fwo = raw.get("FWO_1", b"GATC")
    order = (fwo.decode("ascii", errors="ignore") if isinstance(fwo, bytes) else str(fwo)).upper()
    if len(order) != 4 or any(b not in order for b in "GATC"):
        order = "GATC"
    channels = [raw.get("DATA9"), raw.get("DATA10"), raw.get("DATA11"), raw.get("DATA12")]
    ploc = raw.get("PLOC2") or raw.get("PLOC1")
    pbas = raw.get("PBAS2") or raw.get("PBAS1")
    if not all(channels) or ploc is None or pbas is None:
        raise ValueError(f"missing chromatogram fields in {path.name}")
    bases = (pbas.decode("ascii", errors="ignore") if isinstance(pbas, bytes) else str(pbas)).upper()
    n = min(len(bases), len(ploc))
    trace_len = len(channels[0])
    heights: List[Dict[str, int]] = []
    for i in range(n):
        x = int(ploc[i])
        h = {b: 0 for b in "GATC"}
        if 0 <= x < trace_len:
            for dye, ch in zip(order, channels):
                h[dye] = int(ch[x])
        heights.append(h)
    return bases[:n], heights


def count_internal_inserts(ref_aln: str, query_aln: str, ref_len: int) -> int:
    ref_pos = 0
    inserts = 0
    for rc, qc in zip(ref_aln, query_aln):
        if rc == "-" and qc != "-" and 0 < ref_pos < ref_len:
            inserts += 1
        if rc != "-":
            ref_pos += 1
    return inserts


def trusted_internal_inserts(inserts: int, match_bp: int, ref_len: int, identity: float) -> int:
    if ref_len <= 0 or inserts <= 0:
        return 0
    if match_bp >= int(0.90 * ref_len) and identity >= 0.95:
        return inserts
    return 0


def peak_frac_for_ref_base(heights: Dict[str, int], ref_base: str, orientation: str) -> float:
    rb = ref_base.upper()
    if rb not in "GATC":
        return 0.0
    want = rb if orientation == "as_is" else rb.translate(COMPLEMENT)
    total = sum(heights.get(b, 0) for b in "GATC")
    if total <= 0:
        return 0.0
    return heights.get(want, 0) / total


def align_ab1_support(
    bases: str,
    heights: List[Dict[str, int]],
    ref_seq: str,
    peak_min: float,
) -> Tuple[Dict[int, str], Dict[int, float], float, int, str, int, int]:
    def one(query: str, orientation: str):
        aln = pairwise2.align.localms(
            ref_seq, query, 2.0, -1.0, -7.0, -1.0, one_alignment_only=True
        )
        if not aln:
            return {}, {}, 0.0, 0, 0, 0
        a = aln[0]
        calls: Dict[int, str] = {}
        fracs: Dict[int, float] = {}
        ref_pos = query_pos = 0
        match = aligned = support = 0
        n = len(bases)
        for rc, qc in zip(a.seqA, a.seqB):
            if rc != "-" and qc != "-":
                if 0 <= ref_pos < len(ref_seq) and 0 <= query_pos < n:
                    raw_idx = query_pos if orientation == "as_is" else (n - 1 - query_pos)
                    frac = peak_frac_for_ref_base(heights[raw_idx], rc, orientation)
                    calls[ref_pos] = qc
                    fracs[ref_pos] = frac
                    aligned += 1
                    if qc == rc:
                        match += 1
                        if frac >= peak_min:
                            support += 1
            if rc != "-":
                ref_pos += 1
            if qc != "-":
                query_pos += 1
        ident = (match / aligned) if aligned else 0.0
        raw_ins = count_internal_inserts(a.seqA, a.seqB, len(ref_seq))
        inserts = trusted_internal_inserts(raw_ins, match, len(ref_seq), ident)
        return calls, fracs, ident, aligned, support, inserts

    c1, f1, i1, a1, s1, n1 = one(bases, "as_is")
    c2, f2, i2, a2, s2, n2 = one(str(Seq(bases).reverse_complement()), "RC_of_read")
    if (s2, -n2, a2, i2) > (s1, -n1, a1, i1):
        return c2, f2, i2, a2, "RC_of_read", s2, n2
    return c1, f1, i1, a1, "as_is", s1, n1


def evaluate_union(
    ref_seq: str,
    primer_calls: Dict[str, Dict[int, str]],
    primer_fracs: Dict[str, Dict[int, float]],
    primer_inserts: Dict[str, int],
    peak_min: float,
) -> dict:
    ref_len = len(ref_seq)
    success = covered = mismatch_only = peak_fail_only = 0
    for pos in range(ref_len):
        present = called_ref = ok = False
        for primer, calls in primer_calls.items():
            if pos not in calls:
                continue
            present = True
            if calls[pos] == ref_seq[pos]:
                called_ref = True
                if primer_fracs[primer].get(pos, 0.0) >= peak_min:
                    ok = True
        if not present:
            continue
        covered += 1
        if ok:
            success += 1
        elif called_ref:
            peak_fail_only += 1
        else:
            mismatch_only += 1
    insert_bp = int(sum(primer_inserts.values()))
    return {
        "ref_len": ref_len,
        "success_bp": success,
        "fail_bp": ref_len - success,
        "covered_bp": covered,
        "uncovered_bp": ref_len - covered,
        "mismatch_only_bp": mismatch_only,
        "peak_fail_only_bp": peak_fail_only,
        "internal_insert_bp": insert_bp,
        "success_ratio": round(success / ref_len, 4) if ref_len else 0.0,
        "coverage_ratio": round(covered / ref_len, 4) if ref_len else 0.0,
        "n_reads_used": len(primer_calls),
        "perfect": success == ref_len and insert_bp == 0,
    }


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage-2 AB1 screen with peak gate")
    ap.add_argument("--targets", type=Path, required=True)
    ap.add_argument("--assignments", type=Path, required=True)
    ap.add_argument("--reads-dir", type=Path, required=True)
    ap.add_argument("--worth-csv", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--peak-min", type=float, default=0.75)
    args = ap.parse_args()

    targets = load_targets(args.targets)
    assignments = load_assignments(args.assignments)
    validate_assignments(targets, assignments)
    asg_by_file = index_assignments(assignments)
    worth = load_worth_ab1(args.worth_csv)

    by_clone: Dict[str, dict] = {}
    per_read: List[dict] = []
    file_to_meta: Dict[str, dict] = {}
    missing: List[str] = []
    errors: List[str] = []

    for name in worth:
        path = resolve_read_path(args.reads_dir, name, ".ab1")
        if path is None:
            missing.append(name)
            continue
        asg = asg_by_file.get(path.name) or asg_by_file.get(name) or asg_by_file.get(Path(name).with_suffix(".seq").name)
        if asg is None:
            # try stem match against assignment files
            stem = path.stem.lower()
            asg = next((a for a in assignments if Path(a.file).stem.lower() == stem), None)
        if asg is None:
            errors.append(f"unassigned:{path.name}")
            continue
        ref = targets[asg.target_id]
        try:
            bases, heights = extract_ab1(path)
        except Exception as e:  # noqa: BLE001
            errors.append(f"parse:{path.name}:{e}")
            continue

        calls, fracs, ident, aln_len, ori, support_n, inserts = align_ab1_support(
            bases, heights, ref.sequence, args.peak_min
        )
        match_n = sum(1 for pos, b in calls.items() if b == ref.sequence[pos])
        peak_ok_n = sum(
            1
            for pos, b in calls.items()
            if b == ref.sequence[pos] and fracs.get(pos, 0.0) >= args.peak_min
        )
        key = group_key(asg.target_id, asg.clone_id)
        primer = asg.primer
        slot = by_clone.setdefault(
            key,
            {
                "clone_id": asg.clone_id,
                "target_id": asg.target_id,
                "ref_seq": ref.sequence,
                "primer_calls": {},
                "primer_fracs": {},
                "primer_files": {},
                "primer_support": {},
                "primer_inserts": {},
            },
        )
        if support_n > slot["primer_support"].get(primer, -1):
            slot["primer_calls"][primer] = calls
            slot["primer_fracs"][primer] = fracs
            slot["primer_files"][primer] = path.name
            slot["primer_support"][primer] = support_n
            slot["primer_inserts"][primer] = inserts

        meta = {
            "file": path.name,
            "clone_id": asg.clone_id,
            "target_id": asg.target_id,
            "primer": primer,
            "orientation": ori,
            "aligned_len": aln_len,
            "identity": round(ident, 4),
            "match_bp": match_n,
            "peak_ok_bp": peak_ok_n,
            "support_bp": support_n,
            "internal_insert_bp": inserts,
            "read_len": len(bases),
            "peak_min": args.peak_min,
        }
        per_read.append(meta)
        file_to_meta[path.name] = meta

    clone_rows: List[dict] = []
    clone_status: Dict[str, str] = {}
    for key, slot in sorted(by_clone.items()):
        stats = evaluate_union(
            slot["ref_seq"],
            slot["primer_calls"],
            slot["primer_fracs"],
            slot["primer_inserts"],
            args.peak_min,
        )
        status = "PASS" if stats["perfect"] else "FAIL"
        clone_status[key] = status
        primers = sorted(slot["primer_files"])
        clone_rows.append(
            {
                "target_id": slot["target_id"],
                "clone_id": slot["clone_id"],
                "primers": "|".join(primers),
                "files": "|".join(slot["primer_files"][p] for p in primers),
                **stats,
                "peak_min": args.peak_min,
                "stage2_status": status,
            }
        )

    ab1_names = worth
    target_ids = sorted(targets)
    matrix_rows: List[dict] = []
    long_rows: List[dict] = []
    for tid in target_ids:
        row: dict = {"target_id": tid, "ref_len": len(targets[tid].sequence)}
        for fname in ab1_names:
            meta = file_to_meta.get(fname) or file_to_meta.get(Path(fname).name)
            # also try resolved actual names
            if meta is None:
                for m in per_read:
                    if m["file"] == fname or Path(m["file"]).stem == Path(fname).stem:
                        meta = m
                        break
            if meta is None or meta["target_id"] != tid:
                row[fname] = ""
                continue
            key = group_key(meta["target_id"], meta["clone_id"])
            status = clone_status.get(key, "FAIL")
            row[fname] = status
            long_rows.append(
                {
                    "target_id": tid,
                    "file": meta["file"],
                    "clone_id": meta["clone_id"],
                    "primer": meta["primer"],
                    "clone_union_status": status,
                    "orientation": meta["orientation"],
                    "aligned_len": meta["aligned_len"],
                    "identity": meta["identity"],
                    "match_bp": meta["match_bp"],
                    "peak_ok_bp": meta["peak_ok_bp"],
                    "support_bp": meta["support_bp"],
                    "internal_insert_bp": meta["internal_insert_bp"],
                    "peak_min": args.peak_min,
                }
            )
        gene_clones = [r for r in clone_rows if r["target_id"] == tid]
        n_pass = sum(1 for r in gene_clones if r["stage2_status"] == "PASS")
        row["n_clones"] = len(gene_clones)
        row["n_pass_clones"] = n_pass
        row["target_not_failed"] = "YES" if n_pass > 0 else "NO"
        matrix_rows.append(row)

    pass_keys = {
        group_key(r["target_id"], r["clone_id"])
        for r in clone_rows
        if r["stage2_status"] == "PASS"
    }
    files_pass = [
        {
            "ab1_file": m["file"],
            "clone_id": m["clone_id"],
            "target_id": m["target_id"],
            "primer": m["primer"],
            "clone_union_status": "PASS",
        }
        for m in per_read
        if group_key(m["target_id"], m["clone_id"]) in pass_keys
    ]
    targets_view = []
    for tid in target_ids:
        gene_clones = [r for r in clone_rows if r["target_id"] == tid]
        n_pass = sum(1 for r in gene_clones if r["stage2_status"] == "PASS")
        targets_view.append(
            {
                "target_id": tid,
                "ref_len": len(targets[tid].sequence),
                "n_clones": len(gene_clones),
                "n_pass_clones": n_pass,
                "target_not_failed": "YES" if n_pass > 0 else "NO",
                "best_success_ratio": max((r["success_ratio"] for r in gene_clones), default=0.0),
                "pass_clone_ids": "|".join(
                    r["clone_id"] for r in gene_clones if r["stage2_status"] == "PASS"
                ),
            }
        )

    out = args.out_dir
    write_csv(
        out / "stage2_matrix.csv",
        matrix_rows,
        ["target_id", "ref_len", "n_clones", "n_pass_clones", "target_not_failed", *ab1_names],
    )
    write_csv(
        out / "stage2_matrix_long.csv",
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
            "peak_ok_bp",
            "support_bp",
            "internal_insert_bp",
            "peak_min",
        ],
    )
    write_csv(
        out / "stage2_clone_union.csv",
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
            "peak_fail_only_bp",
            "internal_insert_bp",
            "success_ratio",
            "coverage_ratio",
            "n_reads_used",
            "perfect",
            "peak_min",
            "stage2_status",
        ],
    )
    write_csv(
        out / "stage2_files_pass.csv",
        files_pass,
        ["ab1_file", "clone_id", "target_id", "primer", "clone_union_status"],
    )
    write_csv(
        out / "stage2_targets_not_failed.csv",
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
        out / "stage2_genes_not_failed.csv",
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
        out / "stage2_per_read.csv",
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
            "peak_ok_bp",
            "support_bp",
            "internal_insert_bp",
            "read_len",
            "peak_min",
        ],
    )

    n_pass = sum(1 for r in clone_rows if r["stage2_status"] == "PASS")
    n_ok = sum(1 for r in targets_view if r["target_not_failed"] == "YES")
    print(f"reads_dir={args.reads_dir}")
    print(f"worth={len(worth)} parsed={len(per_read)} missing={len(missing)} errors={len(errors)}")
    print(f"clones={len(clone_rows)} PASS={n_pass} FAIL={len(clone_rows) - n_pass}")
    print(f"targets_not_failed={n_ok}/{len(target_ids)}")
    print(f"peak_min={args.peak_min} out_dir={out}")
    if missing[:3]:
        print(f"missing_e.g.={missing[:3]}")
    if errors[:3]:
        print(f"errors_e.g.={errors[:3]}")


if __name__ == "__main__":
    main()
