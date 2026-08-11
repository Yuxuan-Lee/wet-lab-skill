#!/usr/bin/env python3
"""Stage-1 Sanger screen (.seq only) → gene × sequencing-file matrix.

Reference per gene: flank + CDS + flank (`full_insert_dna`).

Per clone, FW/RV (T7/T7T) are a UNION over reference bases:
  SUCCESS — any read covers the site AND calls the reference base
  FAILURE — uncovered, or covered but no read matches

Clone PASS ⇔ every reference base is SUCCESS.
Gene not-failed ⇔ at least one clone PASS.
Files worth AB1 ⇔ all .seq/.ab1 belonging to PASS clones.

Outputs (always written together):
  stage1_matrix.csv           rows=genes, cols=seq files, cell=PASS|FAIL|empty
  stage1_matrix_long.csv      tidy long form of the same matrix
  stage1_clone_union.csv      clone-level union stats
  stage1_files_worth_ab1.csv  sequencing files for stage-2
  stage1_genes_not_failed.csv genes with ≥1 PASS clone
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import warnings
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from Bio.Seq import Seq

warnings.filterwarnings("ignore", category=DeprecationWarning)
from Bio import pairwise2  # noqa: E402

FILE_RE = re.compile(
    r"^S\d+_(?P<sample>[^_]+)_(?P<primer>T7T?|t7t?)_(?P<well>[^.]+)\.(?P<ext>seq|ab1)$",
    re.IGNORECASE,
)


@dataclass
class RefRecord:
    vendor_id: str
    full_insert_dna: str


def load_refs(path: Path) -> Dict[str, RefRecord]:
    refs: Dict[str, RefRecord] = {}
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            vid = row["vendor_id"].strip()
            dna = row["full_insert_dna"].strip().upper()
            if not vid or not dna:
                raise SystemExit(f"missing vendor_id/full_insert_dna in {path}")
            refs[vid] = RefRecord(vid, dna)
    if not refs:
        raise SystemExit(f"no refs in {path}")
    return refs


def load_aliases(path: Optional[Path]) -> Dict[str, str]:
    if path is None or not path.is_file():
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    aliases = data.get("aliases", data)
    return {str(k): str(v) for k, v in aliases.items()}


def build_her2_prefix_map(refs: Dict[str, RefRecord]) -> Dict[str, str]:
    """h2391-1 → 2391; h2950-1/h2950-2 → 2950 / 29502; solo h6277-2 → 6277."""
    by_task: Dict[str, List[int]] = defaultdict(list)
    for vid in refs:
        m = re.match(r"^h(\d+)-(\d+)$", vid)
        if not m:
            raise SystemExit(
                f"vendor_id {vid!r} is not h<task>-<seq>; "
                "use matching refs or extend naming in stage1_seq_matrix.py"
            )
        by_task[m.group(1)].append(int(m.group(2)))

    prefix_map: Dict[str, str] = {}
    for task, seqs in by_task.items():
        seqs = sorted(set(seqs))
        if len(seqs) == 1:
            prefix_map[task] = f"h{task}-{seqs[0]}"
        else:
            for s in seqs:
                if s == 1:
                    prefix_map[task] = f"h{task}-1"
                else:
                    prefix_map[f"{task}{s}"] = f"h{task}-{s}"
    return prefix_map


def resolve_vendor(
    sample: str,
    prefix_map: Dict[str, str],
    refs: Dict[str, RefRecord],
    aliases: Dict[str, str],
) -> Tuple[Optional[RefRecord], str, str]:
    if "-" not in sample:
        return None, sample, "bad_sample_name"
    prefix, clone_n = sample.rsplit("-", 1)
    if not clone_n.isdigit():
        return None, sample, "bad_clone_suffix"

    if prefix in aliases:
        vid = aliases[prefix]
        if vid not in refs:
            return None, sample, f"alias_missing_ref:{vid}"
        return refs[vid], sample, f"alias:{prefix}->{vid}"

    vid = prefix_map.get(prefix)
    if vid is None:
        return None, sample, "unmapped"
    return refs[vid], sample, "direct"


def read_seq_file(p: Path) -> str:
    lines = [
        x.strip().upper()
        for x in p.read_text(encoding="utf-8", errors="ignore").splitlines()
    ]
    return "".join(x for x in lines if x and not x.startswith(">"))


def count_internal_inserts(ref_aln: str, query_aln: str, ref_len: int) -> int:
    """Query bases inside the reference span that do not map to any ref position.

    Leading/trailing vector sequence (before ref starts or after ref ends) is ignored.
    Internal query-only columns are insertions (or unmapped extras) → counted.
    """
    ref_pos = 0
    inserts = 0
    for rc, qc in zip(ref_aln, query_aln):
        if rc == "-" and qc != "-" and 0 < ref_pos < ref_len:
            inserts += 1
        if rc != "-":
            ref_pos += 1
    return inserts


def trusted_internal_inserts(inserts: int, match_bp: int, ref_len: int, identity: float) -> int:
    """Only trust insert calls from near-complete, high-identity alignments.

    Poor opposite-primer alignments often create gapped artifacts; those must not
    veto a clean forward (or reverse) read. True internal insertions still appear
    on a high-identity full-span alignment.
    """
    if ref_len <= 0 or inserts <= 0:
        return 0
    if match_bp >= int(0.90 * ref_len) and identity >= 0.95:
        return inserts
    return 0


def align_calls(read_seq: str, ref_seq: str) -> Tuple[Dict[int, str], float, int, str, int]:
    """Best orientation local alignment → ref_pos → base call (+ trusted inserts)."""

    def one(q: str) -> Tuple[Dict[int, str], float, int, int]:
        aln = pairwise2.align.localms(
            ref_seq, q, 2.0, -1.0, -7.0, -1.0, one_alignment_only=True
        )
        if not aln:
            return {}, 0.0, 0, 0
        a = aln[0]
        # pairwise2 returns full gapped strings; walk from ref_pos=0.
        # Do NOT use a.start (it double-counts leading gaps in seqA).
        calls: Dict[int, str] = {}
        ref_pos = 0
        match = 0
        aligned = 0
        for rc, qc in zip(a.seqA, a.seqB):
            if rc != "-" and qc != "-" and 0 <= ref_pos < len(ref_seq):
                calls[ref_pos] = qc
                aligned += 1
                if rc == qc:
                    match += 1
            if rc != "-":
                ref_pos += 1
        ident = (match / aligned) if aligned else 0.0
        raw_inserts = count_internal_inserts(a.seqA, a.seqB, len(ref_seq))
        inserts = trusted_internal_inserts(raw_inserts, match, len(ref_seq), ident)
        return calls, ident, aligned, inserts

    c1, i1, a1, n1 = one(read_seq)
    c2, i2, a2, n2 = one(str(Seq(read_seq).reverse_complement()))
    m1 = sum(1 for pos, b in c1.items() if b == ref_seq[pos])
    m2 = sum(1 for pos, b in c2.items() if b == ref_seq[pos])
    # Prefer more matches, fewer internal inserts, then longer alignment.
    if (m2, -n2, a2, i2) > (m1, -n1, a1, i1):
        return c2, i2, a2, "RC_of_read", n2
    return c1, i1, a1, "as_is", n1


def evaluate_union(
    ref_seq: str,
    primer_calls: Dict[str, Dict[int, str]],
    primer_inserts: Dict[str, int],
) -> dict:
    ref_len = len(ref_seq)
    pos_bases: Dict[int, Set[str]] = defaultdict(set)
    covered: Set[int] = set()
    for calls in primer_calls.values():
        for pos, b in calls.items():
            covered.add(pos)
            pos_bases[pos].add(b)

    success = 0
    mismatch_only = 0
    for pos in range(ref_len):
        bases = pos_bases.get(pos, set())
        if ref_seq[pos] in bases:
            success += 1
        elif bases:
            mismatch_only += 1

    uncovered = ref_len - len(covered)
    fail = ref_len - success
    insert_bp = int(sum(primer_inserts.values()))
    return {
        "ref_len": ref_len,
        "success_bp": success,
        "fail_bp": fail,
        "covered_bp": len(covered),
        "uncovered_bp": uncovered,
        "mismatch_only_bp": mismatch_only,
        "internal_insert_bp": insert_bp,
        "success_ratio": round(success / ref_len, 4) if ref_len else 0.0,
        "coverage_ratio": round(len(covered) / ref_len, 4) if ref_len else 0.0,
        "perfect": success == ref_len and insert_bp == 0,
    }


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Stage-1 .seq screen → gene×file matrix + two maintenance views"
    )
    ap.add_argument("--seq-dir", type=Path, required=True, help="folder with *.seq")
    ap.add_argument(
        "--insert-csv",
        type=Path,
        required=True,
        help="pcr_insert_flanks.csv (vendor_id, full_insert_dna)",
    )
    ap.add_argument(
        "--aliases",
        type=Path,
        default=None,
        help="JSON aliases for sample prefix → vendor_id",
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    refs = load_refs(args.insert_csv)
    aliases = load_aliases(args.aliases)
    prefix_map = build_her2_prefix_map(refs)
    gene_ids = sorted(refs)

    seq_files = sorted(args.seq_dir.glob("*.seq"))
    if not seq_files:
        raise SystemExit(f"no .seq files in {args.seq_dir}")

    # clone_id -> slot
    by_clone: Dict[str, dict] = {}
    per_read: List[dict] = []
    unmapped: List[str] = []
    file_to_meta: Dict[str, dict] = {}

    for p in seq_files:
        m = FILE_RE.match(p.name)
        if not m:
            unmapped.append(p.name)
            continue
        sample = m.group("sample")
        primer = m.group("primer").upper()
        ref, clone_id, map_note = resolve_vendor(sample, prefix_map, refs, aliases)
        if ref is None:
            unmapped.append(p.name)
            continue

        seq = read_seq_file(p)
        calls, ident, aln_len, ori, inserts = align_calls(seq, ref.full_insert_dna)
        match_n = sum(1 for pos, b in calls.items() if b == ref.full_insert_dna[pos])
        mism_n = len(calls) - match_n

        slot = by_clone.setdefault(
            clone_id,
            {
                "clone_id": clone_id,
                "vendor_id": ref.vendor_id,
                "map_note": map_note,
                "ref_seq": ref.full_insert_dna,
                "primer_calls": {},
                "primer_files": {},
                "primer_inserts": {},
            },
        )
        prev = slot["primer_calls"].get(primer)
        if prev is None or len(calls) > len(prev):
            slot["primer_calls"][primer] = calls
            slot["primer_files"][primer] = p.name
            slot["primer_inserts"][primer] = inserts

        meta = {
            "file": p.name,
            "clone_id": clone_id,
            "vendor_id": ref.vendor_id,
            "primer": primer,
            "map_note": map_note,
            "orientation": ori,
            "aligned_len": aln_len,
            "identity": round(ident, 4),
            "match_bp": match_n,
            "mismatch_bp": mism_n,
            "internal_insert_bp": inserts,
            "read_len": len(seq),
        }
        per_read.append(meta)
        file_to_meta[p.name] = meta

    clone_rows: List[dict] = []
    clone_status: Dict[str, str] = {}
    for clone_id, slot in sorted(by_clone.items()):
        stats = evaluate_union(slot["ref_seq"], slot["primer_calls"], slot["primer_inserts"])
        status = "PASS" if stats["perfect"] else "FAIL"
        clone_status[clone_id] = status
        primers = sorted(slot["primer_files"])
        clone_rows.append(
            {
                "clone_id": clone_id,
                "vendor_id": slot["vendor_id"],
                "map_note": slot["map_note"],
                "primers": "|".join(primers),
                "files": "|".join(slot["primer_files"][x] for x in primers),
                **stats,
                "stage1_status": status,
            }
        )

    # ---- matrix: rows=genes, cols=seq files ----
    file_names = [p.name for p in seq_files]
    matrix_rows: List[dict] = []
    long_rows: List[dict] = []
    for gene in gene_ids:
        row: dict = {"vendor_id": gene, "ref_len": len(refs[gene].full_insert_dna)}
        for fname in file_names:
            meta = file_to_meta.get(fname)
            if meta is None or meta["vendor_id"] != gene:
                row[fname] = ""
                continue
            status = clone_status[meta["clone_id"]]
            row[fname] = status
            long_rows.append(
                {
                    "vendor_id": gene,
                    "file": fname,
                    "clone_id": meta["clone_id"],
                    "primer": meta["primer"],
                    "clone_union_status": status,
                    "map_note": meta["map_note"],
                    "orientation": meta["orientation"],
                    "aligned_len": meta["aligned_len"],
                    "identity": meta["identity"],
                    "match_bp": meta["match_bp"],
                    "mismatch_bp": meta["mismatch_bp"],
                    "internal_insert_bp": meta["internal_insert_bp"],
                    "read_len": meta["read_len"],
                }
            )
        # gene-level rollup columns (left side helpers)
        gene_clones = [r for r in clone_rows if r["vendor_id"] == gene]
        n_pass = sum(1 for r in gene_clones if r["stage1_status"] == "PASS")
        row["n_clones"] = len(gene_clones)
        row["n_pass_clones"] = n_pass
        row["gene_not_failed"] = "YES" if n_pass > 0 else "NO"
        matrix_rows.append(row)

    # reorder matrix columns: helpers first, then files
    matrix_fields = [
        "vendor_id",
        "ref_len",
        "n_clones",
        "n_pass_clones",
        "gene_not_failed",
        *file_names,
    ]

    # ---- view 1: files worth AB1 ----
    pass_clones = {r["clone_id"] for r in clone_rows if r["stage1_status"] == "PASS"}
    files_worth: List[dict] = []
    for meta in per_read:
        if meta["clone_id"] not in pass_clones:
            continue
        stem = Path(meta["file"]).stem
        ab1_name = stem + ".ab1"
        ab1_path = args.seq_dir / ab1_name
        files_worth.append(
            {
                "seq_file": meta["file"],
                "ab1_file": ab1_name if ab1_path.is_file() else "",
                "ab1_exists": "YES" if ab1_path.is_file() else "NO",
                "clone_id": meta["clone_id"],
                "vendor_id": meta["vendor_id"],
                "primer": meta["primer"],
                "clone_union_status": "PASS",
                "worth_ab1": "YES",
            }
        )

    # ---- view 2: genes not failed ----
    genes_view: List[dict] = []
    for gene in gene_ids:
        gene_clones = [r for r in clone_rows if r["vendor_id"] == gene]
        n_pass = sum(1 for r in gene_clones if r["stage1_status"] == "PASS")
        best = max((r["success_ratio"] for r in gene_clones), default=0.0)
        pass_ids = [r["clone_id"] for r in gene_clones if r["stage1_status"] == "PASS"]
        genes_view.append(
            {
                "vendor_id": gene,
                "ref_len": len(refs[gene].full_insert_dna),
                "n_clones": len(gene_clones),
                "n_pass_clones": n_pass,
                "gene_not_failed": "YES" if n_pass > 0 else "NO",
                "best_success_ratio": best,
                "pass_clone_ids": "|".join(pass_ids),
            }
        )

    out = args.out_dir
    write_csv(out / "stage1_matrix.csv", matrix_rows, matrix_fields)
    write_csv(
        out / "stage1_matrix_long.csv",
        long_rows,
        [
            "vendor_id",
            "file",
            "clone_id",
            "primer",
            "clone_union_status",
            "map_note",
            "orientation",
            "aligned_len",
            "identity",
            "match_bp",
            "mismatch_bp",
            "internal_insert_bp",
            "read_len",
        ],
    )
    write_csv(
        out / "stage1_clone_union.csv",
        clone_rows,
        [
            "clone_id",
            "vendor_id",
            "map_note",
            "primers",
            "files",
            "ref_len",
            "success_bp",
            "fail_bp",
            "covered_bp",
            "uncovered_bp",
            "mismatch_only_bp",
            "internal_insert_bp",
            "success_ratio",
            "coverage_ratio",
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
            "vendor_id",
            "primer",
            "clone_union_status",
            "worth_ab1",
        ],
    )
    write_csv(
        out / "stage1_genes_not_failed.csv",
        genes_view,
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
            "vendor_id",
            "primer",
            "map_note",
            "orientation",
            "aligned_len",
            "identity",
            "match_bp",
            "mismatch_bp",
            "internal_insert_bp",
            "read_len",
        ],
    )

    n_pass_clones = sum(1 for r in clone_rows if r["stage1_status"] == "PASS")
    n_genes_ok = sum(1 for r in genes_view if r["gene_not_failed"] == "YES")
    print(f"seq_dir={args.seq_dir}")
    print(f"genes={len(gene_ids)} clones={len(clone_rows)} seq_files={len(seq_files)}")
    print(f"clone_PASS={n_pass_clones} clone_FAIL={len(clone_rows) - n_pass_clones}")
    print(f"genes_not_failed={n_genes_ok}/{len(gene_ids)}")
    print(f"files_worth_ab1={len(files_worth)}")
    print(f"out_dir={out}")
    if unmapped:
        print(f"unmapped_files={len(unmapped)} e.g. {unmapped[:3]}")


if __name__ == "__main__":
    main()
