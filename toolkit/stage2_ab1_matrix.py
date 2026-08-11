#!/usr/bin/env python3
"""Stage-2 Sanger screen (.ab1) → gene × AB1-file matrix.

Same FW∪RV union logic as stage-1, plus peak gate:
  A reaction supports a reference base only if it calls that base AND
  correct-base channel / sum(G+A+T+C) at PLOC >= peak_min (default 0.75).
  Channel order from FWO_1; signal = trace height at basecall locus.

Only AB1 files in stage1_files_worth_ab1.csv are analyzed.
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
from Bio import SeqIO, pairwise2  # noqa: E402

FILE_RE = re.compile(
    r"^S\d+_(?P<sample>[^_]+)_(?P<primer>T7T?|t7t?)_(?P<well>[^.]+)\.(?P<ext>seq|ab1)$",
    re.IGNORECASE,
)
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


@dataclass
class RefRecord:
    vendor_id: str
    full_insert_dna: str


def load_refs(path: Path) -> Dict[str, RefRecord]:
    refs: Dict[str, RefRecord] = {}
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            refs[row["vendor_id"].strip()] = RefRecord(
                row["vendor_id"].strip(), row["full_insert_dna"].strip().upper()
            )
    return refs


def load_aliases(path: Optional[Path]) -> Dict[str, str]:
    if path is None or not path.is_file():
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    aliases = data.get("aliases", data)
    return {str(k): str(v) for k, v in aliases.items()}


def build_her2_prefix_map(refs: Dict[str, RefRecord]) -> Dict[str, str]:
    by_task: Dict[str, List[int]] = defaultdict(list)
    for vid in refs:
        m = re.match(r"^h(\d+)-(\d+)$", vid)
        if not m:
            raise SystemExit(f"unexpected vendor_id: {vid}")
        by_task[m.group(1)].append(int(m.group(2)))
    prefix_map: Dict[str, str] = {}
    for task, seqs in by_task.items():
        seqs = sorted(set(seqs))
        if len(seqs) == 1:
            prefix_map[task] = f"h{task}-{seqs[0]}"
        else:
            for s in seqs:
                prefix_map[task if s == 1 else f"{task}{s}"] = f"h{task}-{s}"
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


def load_worth_ab1(path: Path) -> List[str]:
    names: List[str] = []
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            ab1 = (row.get("ab1_file") or "").strip()
            if ab1:
                names.append(ab1)
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
    """Query bases inside the reference span that do not map to any ref position."""
    ref_pos = 0
    inserts = 0
    for rc, qc in zip(ref_aln, query_aln):
        if rc == "-" and qc != "-" and 0 < ref_pos < ref_len:
            inserts += 1
        if rc != "-":
            ref_pos += 1
    return inserts


def trusted_internal_inserts(inserts: int, match_bp: int, ref_len: int, identity: float) -> int:
    """Only trust insert calls from near-complete, high-identity alignments."""
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
        ref_pos = 0
        query_pos = 0
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
        raw_inserts = count_internal_inserts(a.seqA, a.seqB, len(ref_seq))
        inserts = trusted_internal_inserts(raw_inserts, match, len(ref_seq), ident)
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
        present = False
        called_ref = False
        ok = False
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
        "perfect": success == ref_len and insert_bp == 0,
    }


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage-2 AB1 screen with peak-fraction gate")
    ap.add_argument("--ab1-dir", type=Path, required=True)
    ap.add_argument("--insert-csv", type=Path, required=True)
    ap.add_argument("--worth-csv", type=Path, required=True)
    ap.add_argument("--aliases", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--peak-min", type=float, default=0.75)
    args = ap.parse_args()

    refs = load_refs(args.insert_csv)
    aliases = load_aliases(args.aliases)
    prefix_map = build_her2_prefix_map(refs)
    gene_ids = sorted(refs)
    worth = load_worth_ab1(args.worth_csv)
    if not worth:
        raise SystemExit(f"no AB1 names in {args.worth_csv}")

    by_clone: Dict[str, dict] = {}
    per_read: List[dict] = []
    file_to_meta: Dict[str, dict] = {}
    missing: List[str] = []
    errors: List[str] = []

    for name in worth:
        p = args.ab1_dir / name
        if not p.is_file():
            missing.append(name)
            continue
        m = FILE_RE.match(name)
        if not m:
            errors.append(f"bad_name:{name}")
            continue
        ref, clone_id, map_note = resolve_vendor(
            m.group("sample"), prefix_map, refs, aliases
        )
        if ref is None:
            errors.append(f"unmapped:{name}")
            continue
        primer = m.group("primer").upper()
        try:
            bases, heights = extract_ab1(p)
        except Exception as e:  # noqa: BLE001
            errors.append(f"parse:{name}:{e}")
            continue

        calls, fracs, ident, aln_len, ori, support_n, inserts = align_ab1_support(
            bases, heights, ref.full_insert_dna, args.peak_min
        )
        match_n = sum(1 for pos, b in calls.items() if b == ref.full_insert_dna[pos])
        peak_ok_n = sum(
            1
            for pos, b in calls.items()
            if b == ref.full_insert_dna[pos] and fracs.get(pos, 0.0) >= args.peak_min
        )

        slot = by_clone.setdefault(
            clone_id,
            {
                "clone_id": clone_id,
                "vendor_id": ref.vendor_id,
                "map_note": map_note,
                "ref_seq": ref.full_insert_dna,
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
            slot["primer_files"][primer] = name
            slot["primer_support"][primer] = support_n
            slot["primer_inserts"][primer] = inserts

        meta = {
            "file": name,
            "clone_id": clone_id,
            "vendor_id": ref.vendor_id,
            "primer": primer,
            "map_note": map_note,
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
        file_to_meta[name] = meta

    clone_rows: List[dict] = []
    clone_status: Dict[str, str] = {}
    for clone_id, slot in sorted(by_clone.items()):
        stats = evaluate_union(
            slot["ref_seq"],
            slot["primer_calls"],
            slot["primer_fracs"],
            slot["primer_inserts"],
            args.peak_min,
        )
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
                "peak_min": args.peak_min,
                "stage2_status": status,
            }
        )

    matrix_rows: List[dict] = []
    long_rows: List[dict] = []
    for gene in gene_ids:
        row: dict = {"vendor_id": gene, "ref_len": len(refs[gene].full_insert_dna)}
        for fname in worth:
            meta = file_to_meta.get(fname)
            if meta is None or meta["vendor_id"] != gene:
                row[fname] = ""
                continue
            status = clone_status.get(meta["clone_id"], "FAIL")
            row[fname] = status
            long_rows.append(
                {
                    "vendor_id": gene,
                    "file": fname,
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
        gene_clones = [r for r in clone_rows if r["vendor_id"] == gene]
        n_pass = sum(1 for r in gene_clones if r["stage2_status"] == "PASS")
        row["n_clones"] = len(gene_clones)
        row["n_pass_clones"] = n_pass
        row["gene_not_failed"] = "YES" if n_pass > 0 else "NO"
        matrix_rows.append(row)

    pass_clones = {r["clone_id"] for r in clone_rows if r["stage2_status"] == "PASS"}
    files_pass = [
        {
            "ab1_file": m["file"],
            "clone_id": m["clone_id"],
            "vendor_id": m["vendor_id"],
            "primer": m["primer"],
            "clone_union_status": "PASS",
        }
        for m in per_read
        if m["clone_id"] in pass_clones
    ]
    genes_view = []
    for gene in gene_ids:
        gene_clones = [r for r in clone_rows if r["vendor_id"] == gene]
        n_pass = sum(1 for r in gene_clones if r["stage2_status"] == "PASS")
        genes_view.append(
            {
                "vendor_id": gene,
                "ref_len": len(refs[gene].full_insert_dna),
                "n_clones": len(gene_clones),
                "n_pass_clones": n_pass,
                "gene_not_failed": "YES" if n_pass > 0 else "NO",
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
        ["vendor_id", "ref_len", "n_clones", "n_pass_clones", "gene_not_failed", *worth],
    )
    write_csv(
        out / "stage2_matrix_long.csv",
        long_rows,
        [
            "vendor_id",
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
            "peak_fail_only_bp",
            "internal_insert_bp",
            "success_ratio",
            "coverage_ratio",
            "perfect",
            "peak_min",
            "stage2_status",
        ],
    )
    write_csv(
        out / "stage2_files_pass.csv",
        files_pass,
        ["ab1_file", "clone_id", "vendor_id", "primer", "clone_union_status"],
    )
    write_csv(
        out / "stage2_genes_not_failed.csv",
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
        out / "stage2_per_read.csv",
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
            "peak_ok_bp",
            "support_bp",
            "internal_insert_bp",
            "read_len",
            "peak_min",
        ],
    )

    n_pass = sum(1 for r in clone_rows if r["stage2_status"] == "PASS")
    n_genes = sum(1 for r in genes_view if r["gene_not_failed"] == "YES")
    print(f"ab1_dir={args.ab1_dir}")
    print(f"worth={len(worth)} parsed={len(per_read)} missing={len(missing)} errors={len(errors)}")
    print(f"clones={len(clone_rows)} PASS={n_pass} FAIL={len(clone_rows) - n_pass}")
    print(f"genes_not_failed={n_genes}/{len(gene_ids)}")
    print(f"peak_min={args.peak_min} out_dir={out}")
    if missing[:3]:
        print(f"missing_e.g.={missing[:3]}")
    if errors[:3]:
        print(f"errors_e.g.={errors[:3]}")


if __name__ == "__main__":
    main()
