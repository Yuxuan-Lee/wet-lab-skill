#!/usr/bin/env python3
"""Export final target × clone/sample matrix for collecting return plasmids.

One sample/clone = all primers assigned to that clone (2, 3, 4, ...), not one file.
"""

from __future__ import annotations

import argparse
import csv
import html
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set


def load_clones(path: Path, status_key: str) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            tid = (row.get("target_id") or row.get("vendor_id") or "").strip()
            cid = (row.get("clone_id") or "").strip()
            key = f"{tid}::{cid}"
            row["_status"] = row[status_key]
            row["_target_id"] = tid
            row["_clone_id"] = cid
            out[key] = row
    return out


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def sample_status(s1: dict | None, s2: dict | None) -> str:
    if s2 is not None and s2["_status"] == "PASS":
        return "CORRECT"
    if s1 is not None and s1["_status"] == "PASS":
        return "SEQ_OK_AB1_FAIL"
    if s1 is not None or s2 is not None:
        return "FAIL"
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description="Final target×sample matrix for plasmid return")
    ap.add_argument("--stage1-clones", type=Path, required=True)
    ap.add_argument("--stage2-clones", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument(
        "--targets",
        type=Path,
        default=None,
        help="optional targets CSV/FASTA to fix row order",
    )
    args = ap.parse_args()

    s1 = load_clones(args.stage1_clones, "stage1_status")
    s2 = load_clones(args.stage2_clones, "stage2_status")

    by_target: Dict[str, List[str]] = defaultdict(list)
    seen: Dict[str, Set[str]] = defaultdict(set)
    for src in (s1, s2):
        for key, row in src.items():
            tid = row["_target_id"]
            cid = row["_clone_id"]
            if cid not in seen[tid]:
                seen[tid].add(cid)
                by_target[tid].append(cid)

    if args.targets and args.targets.is_file():
        # reuse common loader lightly
        from common import load_targets

        target_order = list(load_targets(args.targets).keys())
    else:
        target_order = sorted(by_target)

    max_n = max((len(v) for v in by_target.values()), default=0)
    sample_cols = [f"sample_{i}" for i in range(1, max_n + 1)]

    matrix_rows: List[dict] = []
    long_rows: List[dict] = []
    return_rows: List[dict] = []

    for tid in target_order:
        samples = sorted(by_target.get(tid, []), key=lambda x: (len(x), x))
        row: dict = {
            "target_id": tid,
            "n_samples": len(samples),
            "n_correct": 0,
            "target_has_correct": "NO",
            "correct_sample_ids": "",
        }
        correct_ids: List[str] = []
        for i, col in enumerate(sample_cols):
            if i >= len(samples):
                row[col] = ""
                row[f"{col}_id"] = ""
                continue
            sid = samples[i]
            key = f"{tid}::{sid}"
            st = sample_status(s1.get(key), s2.get(key))
            row[col] = st
            row[f"{col}_id"] = sid
            files = (s2.get(key) or s1.get(key) or {}).get("files", "")
            long_rows.append(
                {
                    "target_id": tid,
                    "sample_id": sid,
                    "replicate": i + 1,
                    "status": st,
                    "stage1_status": s1[key]["_status"] if key in s1 else "",
                    "stage2_status": s2[key]["_status"] if key in s2 else "",
                    "files": files,
                }
            )
            if st == "CORRECT":
                correct_ids.append(sid)
                return_rows.append(
                    {
                        "target_id": tid,
                        "sample_id": sid,
                        "replicate": i + 1,
                        "status": "CORRECT",
                        "action": "collect return plasmid",
                        "files": (s2.get(key) or {}).get("files", ""),
                    }
                )
        row["n_correct"] = len(correct_ids)
        row["target_has_correct"] = "YES" if correct_ids else "NO"
        row["correct_sample_ids"] = "|".join(correct_ids)
        matrix_rows.append(row)

    out = args.out_dir
    wide_fields = [
        "target_id",
        "n_samples",
        "n_correct",
        "target_has_correct",
        "correct_sample_ids",
        *[x for col in sample_cols for x in (col, f"{col}_id")],
    ]
    write_csv(out / "final_target_sample_matrix.csv", matrix_rows, wide_fields)
    # alias for older name
    write_csv(
        out / "final_gene_sample_matrix.csv",
        [
            {
                **{("vendor_id" if k == "target_id" else "gene_has_correct" if k == "target_has_correct" else k): v for k, v in r.items()}
            }
            for r in matrix_rows
        ],
        ["vendor_id" if f == "target_id" else "gene_has_correct" if f == "target_has_correct" else f for f in wide_fields],
    )
    write_csv(
        out / "final_target_sample_long.csv",
        long_rows,
        [
            "target_id",
            "sample_id",
            "replicate",
            "status",
            "stage1_status",
            "stage2_status",
            "files",
        ],
    )
    write_csv(
        out / "final_plasmid_return_list.csv",
        return_rows,
        ["target_id", "sample_id", "replicate", "status", "action", "files"],
    )

    color = {
        "CORRECT": "#c6f6d5",
        "SEQ_OK_AB1_FAIL": "#fefcbf",
        "FAIL": "#fed7d7",
        "": "#f7fafc",
    }
    parts = [
        "<!DOCTYPE html><html><head><meta charset='utf-8'>",
        "<title>Final target × sample matrix</title>",
        "<style>",
        "body{font-family:Segoe UI,Arial,sans-serif;margin:24px;color:#1a202c}",
        "table{border-collapse:collapse;font-size:13px}",
        "th,td{border:1px solid #cbd5e0;padding:6px 8px;text-align:center}",
        "th{background:#edf2f7} td.gene{text-align:left;font-weight:600}",
        ".legend span{display:inline-block;padding:2px 8px;margin-right:8px;border:1px solid #cbd5e0}",
        "</style></head><body>",
        "<h1>Final matrix: target × sample</h1>",
        "<p>One sample/clone merges all assigned primers (2–N reads). "
        "<b>CORRECT</b> = ready to collect return plasmid.</p>",
        "<div class='legend'>",
        f"<span style='background:{color['CORRECT']}'>CORRECT</span>",
        f"<span style='background:{color['SEQ_OK_AB1_FAIL']}'>SEQ_OK_AB1_FAIL</span>",
        f"<span style='background:{color['FAIL']}'>FAIL</span>",
        "</div><table><tr><th>target_id</th><th>n_correct</th>",
    ]
    for i in range(1, max_n + 1):
        parts.append(f"<th>sample_{i}</th>")
    parts.append("</tr>")
    for row in matrix_rows:
        parts.append("<tr>")
        parts.append(f"<td class='gene'>{html.escape(row['target_id'])}</td>")
        parts.append(f"<td>{row['n_correct']}/{row['n_samples']}</td>")
        for i in range(1, max_n + 1):
            st = row.get(f"sample_{i}", "")
            sid = row.get(f"sample_{i}_id", "")
            label = f"{html.escape(sid)}<br><b>{html.escape(st)}</b>" if sid else ""
            parts.append(f"<td style='background:{color.get(st, color[''])}'>{label}</td>")
        parts.append("</tr>")
    parts.append("</table>")
    parts.append("<h2>Return plasmid list (CORRECT)</h2><table>")
    parts.append("<tr><th>target_id</th><th>sample_id</th><th>replicate</th><th>action</th></tr>")
    for r in return_rows:
        parts.append(
            "<tr>"
            f"<td class='gene'>{html.escape(r['target_id'])}</td>"
            f"<td>{html.escape(r['sample_id'])}</td>"
            f"<td>{r['replicate']}</td>"
            f"<td>{html.escape(r['action'])}</td>"
            "</tr>"
        )
    if not return_rows:
        parts.append("<tr><td colspan='4'>No CORRECT samples</td></tr>")
    parts.append("</table></body></html>")
    html_path = out / "final_target_sample_matrix.html"
    html_path.write_text("".join(parts), encoding="utf-8")
    (out / "final_gene_sample_matrix.html").write_text("".join(parts), encoding="utf-8")

    n_correct = len(return_rows)
    n_ok = sum(1 for r in matrix_rows if r["target_has_correct"] == "YES")
    print(f"targets={len(matrix_rows)} targets_with_correct={n_ok}")
    print(f"correct_samples={n_correct}")
    print(f"wrote {out / 'final_target_sample_matrix.csv'}")
    print(f"wrote {out / 'final_plasmid_return_list.csv'}")
    print(f"wrote {html_path}")


if __name__ == "__main__":
    main()
