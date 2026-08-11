#!/usr/bin/env python3
"""Export final gene × sample matrix for plasmid return.

A sample = one clone (T7 + T7T counted as one). Cell marks whether that
sample is correct after stage-2 (default). Also writes a compact HTML table
and a plasmid-return list.
"""

from __future__ import annotations

import argparse
import csv
import html
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set


def load_clones(path: Path, status_key: str) -> Dict[str, dict]:
    """clone_id -> row."""
    out: Dict[str, dict] = {}
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            out[row["clone_id"]] = row
            row["_status"] = row[status_key]
    return out


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def sample_status(s1: dict | None, s2: dict | None) -> str:
    """Final label for plasmid collection."""
    if s2 is not None and s2["_status"] == "PASS":
        return "CORRECT"
    if s1 is not None and s1["_status"] == "PASS":
        return "SEQ_OK_AB1_FAIL"
    if s1 is not None or s2 is not None:
        return "FAIL"
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description="Final gene×sample matrix for plasmid return")
    ap.add_argument("--stage1-clones", type=Path, required=True)
    ap.add_argument("--stage2-clones", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument(
        "--genes-csv",
        type=Path,
        default=None,
        help="optional pcr_insert_flanks.csv to fix gene row order",
    )
    args = ap.parse_args()

    s1 = load_clones(args.stage1_clones, "stage1_status")
    s2 = load_clones(args.stage2_clones, "stage2_status")

    # gene -> ordered sample ids (prefer stage1 universe = all sequenced clones)
    by_gene: Dict[str, List[str]] = defaultdict(list)
    seen: Dict[str, Set[str]] = defaultdict(set)
    for src in (s1, s2):
        for cid, row in src.items():
            g = row["vendor_id"]
            if cid not in seen[g]:
                seen[g].add(cid)
                by_gene[g].append(cid)

    if args.genes_csv and args.genes_csv.is_file():
        with args.genes_csv.open(encoding="utf-8-sig", newline="") as f:
            gene_order = [r["vendor_id"].strip() for r in csv.DictReader(f)]
    else:
        gene_order = sorted(by_gene)

    # normalize sample columns: use replicate index 1..n per gene for a clean grid,
    # but keep real sample_id in long form / return list.
    max_n = max((len(v) for v in by_gene.values()), default=0)
    sample_cols = [f"sample_{i}" for i in range(1, max_n + 1)]

    matrix_rows: List[dict] = []
    long_rows: List[dict] = []
    return_rows: List[dict] = []

    for gene in gene_order:
        samples = sorted(by_gene.get(gene, []), key=lambda x: (len(x), x))
        row: dict = {
            "vendor_id": gene,
            "n_samples": len(samples),
            "n_correct": 0,
            "gene_has_correct": "NO",
            "correct_sample_ids": "",
        }
        correct_ids: List[str] = []
        for i, col in enumerate(sample_cols):
            if i >= len(samples):
                row[col] = ""
                row[f"{col}_id"] = ""
                continue
            sid = samples[i]
            st = sample_status(s1.get(sid), s2.get(sid))
            row[col] = st
            row[f"{col}_id"] = sid
            long_rows.append(
                {
                    "vendor_id": gene,
                    "sample_id": sid,
                    "replicate": i + 1,
                    "status": st,
                    "stage1_status": s1[sid]["_status"] if sid in s1 else "",
                    "stage2_status": s2[sid]["_status"] if sid in s2 else "",
                    "files": (s2.get(sid) or s1.get(sid) or {}).get("files", ""),
                }
            )
            if st == "CORRECT":
                correct_ids.append(sid)
                return_rows.append(
                    {
                        "vendor_id": gene,
                        "sample_id": sid,
                        "replicate": i + 1,
                        "status": "CORRECT",
                        "action": "收集返样质粒",
                        "files": (s2.get(sid) or {}).get("files", ""),
                    }
                )
        row["n_correct"] = len(correct_ids)
        row["gene_has_correct"] = "YES" if correct_ids else "NO"
        row["correct_sample_ids"] = "|".join(correct_ids)
        matrix_rows.append(row)

    out = args.out_dir
    wide_fields = [
        "vendor_id",
        "n_samples",
        "n_correct",
        "gene_has_correct",
        "correct_sample_ids",
        *[x for col in sample_cols for x in (col, f"{col}_id")],
    ]
    write_csv(out / "final_gene_sample_matrix.csv", matrix_rows, wide_fields)
    write_csv(
        out / "final_gene_sample_long.csv",
        long_rows,
        [
            "vendor_id",
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
        ["vendor_id", "sample_id", "replicate", "status", "action", "files"],
    )

    # HTML
    color = {
        "CORRECT": "#c6f6d5",
        "SEQ_OK_AB1_FAIL": "#fefcbf",
        "FAIL": "#fed7d7",
        "": "#f7fafc",
    }
    parts = [
        "<!DOCTYPE html><html><head><meta charset='utf-8'>",
        "<title>Final gene × sample matrix</title>",
        "<style>",
        "body{font-family:Segoe UI,Arial,sans-serif;margin:24px;color:#1a202c}",
        "h1{font-size:20px} h2{font-size:16px;margin-top:28px}",
        "table{border-collapse:collapse;font-size:13px}",
        "th,td{border:1px solid #cbd5e0;padding:6px 8px;text-align:center}",
        "th{background:#edf2f7} td.gene{text-align:left;font-weight:600}",
        ".legend span{display:inline-block;padding:2px 8px;margin-right:8px;border:1px solid #cbd5e0}",
        "</style></head><body>",
        "<h1>最终矩阵：基因 × 样品（正反向合并）</h1>",
        "<p>目的：标注每个想合成的基因，哪些样品测序正确，便于收集返样质粒。"
        "同一样品的 T7/T7T 计为一次测序。</p>",
        "<div class='legend'>",
        f"<span style='background:{color['CORRECT']}'>CORRECT 可返样</span>",
        f"<span style='background:{color['SEQ_OK_AB1_FAIL']}'>SEQ_OK_AB1_FAIL</span>",
        f"<span style='background:{color['FAIL']}'>FAIL</span>",
        "</div>",
        "<h2>基因 × 样品</h2><table><tr><th>vendor_id</th><th>n_correct</th>",
    ]
    for i in range(1, max_n + 1):
        parts.append(f"<th>sample_{i}</th>")
    parts.append("</tr>")
    for row in matrix_rows:
        parts.append("<tr>")
        parts.append(f"<td class='gene'>{html.escape(row['vendor_id'])}</td>")
        parts.append(f"<td>{row['n_correct']}/{row['n_samples']}</td>")
        for i in range(1, max_n + 1):
            st = row.get(f"sample_{i}", "")
            sid = row.get(f"sample_{i}_id", "")
            label = f"{html.escape(sid)}<br><b>{html.escape(st)}</b>" if sid else ""
            parts.append(
                f"<td style='background:{color.get(st, color[''])}'>{label}</td>"
            )
        parts.append("</tr>")
    parts.append("</table>")

    parts.append("<h2>返样质粒清单（CORRECT）</h2><table>")
    parts.append("<tr><th>vendor_id</th><th>sample_id</th><th>replicate</th><th>action</th></tr>")
    for r in return_rows:
        parts.append(
            "<tr>"
            f"<td class='gene'>{html.escape(r['vendor_id'])}</td>"
            f"<td>{html.escape(r['sample_id'])}</td>"
            f"<td>{r['replicate']}</td>"
            f"<td>{html.escape(r['action'])}</td>"
            "</tr>"
        )
    if not return_rows:
        parts.append("<tr><td colspan='4'>无 CORRECT 样品</td></tr>")
    parts.append("</table></body></html>")
    html_path = out / "final_gene_sample_matrix.html"
    html_path.write_text("".join(parts), encoding="utf-8")

    n_correct = len(return_rows)
    n_genes_ok = sum(1 for r in matrix_rows if r["gene_has_correct"] == "YES")
    print(f"genes={len(matrix_rows)} genes_with_correct={n_genes_ok}")
    print(f"correct_samples={n_correct}")
    print(f"wrote {out / 'final_gene_sample_matrix.csv'}")
    print(f"wrote {out / 'final_gene_sample_long.csv'}")
    print(f"wrote {out / 'final_plasmid_return_list.csv'}")
    print(f"wrote {html_path}")


if __name__ == "__main__":
    main()
