#!/usr/bin/env python3
"""Optional helper: draft targets.csv + assignments.csv from a HER2-style batch.

This is ONLY a convenience for one naming convention. The main pipeline does not
require HER2 names — provide your own targets/assignments for general use.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--insert-csv", type=Path, required=True, help="pcr_insert_flanks.csv")
    ap.add_argument("--reads-dir", type=Path, required=True)
    ap.add_argument("--aliases", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    aliases = {}
    if args.aliases and args.aliases.is_file():
        data = json.loads(args.aliases.read_text(encoding="utf-8"))
        aliases = data.get("aliases", data)

    targets = []
    by_task: dict = defaultdict(list)
    with args.insert_csv.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            tid = row["vendor_id"].strip()
            seq = row["full_insert_dna"].strip().upper()
            targets.append({"target_id": tid, "sequence": seq})
            m = re.match(r"^h(\d+)-(\d+)$", tid)
            if m:
                by_task[m.group(1)].append(int(m.group(2)))

    prefix_map = {}
    for task, seqs in by_task.items():
        seqs = sorted(set(seqs))
        if len(seqs) == 1:
            prefix_map[task] = f"h{task}-{seqs[0]}"
        else:
            for s in seqs:
                prefix_map[task if s == 1 else f"{task}{s}"] = f"h{task}-{s}"

    file_re = re.compile(
        r"^S\d+_(?P<sample>[^_]+)_(?P<primer>T7T?|t7t?)_(?P<well>[^.]+)\.(?P<ext>seq|ab1)$",
        re.I,
    )
    asg = []
    for p in sorted(args.reads_dir.glob("*.seq")):
        m = file_re.match(p.name)
        if not m:
            continue
        sample = m.group("sample")
        if "-" not in sample:
            continue
        prefix, _clone_n = sample.rsplit("-", 1)
        tid = aliases.get(prefix) or prefix_map.get(prefix)
        if not tid:
            continue
        asg.append(
            {
                "file": p.name,
                "target_id": tid,
                "clone_id": sample,
                "primer": m.group("primer").upper(),
            }
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "targets.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["target_id", "sequence"])
        w.writeheader()
        w.writerows(targets)
    with (args.out_dir / "assignments.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["file", "target_id", "clone_id", "primer"])
        w.writeheader()
        w.writerows(asg)
    print(f"wrote {args.out_dir / 'targets.csv'} ({len(targets)} targets)")
    print(f"wrote {args.out_dir / 'assignments.csv'} ({len(asg)} rows)")


if __name__ == "__main__":
    main()
