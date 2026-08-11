#!/usr/bin/env python3
"""Shared loaders for general Sanger QC inputs."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass
class Target:
    target_id: str
    sequence: str


@dataclass
class Assignment:
    file: str
    target_id: str
    clone_id: str
    primer: str


def _norm_dna(s: str, *, target_id: str = "") -> str:
    """Strip whitespace only; refuse IUPAC/N/other non-ACGT (never silently delete)."""
    cleaned = re.sub(r"\s+", "", s).upper()
    for i, ch in enumerate(cleaned):
        if ch not in "ACGT":
            where = f"target {target_id}" if target_id else "sequence"
            raise SystemExit(
                f"Invalid target sequence: {where} contains unsupported base "
                f"'{ch}' at position {i + 1}. Only A/C/G/T are allowed "
                f"(whitespace is stripped; N/IUPAC are not silently removed)."
            )
    if not cleaned:
        where = f"target {target_id}" if target_id else "sequence"
        raise SystemExit(f"Invalid target sequence: {where} is empty after stripping whitespace.")
    return cleaned


def load_targets(path: Path) -> Dict[str, Target]:
    """Load target sequences from CSV or FASTA.

    CSV columns (any alias works):
      target_id|id|name|vendor_id
      sequence|seq|dna|full_insert_dna

    FASTA: header id after '>' is target_id.
    """
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"targets file not found: {path}")

    out: Dict[str, Target] = {}
    if path.suffix.lower() in {".fa", ".fasta", ".fna"}:
        tid = None
        chunks: List[str] = []
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if tid and chunks:
                    out[tid] = Target(tid, _norm_dna("".join(chunks), target_id=tid))
                tid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if tid and chunks:
            out[tid] = Target(tid, _norm_dna("".join(chunks), target_id=tid))
    else:
        with path.open(encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f):
                lower = {k.lower().strip(): (v or "").strip() for k, v in row.items() if k}
                tid = (
                    lower.get("target_id")
                    or lower.get("id")
                    or lower.get("name")
                    or lower.get("vendor_id")
                    or ""
                )
                seq = (
                    lower.get("sequence")
                    or lower.get("seq")
                    or lower.get("dna")
                    or lower.get("full_insert_dna")
                    or ""
                )
                if not tid or not seq:
                    raise SystemExit(f"bad targets row in {path}: need target_id + sequence")
                out[tid] = Target(tid, _norm_dna(seq, target_id=tid))

    if not out:
        raise SystemExit(f"no targets loaded from {path}")
    return out


def load_assignments(path: Path) -> List[Assignment]:
    """Load which read files belong to which target/clone.

    Required columns:
      file, target_id, clone_id
    Optional:
      primer  (defaults to 'read')

    `file` is the basename (.seq or .ab1). The same stem may be used for both stages.
    """
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"assignments file not found: {path}")

    rows: List[Assignment] = []
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            lower = {k.lower().strip(): (v or "").strip() for k, v in row.items() if k}
            file_name = lower.get("file") or lower.get("filename") or lower.get("read") or ""
            tid = lower.get("target_id") or lower.get("target") or lower.get("vendor_id") or ""
            clone = lower.get("clone_id") or lower.get("clone") or lower.get("sample_id") or lower.get("sample") or ""
            primer = lower.get("primer") or lower.get("reaction") or "read"
            if not file_name or not tid or not clone:
                raise SystemExit(
                    f"bad assignments row in {path}: need file, target_id, clone_id"
                )
            rows.append(
                Assignment(
                    file=Path(file_name).name,
                    target_id=tid,
                    clone_id=clone,
                    primer=primer.upper(),
                )
            )
    if not rows:
        raise SystemExit(f"no assignments loaded from {path}")
    return rows


def index_assignments(rows: Iterable[Assignment]) -> Dict[str, Assignment]:
    """Map basename -> assignment (last wins if duplicates)."""
    return {a.file: a for a in rows}


def resolve_read_path(reads_dir: Path, name: str, prefer_ext: str) -> Optional[Path]:
    """Find a read file by exact name or by switching .seq/.ab1 stem."""
    reads_dir = Path(reads_dir)
    direct = reads_dir / name
    if direct.is_file():
        return direct
    stem = Path(name).stem
    cand = reads_dir / f"{stem}{prefer_ext}"
    if cand.is_file():
        return cand
    # case-insensitive scan of stem
    want = stem.lower()
    for p in reads_dir.iterdir():
        if p.is_file() and p.stem.lower() == want and p.suffix.lower() == prefer_ext.lower():
            return p
    return None


def validate_assignments(targets: Dict[str, Target], assignments: List[Assignment]) -> None:
    missing = sorted({a.target_id for a in assignments if a.target_id not in targets})
    if missing:
        raise SystemExit(
            "assignments reference unknown target_id(s): "
            + ", ".join(missing)
            + f"\nknown targets: {', '.join(sorted(targets))}"
        )


def group_key(target_id: str, clone_id: str) -> str:
    return f"{target_id}::{clone_id}"


def parse_group_key(key: str) -> Tuple[str, str]:
    tid, clone = key.split("::", 1)
    return tid, clone
