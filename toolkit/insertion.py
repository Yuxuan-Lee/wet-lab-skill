"""Reference-boundary insertion detection (separate from base-position union).

Base correctness remains a per-reference-position union elsewhere.
Insertions are events *between* reference bases and are gated by local
alignment quality around each boundary — never by whole-target coverage.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# --- configurable thresholds (keep centralized) ---
LOCAL_IDENTITY_MIN = 0.95
FLANK_MATCH_MIN_BP = 20
EDGE_MARGIN_BP = 25
LOCAL_WINDOW_BP = 40
AB1_PEAK_MIN_DEFAULT = 0.75
AB1_INSERT_MIN_FRAC = 0.50  # inserted bases: at least some peak signal
AB1_FLANK_MIN_FRAC = 0.75  # flanks should look like clean Sanger peaks

# Stage-1 merge statuses
STATUS_NO_INSERT = "NO_INSERT_EVIDENCE"
STATUS_CANDIDATE = "INSERT_CANDIDATE"
STATUS_STRONG = "STRONG_INSERT"
STATUS_CONFLICT = "INSERT_CONFLICT"

# Stage-2 validation statuses
STATUS_CONFIRMED = "CONFIRMED_INSERT"
STATUS_POSSIBLE = "POSSIBLE_INSERT"


@dataclass
class InsertionCandidate:
    after_ref_pos: int  # insertion after this 1-based? use 0-based: after ref index
    # Convention: after_ref_pos = number of ref bases consumed before the gap
    # (insertion lies between ref[after_ref_pos-1] and ref[after_ref_pos] when after>0)
    inserted_seq: str
    interval_start: int  # inclusive 0-based ref index of ambiguous boundary start
    interval_end: int  # inclusive 0-based; insertion after any pos in [start, end+1) equiv
    left_flank_match_bp: int = 0
    right_flank_match_bp: int = 0
    local_identity: float = 0.0
    distance_to_read_edge: int = 0
    query_start: int = -1  # 0-based index in query of first inserted base
    query_end: int = -1  # exclusive
    high_confidence: bool = False
    primer: str = ""

    @property
    def insert_len(self) -> int:
        return len(self.inserted_seq)

    def summary(self) -> str:
        return (
            f"after={self.after_ref_pos}:"
            f"{self.interval_start}-{self.interval_end}:"
            f"+{self.inserted_seq}(L{self.left_flank_match_bp}/R{self.right_flank_match_bp}/"
            f"id={self.local_identity:.3f}/edge={self.distance_to_read_edge})"
        )


@dataclass
class NoInsertSpan:
    """High-confidence continuous match spanning a ref interval (supports NO_INSERT)."""

    primer: str
    start: int  # inclusive 0-based ref
    end: int  # exclusive 0-based ref
    identity: float = 1.0


@dataclass
class CloneInsertResult:
    status: str
    internal_insert_bp: int
    evidence: str
    events: List[dict] = field(default_factory=list)


def extract_gap_runs(
    ref_aln: str,
    query_aln: str,
) -> List[dict]:
    """Extract REF-gap + QUERY-base runs as raw insertion candidates.

    after_ref_pos: count of non-gap ref bases seen before the gap run
      (0 = before first ref base; ref_len = after last).
    """
    runs: List[dict] = []
    ref_pos = 0
    query_pos = 0
    i = 0
    n = len(ref_aln)
    while i < n:
        rc, qc = ref_aln[i], query_aln[i]
        if rc == "-" and qc != "-":
            start_ref = ref_pos
            start_q = query_pos
            inserted = []
            while i < n and ref_aln[i] == "-" and query_aln[i] != "-":
                inserted.append(query_aln[i])
                query_pos += 1
                i += 1
            runs.append(
                {
                    "after_ref_pos": start_ref,
                    "inserted_seq": "".join(inserted),
                    "query_start": start_q,
                    "query_end": query_pos,
                }
            )
            continue
        if rc != "-":
            ref_pos += 1
        if qc != "-":
            query_pos += 1
        i += 1
    return runs


def _flank_stats(
    ref_aln: str,
    query_aln: str,
    after_ref_pos: int,
    window: int = LOCAL_WINDOW_BP,
) -> Tuple[int, int, float]:
    """Left/right continuous match counts and local identity around boundary."""
    # Map alignment columns to ref positions
    cols_left: List[Tuple[str, str]] = []
    cols_right: List[Tuple[str, str]] = []
    ref_pos = 0
    for rc, qc in zip(ref_aln, query_aln):
        if rc != "-":
            if ref_pos < after_ref_pos:
                cols_left.append((rc, qc))
            else:
                cols_right.append((rc, qc))
            ref_pos += 1

    left = cols_left[-window:] if cols_left else []
    right = cols_right[:window] if cols_right else []

    def match_run(cols: Sequence[Tuple[str, str]], reverse: bool) -> int:
        seq = list(reversed(cols)) if reverse else list(cols)
        n = 0
        for rc, qc in seq:
            if qc == "-" or rc != qc:
                break
            n += 1
        return n

    left_match = match_run(left, reverse=True)
    right_match = match_run(right, reverse=False)

    both = left + right
    aligned = [(rc, qc) for rc, qc in both if rc != "-" and qc != "-"]
    if not aligned:
        local_id = 0.0
    else:
        local_id = sum(1 for rc, qc in aligned if rc == qc) / len(aligned)
    return left_match, right_match, local_id


def homopolymer_interval(ref: str, after_ref_pos: int, inserted_seq: str) -> Tuple[int, int]:
    """Widen boundary when insert is a homopolymer run matching adjacent bases.

    Returns inclusive 0-based interval of ref positions that are equivalent
    placement sites for the insertion (documented as positions around the gap).
    """
    if not ref:
        return (0, 0)
    # Default: single boundary after after_ref_pos-1 .. after_ref_pos
    # interval of "positions" for messaging: after_ref_pos itself as anchor
    lo = max(0, after_ref_pos - 1)
    hi = min(len(ref) - 1, after_ref_pos)

    if not inserted_seq:
        return (lo, hi)

    base = inserted_seq[0]
    if len(set(inserted_seq)) != 1:
        return (max(0, after_ref_pos - 1), min(len(ref) - 1, max(after_ref_pos, after_ref_pos)))

    # Expand through identical homopolymer on either side of the gap
    left = after_ref_pos - 1
    while left >= 0 and ref[left] == base:
        left -= 1
    left += 1
    right = after_ref_pos
    while right < len(ref) and ref[right] == base:
        right += 1
    right -= 1
    if right < left:
        return (lo, hi)
    return (left, right)


def score_insertion(
    ref_seq: str,
    ref_aln: str,
    query_aln: str,
    query_len: int,
    raw: dict,
    primer: str = "",
) -> InsertionCandidate:
    after = int(raw["after_ref_pos"])
    inserted = str(raw["inserted_seq"]).upper()
    q0 = int(raw.get("query_start", -1))
    q1 = int(raw.get("query_end", -1))
    left_m, right_m, local_id = _flank_stats(ref_aln, query_aln, after)
    # distance to nearer read edge using query coordinates of the insert
    if q0 >= 0 and q1 >= 0 and query_len > 0:
        edge = min(q0, max(0, query_len - q1))
    else:
        edge = 0
    iv0, iv1 = homopolymer_interval(ref_seq, after, inserted)
    cand = InsertionCandidate(
        after_ref_pos=after,
        inserted_seq=inserted,
        interval_start=iv0,
        interval_end=iv1,
        left_flank_match_bp=left_m,
        right_flank_match_bp=right_m,
        local_identity=round(local_id, 4),
        distance_to_read_edge=edge,
        query_start=q0,
        query_end=q1,
        primer=primer,
    )
    cand.high_confidence = is_high_confidence(cand, ref_len=len(ref_seq))
    return cand


def is_high_confidence(cand: InsertionCandidate, ref_len: int = 0) -> bool:
    """Local flank gate — intentionally independent of whole-target match_bp."""
    if cand.insert_len <= 0:
        return False
    # Must be internal to the reference (not leading/trailing vector gaps)
    if cand.after_ref_pos <= 0:
        return False
    if ref_len and cand.after_ref_pos >= ref_len:
        return False
    if cand.left_flank_match_bp < FLANK_MATCH_MIN_BP:
        return False
    if cand.right_flank_match_bp < FLANK_MATCH_MIN_BP:
        return False
    if cand.local_identity < LOCAL_IDENTITY_MIN:
        return False
    if cand.distance_to_read_edge < EDGE_MARGIN_BP:
        return False
    return True


def candidates_from_alignment(
    ref_seq: str,
    ref_aln: str,
    query_aln: str,
    query_len: int,
    primer: str = "",
) -> List[InsertionCandidate]:
    out: List[InsertionCandidate] = []
    for raw in extract_gap_runs(ref_aln, query_aln):
        out.append(score_insertion(ref_seq, ref_aln, query_aln, query_len, raw, primer=primer))
    return out


def no_insert_spans_from_alignment(
    ref_aln: str,
    query_aln: str,
    primer: str = "",
    min_run: int = FLANK_MATCH_MIN_BP * 2,
) -> List[NoInsertSpan]:
    """Continuous perfect match runs on the reference (NO_INSERT support)."""
    spans: List[NoInsertSpan] = []
    ref_pos = 0
    run_start: Optional[int] = None
    for rc, qc in zip(ref_aln, query_aln):
        if rc == "-":
            continue
        matched = qc != "-" and rc == qc
        if matched:
            if run_start is None:
                run_start = ref_pos
        else:
            if run_start is not None and ref_pos - run_start >= min_run:
                spans.append(NoInsertSpan(primer=primer, start=run_start, end=ref_pos))
            run_start = None
        ref_pos += 1
    if run_start is not None and ref_pos - run_start >= min_run:
        spans.append(NoInsertSpan(primer=primer, start=run_start, end=ref_pos))
    return spans


def spans_boundary(span: NoInsertSpan, cand: InsertionCandidate) -> bool:
    """True if a NO_INSERT span covers both sides of the insertion boundary."""
    # Boundary sits after after_ref_pos (between after_ref_pos-1 and after_ref_pos)
    # Require span to include at least FLANK_MATCH_MIN_BP on each side conceptually:
    # span.start < after_ref_pos < span.end  with enough margin
    after = cand.after_ref_pos
    if not (span.start < after < span.end):
        # also allow interval overlap for homopolymer
        if not (span.start <= cand.interval_end and span.end > cand.interval_start):
            return False
        if not (span.start < after or after < span.end):
            return False
    left = after - span.start
    right = span.end - after
    return left >= FLANK_MATCH_MIN_BP and right >= FLANK_MATCH_MIN_BP


def _same_event(a: InsertionCandidate, b: InsertionCandidate) -> bool:
    if a.inserted_seq != b.inserted_seq:
        # allow same-length homopolymer of same base
        if not (
            len(a.inserted_seq) == len(b.inserted_seq)
            and len(set(a.inserted_seq)) == 1
            and a.inserted_seq == b.inserted_seq
        ):
            return False
    # overlapping / adjacent intervals
    return not (a.interval_end < b.interval_start - 1 or b.interval_end < a.interval_start - 1)


def merge_clone_insertions(
    primer_candidates: Dict[str, Sequence[InsertionCandidate]],
    primer_no_insert: Optional[Dict[str, Sequence[NoInsertSpan]]] = None,
) -> CloneInsertResult:
    """Merge per-primer insertion evidence by reference boundary."""
    primer_no_insert = primer_no_insert or {}
    hc: List[InsertionCandidate] = []
    for primer, cands in primer_candidates.items():
        for c in cands:
            if c.high_confidence:
                if not c.primer:
                    c.primer = primer
                hc.append(c)

    if not hc:
        return CloneInsertResult(
            status=STATUS_NO_INSERT,
            internal_insert_bp=0,
            evidence="",
            events=[],
        )

    # Cluster HC inserts into events
    clusters: List[List[InsertionCandidate]] = []
    for c in hc:
        placed = False
        for cl in clusters:
            if _same_event(c, cl[0]):
                cl.append(c)
                placed = True
                break
        if not placed:
            clusters.append([c])

    event_statuses: List[str] = []
    events: List[dict] = []
    total_bp = 0
    evidence_bits: List[str] = []

    for cl in clusters:
        primers = sorted({c.primer for c in cl})
        rep = cl[0]
        total_bp += rep.insert_len
        # Conflict if another primer has HC NO_INSERT spanning boundary and no HC insert there
        conflict = False
        for primer, spans in primer_no_insert.items():
            if primer in primers:
                continue
            # If this primer also has an HC insert for same event, not a NO_INSERT witness
            other_hc = primer_candidates.get(primer, [])
            if any(x.high_confidence and _same_event(x, rep) for x in other_hc):
                continue
            for sp in spans:
                if spans_boundary(sp, rep):
                    conflict = True
                    break
            if conflict:
                break

        if conflict:
            st = STATUS_CONFLICT
        elif len(primers) >= 2:
            st = STATUS_STRONG
        else:
            st = STATUS_CANDIDATE

        event_statuses.append(st)
        ev = {
            "status": st,
            "after_ref_pos": rep.after_ref_pos,
            "interval_start": rep.interval_start,
            "interval_end": rep.interval_end,
            "inserted_seq": rep.inserted_seq,
            "insert_len": rep.insert_len,
            "primers": "|".join(primers),
            "n_hc_reads": len(cl),
        }
        events.append(ev)
        evidence_bits.append(
            f"{st}:ins@{rep.interval_start}-{rep.interval_end}+{rep.inserted_seq}[{','.join(primers)}]"
        )

    # Aggregate: conflict > strong > candidate
    if STATUS_CONFLICT in event_statuses:
        status = STATUS_CONFLICT
    elif STATUS_STRONG in event_statuses:
        status = STATUS_STRONG
    else:
        status = STATUS_CANDIDATE

    return CloneInsertResult(
        status=status,
        internal_insert_bp=total_bp,
        evidence=";".join(evidence_bits),
        events=events,
    )


def validate_insert_with_ab1(
    primer_candidates: Dict[str, Sequence[InsertionCandidate]],
    primer_no_insert: Dict[str, Sequence[NoInsertSpan]],
    primer_insert_peak_ok: Dict[str, bool],
    primer_flank_peak_ok: Dict[str, bool],
) -> CloneInsertResult:
    """Promote/demote stage-1 insert evidence using AB1 peak checks.

    If AB1 peak validation is inconclusive, keep POSSIBLE_INSERT (never silent drop).
    """
    base = merge_clone_insertions(primer_candidates, primer_no_insert)
    if base.status == STATUS_NO_INSERT:
        return base

    if base.status == STATUS_CONFLICT:
        return CloneInsertResult(
            status=STATUS_CONFLICT,
            internal_insert_bp=base.internal_insert_bp,
            evidence=base.evidence,
            events=[{**e, "status": STATUS_CONFLICT} for e in base.events],
        )

    any_peak_ok = False
    any_peak_fail = False
    any_unknown = False
    for ev in base.events:
        for p in str(ev["primers"]).split("|"):
            if not p:
                continue
            if p not in primer_insert_peak_ok or p not in primer_flank_peak_ok:
                any_unknown = True
                continue
            if primer_insert_peak_ok[p] and primer_flank_peak_ok[p]:
                any_peak_ok = True
            else:
                any_peak_fail = True

    # Confirmed only when peaks clearly support; otherwise POSSIBLE (never silent clear)
    if any_peak_ok and not any_peak_fail and not any_unknown:
        st = STATUS_CONFIRMED
    else:
        st = STATUS_POSSIBLE

    return CloneInsertResult(
        status=st,
        internal_insert_bp=base.internal_insert_bp,
        evidence=base.evidence.replace(STATUS_CANDIDATE, st).replace(STATUS_STRONG, st),
        events=[{**e, "status": st} for e in base.events],
    )


def ab1_peak_ok_for_insert(
    heights: Sequence[Dict[str, int]],
    query_indices: Sequence[int],
    called_bases: Sequence[str],
    orientation: str,
    peak_min: float = AB1_INSERT_MIN_FRAC,
) -> bool:
    """True if inserted-query bases show non-trivial peak for the called base."""
    complement = str.maketrans("GATC", "CTAG")
    if not query_indices:
        return False
    ok = 0
    for qi, base in zip(query_indices, called_bases):
        if qi < 0 or qi >= len(heights):
            return False
        b = base.upper()
        if b not in "GATC":
            return False
        want = b if orientation == "as_is" else b.translate(complement)
        h = heights[qi]
        total = sum(h.get(x, 0) for x in "GATC")
        if total <= 0:
            return False
        if h.get(want, 0) / total >= peak_min:
            ok += 1
    return ok == len(query_indices)


def ab1_peak_ok_for_flanks(
    heights: Sequence[Dict[str, int]],
    ref_to_query: Dict[int, int],
    ref_seq: str,
    after_ref_pos: int,
    orientation: str,
    flank: int = FLANK_MATCH_MIN_BP,
    peak_min: float = AB1_FLANK_MIN_FRAC,
) -> bool:
    COMPLEMENT = str.maketrans("GATC", "CTAG")
    left_positions = list(range(max(0, after_ref_pos - flank), after_ref_pos))
    right_positions = list(range(after_ref_pos, min(len(ref_seq), after_ref_pos + flank)))
    checked = 0
    ok = 0
    for pos in left_positions + right_positions:
        if pos not in ref_to_query:
            continue
        qi = ref_to_query[pos]
        if qi < 0 or qi >= len(heights):
            continue
        rb = ref_seq[pos]
        want = rb if orientation == "as_is" else rb.translate(COMPLEMENT)
        h = heights[qi]
        total = sum(h.get(x, 0) for x in "GATC")
        if total <= 0:
            continue
        checked += 1
        if h.get(want, 0) / total >= peak_min:
            ok += 1
    if checked < max(10, flank // 2):
        return False
    return (ok / checked) >= 0.80


def candidate_to_dict(c: InsertionCandidate) -> dict:
    return asdict(c)


def format_evidence(result: CloneInsertResult) -> str:
    return result.evidence
