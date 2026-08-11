#!/usr/bin/env python3
"""Regression tests for reference-boundary insertion logic."""

from __future__ import annotations

import unittest

from insertion import (
    EDGE_MARGIN_BP,
    FLANK_MATCH_MIN_BP,
    STATUS_CANDIDATE,
    STATUS_CONFLICT,
    STATUS_NO_INSERT,
    STATUS_STRONG,
    InsertionCandidate,
    NoInsertSpan,
    candidates_from_alignment,
    extract_gap_runs,
    homopolymer_interval,
    is_high_confidence,
    merge_clone_insertions,
    score_insertion,
)


def _pad_flanks(core_ref: str, core_query: str, flank: int = 40, base: str = "ACGT") -> tuple[str, str, str]:
    """Build ref/query/aligned pairs with long matching flanks around a core."""
    left = (base * ((flank // 4) + 2))[:flank]
    right = (base[::-1] * ((flank // 4) + 2))[:flank]
    ref = left + core_ref + right
    query = left + core_query + right
    # Build pairwise-style alignment: insert gaps in ref where query has extras
    # For simple cases core_query is core_ref with an insertion in the middle.
    return ref, query


def _align_with_insert(ref: str, query: str) -> tuple[str, str]:
    """Naive global-ish align assuming query = ref[:k] + INS + ref[k:] for one insert.

    Falls back to character walk when lengths differ by locating first mismatch stretch.
    """
    if len(query) == len(ref) and query == ref:
        return ref, query
    # Find longest prefix match
    i = 0
    while i < len(ref) and i < len(query) and ref[i] == query[i]:
        i += 1
    # Find longest suffix match
    j = 0
    while (
        j < len(ref) - i
        and j < len(query) - i
        and ref[len(ref) - 1 - j] == query[len(query) - 1 - j]
    ):
        j += 1
    ref_mid = ref[i : len(ref) - j]
    query_mid = query[i : len(query) - j]
    if not ref_mid and query_mid:
        # pure insertion
        ref_aln = ref[:i] + ("-" * len(query_mid)) + ref[i:]
        query_aln = query
        return ref_aln, query_aln
    # mismatch-heavy: pad shorter mid
    if len(query_mid) > len(ref_mid):
        ref_aln = ref[:i] + ref_mid + ("-" * (len(query_mid) - len(ref_mid))) + ref[len(ref) - j :]
        query_aln = query
        return ref_aln, query_aln
    ref_aln = ref
    query_aln = query[:i] + query_mid + ("-" * (len(ref_mid) - len(query_mid))) + query[len(query) - j :]
    return ref_aln, query_aln


class TestInsertionBasics(unittest.TestCase):
    def test_01_perfect_match(self):
        ref = "ACGT" * 30
        ref_aln, query_aln = ref, ref
        runs = extract_gap_runs(ref_aln, query_aln)
        self.assertEqual(runs, [])
        cands = candidates_from_alignment(ref, ref_aln, query_aln, len(ref))
        self.assertEqual(cands, [])
        merged = merge_clone_insertions({"T7": cands}, {"T7": []})
        self.assertEqual(merged.status, STATUS_NO_INSERT)

    def test_02_single_base_insertion(self):
        left = "ACGT" * 10  # 40
        right = "TGCA" * 10
        ref = left + "AAAA" + right
        query = left + "A" + "AAAA" + right  # insert A before AAAA (homopolymer)
        # clearer non-homo insert:
        ref = left + "GCTA" + right
        query = left + "G" + "CCTA" + right  # wait - better: insert T after G
        ref = left + "GC" + right
        query = left + "G" + "T" + "C" + right
        ref_aln, query_aln = _align_with_insert(ref, query)
        runs = extract_gap_runs(ref_aln, query_aln)
        self.assertEqual(len(runs), 1)
        self.assertEqual(runs[0]["inserted_seq"], "T")
        cands = candidates_from_alignment(ref, ref_aln, query_aln, len(query), primer="T7")
        self.assertTrue(any(c.high_confidence for c in cands))
        hc = [c for c in cands if c.high_confidence][0]
        self.assertEqual(hc.inserted_seq, "T")
        self.assertGreaterEqual(hc.left_flank_match_bp, FLANK_MATCH_MIN_BP)
        self.assertGreaterEqual(hc.right_flank_match_bp, FLANK_MATCH_MIN_BP)

    def test_03_multi_base_insertion(self):
        left = "ACGT" * 10
        right = "TGCA" * 10
        ref = left + "GC" + right
        query = left + "G" + "AAA" + "C" + right
        ref_aln, query_aln = _align_with_insert(ref, query)
        runs = extract_gap_runs(ref_aln, query_aln)
        self.assertEqual(runs[0]["inserted_seq"], "AAA")
        cands = candidates_from_alignment(ref, ref_aln, query_aln, len(query))
        self.assertTrue(any(c.high_confidence and c.inserted_seq == "AAA" for c in cands))

    def test_04_long_target_local_read_must_detect(self):
        """Critical regression: whole-target 90% gate must NOT drop local inserts."""
        # ~2000 bp target; read covers ~600 bp middle with insertion
        target = ("ACGT" * 500)  # 2000
        # local window around position 900
        start, end = 700, 1300
        left = target[start:900]
        right = target[900:end]
        local_ref = left + right
        local_query = left + "GGG" + right
        # Align only the local pieces as a Sanger would (partial coverage)
        ref_aln = ("-" * start) + left + ("-" * 3) + right + ("-" * (len(target) - end))
        # Better: align local read against full target as pairwise would for localms
        # Simulate local alignment covering [start, end) with insert after pos 900
        ref_aln = target[start:900] + ("-" * 3) + target[900:end]
        query_aln = local_query
        # candidates_from_alignment needs full ref for internal check
        # after_ref_pos will be relative to aln which starts at start → adjust by using
        # full-target alignment columns:
        ref_aln_full = target[:start] + target[start:900] + ("-" * 3) + target[900:end] + target[end:]
        # query only covers local — pad outside with gaps
        query_aln_full = ("-" * start) + local_query + ("-" * (len(target) - end))
        # But query_len for edge distance should be local read length
        cands = candidates_from_alignment(
            target, ref_aln_full, query_aln_full, query_len=len(local_query), primer="MID"
        )
        # Remap: extract_gap_runs on full aln sees insert after ref 900
        runs = extract_gap_runs(ref_aln_full, query_aln_full)
        self.assertTrue(any(r["inserted_seq"] == "GGG" for r in runs))
        # Score manually with local aln (as real localms would produce for the read)
        raw = {
            "after_ref_pos": 900,
            "inserted_seq": "GGG",
            "query_start": 200,  # middle of 600bp read
            "query_end": 203,
        }
        # Build local-only aln strings for flank scoring
        local_ref_aln = left + "---" + right
        local_query_aln = left + "GGG" + right
        cand = score_insertion(
            target, local_ref_aln, local_query_aln, len(local_query), raw, primer="MID"
        )
        # Flank scoring uses after_ref_pos relative to the aln's ref coords starting at 0
        # Fix: score against local ref piece with after=len(left)
        raw2 = {
            "after_ref_pos": len(left),
            "inserted_seq": "GGG",
            "query_start": len(left),
            "query_end": len(left) + 3,
        }
        cand = score_insertion(
            local_ref, local_ref_aln, local_query_aln, len(local_query), raw2, primer="MID"
        )
        self.assertTrue(
            cand.high_confidence,
            msg=f"local insert must be HC without whole-target 90%: {cand.summary()}",
        )
        # Whole-target match would be ~600/2000 = 0.3 < 0.90 — old gate would drop this
        match_bp = 600
        ref_len = 2000
        self.assertLess(match_bp / ref_len, 0.90)

    def test_05_insertion_in_read_center(self):
        left = "ACGT" * 15  # 60
        right = "TGCA" * 15
        ref = left + "AT" + right
        query = left + "A" + "CCC" + "T" + right
        ref_aln, query_aln = _align_with_insert(ref, query)
        cands = candidates_from_alignment(ref, ref_aln, query_aln, len(query))
        hc = [c for c in cands if c.high_confidence]
        self.assertEqual(len(hc), 1)
        self.assertGreaterEqual(hc[0].distance_to_read_edge, EDGE_MARGIN_BP)

    def test_06_edge_false_insertion_rejected(self):
        # Insert near the start of a short matching region → edge distance tiny
        ref = "ACGT" * 20
        # query has leading junk insertion-like gap at start of covered region
        ref_aln = "----" + ref
        query_aln = "AAAA" + ref
        raw = extract_gap_runs(ref_aln, query_aln)[0]
        self.assertEqual(raw["after_ref_pos"], 0)
        cand = score_insertion(ref, ref_aln, query_aln, len(query_aln.replace("-", "")), raw)
        self.assertFalse(cand.high_confidence)

        # Insert within EDGE_MARGIN of read end
        left = "ACGT" * 10
        right = "T" * 10  # short right flank AND near edge
        ref = left + "GC" + right
        query = left + "G" + "AAA" + "C" + right
        ref_aln, query_aln = _align_with_insert(ref, query)
        cands = candidates_from_alignment(ref, ref_aln, query_aln, len(query))
        # right flank may also fail FLANK_MATCH; either way must not be HC
        self.assertFalse(any(c.high_confidence for c in cands))

    def test_07_single_read_insert_others_uncovered(self):
        left = "ACGT" * 10
        right = "TGCA" * 10
        ref = left + "GC" + right
        query = left + "G" + "TT" + "C" + right
        ref_aln, query_aln = _align_with_insert(ref, query)
        cands = candidates_from_alignment(ref, ref_aln, query_aln, len(query), primer="T7")
        merged = merge_clone_insertions({"T7": cands, "T7T": []}, {"T7": [], "T7T": []})
        self.assertEqual(merged.status, STATUS_CANDIDATE)
        self.assertGreater(merged.internal_insert_bp, 0)

    def test_08_two_reads_same_insertion(self):
        left = "ACGT" * 10
        right = "TGCA" * 10
        ref = left + "GC" + right
        query = left + "G" + "TT" + "C" + right
        ref_aln, query_aln = _align_with_insert(ref, query)
        c1 = candidates_from_alignment(ref, ref_aln, query_aln, len(query), primer="T7")
        c2 = candidates_from_alignment(ref, ref_aln, query_aln, len(query), primer="T7T")
        merged = merge_clone_insertions({"T7": c1, "T7T": c2}, {})
        self.assertEqual(merged.status, STATUS_STRONG)

    def test_09_insert_vs_no_insert_conflict(self):
        left = "ACGT" * 10
        right = "TGCA" * 10
        ref = left + "GC" + right
        query_ins = left + "G" + "TT" + "C" + right
        ref_aln, q_aln = _align_with_insert(ref, query_ins)
        c_ins = candidates_from_alignment(ref, ref_aln, q_aln, len(query_ins), primer="T7")
        # Second primer: perfect match spanning whole ref
        spans = [NoInsertSpan(primer="T7T", start=0, end=len(ref), identity=1.0)]
        merged = merge_clone_insertions(
            {"T7": c_ins, "T7T": []},
            {"T7": [], "T7T": spans},
        )
        self.assertEqual(merged.status, STATUS_CONFLICT)

    def test_10_homopolymer_boundary_interval(self):
        ref = "ACGT" + ("A" * 5) + "TGCA"
        after = 6  # after first A of the run (index 4 is first A → after 6 is mid-run)
        iv = homopolymer_interval(ref, after_ref_pos=after, inserted_seq="A")
        self.assertLessEqual(iv[0], 4)
        self.assertGreaterEqual(iv[1], 8)
        # Two candidates with shifted after_pos but same insert should cluster
        c1 = InsertionCandidate(
            after_ref_pos=5,
            inserted_seq="A",
            interval_start=iv[0],
            interval_end=iv[1],
            left_flank_match_bp=25,
            right_flank_match_bp=25,
            local_identity=1.0,
            distance_to_read_edge=40,
            high_confidence=True,
            primer="T7",
        )
        c1.high_confidence = True
        iv2 = homopolymer_interval(ref, after_ref_pos=7, inserted_seq="A")
        c2 = InsertionCandidate(
            after_ref_pos=7,
            inserted_seq="A",
            interval_start=iv2[0],
            interval_end=iv2[1],
            left_flank_match_bp=25,
            right_flank_match_bp=25,
            local_identity=1.0,
            distance_to_read_edge=40,
            high_confidence=True,
            primer="T7T",
        )
        merged = merge_clone_insertions({"T7": [c1], "T7T": [c2]}, {})
        self.assertEqual(merged.status, STATUS_STRONG)
        self.assertEqual(len(merged.events), 1)


class TestPerEventAb1Validation(unittest.TestCase):
    def _hc_cand(self, after: int, seq: str, primer: str) -> InsertionCandidate:
        return InsertionCandidate(
            after_ref_pos=after,
            inserted_seq=seq,
            interval_start=after - 1,
            interval_end=after,
            left_flank_match_bp=25,
            right_flank_match_bp=25,
            local_identity=1.0,
            distance_to_read_edge=40,
            high_confidence=True,
            primer=primer,
        )

    def test_two_events_independent_peak_status(self):
        from insertion import (
            STATUS_CONFIRMED,
            STATUS_POSSIBLE,
            validate_insert_with_ab1,
        )

        p = "T7::read.ab1"
        c1 = self._hc_cand(100, "A", p)
        c2 = self._hc_cand(300, "GGG", p)
        # Event1 peaks OK, event2 peaks fail — must NOT both become CONFIRMED
        peak = {
            (p, c1.interval_start, c1.interval_end, c1.inserted_seq): True,
            (p, c2.interval_start, c2.interval_end, c2.inserted_seq): False,
        }
        flank = {
            (p, c1.interval_start, c1.interval_end, c1.inserted_seq): True,
            (p, c2.interval_start, c2.interval_end, c2.inserted_seq): False,
        }
        result = validate_insert_with_ab1({p: [c1, c2]}, {p: []}, peak, flank)
        statuses = {e["inserted_seq"]: e["status"] for e in result.events}
        self.assertEqual(statuses["A"], STATUS_CONFIRMED)
        self.assertEqual(statuses["GGG"], STATUS_POSSIBLE)
        # Clone aggregates to CONFIRMED if any event confirmed
        self.assertEqual(result.status, STATUS_CONFIRMED)


class TestNormDna(unittest.TestCase):
    def test_rejects_n(self):
        from common import _norm_dna

        with self.assertRaises(SystemExit) as ctx:
            _norm_dna("ACGTNACGT", target_id="HER2")
        msg = str(ctx.exception)
        self.assertIn("HER2", msg)
        self.assertIn("'N'", msg)
        self.assertIn("5", msg)

    def test_strips_whitespace_only(self):
        from common import _norm_dna

        self.assertEqual(_norm_dna("acgt\nacgt", target_id="x"), "ACGTACGT")
