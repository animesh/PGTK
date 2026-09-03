#!/usr/bin/env python3
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_finding_igv_reviews as review
import prepare_event_igv_tracks as tracks


class FakeRead:
    def __init__(self, name, start, flag, cigar, sequence, blocks):
        self.query_name = name
        self.reference_name = '19'
        self.reference_start = start
        self.flag = flag
        self.cigarstring = cigar
        self.query_sequence = sequence
        self._blocks = blocks
        self.is_unmapped = False
        self.is_secondary = False
        self.is_supplementary = False
        self.is_duplicate = False
        self.is_qcfail = False
        self.cigartuples = []

    def get_blocks(self):
        return self._blocks


splice_span = FakeRead('read1', 100, 0, '10M1000N10M', 'A' * 20, [(100, 110), (1110, 1120)])
aligned = FakeRead('read1', 100, 0, '20M', 'A' * 20, [(100, 120)])
assert not review.read_observes_interval(splice_span, 500, 501)
assert review.read_observes_interval(aligned, 105, 106)
assert review.alignment_key_values(splice_span, '19')[0] != review.alignment_key_values(aligned, '19')[0]
assert tracks.alignment_key(splice_span, '19') == review.alignment_key_values(splice_span, '19')[0]
assert tracks.alignment_key(aligned, '19') == review.alignment_key_values(aligned, '19')[0]
aligned.cigartuples = [(0, 10), (2, 1), (0, 9)]
assert review.browser_safe_display(aligned, 'A', 'G')
print('IGV event identity regression tests: PASS')
