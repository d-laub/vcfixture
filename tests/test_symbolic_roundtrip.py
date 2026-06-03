from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings

from vcfixture import strategies as S

cyvcf2 = pytest.importorskip("cyvcf2")


@settings(max_examples=50, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(S.symbolic_documents())
def test_symbolic_alts_round_trip_through_cyvcf2(doc):
    truth = doc.truth()
    d = tempfile.mkdtemp()
    path = doc.write(Path(d) / "x.vcf.gz", bgzip=True, index=True)
    vf = cyvcf2.VCF(str(path))
    for ri, variant in enumerate(vf):
        for ai, alt in enumerate(variant.ALT):
            at = truth.alts_truth[ri][ai]
            assert alt == doc.records[ri].alts[ai].render()
            if at.kind == "SYMBOLIC":
                assert at.is_sequence is False
                got = variant.INFO.get("SVLEN")
                got_i = got[ai] if isinstance(got, (list, tuple)) else got
                assert abs(int(got_i)) == at.svlen
