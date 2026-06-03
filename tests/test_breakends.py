from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Bnd
from vcfixture.build import VcfBuilder

cyvcf2 = pytest.importorskip("cyvcf2")

CASES = ["T[chr2:5[", "]chr2:5]T", "[chr2:5[T", "T]chr2:5]", ".TGCA", "TGCA."]


@pytest.mark.parametrize("bnd", CASES)
def test_breakend_flagged_and_round_trips(bnd):
    doc = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 1000), ("chr2", 1000)])
        .fmt("GT")
        .info("MATEID", Number.A, Type.STRING)
        .record(
            "chr1",
            100,
            ref="T",
            alt=[Bnd.parse(bnd)],
            gt=["0/1"],
            info={"MATEID": ["mate1"]},
        )
        .build()
    )
    t = doc.truth()
    assert t.alts_truth[0][0].kind == "BND"
    assert t.alts_truth[0][0].is_sequence is False
    assert t.variant_class == ["BND"]

    d = tempfile.mkdtemp()
    path = doc.write(Path(d) / "b.vcf.gz", bgzip=True, index=True)
    vf = cyvcf2.VCF(str(path))
    variant = next(iter(vf))
    assert variant.POS == 100
    assert variant.REF == "T"
    assert variant.ALT[0] == bnd
    assert variant.INFO.get("MATEID") == "mate1"
