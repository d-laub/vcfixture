import os
import tempfile

import pysam
import pytest

from vcfixture.build import VcfBuilder
from vcfixture.reference import Reference

cyvcf2 = pytest.importorskip("cyvcf2")


def _ref(tmp):
    fa = os.path.join(tmp, "ref.fa")
    with open(fa, "w") as f:
        f.write(">chr1\n" + "ACGTACGTAC" * 10 + "\n")
    pysam.faidx(fa)
    return Reference(fa), fa


def test_reference_anchored_doc_round_trips():
    tmp = tempfile.mkdtemp()
    ref, _ = _ref(tmp)
    b = VcfBuilder(samples=["s1"], contigs=[("chr1", 100)]).fmt("GT")
    specs = [(5, "SNP"), (15, "DEL"), (25, "INS"), (35, "MNP")]
    expected = []
    for pos0, klass in specs:
        rref, alts = ref.draw_ref_alt("chr1", pos0, klass=klass)
        b.record("chr1", pos0 + 1, ref=rref, alt=alts, gt=["0|1"])  # 1-based POS
        expected.append((pos0 + 1, rref))
    doc = b.build()
    path = doc.write(os.path.join(tmp, "r.vcf.gz"), bgzip=True, index=True)
    vf = cyvcf2.VCF(str(path))
    got = [(v.POS, v.REF) for v in vf]
    assert got == expected
    for pos1, rref in expected:
        assert ref.seq("chr1", pos1 - 1, len(rref)) == rref
