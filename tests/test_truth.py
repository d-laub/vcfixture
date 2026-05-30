import numpy as np

from vcfixture.genotype import Genotype
from vcfixture.model import ContigDef, Record, VcfDocument
from vcfixture.truth import derive_truth


def _doc():
    rec = Record(
        chrom="chr1",
        pos=81262,
        ids=None,
        ref="GAT",
        alts=("A",),
        qual=None,
        filters=None,
        info={"AF": [0.5]},
        fmt_keys=("GT", "DS"),
        samples=(
            {"GT": Genotype.parse("0|1"), "DS": [1.0]},
            {"GT": Genotype.parse("./."), "DS": [None]},
        ),
    )
    return VcfDocument(
        "VCFv4.5", (), (), (), (ContigDef("chr1", 1000),), ("s1", "s2"), (rec,)
    )


def test_genotype_matrix_with_missing_sentinel():
    t = derive_truth(_doc())
    assert t.genotypes.shape == (1, 2, 2)
    np.testing.assert_array_equal(t.genotypes[0, 0], [0, 1])
    np.testing.assert_array_equal(t.genotypes[0, 1], [-1, -1])


def test_phasing_matrix():
    t = derive_truth(_doc())
    assert t.phasing.shape == (1, 2)
    assert bool(t.phasing[0, 0]) is True
    assert bool(t.phasing[0, 1]) is False


def test_pos_ref_alt_and_class():
    t = derive_truth(_doc())
    np.testing.assert_array_equal(t.pos, [81262])
    assert t.ref == ["GAT"]
    assert t.alts == [["A"]]
    # GAT->A changes the first base (G->A), so it is a complex substitution
    # (DELINS), not a clean anchored deletion. Matches classify() in variants.py.
    assert t.variant_class == ["DELINS"]


def test_info_and_format_echoed():
    t = derive_truth(_doc())
    assert t.info[0]["AF"] == [0.5]
    assert t.format[0][0]["DS"] == [1.0]
    assert t.format[0][1]["DS"] == [None]
