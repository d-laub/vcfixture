"""Reproduce genoray's hand-authored biallelic fixture via the builder and
assert the derived truth matches the numpy literals currently hand-coded in
genoray/tests/test_vcf.py."""

import numpy as np

from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Seq
from vcfixture.build import VcfBuilder


def _genoray_biallelic():
    return (
        VcfBuilder(
            samples=["sample1", "sample2"],
            contigs=[("chr1", None), ("chr2", None), ("chr3", None)],
        )
        .fmt("GT")
        .fmt("DS", Number.A, Type.FLOAT)
        .record(
            "chr1",
            81262,
            ref="GAT",
            alt=[Seq("A")],
            gt=["0|1", "1|1"],
            DS=[[1.0], [2.0]],
        )
        .record(
            "chr1",
            81262,
            ref="G",
            alt=[Seq("A")],
            gt=["./.", "0/1"],
            DS=[[None], [1.0]],
        )
        .record(
            "chr1",
            81265,
            ref="T",
            alt=[Seq("C")],
            gt=["1|0", "./."],
            DS=[[0.9], [None]],
        )
        .build()
    )


def test_truth_matches_handcoded_genos_for_region():
    truth = _genoray_biallelic().truth()
    np.testing.assert_array_equal(truth.genotypes[0], [[0, 1], [1, 1]])
    np.testing.assert_array_equal(truth.genotypes[1], [[-1, -1], [0, 1]])
    np.testing.assert_array_equal(truth.phasing[0], [True, True])
    np.testing.assert_array_equal(truth.phasing[1], [False, False])
    assert truth.format[0][0]["DS"] == [1.0]
    assert truth.format[1][0]["DS"] == [None]


def test_renders_and_round_trips(tmp_path):
    import pysam

    path = _genoray_biallelic().write(tmp_path / "g.vcf.gz", bgzip=True, index=True)
    vf = pysam.VariantFile(str(path))
    rows = list(vf.fetch("chr1", 81260, 81266))
    assert len(rows) == 3
