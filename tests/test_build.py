import numpy as np
import pytest

from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Seq
from vcfixture.build import VcfBuilder


def test_build_biallelic_with_dosage():
    doc = (
        VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 1000)])
        .info("AF", Number.A, Type.FLOAT)
        .fmt("GT")
        .fmt("DS", Number.A, Type.FLOAT)
        .record(
            "chr1",
            81262,
            ref="GAT",
            alt=[Seq("A")],
            gt=["0|1", "1|1"],
            info={"AF": [0.5]},
            DS=[[1.0], [2.0]],
        )
        .build()
    )
    t = doc.truth()
    np.testing.assert_array_equal(t.genotypes[0], [[0, 1], [1, 1]])
    assert t.format[0][0]["DS"] == [1.0]
    assert "##fileformat=VCFv4.5" in doc.render()


def test_reserved_field_by_name():
    doc = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", None)])
        .fmt("GT")
        .record("chr1", 1, ref="A", alt=[Seq("T")], gt=["0/1"])
        .build()
    )
    assert doc.format_defs[0].id == "GT"


def test_undefined_format_field_raises():
    b = VcfBuilder(samples=["s1"], contigs=[("chr1", None)]).fmt("GT")
    with pytest.raises(ValueError, match="not declared"):
        b.record("chr1", 1, ref="A", alt=[Seq("T")], gt=["0/1"], DS=[[1.0]])


def test_gt_index_out_of_range_raises():
    b = VcfBuilder(samples=["s1"], contigs=[("chr1", None)]).fmt("GT")
    with pytest.raises(ValueError, match="allele index"):
        b.record("chr1", 1, ref="A", alt=[Seq("T")], gt=["0/5"])


def test_cardinality_mismatch_raises():
    b = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", None)])
        .fmt("GT")
        .fmt("AD", Number.R, Type.INTEGER)
    )
    with pytest.raises(ValueError, match="cardinality"):
        b.record("chr1", 1, ref="A", alt=[Seq("T")], gt=["0/1"], AD=[[5]])


def test_builder_version_sets_header():
    from vcfixture import VcfBuilder, VcfVersion

    doc = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 1000)], version=VcfVersion.V4_2)
        .fmt("GT")
        .build()
    )
    assert doc.render().startswith("##fileformat=VCFv4.2\n")


def test_builder_rejects_svclaim_before_4_4():
    from vcfixture import VcfBuilder, VcfVersion

    b = VcfBuilder(samples=["s1"], contigs=[("chr1", 1000)], version=VcfVersion.V4_3)
    with pytest.raises(ValueError, match="introduced in VCFv4.4"):
        b.info("SVCLAIM")


def test_svlen_number_a_count_enforced_at_4_4():
    from vcfixture import Sym, VcfBuilder, VcfVersion

    b = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 100000)], version=VcfVersion.V4_4)
        .fmt("GT")
        .info("SVLEN")
    )
    with pytest.raises(ValueError, match="cardinality"):
        # one ALT but two SVLEN values: Number=A requires exactly n_alt
        b.record(
            "chr1", 10, ref="A", alt=[Sym("INS")], gt=["0/1"], info={"SVLEN": [50, 60]}
        )


def test_svlen_any_count_allowed_at_4_3():
    from vcfixture import Sym, VcfBuilder, VcfVersion

    b = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 100000)], version=VcfVersion.V4_3)
        .fmt("GT")
        .info("SVLEN")
    )
    # Number=. at 4.3: no cardinality enforcement; single value for single ALT is fine.
    b.record("chr1", 10, ref="A", alt=[Sym("INS")], gt=["0/1"], info={"SVLEN": [30]})
    assert b.build().render().startswith("##fileformat=VCFv4.3\n")


def test_symbolic_del_no_svclaim_required_before_4_4():
    from vcfixture import Sym, VcfBuilder, VcfVersion

    # DEL requires SVCLAIM at >= 4.4, but SVCLAIM does not exist pre-4.4, so the
    # requirement must not apply there.
    b = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 100000)], version=VcfVersion.V4_3)
        .fmt("GT")
        .info("SVLEN")
    )
    b.record("chr1", 10, ref="A", alt=[Sym("DEL")], gt=["0/1"], info={"SVLEN": [-200]})
    assert "DEL" in b.build().render()
