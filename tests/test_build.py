import numpy as np
import pytest
from vcforge.build import VcfBuilder
from vcforge._spec.number import Number
from vcforge._spec.types import Type

def test_build_biallelic_with_dosage():
    doc = (VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 1000)])
           .info("AF", Number.A, Type.FLOAT)
           .fmt("GT")
           .fmt("DS", Number.A, Type.FLOAT)
           .record("chr1", 81262, ref="GAT", alt=["A"],
                   gt=["0|1", "1|1"], info={"AF": [0.5]},
                   DS=[[1.0], [2.0]])
           .build())
    t = doc.truth()
    np.testing.assert_array_equal(t.genotypes[0], [[0, 1], [1, 1]])
    assert t.format[0][0]["DS"] == [1.0]
    assert "##fileformat=VCFv4.5" in doc.render()

def test_reserved_field_by_name():
    doc = (VcfBuilder(samples=["s1"], contigs=[("chr1", None)])
           .fmt("GT")
           .record("chr1", 1, ref="A", alt=["T"], gt=["0/1"])
           .build())
    assert doc.format_defs[0].id == "GT"

def test_undefined_format_field_raises():
    b = (VcfBuilder(samples=["s1"], contigs=[("chr1", None)]).fmt("GT"))
    with pytest.raises(ValueError, match="not declared"):
        b.record("chr1", 1, ref="A", alt=["T"], gt=["0/1"], DS=[[1.0]])

def test_gt_index_out_of_range_raises():
    b = (VcfBuilder(samples=["s1"], contigs=[("chr1", None)]).fmt("GT"))
    with pytest.raises(ValueError, match="allele index"):
        b.record("chr1", 1, ref="A", alt=["T"], gt=["0/5"])

def test_cardinality_mismatch_raises():
    b = (VcfBuilder(samples=["s1"], contigs=[("chr1", None)])
         .fmt("GT").fmt("AD", Number.R, Type.INTEGER))
    with pytest.raises(ValueError, match="cardinality"):
        b.record("chr1", 1, ref="A", alt=["T"], gt=["0/1"], AD=[[5]])
