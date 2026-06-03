from __future__ import annotations

import numpy as np

from vcfixture import VcfVersion
from vcfixture.allele import Seq, Star, Sym, Unspecified
from vcfixture.genotype import Genotype
from vcfixture.model import ContigDef, Record, VcfDocument


def _doc(alts: tuple, info: dict) -> VcfDocument:
    rec = Record(
        chrom="chr1",
        pos=100,
        ids=None,
        ref="G",
        alts=alts,
        qual=None,
        filters=None,
        info=info,
        fmt_keys=("GT",),
        samples=({"GT": Genotype.parse("0/1")},),
    )
    return VcfDocument(
        version=VcfVersion.V4_5,
        info_defs=(),
        format_defs=(),
        filter_defs=(),
        contigs=(ContigDef("chr1", 1000),),
        samples=("s1",),
        records=(rec,),
    )


def test_symbolic_del_geometry_and_flag():
    t = _doc((Sym.deletion(),), {"SVLEN": [50]}).truth()
    at = t.alts_truth[0][0]
    assert at.kind == "SYMBOLIC" and at.is_sequence is False
    assert at.sv_type == "DEL" and at.svlen == 50 and at.sv_end == 150
    assert t.variant_class == ["SV_DEL"]


def test_insertion_has_no_end():
    t = _doc((Sym.insertion(),), {"SVLEN": [30]}).truth()
    at = t.alts_truth[0][0]
    assert at.svlen == 30 and at.sv_end is None


def test_negative_svlen_normalized_absolute():
    t = _doc((Sym.deletion(),), {"SVLEN": [-50]}).truth()
    assert t.alts_truth[0][0].svlen == 50


def test_multiallelic_per_allele_svlen_indexing():
    t = _doc((Sym.deletion(), Sym.duplication()), {"SVLEN": [50, 80]}).truth()
    a0, a1 = t.alts_truth[0]
    assert a0.svlen == 50 and a0.sv_end == 150
    assert a1.svlen == 80 and a1.sv_end == 180


def test_scalar_svlen_normalized_for_single_alt():
    t = _doc((Sym.deletion(),), {"SVLEN": 50}).truth()
    assert t.alts_truth[0][0].svlen == 50


def test_unspecified_and_spanning_and_mixed_mask():
    t = _doc((Seq("A"), Star(), Unspecified()), {}).truth()
    kinds = [a.kind for a in t.alts_truth[0]]
    assert kinds == ["SNP", "SPANNING_DEL", "UNSPECIFIED"]
    np.testing.assert_array_equal(t.is_sequence_mask[0], [True, False, False])
    assert t.variant_class == ["MULTIALLELIC"]
