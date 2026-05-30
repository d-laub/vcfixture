from vcfixture._spec.reserved import reserved
from vcfixture.genotype import Genotype
from vcfixture.model import ContigDef, Record, VcfDocument


def _doc():
    rec = Record(
        chrom="chr1",
        pos=81262,
        ids=None,
        ref="GAT",
        alts=("A",),
        qual=None,
        filters=None,
        info={},
        fmt_keys=("GT", "DS"),
        samples=(
            {"GT": Genotype.parse("0|1"), "DS": [1.0]},
            {"GT": Genotype.parse("1|1"), "DS": [2.0]},
        ),
    )
    return VcfDocument(
        fileformat="VCFv4.5",
        info_defs=(),
        format_defs=(reserved("GT", "FORMAT"),),
        filter_defs=(),
        contigs=(ContigDef("chr1", 1000),),
        samples=("s1", "s2"),
        records=(rec,),
    )


def test_document_is_frozen_and_holds_records():
    doc = _doc()
    assert doc.samples == ("s1", "s2")
    assert doc.records[0].alts == ("A",)
    assert doc.max_ploidy() == 2


def test_ploidy_varies_uses_max():
    rec = Record(
        "chr1",
        1,
        None,
        "A",
        ("T",),
        None,
        None,
        {},
        ("GT",),
        ({"GT": Genotype.parse("0/0/1")},),
    )
    doc = VcfDocument(
        "VCFv4.5", (), (), (), (ContigDef("chr1", None),), ("s1",), (rec,)
    )
    assert doc.max_ploidy() == 3
