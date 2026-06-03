from hypothesis.vendor.pretty import pretty

from vcfixture import VcfVersion
from vcfixture._spec.fielddef import FieldDef
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Seq
from vcfixture.genotype import Genotype
from vcfixture.model import ContigDef, Record, VcfDocument
from vcfixture.reference import ReferenceSpec, RepeatFeature


def test_number_singletons_not_in_dataclass_fields():
    # ClassVar singletons must not leak into the dataclass field set, or any
    # field-walking pretty-printer (Hypothesis) recurses into them forever.
    assert set(Number.__dataclass_fields__) == {"kind", "count"}


def test_number_pretty_does_not_recurse():
    # Hypothesis's vendored pretty-printer is what produces "Falsifying example"
    # output. It must not explode Number into its singletons.
    out = pretty(Number.G)
    assert "ONE=Number" not in out
    assert len(out) < 200


def test_number_repr_compact():
    assert repr(Number.G) == "Number(G)"
    assert repr(Number.A) == "Number(A)"
    assert repr(Number.R) == "Number(R)"
    assert repr(Number.DOT) == "Number(.)"
    assert repr(Number.ONE) == "Number(1)"
    assert repr(Number.fixed(2)) == "Number(2)"
    assert repr(Number.FLAG) == "Number(FLAG)"


def test_number_pretty_uses_compact_repr():
    # _repr_pretty_ must route Hypothesis's printer through __repr__.
    assert pretty(Number.G) == "Number(G)"


def test_fielddef_repr_compact():
    gt = FieldDef(
        id="GT",
        number=Number.ONE,
        type=Type.STRING,
        description="Genotype",
        kind="FORMAT",
    )
    assert repr(gt) == "FieldDef(GT FORMAT Number=1 Type=String)"

    dp = FieldDef(
        id="DP", number=Number.ONE, type=Type.INTEGER, description="Depth", kind="INFO"
    )
    assert repr(dp) == "FieldDef(DP INFO Number=1 Type=Integer)"


def test_genotype_repr_compact():
    assert repr(Genotype((0, 1), (True,))) == "Genotype(0|1)"
    assert repr(Genotype((0, 1), (False,))) == "Genotype(0/1)"
    assert repr(Genotype((None, None), (False,))) == "Genotype(./.)"


def test_contigdef_repr_compact():
    assert repr(ContigDef(id="chr1", length=200)) == "ContigDef(chr1:200)"
    assert repr(ContigDef(id="chr1")) == "ContigDef(chr1)"


_DEFAULT_ALTS = (Seq("T"), Seq("G"))


def _make_record(alts=_DEFAULT_ALTS, labels=frozenset()):
    return Record(
        chrom="chr1",
        pos=5,
        ids=None,
        ref="A",
        alts=alts,
        qual=None,
        filters=None,
        info={},
        fmt_keys=("GT",),
        samples=(
            {"GT": Genotype((0, 1), (True,))},
            {"GT": Genotype((1, 1), (True,))},
        ),
        labels=labels,
    )


def test_record_repr_compact():
    rec = _make_record(labels=frozenset({"multiallelic", "snp"}))
    assert repr(rec) == "Record(chr1:5 A>T,G ×2 [multiallelic,snp])"


def test_record_repr_no_labels():
    rec = _make_record(labels=frozenset())
    assert repr(rec) == "Record(chr1:5 A>T,G ×2)"


def test_record_repr_empty_alts():
    rec = _make_record(alts=(), labels=frozenset())
    assert repr(rec) == "Record(chr1:5 A>. ×2)"


def test_vcfdocument_repr_compact():
    gt = FieldDef(
        id="GT",
        number=Number.ONE,
        type=Type.STRING,
        description="Genotype",
        kind="FORMAT",
    )
    doc = VcfDocument(
        version=VcfVersion.V4_5,
        info_defs=(),
        format_defs=(gt,),
        filter_defs=(),
        contigs=(),
        samples=("s0", "s1"),
        records=(_make_record(), _make_record()),
    )
    assert repr(doc) == "VcfDocument(VCFv4.5 samples=2 records=2 info=0 format=1)"


def test_repeatfeature_repr_compact():
    rf = RepeatFeature(contig="chr1", pos0=10, motif="AT", count=3)
    assert repr(rf) == "RepeatFeature(chr1@10 AT×3)"


def test_referencespec_repr_compact():
    spec = ReferenceSpec(
        contigs=(("chr1", "A" * 200),),
        repeats=(RepeatFeature(contig="chr1", pos0=10, motif="AT", count=3),),
    )
    assert repr(spec) == "ReferenceSpec(contigs=[chr1:200bp], repeats=1)"


def test_referencespec_repr_pretty_no_full_sequence():
    spec = ReferenceSpec(contigs=(("chr1", "ACGT" * 50),), repeats=())
    out = pretty(spec)
    assert out == "ReferenceSpec(contigs=[chr1:200bp], repeats=0)"
    assert "ACGTACGT" not in out


def test_full_document_pretty_is_bounded():
    # The original bug: a VcfDocument printed ~2k lines via Hypothesis's
    # pretty-printer. With compact reprs + _repr_pretty_, the whole document
    # collapses to a single short line.
    gt = FieldDef(
        id="GT",
        number=Number.ONE,
        type=Type.STRING,
        description="Genotype",
        kind="FORMAT",
    )
    doc = VcfDocument(
        version=VcfVersion.V4_5,
        info_defs=(),
        format_defs=(gt,),
        filter_defs=(),
        contigs=(ContigDef("chr1", 200),),
        samples=("s0", "s1"),
        records=tuple(_make_record() for _ in range(5)),
    )
    out = pretty(doc)
    assert out == "VcfDocument(VCFv4.5 samples=2 records=5 info=0 format=1)"
    assert "ONE=Number" not in out
    assert len(out) < 200
