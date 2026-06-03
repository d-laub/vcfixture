from vcfixture import VcfVersion
from vcfixture._spec.fielddef import FieldDef
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Seq
from vcfixture.genotype import Genotype
from vcfixture.model import ContigDef, Record, VcfDocument
from vcfixture.serialize import render_document


def _doc():
    af = FieldDef("AF", Number.A, Type.FLOAT, "Allele frequency", "INFO")
    db = FieldDef("DB", Number.FLAG, Type.FLAG, "dbSNP", "INFO")
    gt = FieldDef("GT", Number.ONE, Type.STRING, "Genotype", "FORMAT")
    ds = FieldDef("DS", Number.A, Type.FLOAT, "Dosage", "FORMAT")
    rec = Record(
        chrom="chr1",
        pos=81262,
        ids=None,
        ref="GAT",
        alts=(Seq("A"),),
        qual=None,
        filters=None,
        info={"AF": [0.5], "DB": True},
        fmt_keys=("GT", "DS"),
        samples=(
            {"GT": Genotype.parse("0|1"), "DS": [1.0]},
            {"GT": Genotype.parse("./."), "DS": [None]},
        ),
    )
    return VcfDocument(
        VcfVersion.V4_5,
        (af, db),
        (gt, ds),
        (),
        (ContigDef("chr1", 1000),),
        ("s1", "s2"),
        (rec,),
    )


def test_header_first_line_is_fileformat():
    text = render_document(_doc())
    assert text.splitlines()[0] == "##fileformat=VCFv4.5"


def test_chrom_header_includes_samples():
    text = render_document(_doc())
    chrom = [line for line in text.splitlines() if line.startswith("#CHROM")][0]
    assert chrom.split("\t")[-2:] == ["s1", "s2"]
    assert chrom.split("\t")[8] == "FORMAT"


def test_record_line_fields():
    text = render_document(_doc())
    data = [line for line in text.splitlines() if not line.startswith("#")][0]
    cols = data.split("\t")
    assert cols[0:5] == ["chr1", "81262", ".", "GAT", "A"]
    assert cols[5] == "."
    assert cols[6] == "."
    assert cols[7] == "AF=0.5;DB"
    assert cols[8] == "GT:DS"
    assert cols[9] == "0|1:1.0"
    assert cols[10] == "./.:."


def test_non_finite_floats_render_as_missing():
    # VCF has no nan/inf literal; non-finite floats must serialize to ".".
    from vcfixture.serialize import _fmt_scalar, _fmt_value

    assert _fmt_scalar(float("nan")) == "."
    assert _fmt_scalar(float("inf")) == "."
    assert _fmt_value([float("nan"), 0.5]) == ".,0.5"


def test_percent_encoding_of_reserved_chars():
    from vcfixture.serialize import _encode, _fmt_scalar

    assert _encode("a;b") == "a%3Bb"
    assert _encode("a:b,c=d") == "a%3Ab%2Cc%3Dd"
    assert _encode("100%") == "100%25"
    assert _encode("x\ty\n") == "x%09y%0A"
    assert _fmt_scalar("a;b") == "a%3Bb"
    assert _fmt_scalar(5) == "5"
    assert _fmt_scalar(0.5) == "0.5"


def test_pass_filter_renders_PASS():
    doc = _doc()
    rec = doc.records[0]
    rec2 = type(rec)(**{**rec.__dict__, "filters": ()})
    doc2 = type(doc)(**{**doc.__dict__, "records": (rec2,)})
    lines = render_document(doc2).splitlines()
    data = [line for line in lines if not line.startswith("#")][0]
    assert data.split("\t")[6] == "PASS"
