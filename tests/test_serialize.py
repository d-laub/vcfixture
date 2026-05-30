from vcforge.model import Record, ContigDef, VcfDocument
from vcforge.genotype import Genotype
from vcforge.serialize import render_document
from vcforge._spec.fielddef import FieldDef
from vcforge._spec.number import Number
from vcforge._spec.types import Type

def _doc():
    af = FieldDef("AF", Number.A, Type.FLOAT, "Allele frequency", "INFO")
    db = FieldDef("DB", Number.FLAG, Type.FLAG, "dbSNP", "INFO")
    gt = FieldDef("GT", Number.ONE, Type.STRING, "Genotype", "FORMAT")
    ds = FieldDef("DS", Number.A, Type.FLOAT, "Dosage", "FORMAT")
    rec = Record(
        chrom="chr1", pos=81262, ids=None, ref="GAT", alts=("A",),
        qual=None, filters=None,
        info={"AF": [0.5], "DB": True},
        fmt_keys=("GT", "DS"),
        samples=({"GT": Genotype.parse("0|1"), "DS": [1.0]},
                 {"GT": Genotype.parse("./."), "DS": [None]}),
    )
    return VcfDocument("VCFv4.5", (af, db), (gt, ds), (),
                       (ContigDef("chr1", 1000),), ("s1", "s2"), (rec,))

def test_header_first_line_is_fileformat():
    text = render_document(_doc())
    assert text.splitlines()[0] == "##fileformat=VCFv4.5"

def test_chrom_header_includes_samples():
    text = render_document(_doc())
    chrom = [l for l in text.splitlines() if l.startswith("#CHROM")][0]
    assert chrom.split("\t")[-2:] == ["s1", "s2"]
    assert chrom.split("\t")[8] == "FORMAT"

def test_record_line_fields():
    text = render_document(_doc())
    data = [l for l in text.splitlines() if not l.startswith("#")][0]
    cols = data.split("\t")
    assert cols[0:5] == ["chr1", "81262", ".", "GAT", "A"]
    assert cols[5] == "."
    assert cols[6] == "."
    assert cols[7] == "AF=0.5;DB"
    assert cols[8] == "GT:DS"
    assert cols[9] == "0|1:1.0"
    assert cols[10] == "./.:."

def test_pass_filter_renders_PASS():
    doc = _doc()
    rec = doc.records[0]
    rec2 = type(rec)(**{**rec.__dict__, "filters": ()})
    doc2 = type(doc)(**{**doc.__dict__, "records": (rec2,)})
    data = [l for l in render_document(doc2).splitlines() if not l.startswith("#")][0]
    assert data.split("\t")[6] == "PASS"
