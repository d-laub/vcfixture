from hypothesis import given, settings, HealthCheck
from vcforge import strategies as S
from vcforge.model import VcfDocument
from vcforge._spec.number import NumberKind

def test_all_number_type_combos_table_is_exhaustive():
    combos = S.ALL_NUMBER_TYPE_COMBOS
    kinds = {n.kind for (n, t, kind) in combos}
    assert {NumberKind.FIXED, NumberKind.A, NumberKind.R,
            NumberKind.G, NumberKind.DOT, NumberKind.FLAG} <= kinds
    for (n, t, kind) in combos:
        if t.value == "Flag":
            assert kind == "INFO"

def test_all_variant_classes_present():
    assert set(S.ALL_VARIANT_CLASSES) >= {
        "SNP", "MNP", "INS", "DEL", "DELINS", "SPANNING_DEL"}

@settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
@given(S.documents())
def test_documents_are_well_formed(doc: VcfDocument):
    assert isinstance(doc, VcfDocument)
    t = doc.truth()
    assert t.genotypes.shape[0] == len(doc.records)
    assert doc.render().startswith("##fileformat=")

@settings(max_examples=80, suppress_health_check=[HealthCheck.too_slow])
@given(S.documents(max_alt=3))
def test_documents_can_be_multiallelic(doc):
    from vcforge.genotype import Genotype
    for rec in doc.records:
        n_alt = len(rec.alts)
        for s in rec.samples:
            gt = s["GT"]
            for a in gt.alleles:
                assert a is None or a <= n_alt
