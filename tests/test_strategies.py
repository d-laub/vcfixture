from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from vcfixture import strategies as S
from vcfixture._spec.number import NumberKind
from vcfixture.model import VcfDocument
from vcfixture.reference import ReferenceSpec


def test_all_number_type_combos_table_is_exhaustive():
    combos = S.ALL_NUMBER_TYPE_COMBOS
    kinds = {n.kind for (n, t, kind) in combos}
    assert {
        NumberKind.FIXED,
        NumberKind.A,
        NumberKind.R,
        NumberKind.G,
        NumberKind.DOT,
        NumberKind.FLAG,
    } <= kinds
    for _n, t, kind in combos:
        if t.value == "Flag":
            assert kind == "INFO"


def test_all_variant_classes_present():
    assert set(S.ALL_VARIANT_CLASSES) >= {
        "SNP",
        "MNP",
        "INS",
        "DEL",
        "DELINS",
        "SPANNING_DEL",
    }


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
    for rec in doc.records:
        n_alt = len(rec.alts)
        for s in rec.samples:
            gt = s["GT"]
            for a in gt.alleles:
                assert a is None or a <= n_alt


@settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
@given(S.references())
def test_references_are_well_formed(spec: ReferenceSpec):
    assert isinstance(spec, ReferenceSpec)
    assert len(spec.contigs) >= 1
    for _cid, seq in spec.contigs:
        assert len(seq) >= 1 and set(seq) <= set("ACGT")
    # planted repeats actually appear at their advertised loci
    for rf in spec.repeats:
        assert spec.seq(rf.contig, rf.pos0, rf.length) == rf.motif * rf.count


@settings(
    max_examples=40,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much],
)
@given(st.data())
def test_reference_consistent_and_labeled(data):
    spec = data.draw(S.references(max_repeats=3))
    doc = data.draw(
        S.documents(
            reference=spec,
            violations=frozenset({"multiallelic", "non_atomic", "non_left_aligned"}),
        )
    )
    truth = doc.truth()
    # Every REF matches the reference sequence at its position.
    for rec in doc.records:
        assert spec.seq(rec.chrom, rec.pos - 1, len(rec.ref)) == rec.ref
    # Records are position-sorted per contig (norm/consensus/gvl require this).
    last: dict[str, int] = {}
    for rec in doc.records:
        assert rec.pos >= last.get(rec.chrom, 0)
        last[rec.chrom] = rec.pos
    # truth lines up with the document.
    assert truth.genotypes.shape[0] == len(doc.records)
    assert truth.labels == [r.labels for r in doc.records]
    # Provenance labels only use the known vocabulary.
    allowed = {"multiallelic", "non_atomic", "off_anchor", "tandem_repeat"}
    for lbls in truth.labels:
        assert lbls <= allowed


@settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow])
@given(S.documents())
def test_documents_back_compat_unlabeled(doc):
    # Reference-free documents still work and carry no labels.
    assert all(r.labels == frozenset() for r in doc.records)


@settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow])
@given(S.reference_and_documents())
def test_reference_and_documents_tuple(triple):
    spec, doc, truth = triple
    assert isinstance(spec, ReferenceSpec)
    assert truth.genotypes.shape[0] == len(doc.records)
    # Default (no violations) -> reference-consistent, unlabeled.
    for rec in doc.records:
        assert spec.seq(rec.chrom, rec.pos - 1, len(rec.ref)) == rec.ref
        assert rec.labels == frozenset()
