from vcfixture import VcfBuilder
from vcfixture.allele import Seq


def _builder() -> VcfBuilder:
    return VcfBuilder(samples=["s0"], contigs=[("chr1", 1000)]).fmt("GT")


def test_labels_default_empty():
    b = _builder().record("chr1", 10, ref="A", alt=[Seq("C")], gt=["0|1"])
    doc = b.build()
    assert doc.records[0].labels == frozenset()
    assert doc.truth().labels == [frozenset()]


def test_labels_carried_record_to_truth():
    b = _builder().record(
        "chr1", 10, ref="A", alt=[Seq("C")], gt=["0|1"], labels=["off_anchor", "x"]
    )
    doc = b.build()
    assert doc.records[0].labels == frozenset({"off_anchor", "x"})
    assert doc.truth().labels == [frozenset({"off_anchor", "x"})]


def test_labels_not_serialized():
    b = _builder().record(
        "chr1", 10, ref="A", alt=[Seq("C")], gt=["0|1"], labels=["off_anchor"]
    )
    text = b.build().render()
    assert "off_anchor" not in text
