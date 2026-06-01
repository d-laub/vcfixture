from pathlib import Path

import pysam

from vcfixture.reference import Reference, ReferenceBuilder, RepeatFeature


def _make_fasta(tmp_path: Path) -> Path:
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "ACGTACGTAC" * 5 + "\n")
    pysam.faidx(str(fa))
    return fa


def test_base_and_seq(tmp_path):
    ref = Reference(_make_fasta(tmp_path))
    assert ref.base("chr1", 0) == "A"
    assert ref.seq("chr1", 0, 4) == "ACGT"


def test_draw_snp_ref_matches_sequence(tmp_path):
    ref = Reference(_make_fasta(tmp_path))
    rec_ref, alts = ref.draw_ref_alt("chr1", pos0=0, klass="SNP", alt_index=1)
    assert rec_ref == "A"
    assert alts[0] != "A" and len(alts[0]) == 1


def test_draw_deletion_ref_starts_with_sequence(tmp_path):
    ref = Reference(_make_fasta(tmp_path))
    rec_ref, alts = ref.draw_ref_alt("chr1", pos0=0, klass="DEL", del_len=2)
    assert rec_ref == ref.seq("chr1", 0, 3)
    assert alts == [rec_ref[0]]


def test_builder_random_fill_is_seeded_and_acgt():
    a = ReferenceBuilder(seed=0).add_contig("chr1", 200).build()
    b = ReferenceBuilder(seed=0).add_contig("chr1", 200).build()
    assert a.contigs == b.contigs  # deterministic
    seq = a.seq("chr1", 0, 200)
    assert len(seq) == 200 and set(seq) <= set("ACGT")


def test_set_base_and_set_seq_overwrite():
    spec = (
        ReferenceBuilder(seed=1)
        .add_contig("chr1", 100)
        .set_base("chr1", 10, "A")
        .set_seq("chr1", 20, "GATTACA")
        .build()
    )
    assert spec.base("chr1", 10) == "A"
    assert spec.seq("chr1", 20, 7) == "GATTACA"


def test_tandem_repeat_writes_and_records_feature():
    spec = (
        ReferenceBuilder(seed=2)
        .add_contig("chr1", 100)
        .tandem_repeat("chr1", 30, "AG", 5)
        .build()
    )
    assert spec.seq("chr1", 30, 10) == "AGAGAGAGAG"
    assert spec.repeats == (RepeatFeature("chr1", 30, "AG", 5),)
    assert spec.repeats[0].length == 10


def test_set_seq_out_of_bounds_raises():
    import pytest

    rb = ReferenceBuilder(seed=0).add_contig("chr1", 10)
    with pytest.raises(ValueError):
        rb.set_seq("chr1", 8, "ACGT")  # runs past length 10
