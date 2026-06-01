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


def test_referencespec_write_roundtrips(tmp_path):
    spec = (
        ReferenceBuilder(seed=3)
        .add_contig("chr1", 300)
        .add_contig("chr2", 150)
        .set_seq("chr1", 50, "GATTACA")
        .build()
    )
    out = spec.write(tmp_path / "ref.fa.bgz")
    assert out.exists()
    assert (out.parent / (out.name + ".fai")).exists()
    with pysam.FastaFile(str(out)) as fa:
        assert fa.fetch("chr1", 50, 57).upper() == "GATTACA"
        assert fa.fetch("chr1", 0, 300).upper() == spec.seq("chr1", 0, 300)
        assert fa.references == ("chr1", "chr2") or set(fa.references) == {
            "chr1",
            "chr2",
        }


def test_referencespec_write_plain(tmp_path):
    spec = ReferenceBuilder(seed=4).add_contig("chr1", 80).build()
    out = spec.write(tmp_path / "ref.fa", bgzip=False)
    assert out.exists()
    with pysam.FastaFile(str(out)) as fa:
        assert fa.fetch("chr1", 0, 80).upper() == spec.seq("chr1", 0, 80)


def test_referencespec_draw_ref_alt_matches_sequence():
    spec = (
        ReferenceBuilder(seed=5)
        .add_contig("chr1", 100)
        .set_seq("chr1", 10, "ACGT")
        .build()
    )
    ref, alts = spec.draw_ref_alt("chr1", pos0=10, klass="SNP", alt_index=1)
    assert ref == "A" and len(alts) == 1 and alts[0] != "A"
    dref, dalts = spec.draw_ref_alt("chr1", pos0=10, klass="DEL", del_len=2)
    assert dref == spec.seq("chr1", 10, 3) and dalts == [dref[0]]
