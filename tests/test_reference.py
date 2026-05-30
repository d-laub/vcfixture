from pathlib import Path

import pysam

from vcfixture.reference import Reference


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
