from __future__ import annotations
import pysam

_BASES = "ACGT"

class Reference:
    """Thin pysam.FastaFile wrapper that draws spec-correct REF/ALT.

    Deterministic given its arguments (alt_index / del_len / ins_seq);
    Hypothesis supplies randomness by choosing those arguments.
    """

    def __init__(self, fasta_path):
        self._fa = pysam.FastaFile(str(fasta_path))

    def base(self, contig: str, pos0: int) -> str:
        return self._fa.fetch(contig, pos0, pos0 + 1).upper()

    def seq(self, contig: str, start0: int, length: int) -> str:
        return self._fa.fetch(contig, start0, start0 + length).upper()

    def draw_ref_alt(self, contig, pos0, klass, *, alt_index=1,
                     del_len=1, ins_seq="T", mnp_len=2):
        if klass == "SNP":
            r = self.base(contig, pos0)
            alt = _BASES[(_BASES.index(r) + alt_index) % 4]
            return r, [alt]
        if klass == "MNP":
            r = self.seq(contig, pos0, mnp_len)
            alt = "".join(_BASES[(_BASES.index(b) + alt_index) % 4] for b in r)
            return r, [alt]
        if klass == "INS":
            anchor = self.base(contig, pos0)
            return anchor, [anchor + ins_seq]
        if klass == "DEL":
            r = self.seq(contig, pos0, del_len + 1)
            return r, [r[0]]
        if klass == "DELINS":
            r = self.seq(contig, pos0, mnp_len)
            return r, [ins_seq]
        if klass == "SPANNING_DEL":
            return self.base(contig, pos0), ["*"]
        raise ValueError(f"unknown class {klass!r}")
