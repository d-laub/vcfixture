from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pysam

from ._typing import StrPath

_BASES = "ACGT"
_BASES_ARR = np.frombuffer(b"ACGT", dtype="S1")


class Reference:
    """Thin pysam.FastaFile wrapper that draws spec-correct REF/ALT.

    Deterministic given its arguments (alt_index / del_len / ins_seq);
    Hypothesis supplies randomness by choosing those arguments.
    """

    def __init__(self, fasta_path: StrPath):
        self._fa = pysam.FastaFile(str(fasta_path))

    def base(self, contig: str, pos0: int) -> str:
        return self._fa.fetch(contig, pos0, pos0 + 1).upper()

    def seq(self, contig: str, start0: int, length: int) -> str:
        return self._fa.fetch(contig, start0, start0 + length).upper()

    def draw_ref_alt(
        self,
        contig: str,
        pos0: int,
        klass: str,
        *,
        alt_index: int = 1,
        del_len: int = 1,
        ins_seq: str = "T",
        mnp_len: int = 2,
    ) -> tuple[str, list[str]]:
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


@dataclass(frozen=True)
class RepeatFeature:
    """A tandem repeat planted into a reference, for provenance."""

    contig: str
    pos0: int  # 0-based start of the repeat run
    motif: str
    count: int

    @property
    def length(self) -> int:
        return len(self.motif) * self.count


@dataclass(frozen=True)
class ReferenceSpec:
    """Immutable in-memory reference: contig sequences + planted repeats."""

    contigs: tuple[tuple[str, str], ...]  # (id, sequence)
    repeats: tuple[RepeatFeature, ...] = ()

    def _seq_for(self, contig: str) -> str:
        for cid, seq in self.contigs:
            if cid == contig:
                return seq
        raise KeyError(contig)

    def length(self, contig: str) -> int:
        return len(self._seq_for(contig))

    def base(self, contig: str, pos0: int) -> str:
        return self._seq_for(contig)[pos0]

    def seq(self, contig: str, start0: int, length: int) -> str:
        return self._seq_for(contig)[start0 : start0 + length]


class ReferenceBuilder:
    """Mutable builder for a synthetic reference.

    Random-fills contigs (seeded), supports single-base / multi-nucleotide /
    tandem-repeat overwrites, then ``build()`` freezes a ``ReferenceSpec``.
    """

    def __init__(self, seed: int = 0):
        self._rng = np.random.default_rng(seed)
        self._seqs: dict[str, np.ndarray] = {}
        self._order: list[str] = []
        self._repeats: list[RepeatFeature] = []

    def add_contig(self, id: str, length: int) -> ReferenceBuilder:
        if id in self._seqs:
            raise ValueError(f"contig {id!r} already added")
        self._seqs[id] = self._rng.choice(_BASES_ARR, size=length)
        self._order.append(id)
        return self

    def set_base(self, contig: str, pos0: int, base: str) -> ReferenceBuilder:
        if len(base) != 1:
            raise ValueError(f"set_base expects one base, got {base!r}")
        self._seqs[contig][pos0] = base.encode()
        return self

    def set_seq(self, contig: str, pos0: int, seq: str) -> ReferenceBuilder:
        arr = self._seqs[contig]
        if pos0 < 0 or pos0 + len(seq) > arr.size:
            raise ValueError(
                f"set_seq {contig}:{pos0}+{len(seq)} runs past length {arr.size}"
            )
        arr[pos0 : pos0 + len(seq)] = np.frombuffer(seq.encode(), dtype="S1")
        return self

    def tandem_repeat(
        self, contig: str, pos0: int, motif: str, n: int
    ) -> ReferenceBuilder:
        self.set_seq(contig, pos0, motif * n)
        self._repeats.append(RepeatFeature(contig, pos0, motif, n))
        return self

    def build(self) -> ReferenceSpec:
        contigs = tuple(
            (cid, self._seqs[cid].tobytes().decode()) for cid in self._order
        )
        return ReferenceSpec(contigs=contigs, repeats=tuple(self._repeats))
