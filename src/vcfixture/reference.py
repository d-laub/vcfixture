from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam

from ._typing import StrPath

_BASES = "ACGT"
_BASES_ARR = np.frombuffer(b"ACGT", dtype="S1")


def _draw_ref_alt(
    base_fn: Callable[[str, int], str],
    seq_fn: Callable[[str, int, int], str],
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
        r = base_fn(contig, pos0)
        alt = _BASES[(_BASES.index(r) + alt_index) % 4]
        return r, [alt]
    if klass == "MNP":
        r = seq_fn(contig, pos0, mnp_len)
        alt = "".join(_BASES[(_BASES.index(b) + alt_index) % 4] for b in r)
        return r, [alt]
    if klass == "INS":
        anchor = base_fn(contig, pos0)
        return anchor, [anchor + ins_seq]
    if klass == "DEL":
        r = seq_fn(contig, pos0, del_len + 1)
        return r, [r[0]]
    if klass == "DELINS":
        r = seq_fn(contig, pos0, mnp_len)
        return r, [ins_seq]
    if klass == "SPANNING_DEL":
        return base_fn(contig, pos0), ["*"]
    raise ValueError(f"unknown class {klass!r}")


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
        return _draw_ref_alt(
            self.base,
            self.seq,
            contig,
            pos0,
            klass,
            alt_index=alt_index,
            del_len=del_len,
            ins_seq=ins_seq,
            mnp_len=mnp_len,
        )


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

    def write(self, path: StrPath, *, bgzip: bool = True, index: bool = True) -> Path:
        """Write a 60-col FASTA; bgzip + faidx it via pysam. Returns the path."""
        path = Path(path)
        text_lines: list[str] = []
        for cid, seq in self.contigs:
            text_lines.append(f">{cid}")
            text_lines.extend(seq[i : i + 60] for i in range(0, len(seq), 60))
        fasta_text = "\n".join(text_lines) + "\n"

        if bgzip:
            plain = path.with_name(path.name + ".tmp.fa")
            plain.write_text(fasta_text)
            pysam.tabix_compress(str(plain), str(path), force=True)
            plain.unlink()
        else:
            path.write_text(fasta_text)

        if index:
            pysam.faidx(str(path))  # writes <path>.fai (+ .gzi when bgzipped)
        return path

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
        return _draw_ref_alt(
            self.base,
            self.seq,
            contig,
            pos0,
            klass,
            alt_index=alt_index,
            del_len=del_len,
            ins_seq=ins_seq,
            mnp_len=mnp_len,
        )


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
