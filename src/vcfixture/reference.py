from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam

from ._repr import CompactRepr, override
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
    """Reference backed by an existing FASTA file.

    Wraps a ``pysam.FastaFile`` to provide single-base and subsequence
    lookup, and to derive spec-correct REF/ALT allele pairs anchored to the
    real reference sequence.  Use ``ReferenceBuilder`` / ``ReferenceSpec``
    instead when you need a fully synthetic reference without an on-disk FASTA.

    Attributes:
        _fa: The underlying ``pysam.FastaFile`` handle.
    """

    def __init__(self, fasta_path: StrPath):
        """Open a FASTA file for reference-aware allele generation.

        Args:
            fasta_path: Path to the FASTA (or bgzipped FASTA) file.  A
                companion ``.fai`` (and ``.gzi`` when bgzipped) index must
                exist alongside it.
        """
        self._fa = pysam.FastaFile(str(fasta_path))

    def base(self, contig: str, pos0: int) -> str:
        """Return the single reference base at 0-based ``pos0`` on ``contig``.

        Args:
            contig: Contig/chromosome name as it appears in the FASTA header.
            pos0: 0-based position of the base to fetch.

        Returns:
            A single uppercase nucleotide character (``"A"``, ``"C"``,
            ``"G"``, or ``"T"``).
        """
        return self._fa.fetch(contig, pos0, pos0 + 1).upper()

    def seq(self, contig: str, start0: int, length: int) -> str:
        """Return ``length`` reference bases from 0-based ``start0`` on ``contig``.

        Args:
            contig: Contig/chromosome name as it appears in the FASTA header.
            start0: 0-based start position of the subsequence.
            length: Number of bases to fetch.

        Returns:
            An uppercase nucleotide string of exactly ``length`` characters.
        """
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
        """Draw a REF allele matching the reference plus a realistic ALT allele.

        The result is deterministic given its arguments; Hypothesis supplies
        randomness by choosing those arguments via its strategies.

        Args:
            contig: Contig/chromosome name.
            pos0: 0-based position of the variant anchor base.
            klass: Variant class — one of ``"SNP"``, ``"MNP"``, ``"INS"``,
                ``"DEL"``, ``"DELINS"``, or ``"SPANNING_DEL"``.
            alt_index: Cyclic offset into ``"ACGT"`` used to choose the ALT
                base for SNPs and MNPs (default ``1``).
            del_len: Number of deleted bases (not counting the anchor) for
                ``"DEL"`` variants (default ``1``).
            ins_seq: Sequence inserted after the anchor for ``"INS"`` and
                ``"DELINS"`` variants (default ``"T"``).
            mnp_len: Length of the REF/ALT run for ``"MNP"`` and ``"DELINS"``
                variants (default ``2``).

        Returns:
            A ``(ref, alts)`` tuple where ``ref`` is the REF allele string and
            ``alts`` is a list containing the single ALT allele string.

        Raises:
            ValueError: ``klass`` is not one of the recognized variant classes.
        """
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
class RepeatFeature(CompactRepr):
    """A tandem repeat planted into a synthetic reference, for provenance.

    Records the location and composition of a repeat written by
    ``ReferenceBuilder.tandem_repeat`` so callers can generate
    tandem-repeat expansion/contraction variants that are anchored to a
    known repeat locus.

    Attributes:
        contig: Contig name where the repeat was planted.
        pos0: 0-based start position of the repeat run.
        motif: Repeated unit sequence (e.g. ``"CAG"``).
        count: Number of times ``motif`` is repeated.
    """

    contig: str
    pos0: int  # 0-based start of the repeat run
    motif: str
    count: int

    @property
    def length(self) -> int:
        """Total length of the repeat run in bases (``len(motif) * count``)."""
        return len(self.motif) * self.count

    @override
    def __repr__(self) -> str:
        return f"RepeatFeature({self.contig}@{self.pos0} {self.motif}×{self.count})"


@dataclass(frozen=True)
class ReferenceSpec(CompactRepr):
    """Immutable in-memory reference: contig sequences and planted repeats.

    Produced by ``ReferenceBuilder.build()``.  Provides the same
    ``base``/``seq``/``draw_ref_alt`` interface as ``Reference`` but operates
    entirely in memory, without requiring an on-disk FASTA file until
    ``write()`` is called.

    Attributes:
        contigs: Ordered ``(id, sequence)`` pairs for each contig.
        repeats: ``RepeatFeature`` records for every tandem repeat planted via
            ``ReferenceBuilder.tandem_repeat``.
    """

    contigs: tuple[tuple[str, str], ...]  # (id, sequence)
    repeats: tuple[RepeatFeature, ...] = ()

    @override
    def __repr__(self) -> str:
        contigs = ", ".join(f"{cid}:{len(seq)}bp" for cid, seq in self.contigs)
        return f"ReferenceSpec(contigs=[{contigs}], repeats={len(self.repeats)})"

    def _seq_for(self, contig: str) -> str:
        for cid, seq in self.contigs:
            if cid == contig:
                return seq
        raise KeyError(contig)

    def length(self, contig: str) -> int:
        """Return the total length in bases of ``contig``.

        Args:
            contig: Contig name; must match one of the IDs in ``contigs``.

        Returns:
            Length of the contig sequence in bases.

        Raises:
            KeyError: ``contig`` is not present in this reference.
        """
        return len(self._seq_for(contig))

    def base(self, contig: str, pos0: int) -> str:
        """Return the single reference base at 0-based ``pos0`` on ``contig``.

        Args:
            contig: Contig name; must match one of the IDs in ``contigs``.
            pos0: 0-based position of the base to fetch.

        Returns:
            A single uppercase nucleotide character.

        Raises:
            KeyError: ``contig`` is not present in this reference.
        """
        return self._seq_for(contig)[pos0]

    def seq(self, contig: str, start0: int, length: int) -> str:
        """Return ``length`` reference bases from 0-based ``start0`` on ``contig``.

        Args:
            contig: Contig name; must match one of the IDs in ``contigs``.
            start0: 0-based start position of the subsequence.
            length: Number of bases to fetch.

        Returns:
            An uppercase nucleotide string of exactly ``length`` characters.

        Raises:
            KeyError: ``contig`` is not present in this reference.
        """
        return self._seq_for(contig)[start0 : start0 + length]

    def write(self, path: StrPath, *, bgzip: bool = True, index: bool = True) -> Path:
        """Write a 60-column FASTA to ``path``; optionally bgzip and index it.

        Args:
            path: Destination path for the written file.
            bgzip: If ``True`` (default), compress the output with bgzip via
                ``pysam.tabix_compress``.  The path is used as-is; append
                ``".gz"`` yourself if you want the conventional extension.
            index: If ``True`` (the default), write a ``.fai`` index
                alongside the output file via ``pysam.faidx``; a ``.gzi``
                block index is additionally written when ``bgzip=True``.

        Returns:
            The path that was written (exactly as given; no extension is
            appended automatically).
        """
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
        """Draw a REF allele matching the reference plus a realistic ALT allele.

        The result is deterministic given its arguments; Hypothesis supplies
        randomness by choosing those arguments via its strategies.

        Args:
            contig: Contig name; must match one of the IDs in ``contigs``.
            pos0: 0-based position of the variant anchor base.
            klass: Variant class — one of ``"SNP"``, ``"MNP"``, ``"INS"``,
                ``"DEL"``, ``"DELINS"``, or ``"SPANNING_DEL"``.
            alt_index: Cyclic offset into ``"ACGT"`` used to choose the ALT
                base for SNPs and MNPs (default ``1``).
            del_len: Number of deleted bases (not counting the anchor) for
                ``"DEL"`` variants (default ``1``).
            ins_seq: Sequence inserted after the anchor for ``"INS"`` and
                ``"DELINS"`` variants (default ``"T"``).
            mnp_len: Length of the REF/ALT run for ``"MNP"`` and ``"DELINS"``
                variants (default ``2``).

        Returns:
            A ``(ref, alts)`` tuple where ``ref`` is the REF allele string and
            ``alts`` is a list containing the single ALT allele string.

        Raises:
            ValueError: ``klass`` is not one of the recognized variant classes.
            KeyError: ``contig`` is not present in this reference.
        """
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
    """Mutable builder for a fully synthetic reference genome.

    Random-fills each declared contig (seeded for reproducibility), supports
    single-base, multi-nucleotide, and tandem-repeat overwrites, then
    ``build()`` freezes the result into an immutable ``ReferenceSpec``.

    Use this class when you need a reference without an existing FASTA file.
    For an existing FASTA, use ``Reference`` instead.
    """

    def __init__(self, seed: int = 0):
        """Create a builder with a seeded random-number generator.

        Args:
            seed: Seed for ``numpy.random.default_rng``; controls the
                random bases used to fill newly added contigs (default ``0``).
        """
        self._rng = np.random.default_rng(seed)
        self._seqs: dict[str, np.ndarray] = {}
        self._order: list[str] = []
        self._repeats: list[RepeatFeature] = []

    def add_contig(self, id: str, length: int) -> ReferenceBuilder:
        """Add a new contig filled with random bases.

        The bases are drawn from ``"ACGT"`` using the builder's seeded RNG,
        so results are reproducible for the same ``seed`` and call order.

        Args:
            id: Contig identifier (e.g. ``"chr1"``).
            length: Length of the contig in bases.

        Returns:
            The builder, for chaining.

        Raises:
            ValueError: A contig with ``id`` has already been added.
        """
        if id in self._seqs:
            raise ValueError(f"contig {id!r} already added")
        self._seqs[id] = self._rng.choice(_BASES_ARR, size=length)
        self._order.append(id)
        return self

    def set_base(self, contig: str, pos0: int, base: str) -> ReferenceBuilder:
        """Overwrite a single base at 0-based ``pos0`` on ``contig``.

        Args:
            contig: Contig name; must have been added via ``add_contig``.
            pos0: 0-based position to overwrite.
            base: Replacement base — must be a single character.

        Returns:
            The builder, for chaining.

        Raises:
            ValueError: ``base`` is not exactly one character.
        """
        if len(base) != 1:
            raise ValueError(f"set_base expects one base, got {base!r}")
        self._seqs[contig][pos0] = base.encode()
        return self

    def set_seq(self, contig: str, pos0: int, seq: str) -> ReferenceBuilder:
        """Overwrite a run of bases starting at 0-based ``pos0`` on ``contig``.

        Args:
            contig: Contig name; must have been added via ``add_contig``.
            pos0: 0-based start position of the region to overwrite.
            seq: Replacement sequence; must fit within the contig boundaries.

        Returns:
            The builder, for chaining.

        Raises:
            ValueError: The range ``[pos0, pos0 + len(seq))`` extends past the
                end of the contig.
        """
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
        """Plant a tandem repeat and record it for downstream provenance.

        Writes ``motif * n`` into the contig starting at ``pos0`` and appends
        a ``RepeatFeature`` to the spec, making the repeat locus addressable
        for tandem-repeat variant generation.

        Args:
            contig: Contig name; must have been added via ``add_contig``.
            pos0: 0-based start position of the repeat run.
            motif: Repeated unit sequence (e.g. ``"CAG"``).
            n: Number of times ``motif`` is repeated.

        Returns:
            The builder, for chaining.

        Raises:
            ValueError: The repeat run extends past the end of the contig.
        """
        self.set_seq(contig, pos0, motif * n)
        self._repeats.append(RepeatFeature(contig, pos0, motif, n))
        return self

    def build(self) -> ReferenceSpec:
        """Freeze the builder into an immutable ``ReferenceSpec``.

        Returns:
            A ``ReferenceSpec`` containing all contigs (in the order they were
            added) and all ``RepeatFeature`` records planted so far.
        """
        contigs = tuple(
            (cid, self._seqs[cid].tobytes().decode()) for cid in self._order
        )
        return ReferenceSpec(contigs=contigs, repeats=tuple(self._repeats))
