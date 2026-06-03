from __future__ import annotations

import re
from dataclasses import dataclass

from ._repr import CompactRepr, override

_TOKEN_RE = re.compile(r"([|/])")


@dataclass(frozen=True)
class Genotype(CompactRepr):
    """An immutable VCF genotype call (the ``GT`` field).

    A ``Genotype`` stores the allele indices and the phasing separator
    between each consecutive pair.  Missing alleles are represented as
    ``None``.

    Attributes:
        alleles: Allele indices in call order; ``None`` denotes a missing
            allele (``.`` in VCF text).
        phased: One boolean per separator between consecutive alleles;
            ``True`` means the separator was ``|`` (phased), ``False``
            means ``/`` (unphased).  Length is ``ploidy - 1``.
    """

    alleles: tuple[int | None, ...]  # None == missing allele
    phased: tuple[bool, ...]  # one per separator; len == len(alleles) - 1

    @classmethod
    def parse(cls, s: str) -> Genotype:
        """Parse a VCF genotype string into a Genotype.

        Args:
            s: The GT field text (e.g. ``"0|1"``, ``"./."``).  ``|``
                separates phased alleles, ``/`` unphased, and ``.``
                denotes a missing allele.

        Returns:
            The parsed genotype.
        """
        parts = _TOKEN_RE.split(s)  # "0|1" -> ["0","|","1"]
        alleles: list[int | None] = []
        phased: list[bool] = []
        for i, tok in enumerate(parts):
            if i % 2 == 0:
                alleles.append(None if tok == "." else int(tok))
            else:
                phased.append(tok == "|")
        return cls(tuple(alleles), tuple(phased))

    @property
    def ploidy(self) -> int:
        """Number of alleles in this genotype call."""
        return len(self.alleles)

    @property
    def is_phased(self) -> bool:
        """Return ``True`` if all separators in this genotype are phased (``|``).

        A haploid genotype with no separators returns ``False``.
        """
        return len(self.phased) > 0 and all(self.phased)

    def render(self) -> str:
        """Render this genotype back to VCF GT field text.

        Returns:
            The genotype string (e.g. ``"0|1"``, ``"./."``), with ``|``
            for phased separators, ``/`` for unphased, and ``.`` for
            missing alleles.
        """
        out = ["." if a is None else str(a) for a in self.alleles]
        seps = ["|" if p else "/" for p in self.phased]
        chars: list[str] = [out[0]]
        for sep, allele in zip(seps, out[1:], strict=True):
            chars.append(sep)
            chars.append(allele)
        return "".join(chars)

    @override
    def __repr__(self) -> str:
        return f"Genotype({self.render()})"
