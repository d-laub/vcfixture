from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, ClassVar

from .._repr import CompactRepr, override


class NumberKind(str, Enum):
    """Discriminant tag for ``Number`` cardinality variants.

    Each member corresponds to one of the VCF spec's ``Number=`` tokens:
    ``FIXED`` for a literal integer count, ``A``/``R``/``G`` for
    allele-relative counts, ``DOT`` for variable-length, and ``FLAG`` for
    the zero-value presence/absence sentinel.
    """

    FIXED = "FIXED"
    """An exact, non-negative integer count supplied at construction time."""

    A = "A"
    """One value per ALT allele (``Number=A`` in the VCF header)."""

    R = "R"
    """One value per allele including REF (``Number=R``; ALT count + 1)."""

    G = "G"
    """One value per possible genotype (``Number=G``; depends on ploidy)."""

    DOT = "."
    """Variable or unknown count (``Number=.``); resolved count is ``None``."""

    FLAG = "FLAG"
    """Presence/absence flag; always has count zero (``Number=0``)."""


@dataclass(frozen=True)
class Number(CompactRepr):
    """VCF ``Number=`` cardinality descriptor.

    ``Number`` encodes how many values a VCF INFO or FORMAT field carries
    for each record.  It is either a concrete count (``ONE``, ``FIXED(n)``)
    or an allele-relative / variable symbolic count (``A``, ``R``, ``G``,
    ``DOT``, ``FLAG``).

    Prefer the pre-built singletons for common cases:

    - ``Number.ONE``  — exactly one value (``Number=1``)
    - ``Number.A``    — one value per ALT allele
    - ``Number.R``    — one value per allele (REF + ALTs)
    - ``Number.G``    — one value per possible genotype
    - ``Number.DOT``  — variable / unknown count
    - ``Number.FLAG`` — zero-value presence flag

    Use ``Number.fixed(n)`` to construct arbitrary fixed counts.
    """

    kind: NumberKind
    count: int | None = None  # set only for FIXED

    # Canonical singletons, assigned after the class body. Declared under
    # TYPE_CHECKING so type checkers know these attributes exist on Number,
    # while keeping them out of __dataclass_fields__ at runtime (otherwise
    # field-walking pretty-printers recurse into them — see tests/test_repr.py).
    if TYPE_CHECKING:
        ONE: ClassVar[Number]
        A: ClassVar[Number]
        R: ClassVar[Number]
        G: ClassVar[Number]
        DOT: ClassVar[Number]
        FLAG: ClassVar[Number]

    @classmethod
    def fixed(cls, n: int) -> Number:
        """Construct a ``FIXED(n)`` Number with an explicit count.

        Args:
            n: The exact value count; must be non-negative.

        Returns:
            A ``Number`` whose cardinality is always ``n``.

        Raises:
            ValueError: If ``n`` is negative.
        """
        if n < 0:
            raise ValueError("fixed Number must be >= 0")
        return cls(NumberKind.FIXED, n)

    def header_str(self) -> str:
        """Return the VCF header token for this Number.

        Returns:
            The token written after ``Number=`` in a VCF meta-information
            line: the literal integer for ``FIXED`` counts (``"1"``,
            ``"2"``, …), ``"0"`` for ``FLAG``, and the symbolic letter
            (``"A"``, ``"R"``, ``"G"``, ``"."``) for the rest.
        """
        if self.kind is NumberKind.FIXED:
            return str(self.count)
        if self.kind is NumberKind.FLAG:
            return "0"
        return self.kind.value

    @override
    def __repr__(self) -> str:
        if self.kind is NumberKind.FIXED:
            tok = str(self.count)
        elif self.kind is NumberKind.FLAG:
            tok = "FLAG"
        else:
            tok = self.kind.value
        return f"Number({tok})"

    def cardinality(self, n_alt: int, ploidy: int) -> int | None:
        """Resolve this Number to a concrete value count for one record.

        Args:
            n_alt: Number of ALT alleles in the record.
            ploidy: Sample ploidy (used for ``Number=G``).

        Returns:
            The required value count, or ``None`` when the count is
            unbounded (``Number=.``) and cannot be resolved statically.
        """
        k = self.kind
        if k is NumberKind.FIXED:
            return self.count
        if k is NumberKind.FLAG:
            return 0
        if k is NumberKind.A:
            return n_alt
        if k is NumberKind.R:
            return n_alt + 1
        if k is NumberKind.G:
            n_alleles = n_alt + 1
            return math.comb(n_alleles + ploidy - 1, ploidy)
        return None  # DOT / variable


Number.ONE = Number(NumberKind.FIXED, 1)
Number.A = Number(NumberKind.A)
Number.R = Number(NumberKind.R)
Number.G = Number(NumberKind.G)
Number.DOT = Number(NumberKind.DOT)
Number.FLAG = Number(NumberKind.FLAG)
