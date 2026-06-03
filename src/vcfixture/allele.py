from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TypeAlias

from ._repr import CompactRepr, override

_SEQ_RE = re.compile(r"^[ACGTNacgtn]+$")
_SV_FIRST_TYPES = frozenset({"DEL", "INS", "DUP", "INV", "CNV"})
# t[p[ , t]p] , [p[t , ]p]t  — both brackets identical (\2 backref); mate is chr:pos
_BND_PAIRED_RE = re.compile(r"^([ACGTNacgtn]*)([\[\]])([^\[\]]+:\d+)\2([ACGTNacgtn]*)$")
_BND_SINGLE_RE = re.compile(r"^(?:\.[ACGTNacgtn]+|[ACGTNacgtn]+\.)$")


@dataclass(frozen=True)
class SequenceAllele(CompactRepr):
    """A literal nucleotide ALT (or REF) allele, e.g. ``A``, ``GAT``.

    Attributes:
        bases: The allele sequence; must match ``[ACGTNacgtn]+``.
    """

    bases: str

    def __post_init__(self) -> None:
        if not _SEQ_RE.match(self.bases):
            raise ValueError(
                f"sequence allele bases must be [ACGTN]+, got {self.bases!r}"
            )

    def render(self) -> str:
        return self.bases

    @override
    def __repr__(self) -> str:
        return f"Seq({self.bases})"


@dataclass(frozen=True)
class SpanningDeletion(CompactRepr):
    """The spanning-deletion ALT ``*`` (alias ``Star``)."""

    def render(self) -> str:
        return "*"

    @override
    def __repr__(self) -> str:
        return "Star()"


@dataclass(frozen=True)
class SymbolicAllele(CompactRepr):
    """A symbolic structural-variant ALT, e.g. ``<DEL>`` or ``<DUP:TANDEM>``
    (alias ``Sym``).

    Attributes:
        first_type: SV type; must be one of ``DEL``, ``INS``, ``DUP``, ``INV``, ``CNV``.
        subtypes: Optional subtype tokens joined by ``:`` when rendered.
    """

    first_type: str
    subtypes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.first_type not in _SV_FIRST_TYPES:
            raise ValueError(
                f"symbolic SV first type must be one of {sorted(_SV_FIRST_TYPES)}, "
                f"got {self.first_type!r}"
            )

    @property
    def type_str(self) -> str:
        """The inner ``<...>`` token, e.g. ``DEL`` or ``DUP:TANDEM``."""
        return ":".join((self.first_type, *self.subtypes))

    def render(self) -> str:
        return f"<{self.type_str}>"

    @override
    def __repr__(self) -> str:
        return f"Sym(<{self.type_str}>)"

    @classmethod
    def deletion(cls, *subtypes: str) -> SymbolicAllele:
        """Construct a ``<DEL[:subtypes...]>`` symbolic allele."""
        return cls("DEL", subtypes)

    @classmethod
    def insertion(cls, *subtypes: str) -> SymbolicAllele:
        """Construct an ``<INS[:subtypes...]>`` symbolic allele."""
        return cls("INS", subtypes)

    @classmethod
    def duplication(cls, *subtypes: str) -> SymbolicAllele:
        """Construct a ``<DUP[:subtypes...]>`` symbolic allele."""
        return cls("DUP", subtypes)

    @classmethod
    def inversion(cls, *subtypes: str) -> SymbolicAllele:
        """Construct an ``<INV[:subtypes...]>`` symbolic allele."""
        return cls("INV", subtypes)

    @classmethod
    def cnv(cls, *subtypes: str) -> SymbolicAllele:
        """Construct a ``<CNV[:subtypes...]>`` symbolic allele."""
        return cls("CNV", subtypes)


@dataclass(frozen=True)
class UnspecifiedAllele(CompactRepr):
    """The unspecified/non-ref ALT ``<*>`` (alias ``Unspecified``)."""

    def render(self) -> str:
        return "<*>"

    @override
    def __repr__(self) -> str:
        return "Unspecified()"


@dataclass(frozen=True)
class BreakendAllele(CompactRepr):
    """A breakend ALT replacement string, e.g. ``G[chr2:321[`` (alias ``Bnd``).

    Attributes:
        raw: The verbatim ALT text.
        single: ``True`` for single-breakend forms (``.t`` / ``t.``);
            ``False`` for paired forms (``t[p[``, ``t]p]``, ``[p[t``, ``]p]t``).
    """

    raw: str
    single: bool = False

    def render(self) -> str:
        return self.raw

    @override
    def __repr__(self) -> str:
        return f"Bnd({self.raw})"

    @classmethod
    def parse(cls, s: str) -> BreakendAllele:
        """Parse a breakend replacement string into a ``BreakendAllele``.

        Args:
            s: The ALT text — paired forms (``t[p[``, ``t]p]``, ``[p[t``,
                ``]p]t``) or single-breakend forms (``.t`` / ``t.``).

        Returns:
            The parsed breakend, with ``single=True`` for single-breakend forms.

        Raises:
            ValueError: ``s`` is not a valid breakend replacement string.
        """
        if _BND_SINGLE_RE.match(s):
            return cls(raw=s, single=True)
        if _BND_PAIRED_RE.match(s):
            return cls(raw=s, single=False)
        raise ValueError(f"not a valid breakend replacement string: {s!r}")


# Union of all five allele types.  Anywhere vcfixture accepts or returns an
# allele value, this is the annotated type.
Allele: TypeAlias = (
    SequenceAllele
    | SpanningDeletion
    | SymbolicAllele
    | UnspecifiedAllele
    | BreakendAllele
)

# Ergonomic aliases for terse fixtures: Seq / Sym / Star / Unspecified / Bnd.
Seq = SequenceAllele
Sym = SymbolicAllele
Star = SpanningDeletion
Unspecified = UnspecifiedAllele
Bnd = BreakendAllele


def classify_allele(alt: str) -> Allele:
    """Parse a raw ALT string into a typed Allele (syntactic dispatch)."""
    if alt == "*":
        return SpanningDeletion()
    if alt == "<*>":
        return UnspecifiedAllele()
    if alt.startswith("<") and alt.endswith(">"):
        parts = alt[1:-1].split(":")
        return SymbolicAllele(parts[0], tuple(parts[1:]))
    if "[" in alt or "]" in alt:
        return BreakendAllele.parse(alt)
    if len(alt) > 1 and (alt.startswith(".") or alt.endswith(".")):
        return BreakendAllele.parse(alt)
    return SequenceAllele(alt)
