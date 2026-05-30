from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import ClassVar


class NumberKind(str, Enum):
    FIXED = "FIXED"
    A = "A"
    R = "R"
    G = "G"
    DOT = "."
    FLAG = "FLAG"


@dataclass(frozen=True)
class Number:
    kind: NumberKind
    count: int | None = None  # set only for FIXED

    # Canonical singletons, assigned after the class body. Declared as ClassVar
    # so type checkers know these attributes exist on Number.
    ONE: ClassVar[Number]
    A: ClassVar[Number]
    R: ClassVar[Number]
    G: ClassVar[Number]
    DOT: ClassVar[Number]
    FLAG: ClassVar[Number]

    @classmethod
    def fixed(cls, n: int) -> Number:
        if n < 0:
            raise ValueError("fixed Number must be >= 0")
        return cls(NumberKind.FIXED, n)

    def header_str(self) -> str:
        if self.kind is NumberKind.FIXED:
            return str(self.count)
        if self.kind is NumberKind.FLAG:
            return "0"
        return self.kind.value

    def cardinality(self, n_alt: int, ploidy: int) -> int | None:
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
