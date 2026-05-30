from __future__ import annotations
import re
from dataclasses import dataclass

_TOKEN_RE = re.compile(r"([|/])")

@dataclass(frozen=True)
class Genotype:
    alleles: tuple[int | None, ...]   # None == missing allele
    phased: tuple[bool, ...]          # one per separator; len == len(alleles) - 1

    @classmethod
    def parse(cls, s: str) -> "Genotype":
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
        return len(self.alleles)

    @property
    def is_phased(self) -> bool:
        return len(self.phased) > 0 and all(self.phased)

    def render(self) -> str:
        out = ["." if a is None else str(a) for a in self.alleles]
        seps = ["|" if p else "/" for p in self.phased]
        chars: list[str] = [out[0]]
        for sep, allele in zip(seps, out[1:]):
            chars.append(sep)
            chars.append(allele)
        return "".join(chars)
