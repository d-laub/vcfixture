from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from ._spec.fielddef import FieldDef
from ._typing import StrPath
from .genotype import Genotype

if TYPE_CHECKING:
    from .truth import GroundTruth


@dataclass(frozen=True)
class ContigDef:
    id: str
    length: int | None = None

    def header_line(self) -> str:
        if self.length is None:
            return f"##contig=<ID={self.id}>"
        return f"##contig=<ID={self.id},length={self.length}>"


@dataclass(frozen=True)
class Record:
    chrom: str
    pos: int  # 1-based
    ids: tuple[str, ...] | None  # None -> "."
    ref: str
    alts: tuple[str, ...]  # may contain "*"
    qual: float | None
    filters: tuple[str, ...] | None  # None -> "."; () -> "PASS"
    info: Mapping[str, Any]  # id -> value(s); Flag -> True
    fmt_keys: tuple[str, ...]  # FORMAT column order
    samples: tuple[Mapping[str, Any], ...]  # per-sample: key -> value(s)/Genotype

    @property
    def n_alt(self) -> int:
        return len(self.alts)


@dataclass(frozen=True)
class VcfDocument:
    fileformat: str
    info_defs: tuple[FieldDef, ...]
    format_defs: tuple[FieldDef, ...]
    filter_defs: tuple[tuple[str, str], ...]  # (id, description)
    contigs: tuple[ContigDef, ...]
    samples: tuple[str, ...]
    records: tuple[Record, ...]

    def max_ploidy(self) -> int:
        p = 1
        for rec in self.records:
            for s in rec.samples:
                gt = s.get("GT")
                if isinstance(gt, Genotype):
                    p = max(p, gt.ploidy)
        return p

    def render(self) -> str:
        from .serialize import render_document

        return render_document(self)

    def truth(self) -> GroundTruth:
        from .truth import derive_truth

        return derive_truth(self)

    def write(self, path: StrPath, *, bgzip: bool = False, index: bool = False) -> Path:
        from . import io

        if bgzip:
            return io.write_bgzip(self, path, index=index)
        return io.write_text(self, path)
