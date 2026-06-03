from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from ._repr import CompactRepr, override
from ._spec.fielddef import FieldDef
from ._spec.version import VcfVersion
from ._typing import StrPath
from .allele import Allele
from .genotype import Genotype

if TYPE_CHECKING:
    from .truth import GroundTruth


@dataclass(frozen=True)
class ContigDef(CompactRepr):
    id: str
    length: int | None = None

    @override
    def __repr__(self) -> str:
        if self.length is None:
            return f"ContigDef({self.id})"
        return f"ContigDef({self.id}:{self.length})"

    def header_line(self) -> str:
        if self.length is None:
            return f"##contig=<ID={self.id}>"
        return f"##contig=<ID={self.id},length={self.length}>"


@dataclass(frozen=True)
class AltDef(CompactRepr):
    id: str
    description: str

    @override
    def __repr__(self) -> str:
        return f"AltDef({self.id})"

    def header_line(self) -> str:
        return f'##ALT=<ID={self.id},Description="{self.description}">'


@dataclass(frozen=True)
class Record(CompactRepr):
    chrom: str
    pos: int  # 1-based
    ids: tuple[str, ...] | None  # None -> "."
    ref: str
    alts: tuple[Allele, ...]  # typed; may include SpanningDeletion()/symbolic/...
    qual: float | None
    filters: tuple[str, ...] | None  # None -> "."; () -> "PASS"
    info: Mapping[str, Any]  # id -> value(s); Flag -> True
    fmt_keys: tuple[str, ...]  # FORMAT column order
    samples: tuple[Mapping[str, Any], ...]  # per-sample: key -> value(s)/Genotype
    labels: frozenset[str] = frozenset()

    @override
    def __repr__(self) -> str:
        alts = ",".join(a.render() for a in self.alts) if self.alts else "."
        out = f"Record({self.chrom}:{self.pos} {self.ref}>{alts} ×{len(self.samples)}"
        if self.labels:
            out += f" [{','.join(sorted(self.labels))}]"
        return out + ")"

    @property
    def n_alt(self) -> int:
        return len(self.alts)


@dataclass(frozen=True)
class VcfDocument(CompactRepr):
    version: VcfVersion
    info_defs: tuple[FieldDef, ...]
    format_defs: tuple[FieldDef, ...]
    filter_defs: tuple[tuple[str, str], ...]  # (id, description)
    contigs: tuple[ContigDef, ...]
    samples: tuple[str, ...]
    records: tuple[Record, ...]
    alt_defs: tuple[AltDef, ...] = ()

    @override
    def __repr__(self) -> str:
        return (
            f"VcfDocument({self.version.value} samples={len(self.samples)} "
            f"records={len(self.records)} info={len(self.info_defs)} "
            f"format={len(self.format_defs)})"
        )

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
