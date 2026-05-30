from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from ._spec.fielddef import FieldDef
from ._spec.number import Number, NumberKind
from ._spec.reserved import reserved
from ._spec.types import Type
from ._typing import StrPath
from .genotype import Genotype
from .model import ContigDef, Record, VcfDocument

if TYPE_CHECKING:
    from .truth import GroundTruth


class VcfBuilder:
    def __init__(
        self,
        samples: Iterable[str],
        contigs: Iterable[tuple[str, int | None]],
        fileformat: str = "VCFv4.5",
    ):
        self._samples = tuple(samples)
        self._contigs = tuple(ContigDef(c[0], c[1]) for c in contigs)
        self._fileformat = fileformat
        self._info_defs: dict[str, FieldDef] = {}
        self._format_defs: dict[str, FieldDef] = {}
        self._filter_defs: list[tuple[str, str]] = []
        self._records: list[Record] = []

    def info(
        self,
        id: str,
        number: Number | None = None,
        type: Type | None = None,
        description: str | None = None,
    ) -> VcfBuilder:
        self._info_defs[id] = self._make_def(id, number, type, description, "INFO")
        return self

    def fmt(
        self,
        id: str,
        number: Number | None = None,
        type: Type | None = None,
        description: str | None = None,
    ) -> VcfBuilder:
        self._format_defs[id] = self._make_def(id, number, type, description, "FORMAT")
        return self

    def filter(self, id: str, description: str) -> VcfBuilder:
        self._filter_defs.append((id, description))
        return self

    @staticmethod
    def _make_def(
        id: str,
        number: Number | None,
        type: Type | None,
        description: str | None,
        kind: str,
    ) -> FieldDef:
        if number is None or type is None:
            try:
                return reserved(id, kind)
            except KeyError:
                raise ValueError(
                    f"{kind} field {id!r} is not a known reserved field; "
                    f"pass number= and type= to declare it explicitly"
                ) from None
        return FieldDef(id, number, type, description or id, kind)

    def record(
        self,
        chrom: str,
        pos: int,
        *,
        ref: str,
        alt: Sequence[str],
        ids: Iterable[str] | None = None,
        qual: float | None = None,
        filter: Iterable[str] | None = None,
        gt: Sequence[str] | None = None,
        info: Mapping[str, object] | None = None,
        **fmt_fields: Sequence[object],
    ) -> VcfBuilder:
        alts = tuple(alt)
        n_alt = len(alts)

        fmt_keys: list[str] = []
        samples: list[dict[str, Any]] = []
        for _ in self._samples:
            samples.append({})

        if gt is not None:
            if "GT" not in self._format_defs:
                raise ValueError("GT not declared; call .fmt('GT')")
            fmt_keys.append("GT")
            for si, s in enumerate(gt):
                geno = Genotype.parse(s)
                for a in geno.alleles:
                    if a is not None and a > n_alt:
                        raise ValueError(
                            f"allele index {a} out of range (n_alt={n_alt})"
                        )
                samples[si]["GT"] = geno

        ploidy = max((s["GT"].ploidy for s in samples if "GT" in s), default=2)

        for key, per_sample in fmt_fields.items():
            if key not in self._format_defs:
                raise ValueError(f"FORMAT field {key!r} not declared")
            fdef = self._format_defs[key]
            fmt_keys.append(key)
            card = fdef.number.cardinality(n_alt, ploidy)
            for si, val in enumerate(per_sample):
                if (
                    card is not None
                    and isinstance(val, (list, tuple))
                    and len(val) != card
                ):
                    raise ValueError(
                        f"{key} cardinality mismatch: expected {card}, got {len(val)}"
                    )
                samples[si][key] = val

        info_dict: dict[str, Any] = {}
        if info:
            for key, val in info.items():
                if key not in self._info_defs:
                    raise ValueError(f"INFO field {key!r} not declared")
                fdef = self._info_defs[key]
                card = fdef.number.cardinality(n_alt, ploidy)
                if (
                    card is not None
                    and fdef.number.kind is not NumberKind.FLAG
                    and isinstance(val, (list, tuple))
                    and len(val) != card
                ):
                    raise ValueError(
                        f"{key} cardinality mismatch: expected {card}, got {len(val)}"
                    )
                info_dict[key] = val

        self._records.append(
            Record(
                chrom=chrom,
                pos=pos,
                ids=tuple(ids) if ids else None,
                ref=ref,
                alts=alts,
                qual=qual,
                filters=tuple(filter) if filter is not None else None,
                info=info_dict,
                fmt_keys=tuple(fmt_keys),
                samples=tuple(samples),
            )
        )
        return self

    def build(self) -> VcfDocument:
        return VcfDocument(
            fileformat=self._fileformat,
            info_defs=tuple(self._info_defs.values()),
            format_defs=tuple(self._format_defs.values()),
            filter_defs=tuple(self._filter_defs),
            contigs=self._contigs,
            samples=self._samples,
            records=tuple(self._records),
        )

    def render(self) -> str:
        return self.build().render()

    def truth(self) -> GroundTruth:
        return self.build().truth()

    def write(self, path: StrPath, *, bgzip: bool = False, index: bool = False) -> Path:
        return self.build().write(path, bgzip=bgzip, index=index)
