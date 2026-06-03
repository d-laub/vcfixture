from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from ._spec.fielddef import FieldDef
from ._spec.number import Number, NumberKind
from ._spec.reserved import reserved
from ._spec.types import Type
from ._spec.version import LATEST, VcfVersion
from ._typing import StrPath
from .allele import (
    Allele,
    BreakendAllele,
    SpanningDeletion,
    SymbolicAllele,
    UnspecifiedAllele,
)
from .genotype import Genotype
from .model import AltDef, ContigDef, Record, VcfDocument

if TYPE_CHECKING:
    from .truth import GroundTruth

_SVCLAIM_RULES: dict[str, frozenset[str]] = {
    "DEL": frozenset({"D", "J", "DJ"}),
    "DUP": frozenset({"D", "J", "DJ"}),
    "CNV": frozenset({"D"}),
    "INS": frozenset({"J"}),
    "INV": frozenset({"J"}),
}
_SVCLAIM_REQUIRED = frozenset({"DEL", "DUP"})
_CN_SVLEN_TYPES = frozenset({"CNV", "DEL", "DUP"})


def _per_allele(value: object, i: int) -> object:
    """Resolve the i-th per-allele entry of a Number=A info value.

    A scalar is treated as the single (index-0) value; a short/missing list
    yields None.
    """
    if isinstance(value, (list, tuple)):
        return value[i] if i < len(value) else None  # type: ignore[arg-type]
    return value if i == 0 else None


class VcfBuilder:
    def __init__(
        self,
        samples: Iterable[str],
        contigs: Iterable[tuple[str, int | None]],
        version: VcfVersion = LATEST,
    ):
        self._samples = tuple(samples)
        self._contigs = tuple(ContigDef(c[0], c[1]) for c in contigs)
        self._version = version
        self._info_defs: dict[str, FieldDef] = {}
        self._format_defs: dict[str, FieldDef] = {}
        self._filter_defs: list[tuple[str, str]] = []
        self._alt_defs: dict[str, str] = {}
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

    def alt(self, id: str, description: str) -> VcfBuilder:
        self._alt_defs[id] = description
        return self

    def _make_def(
        self,
        id: str,
        number: Number | None,
        type: Type | None,
        description: str | None,
        kind: str,
    ) -> FieldDef:
        if number is None or type is None:
            try:
                return reserved(id, kind, self._version)
            except KeyError:
                raise ValueError(
                    f"{kind} field {id!r} is not a known reserved field; "
                    f"pass number= and type= to declare it explicitly"
                ) from None
        return FieldDef(id, number, type, description or id, kind)

    def _validate_alleles(
        self,
        ref: str,
        alts: tuple[Allele, ...],
        info: Mapping[str, object] | None,
    ) -> None:
        info = info or {}
        svlen = info.get("SVLEN")
        svclaim = info.get("SVCLAIM")
        needs_padding = any(
            isinstance(a, (SymbolicAllele, BreakendAllele)) for a in alts
        )
        if needs_padding and len(ref) != 1:
            raise ValueError(
                "symbolic/breakend ALT requires a single preceding REF padding base, "
                f"got REF={ref!r}"
            )
        for i, a in enumerate(alts):
            sv = _per_allele(svlen, i)
            cl = _per_allele(svclaim, i)
            if isinstance(a, SymbolicAllele):
                if sv is None:
                    raise ValueError(f"SVLEN required for symbolic allele {a.render()}")
                allowed = _SVCLAIM_RULES[a.first_type]
                if cl is not None and cl not in allowed:
                    raise ValueError(
                        f"SVCLAIM {cl!r} invalid for {a.render()}; "
                        f"allowed {sorted(allowed)}"
                    )
                if (
                    self._version >= VcfVersion.V4_4
                    and a.first_type in _SVCLAIM_REQUIRED
                    and cl is None
                ):
                    raise ValueError(f"SVCLAIM required for {a.render()} (D/J/DJ)")
            elif isinstance(a, (BreakendAllele, UnspecifiedAllele, SpanningDeletion)):
                if sv is not None:
                    raise ValueError(f"SVLEN must be missing for {a.render()}")

    def record(
        self,
        chrom: str,
        pos: int,
        *,
        ref: str,
        alt: Sequence[Allele],
        ids: Iterable[str] | None = None,
        qual: float | None = None,
        filter: Iterable[str] | None = None,
        gt: Sequence[str] | None = None,
        info: Mapping[str, object] | None = None,
        labels: Iterable[str] | None = None,
        **fmt_fields: Sequence[object],
    ) -> VcfBuilder:
        alts = tuple(alt)
        n_alt = len(alts)
        self._validate_alleles(ref, alts, info)

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

        if "CN" in fmt_keys:
            svlen_val = (info or {}).get("SVLEN")
            cn_svlens = {
                _per_allele(svlen_val, i)
                for i, a in enumerate(alts)
                if isinstance(a, SymbolicAllele) and a.first_type in _CN_SVLEN_TYPES
            }
            if len(cn_svlens) > 1:
                raise ValueError(
                    "FORMAT CN requires equal SVLEN across <CNV>/<DEL>/<DUP> alleles"
                )

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
                labels=frozenset(labels) if labels else frozenset(),
            )
        )
        return self

    def build(self) -> VcfDocument:
        # Auto-describe each symbolic ALT type; explicit .alt() overrides then win.
        alt_ids: dict[str, str] = {}
        for rec in self._records:
            for a in rec.alts:
                if isinstance(a, SymbolicAllele):
                    alt_ids.setdefault(a.type_str, f"{a.type_str} structural variant")
        alt_ids.update(self._alt_defs)
        alt_defs = tuple(AltDef(i, d) for i, d in alt_ids.items())
        return VcfDocument(
            version=self._version,
            info_defs=tuple(self._info_defs.values()),
            format_defs=tuple(self._format_defs.values()),
            filter_defs=tuple(self._filter_defs),
            contigs=self._contigs,
            samples=self._samples,
            records=tuple(self._records),
            alt_defs=alt_defs,
        )

    def render(self) -> str:
        return self.build().render()

    def truth(self) -> GroundTruth:
        return self.build().truth()

    def write(self, path: StrPath, *, bgzip: bool = False, index: bool = False) -> Path:
        return self.build().write(path, bgzip=bgzip, index=index)
