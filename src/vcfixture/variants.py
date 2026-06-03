from __future__ import annotations

from typing_extensions import assert_never

from .allele import (
    Allele,
    BreakendAllele,
    SequenceAllele,
    SpanningDeletion,
    SymbolicAllele,
    UnspecifiedAllele,
)


def snp(ref_base: str, alt_base: str) -> tuple[str, str]:
    return ref_base, alt_base


def mnp(ref_seq: str, alt_seq: str) -> tuple[str, str]:
    return ref_seq, alt_seq


def insertion(anchor: str, inserted: str) -> tuple[str, str]:
    return anchor, anchor + inserted


def deletion(anchor: str, deleted: str) -> tuple[str, str]:
    return anchor + deleted, anchor


def delins(ref_seq: str, alt_seq: str) -> tuple[str, str]:
    return ref_seq, alt_seq


def spanning_deletion(ref_base: str) -> tuple[str, str]:
    return ref_base, "*"


def classify(ref: str, alt: str) -> str:
    if alt == "*":
        return "SPANNING_DEL"
    lr, la = len(ref), len(alt)
    if lr == 1 and la == 1:
        return "SNP"
    if lr == la:
        return "MNP"
    if la > lr and alt.startswith(ref):
        return "INS"
    if lr > la and ref.startswith(alt):
        return "DEL"
    return "DELINS"


def record_class(ref: str, alts: tuple[Allele, ...]) -> str:
    # Covers len>1 (true multiallelic) and len==0 (monomorphic/REF-only, out of
    # v1 scope); single-ALT records fall through to per-allele dispatch below.
    if len(alts) != 1:
        return "MULTIALLELIC"
    a = alts[0]
    if isinstance(a, SequenceAllele):
        return classify(ref, a.bases)
    if isinstance(a, SpanningDeletion):
        return "SPANNING_DEL"
    if isinstance(a, UnspecifiedAllele):
        return "UNSPECIFIED"
    if isinstance(a, BreakendAllele):
        return "BND"
    if isinstance(a, SymbolicAllele):
        return a.first_type if a.first_type == "CNV" else f"SV_{a.first_type}"
    assert_never(a)
