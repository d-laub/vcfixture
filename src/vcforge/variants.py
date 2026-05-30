from __future__ import annotations

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

def record_class(ref: str, alts: tuple[str, ...]) -> str:
    if len(alts) > 1:
        return "MULTIALLELIC"
    return classify(ref, alts[0])
