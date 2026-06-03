from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from typing_extensions import assert_never

from .allele import (
    Allele,
    BreakendAllele,
    SequenceAllele,
    SpanningDeletion,
    SymbolicAllele,
    UnspecifiedAllele,
)
from .genotype import Genotype
from .model import VcfDocument
from .variants import classify, record_class


@dataclass(frozen=True)
class AlleleTruth:
    kind: str  # SNP|MNP|INS|DEL|DELINS|SPANNING_DEL|SYMBOLIC|UNSPECIFIED|BND
    is_sequence: bool  # True iff literal DNA a tool may splice
    sv_type: str | None  # e.g. "DEL"/"DUP:TANDEM" for symbolic; else None
    svlen: int | None  # resolved per-allele (absolute); None where undefined
    sv_end: int | None  # 1-based inclusive end = POS + svlen for DEL/DUP/INV/CNV


# Symbolic SV types that have a reference span (=> computed sv_end). Deliberately a
# subset of allele._SV_FIRST_TYPES, excluding INS (insertion has length but no span).
_SV_SPANNING = frozenset({"DEL", "DUP", "INV", "CNV"})


def _allele_truth(ref: str, pos: int, allele: Allele, svlen_val: object) -> AlleleTruth:
    if isinstance(allele, SequenceAllele):
        return AlleleTruth(classify(ref, allele.bases), True, None, None, None)
    if isinstance(allele, SpanningDeletion):
        return AlleleTruth("SPANNING_DEL", False, None, None, None)
    if isinstance(allele, UnspecifiedAllele):
        return AlleleTruth("UNSPECIFIED", False, None, None, None)
    if isinstance(allele, BreakendAllele):
        return AlleleTruth("BND", False, None, None, None)
    if isinstance(allele, SymbolicAllele):
        svlen = abs(int(svlen_val)) if isinstance(svlen_val, (int, float)) else None
        end = (
            pos + svlen
            if svlen is not None and allele.first_type in _SV_SPANNING
            else None
        )
        return AlleleTruth("SYMBOLIC", False, allele.type_str, svlen, end)
    assert_never(allele)


@dataclass(frozen=True)
class GroundTruth:
    """Decoded ground truth derived from a ``VcfDocument``.

    Produced by ``VcfDocument.truth()`` (or ``derive_truth``).  Every field is
    parallel to the record list of the originating document, so
    ``pos[i]``, ``ref[i]``, ``genotypes[i]``, etc. all describe the same record.

    Attributes:
        samples: Sample names in column order.
        contigs: Contig id for each record (parallel to ``pos``).
        pos: ``(records,)`` int64 array of 1-based positions.
        ref: REF allele string per record.
        alts: ALT allele strings per record (rendered to text).
        variant_class: Variant-class label per record
            (e.g. ``"SNP"``, ``"MNP"``, ``"INS"``, ``"DEL"``, ``"DELINS"``,
            ``"SPANNING_DEL"``).
        genotypes: ``(records, samples, ploidy)`` int32 array; ``-1`` encodes a
            missing allele call.
        phasing: ``(records, samples)`` bool array; ``True`` when the genotype is
            fully phased (``|``-separated in the VCF).
        info: Per-record mapping from INFO field id to decoded value(s).
        format: Per-record, per-sample mapping from FORMAT field id to
            decoded value(s).  GT is excluded (use ``genotypes``/``phasing``).
        labels: Per-record set of FILTER labels.
        alts_truth: Per-record, per-ALT ``AlleleTruth`` objects carrying
            allele kind, sequence flag, SV type, SVLEN, and END.
        is_sequence_mask: Per-record bool array (length = number of ALTs)
            indicating which ALTs are literal DNA sequences a tool may splice.
    """

    samples: tuple[str, ...]
    contigs: tuple[str, ...]
    pos: np.ndarray  # (records,) int64, 1-based
    ref: list[str]
    alts: list[list[str]]
    variant_class: list[str]
    genotypes: np.ndarray  # (records, samples, ploidy) int32, -1 missing
    phasing: np.ndarray  # (records, samples) bool (fully phased)
    info: list[dict[str, object]]  # per record: id -> decoded value(s)
    format: list[list[dict[str, object]]]  # per record, per sample: id -> value(s)
    labels: list[frozenset[str]]  # per record
    alts_truth: list[list[AlleleTruth]]  # per record, per ALT
    is_sequence_mask: list[np.ndarray]  # per record: bool array over ALTs


def derive_truth(doc: VcfDocument) -> GroundTruth:
    n_rec = len(doc.records)
    n_smp = len(doc.samples)
    ploidy = doc.max_ploidy()

    genos = np.full((n_rec, n_smp, ploidy), -1, dtype=np.int32)
    phasing = np.zeros((n_rec, n_smp), dtype=bool)
    pos = np.zeros(n_rec, dtype=np.int64)
    ref: list[str] = []
    alts: list[list[str]] = []
    vclass: list[str] = []
    info: list[dict[str, object]] = []
    fmt: list[list[dict[str, object]]] = []
    labels: list[frozenset[str]] = []
    alts_truth: list[list[AlleleTruth]] = []
    seq_mask: list[np.ndarray] = []

    for ri, rec in enumerate(doc.records):
        pos[ri] = rec.pos
        ref.append(rec.ref)
        alts.append([a.render() for a in rec.alts])
        vclass.append(record_class(rec.ref, rec.alts))
        info.append(dict(rec.info))
        per_sample: list[dict[str, object]] = []
        for si, sample in enumerate(rec.samples):
            gt = sample.get("GT")
            if isinstance(gt, Genotype):
                for ai, allele in enumerate(gt.alleles):
                    genos[ri, si, ai] = -1 if allele is None else allele
                phasing[ri, si] = gt.is_phased
            per_sample.append({k: v for k, v in sample.items() if k != "GT"})
        fmt.append(per_sample)
        labels.append(rec.labels)
        svlen_list = rec.info.get("SVLEN")
        if isinstance(svlen_list, (int, float)):
            svlen_list = [svlen_list]
        per_alt: list[AlleleTruth] = []
        for ai, allele in enumerate(rec.alts):
            sv = (
                svlen_list[ai]  # type: ignore[index]
                if isinstance(svlen_list, (list, tuple)) and ai < len(svlen_list)
                else None
            )
            per_alt.append(_allele_truth(rec.ref, rec.pos, allele, sv))
        alts_truth.append(per_alt)
        seq_mask.append(np.array([a.is_sequence for a in per_alt], dtype=bool))

    return GroundTruth(
        samples=doc.samples,
        contigs=tuple(c.id for c in doc.contigs),
        pos=pos,
        ref=ref,
        alts=alts,
        variant_class=vclass,
        genotypes=genos,
        phasing=phasing,
        info=info,
        format=fmt,
        labels=labels,
        alts_truth=alts_truth,
        is_sequence_mask=seq_mask,
    )
