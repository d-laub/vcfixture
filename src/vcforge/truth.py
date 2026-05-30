from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .genotype import Genotype
from .model import VcfDocument
from .variants import record_class

@dataclass(frozen=True)
class GroundTruth:
    samples: tuple[str, ...]
    contigs: tuple[str, ...]
    pos: np.ndarray                 # (records,) int64, 1-based
    ref: list[str]
    alts: list[list[str]]
    variant_class: list[str]
    genotypes: np.ndarray           # (records, samples, ploidy) int32, -1 missing
    phasing: np.ndarray             # (records, samples) bool (fully phased)
    info: list[dict]                # per record: id -> decoded value(s)
    format: list[list[dict]]        # per record, per sample: id -> value(s)

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
    info: list[dict] = []
    fmt: list[list[dict]] = []

    for ri, rec in enumerate(doc.records):
        pos[ri] = rec.pos
        ref.append(rec.ref)
        alts.append(list(rec.alts))
        vclass.append(record_class(rec.ref, rec.alts))
        info.append(dict(rec.info))
        per_sample: list[dict] = []
        for si, sample in enumerate(rec.samples):
            gt = sample.get("GT")
            if isinstance(gt, Genotype):
                for ai, allele in enumerate(gt.alleles):
                    genos[ri, si, ai] = -1 if allele is None else allele
                phasing[ri, si] = gt.is_phased
            per_sample.append({k: v for k, v in sample.items() if k != "GT"})
        fmt.append(per_sample)

    return GroundTruth(
        samples=doc.samples,
        contigs=tuple(c.id for c in doc.contigs),
        pos=pos, ref=ref, alts=alts, variant_class=vclass,
        genotypes=genos, phasing=phasing, info=info, format=fmt,
    )
