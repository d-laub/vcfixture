from __future__ import annotations


def genotype_ordering(ploidy: int, n_alleles: int) -> list[tuple[int, ...]]:
    """Ordered genotypes per VCF 4.5 'Number=G' ordering."""
    if ploidy < 1:
        raise ValueError("ploidy must be >= 1")

    def rec(p: int) -> list[tuple[int, ...]]:
        if p == 1:
            return [(a,) for a in range(n_alleles)]
        out: list[tuple[int, ...]] = []
        for a in range(n_alleles):
            for prefix in rec(p - 1):
                if prefix[-1] <= a:
                    out.append(prefix + (a,))
        return out

    return rec(ploidy)
