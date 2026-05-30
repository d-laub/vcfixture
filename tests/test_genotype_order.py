from vcfixture._spec.genotype_order import genotype_ordering


def test_diploid_biallelic():
    assert genotype_ordering(ploidy=2, n_alleles=2) == [(0, 0), (0, 1), (1, 1)]


def test_diploid_triallelic():
    assert genotype_ordering(ploidy=2, n_alleles=3) == [
        (0, 0),
        (0, 1),
        (1, 1),
        (0, 2),
        (1, 2),
        (2, 2),
    ]


def test_triploid_triallelic_count():
    out = genotype_ordering(ploidy=3, n_alleles=3)
    assert len(out) == 10
    assert out[0] == (0, 0, 0)
    assert out[-1] == (2, 2, 2)


def test_haploid():
    assert genotype_ordering(ploidy=1, n_alleles=3) == [(0,), (1,), (2,)]
