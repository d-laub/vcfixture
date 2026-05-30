from vcfixture.genotype import Genotype


def test_parse_phased_diploid():
    g = Genotype.parse("0|1")
    assert g.alleles == (0, 1)
    assert g.phased == (True,)
    assert g.is_phased is True
    assert g.ploidy == 2


def test_parse_unphased_with_missing():
    g = Genotype.parse("./1")
    assert g.alleles == (None, 1)
    assert g.phased == (False,)
    assert g.is_phased is False


def test_parse_haploid():
    g = Genotype.parse("1")
    assert g.alleles == (1,)
    assert g.phased == ()


def test_round_trip_render():
    for s in ["0|1", "1/1", "./.", "0/0/1", "1|0|2"]:
        assert Genotype.parse(s).render() == s


def test_render_missing_genotype():
    assert Genotype.parse("./.").render() == "./."
