import pytest
from vcforge._spec.number import Number, NumberKind

def test_header_strings():
    assert Number.ONE.header_str() == "1"
    assert Number.fixed(2).header_str() == "2"
    assert Number.A.header_str() == "A"
    assert Number.R.header_str() == "R"
    assert Number.G.header_str() == "G"
    assert Number.DOT.header_str() == "."
    assert Number.FLAG.header_str() == "0"

def test_cardinality_fixed_a_r():
    assert Number.fixed(3).cardinality(n_alt=5, ploidy=2) == 3
    assert Number.A.cardinality(n_alt=2, ploidy=2) == 2
    assert Number.R.cardinality(n_alt=2, ploidy=2) == 3
    assert Number.FLAG.cardinality(n_alt=2, ploidy=2) == 0

def test_cardinality_g_matches_spec_examples():
    assert Number.G.cardinality(n_alt=1, ploidy=2) == 3
    assert Number.G.cardinality(n_alt=2, ploidy=2) == 6
    assert Number.G.cardinality(n_alt=2, ploidy=3) == 10

def test_cardinality_dot_is_variable():
    assert Number.DOT.cardinality(n_alt=2, ploidy=2) is None

def test_equality_and_hash():
    assert Number.fixed(1) == Number.ONE
    assert {Number.A, Number.A} == {Number.A}
