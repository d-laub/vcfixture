from hypothesis.vendor.pretty import pretty

from vcfixture._spec.number import Number


def test_number_singletons_not_in_dataclass_fields():
    # ClassVar singletons must not leak into the dataclass field set, or any
    # field-walking pretty-printer (Hypothesis) recurses into them forever.
    assert set(Number.__dataclass_fields__) == {"kind", "count"}


def test_number_pretty_does_not_recurse():
    # Hypothesis's vendored pretty-printer is what produces "Falsifying example"
    # output. It must not explode Number into its singletons.
    out = pretty(Number.G)
    assert "ONE=Number" not in out
    assert len(out) < 200


def test_number_repr_compact():
    assert repr(Number.G) == "Number(G)"
    assert repr(Number.A) == "Number(A)"
    assert repr(Number.R) == "Number(R)"
    assert repr(Number.DOT) == "Number(.)"
    assert repr(Number.ONE) == "Number(1)"
    assert repr(Number.fixed(2)) == "Number(2)"
    assert repr(Number.FLAG) == "Number(FLAG)"


def test_number_pretty_uses_compact_repr():
    # _repr_pretty_ must route Hypothesis's printer through __repr__.
    assert pretty(Number.G) == "Number(G)"
