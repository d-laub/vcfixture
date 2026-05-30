from hypothesis import given, settings
from hypothesis import strategies as st

from vcfixture import strategies as S
from vcfixture._spec.fielddef import FieldDef
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type

_SAFE = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")


@settings(max_examples=50)
@given(st.data())
def test_field_value_cardinality_and_type(data):
    fd = FieldDef("XR", Number.R, Type.FLOAT, "x", "INFO")
    v = data.draw(S.field_value(fd, n_alt=1, ploidy=2))
    assert isinstance(v, list) and len(v) == 2
    import numpy as np

    for x in v:
        assert isinstance(x, float)
        assert np.float32(x) == x


def test_field_value_flag_is_true():
    fd = FieldDef("XF", Number.FLAG, Type.FLAG, "x", "INFO")
    import hypothesis

    @hypothesis.given(hypothesis.strategies.data())
    def inner(data):
        assert data.draw(S.field_value(fd, n_alt=2, ploidy=2)) is True

    inner()


@settings(max_examples=30)
@given(st.data())
def test_field_value_G_count_multiallelic(data):
    fd = FieldDef("PL", Number.G, Type.INTEGER, "x", "FORMAT")
    v = data.draw(S.field_value(fd, n_alt=2, ploidy=2))
    assert len(v) == 6
    assert all(isinstance(x, int) for x in v)


@settings(max_examples=30)
@given(st.data())
def test_field_value_string_is_safe_alphabet(data):
    fd = FieldDef("XS", Number.ONE, Type.STRING, "x", "INFO")
    v = data.draw(S.field_value(fd, n_alt=1, ploidy=2))
    assert len(v) == 1
    assert set(v[0]) <= _SAFE and len(v[0]) >= 1
