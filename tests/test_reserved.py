import pytest
from vcforge._spec.reserved import reserved
from vcforge._spec.number import Number
from vcforge._spec.types import Type

def test_reserved_info_af():
    f = reserved("AF", "INFO")
    assert f.number == Number.A and f.type == Type.FLOAT and f.kind == "INFO"

def test_reserved_format_gt():
    f = reserved("GT", "FORMAT")
    assert f.number == Number.ONE and f.type == Type.STRING

def test_reserved_format_pl_is_g_integer():
    f = reserved("PL", "FORMAT")
    assert f.number == Number.G and f.type == Type.INTEGER

def test_reserved_flag():
    f = reserved("DB", "INFO")
    assert f.type == Type.FLAG and f.number == Number.FLAG

def test_unknown_reserved_raises():
    with pytest.raises(KeyError):
        reserved("NOPE", "INFO")
