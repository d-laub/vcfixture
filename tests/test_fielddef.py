import pytest

from vcfixture._spec.fielddef import FieldDef
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type


def test_valid_info_field():
    f = FieldDef("AF", Number.A, Type.FLOAT, "Allele frequency", "INFO")
    assert f.id == "AF"
    assert f.header_line() == (
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">'
    )


def test_flag_must_be_number_zero():
    with pytest.raises(ValueError, match="Flag"):
        FieldDef("DB", Number.ONE, Type.FLAG, "x", "INFO")


def test_flag_must_be_info():
    with pytest.raises(ValueError, match="Flag"):
        FieldDef("DB", Number.FLAG, Type.FLAG, "x", "FORMAT")


def test_format_cannot_be_flag():
    with pytest.raises(ValueError, match="Flag"):
        FieldDef("X", Number.ONE, Type.FLAG, "x", "FORMAT")


def test_id_must_match_regex():
    with pytest.raises(ValueError, match="ID"):
        FieldDef("1bad", Number.ONE, Type.STRING, "x", "INFO")
    FieldDef("1000G", Number.DOT, Type.STRING, "x", "INFO")


def test_kind_must_be_valid():
    with pytest.raises(ValueError, match="kind"):
        FieldDef("X", Number.ONE, Type.STRING, "x", "BOGUS")
