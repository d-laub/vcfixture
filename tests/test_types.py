from vcfixture._spec.types import Type


def test_type_header_strings():
    assert Type.INTEGER.value == "Integer"
    assert Type.FLOAT.value == "Float"
    assert Type.FLAG.value == "Flag"
    assert Type.CHARACTER.value == "Character"
    assert Type.STRING.value == "String"


def test_format_allowed_types_excludes_flag():
    assert Type.FLAG not in Type.format_allowed()
    assert Type.INTEGER in Type.format_allowed()
    assert set(Type.info_allowed()) == set(Type)
