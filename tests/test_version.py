from vcfixture import VcfVersion
from vcfixture._spec.version import LATEST


def test_value_is_fileformat_string():
    assert VcfVersion.V4_3.value == "VCFv4.3"
    assert VcfVersion.V4_5.value == "VCFv4.5"


def test_ordering():
    assert VcfVersion.V4_1 < VcfVersion.V4_4
    assert VcfVersion.V4_4 <= VcfVersion.V4_4
    assert VcfVersion.V4_5 > VcfVersion.V4_3
    assert VcfVersion.V4_4 >= VcfVersion.V4_4


def test_sorted_is_chronological():
    assert sorted(VcfVersion) == [
        VcfVersion.V4_1,
        VcfVersion.V4_2,
        VcfVersion.V4_3,
        VcfVersion.V4_4,
        VcfVersion.V4_5,
    ]


def test_latest_is_v4_5():
    assert LATEST is VcfVersion.V4_5
