from __future__ import annotations

from .fielddef import FieldDef
from .number import Number
from .types import Type
from .version import LATEST, VcfVersion

_INFO = {
    "AA": FieldDef("AA", Number.ONE, Type.STRING, "Ancestral allele", "INFO"),
    "AC": FieldDef("AC", Number.A, Type.INTEGER, "Allele count", "INFO"),
    "AF": FieldDef("AF", Number.A, Type.FLOAT, "Allele frequency", "INFO"),
    "AN": FieldDef("AN", Number.ONE, Type.INTEGER, "Total allele number", "INFO"),
    "DP": FieldDef("DP", Number.ONE, Type.INTEGER, "Combined depth", "INFO"),
    "DB": FieldDef("DB", Number.FLAG, Type.FLAG, "dbSNP membership", "INFO"),
    "H2": FieldDef("H2", Number.FLAG, Type.FLAG, "HapMap2 membership", "INFO"),
    "END": FieldDef(
        "END", Number.ONE, Type.INTEGER, "End position (deprecated)", "INFO"
    ),
    "SVTYPE": FieldDef(
        "SVTYPE", Number.ONE, Type.STRING, "Type of structural variant", "INFO"
    ),
    "SVLEN": FieldDef(
        "SVLEN", Number.A, Type.INTEGER, "Length of structural variant", "INFO"
    ),
    "SVCLAIM": FieldDef(
        "SVCLAIM", Number.A, Type.STRING, "Structural variant claim", "INFO"
    ),
    "CIPOS": FieldDef(
        "CIPOS", Number.DOT, Type.INTEGER, "Confidence interval around POS", "INFO"
    ),
    "CIEND": FieldDef(
        "CIEND", Number.DOT, Type.INTEGER, "Confidence interval around END", "INFO"
    ),
    "CILEN": FieldDef(
        "CILEN", Number.DOT, Type.INTEGER, "Confidence interval around SVLEN", "INFO"
    ),
    "MATEID": FieldDef("MATEID", Number.A, Type.STRING, "ID of mate breakend", "INFO"),
    "PARID": FieldDef("PARID", Number.A, Type.STRING, "ID of partner breakend", "INFO"),
    "IMPRECISE": FieldDef(
        "IMPRECISE", Number.FLAG, Type.FLAG, "Imprecise structural variant", "INFO"
    ),
}
_FORMAT = {
    "GT": FieldDef("GT", Number.ONE, Type.STRING, "Genotype", "FORMAT"),
    "GQ": FieldDef("GQ", Number.ONE, Type.INTEGER, "Genotype quality", "FORMAT"),
    "DP": FieldDef("DP", Number.ONE, Type.INTEGER, "Read depth", "FORMAT"),
    "AD": FieldDef("AD", Number.R, Type.INTEGER, "Allelic depths", "FORMAT"),
    "PL": FieldDef(
        "PL", Number.G, Type.INTEGER, "Phred genotype likelihoods", "FORMAT"
    ),
    "GL": FieldDef("GL", Number.G, Type.FLOAT, "Log10 genotype likelihoods", "FORMAT"),
    "PS": FieldDef("PS", Number.ONE, Type.INTEGER, "Phase set", "FORMAT"),
    "CN": FieldDef("CN", Number.ONE, Type.FLOAT, "Copy number", "FORMAT"),
    "LEN": FieldDef(
        "LEN", Number.ONE, Type.INTEGER, "Length of <*> reference block", "FORMAT"
    ),
}


# Version each reserved field was introduced. Fields absent from these maps
# exist since 4.1. Verified against docs/reference/VCFv4.*.tex: among the current
# registry only SVCLAIM (INFO) and LEN (FORMAT, the <*> reference block) are new
# in 4.4; everything else is defined in the 4.1 spec.
_SINCE_INFO: dict[str, VcfVersion] = {"SVCLAIM": VcfVersion.V4_4}
_SINCE_FORMAT: dict[str, VcfVersion] = {"LEN": VcfVersion.V4_4}

# SVLEN's definition changed at the 4.3 -> 4.4 boundary (the one breaking change):
# Number=. + "difference in length" (signed) became Number=A + "length" (unsigned).
# The _INFO entry above holds the >= 4.4 form; this is the <= 4.3 form.
_SVLEN_PRE_4_4 = FieldDef(
    "SVLEN",
    Number.DOT,
    Type.INTEGER,
    "Difference in length between REF and ALT alleles",
    "INFO",
)


def reserved(id: str, kind: str, version: VcfVersion = LATEST) -> FieldDef:
    table = _INFO if kind == "INFO" else _FORMAT
    fd = table[id]  # KeyError => genuinely unknown reserved id
    since = (_SINCE_INFO if kind == "INFO" else _SINCE_FORMAT).get(id, VcfVersion.V4_1)
    if version < since:
        raise ValueError(
            f"{kind} field {id!r} was introduced in {since.value}; "
            f"not available in {version.value}"
        )
    if id == "SVLEN" and kind == "INFO" and version < VcfVersion.V4_4:
        return _SVLEN_PRE_4_4
    return fd
