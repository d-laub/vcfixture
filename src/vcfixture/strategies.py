from __future__ import annotations

from hypothesis import strategies as st
from hypothesis.strategies import DrawFn

from ._spec.fielddef import FieldDef
from ._spec.number import Number
from ._spec.types import Type
from .build import VcfBuilder
from .model import VcfDocument
from .variants import deletion, delins, insertion, mnp, snp, spanning_deletion


def _build_number_type_combos():
    numbers = [Number.ONE, Number.fixed(2), Number.A, Number.R, Number.G, Number.DOT]
    combos = []
    for kind in ("INFO", "FORMAT"):
        allowed = Type.info_allowed() if kind == "INFO" else Type.format_allowed()
        for n in numbers:
            for t in allowed:
                if t is Type.FLAG:
                    continue
                combos.append((n, t, kind))
        if kind == "INFO":
            combos.append((Number.FLAG, Type.FLAG, "INFO"))
    return combos


ALL_NUMBER_TYPE_COMBOS = _build_number_type_combos()
ALL_VARIANT_CLASSES = ["SNP", "MNP", "INS", "DEL", "DELINS", "SPANNING_DEL"]

_BASES = "ACGT"


@st.composite
def _ref_alt(draw: DrawFn, klass: str) -> tuple[str, str]:
    b = draw(st.sampled_from(_BASES))
    b2 = draw(st.sampled_from(_BASES))
    if klass == "SNP":
        alt = _BASES[(_BASES.index(b) + 1 + draw(st.integers(0, 2))) % 4]
        return snp(b, alt)
    if klass == "MNP":
        return mnp(
            b + b2,
            _BASES[(_BASES.index(b) + 1) % 4] + _BASES[(_BASES.index(b2) + 1) % 4],
        )
    if klass == "INS":
        return insertion(b, draw(st.text(_BASES, min_size=1, max_size=3)))
    if klass == "DEL":
        return deletion(b, draw(st.text(_BASES, min_size=1, max_size=3)))
    if klass == "DELINS":
        return delins(b + b2, draw(st.text(_BASES, min_size=1, max_size=3)))
    return spanning_deletion(b)


@st.composite
def genotypes(draw: DrawFn, ploidy: int, n_alt: int, missing_rate: float = 0.1) -> str:
    alleles = []
    for _ in range(ploidy):
        if draw(st.floats(0, 1)) < missing_rate:
            alleles.append(".")
        else:
            alleles.append(str(draw(st.integers(0, n_alt))))
    phased = draw(st.booleans())
    sep = "|" if phased else "/"
    return sep.join(alleles)


_SAFE_ALNUM = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


@st.composite
def _scalar_value(draw: DrawFn, typ: Type) -> int | float | str:
    if typ is Type.INTEGER:
        return draw(st.integers(min_value=-1000, max_value=1000))
    if typ is Type.FLOAT:
        return draw(
            st.floats(
                min_value=-1e6,
                max_value=1e6,
                allow_nan=False,
                allow_infinity=False,
                width=32,
            )
        )
    if typ is Type.CHARACTER:
        return draw(st.sampled_from(_SAFE_ALNUM))
    return draw(st.text(alphabet=_SAFE_ALNUM, min_size=1, max_size=6))


@st.composite
def field_value(
    draw: DrawFn, fielddef: FieldDef, n_alt: int, ploidy: int
) -> bool | list[int | float | str]:
    """A spec-valid value for `fielddef` at a record with n_alt/ploidy.
    Flag -> True. Otherwise a list of `cardinality` scalars (Number=. picks a
    small random count). Safe alphabets / float32-exact floats."""
    if fielddef.type is Type.FLAG:
        return True
    card = fielddef.number.cardinality(n_alt, ploidy)
    if card is None:
        card = draw(st.integers(min_value=1, max_value=3))
    return [draw(_scalar_value(fielddef.type)) for _ in range(card)]


def _matrix_field_defs() -> tuple[list[FieldDef], list[FieldDef]]:
    """One INFO and one FORMAT FieldDef per classic combo (Flag only as INFO)."""
    info_defs: list[FieldDef] = []
    format_defs: list[FieldDef] = []
    numbers = [
        ("1", Number.ONE),
        ("2", Number.fixed(2)),
        ("A", Number.A),
        ("R", Number.R),
        ("G", Number.G),
        ("D", Number.DOT),
    ]
    types = [
        ("i", Type.INTEGER),
        ("f", Type.FLOAT),
        ("c", Type.CHARACTER),
        ("s", Type.STRING),
    ]
    for nk, num in numbers:
        for tk, typ in types:
            info_defs.append(FieldDef(f"I{nk}{tk}", num, typ, "x", "INFO"))
            format_defs.append(FieldDef(f"F{nk}{tk}", num, typ, "x", "FORMAT"))
    info_defs.append(FieldDef("IFLAG", Number.FLAG, Type.FLAG, "x", "INFO"))
    return info_defs, format_defs


MATRIX_INFO_DEFS, MATRIX_FORMAT_DEFS = _matrix_field_defs()


@st.composite
def documents_with_fields(
    draw: DrawFn, max_samples: int = 3, max_records: int = 3, max_alt: int = 3
) -> VcfDocument:
    n_samples = draw(st.integers(1, max_samples))
    samples = [f"s{i}" for i in range(n_samples)]
    ploidy = draw(st.integers(1, 2))
    b = VcfBuilder(samples=samples, contigs=[("chr1", 100000)])
    b.fmt("GT")
    for fd in MATRIX_INFO_DEFS:
        b.info(fd.id, fd.number, fd.type)
    for fd in MATRIX_FORMAT_DEFS:
        b.fmt(fd.id, fd.number, fd.type)

    n_rec = draw(st.integers(1, max_records))
    pos = 1000
    for _ in range(n_rec):
        n_alt = draw(st.integers(1, max_alt))
        alts = []
        for j in range(n_alt):
            klass = draw(st.sampled_from(ALL_VARIANT_CLASSES))
            if klass == "SPANNING_DEL" and j != n_alt - 1:
                klass = "SNP"
            _, alt = draw(_ref_alt(klass))
            alts.append(alt)
        ref = draw(st.sampled_from(_BASES))
        gts = [draw(genotypes(ploidy, n_alt=n_alt)) for _ in samples]
        info = {}
        for fd in MATRIX_INFO_DEFS:
            info[fd.id] = draw(field_value(fd, n_alt=n_alt, ploidy=ploidy))
        fmt = {}
        for fd in MATRIX_FORMAT_DEFS:
            fmt[fd.id] = [
                draw(field_value(fd, n_alt=n_alt, ploidy=ploidy)) for _ in samples
            ]
        b.record(
            "chr1",
            pos,
            ref=ref,
            alt=alts,
            gt=gts,
            info=info,
            # FORMAT field IDs never collide with record()'s named keyword
            # params (ids/qual/filter), but the checker can't know that when
            # unpacking a str-keyed dict, so it flags a spurious type clash.
            **fmt,  # pyrefly: ignore[bad-argument-type]
        )
        pos += draw(st.integers(1, 50))
    return b.build()


@st.composite
def documents(
    draw: DrawFn, max_samples: int = 3, max_records: int = 4, max_alt: int = 1
) -> VcfDocument:
    n_samples = draw(st.integers(1, max_samples))
    samples = [f"s{i}" for i in range(n_samples)]
    ploidy = draw(st.integers(1, 2))
    b = VcfBuilder(samples=samples, contigs=[("chr1", 100000)]).fmt("GT")

    n_rec = draw(st.integers(1, max_records))
    pos = 1000
    for _ in range(n_rec):
        n_alt = draw(st.integers(1, max_alt))
        alts = []
        for j in range(n_alt):
            klass = draw(st.sampled_from(ALL_VARIANT_CLASSES))
            if klass == "SPANNING_DEL" and j != n_alt - 1:
                klass = "SNP"
            _, alt = draw(_ref_alt(klass))
            alts.append(alt)
        ref = draw(st.sampled_from(_BASES))
        gts = [draw(genotypes(ploidy, n_alt=n_alt)) for _ in samples]
        b.record("chr1", pos, ref=ref, alt=alts, gt=gts)
        pos += draw(st.integers(1, 50))
    return b.build()
