from __future__ import annotations

from hypothesis import strategies as st
from hypothesis.strategies import DrawFn

from ._spec.fielddef import FieldDef
from ._spec.number import Number
from ._spec.types import Type
from ._spec.version import LATEST, VcfVersion
from .allele import _SV_FIRST_TYPES, SymbolicAllele, UnspecifiedAllele, classify_allele
from .build import _SVCLAIM_RULES, VcfBuilder
from .model import VcfDocument
from .reference import ReferenceBuilder, ReferenceSpec
from .truth import GroundTruth
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
    """Draw a VCF GT string for a single sample.

    Constructs the genotype directly from drawn parameters — no rejection
    sampling — so every drawn value is a valid, spec-conformant GT token.

    Args:
        ploidy: Number of allele slots (e.g. 1 for haploid, 2 for diploid).
        n_alt: Number of ALT alleles; allele indices are drawn from ``[0, n_alt]``.
        missing_rate: Probability that any individual allele slot is ``"."``.
    """
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
    """Draw a spec-valid value for a single INFO or FORMAT field.

    Constructs the value directly — no rejection sampling.  Flag fields
    return ``True``; all other fields return a list of scalars whose length
    matches the resolved cardinality.  ``Number=.`` fields pick a small
    random count (1–3).  String/Character values use a safe alphanumeric
    alphabet; floats are float32-exact and finite.

    Args:
        fielddef: Field definition carrying the ``Number`` and ``Type``.
        n_alt: Number of ALT alleles for the current record (used to resolve
            ``Number.A`` and ``Number.R`` cardinalities).
        ploidy: Sample ploidy for the current record (used to resolve
            ``Number.G`` cardinality).
    """
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
    """Draw a ``VcfDocument`` exercising the full Number×Type field matrix.

    Declares every combination of ``Number`` (1, 2, A, R, G, .) and ``Type``
    (Integer, Float, Character, String) as both INFO and FORMAT fields, plus a
    FLAG INFO field.  Every record receives a drawn value for each field,
    providing exhaustive coverage of INFO/FORMAT decoding.  Records are
    reference-free (arbitrary single REF base) and not position-sorted.

    Args:
        max_samples: Upper bound on the number of samples drawn.
        max_records: Upper bound on the number of variant records drawn.
        max_alt: Upper bound on the number of ALT alleles per record.
    """
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
            alt=[classify_allele(a) for a in alts],
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
def references(
    draw: DrawFn,
    *,
    max_contigs: int = 2,
    max_contig_len: int = 2000,
    max_repeats: int = 3,
) -> ReferenceSpec:
    """Draw a small ``ReferenceSpec`` with optional planted tandem repeats.

    Generates a synthetic reference genome with 1–``max_contigs`` named contigs
    and up to ``max_repeats`` non-overlapping tandem repeats.  Repeat
    positions are tracked per contig to guarantee they never overlap.  The
    resulting ``ReferenceSpec`` is suitable for passing to ``documents`` or
    ``reference_and_documents``.

    Args:
        max_contigs: Upper bound on the number of contigs generated.
        max_contig_len: Upper bound on per-contig length in base pairs.
        max_repeats: Upper bound on the number of tandem repeats planted
            across all contigs.
    """
    seed = draw(st.integers(min_value=0, max_value=2**32 - 1))
    rb = ReferenceBuilder(seed=seed)

    n_contigs = draw(st.integers(1, max_contigs))
    lengths: dict[str, int] = {}
    for i in range(n_contigs):
        cid = f"chr{i + 1}"
        length = draw(st.integers(200, max_contig_len))
        rb.add_contig(cid, length)
        lengths[cid] = length

    # Plant repeats with a per-contig cursor so they never overlap.
    cursor = {cid: 50 for cid in lengths}
    n_rep = draw(st.integers(0, max_repeats))
    for _ in range(n_rep):
        cid = draw(st.sampled_from(list(lengths)))
        motif = draw(st.text(_BASES, min_size=1, max_size=3))
        count = draw(st.integers(3, 6))
        rlen = len(motif) * count
        pos0 = cursor[cid]
        if pos0 + rlen + 20 > lengths[cid]:
            continue
        rb.tandem_repeat(cid, pos0, motif, count)
        cursor[cid] = pos0 + rlen + draw(st.integers(20, 60))

    return rb.build()


_DEFAULT_LABELS = {
    "multiallelic": "multiallelic",
    "non_atomic": "non_atomic",
    "off_anchor": "off_anchor",
    "tandem_repeat": "tandem_repeat",
}


@st.composite
def _reference_documents(
    draw: DrawFn,
    reference: ReferenceSpec,
    violations: frozenset[str],
    label_overrides: dict[str, str] | None,
    max_samples: int,
    max_records: int,
    version: VcfVersion,
) -> VcfDocument:
    def lbl(key: str) -> str:
        base = _DEFAULT_LABELS[key]
        return (label_overrides or {}).get(base, base)

    n_samples = draw(st.integers(1, max_samples))
    samples = [f"s{i}" for i in range(n_samples)]
    ploidy = draw(st.integers(1, 2))

    # Prefer a contig that carries a repeat when non_left_aligned is requested.
    repeat_contigs = sorted({rf.contig for rf in reference.repeats})
    if "non_left_aligned" in violations and repeat_contigs:
        contig = draw(st.sampled_from(repeat_contigs))
    else:
        contig = draw(st.sampled_from([cid for cid, _ in reference.contigs]))
    clen = reference.length(contig)
    contig_repeats = [rf for rf in reference.repeats if rf.contig == contig]

    b = VcfBuilder(
        samples=samples,
        contigs=[(cid, reference.length(cid)) for cid, _ in reference.contigs],
        version=version,
    ).fmt("GT")

    enabled = [
        v for v in ("multiallelic", "non_atomic", "non_left_aligned") if v in violations
    ]

    n_rec = draw(st.integers(1, max_records))
    cursor = 10  # low start so planted repeats (pos0 >= 50) are reachable ahead
    for _ in range(n_rec):
        if cursor + 30 >= clen:
            break
        # Decide what kind of record to emit. Anchor `a` defaults to the cursor;
        # only a forward jump (off_anchor) may move it ahead, never behind, so
        # records stay position-sorted.
        a = cursor
        kind = (
            draw(st.sampled_from(["canonical", *enabled])) if enabled else "canonical"
        )

        # off_anchor needs a planted repeat strictly AHEAD of the cursor; if none
        # is available, fall back to a canonical record.
        usable_repeats = [rf for rf in contig_repeats if rf.pos0 >= cursor]
        if kind == "non_left_aligned" and not usable_repeats:
            kind = "canonical"

        labels: set[str] = set()
        if kind == "non_left_aligned":
            rf = draw(st.sampled_from(usable_repeats))
            mlen = len(rf.motif)
            # Anchor inside the repeat (copy index k>=1): a non-left-aligned
            # representation of deleting one motif unit. a = rf.pos0 + mlen*k - 1
            # is always >= cursor (rf.pos0 >= cursor), so order is preserved.
            k = draw(st.integers(1, max(1, rf.count - 1)))
            a = rf.pos0 + mlen * k - 1
            ref = reference.seq(contig, a, 1 + mlen)
            alts = [ref[0]]
            labels = {lbl("off_anchor"), lbl("tandem_repeat")}
        elif kind == "multiallelic":
            r = reference.base(contig, a)
            others = [x for x in "ACGT" if x != r]
            i = draw(st.integers(0, len(others) - 1))
            alt1 = others[i]
            alt2 = others[(i + 1) % len(others)]
            ref, alts = r, [alt1, alt2]
            labels = {lbl("multiallelic")}
        elif kind == "non_atomic":
            ref, alts = reference.draw_ref_alt(contig, a, "MNP", mnp_len=2)
            labels = {lbl("non_atomic")}
        else:  # canonical: SNP, or a left-aligned 1bp DEL when context allows
            want_del = draw(st.booleans())
            if (
                want_del
                and a + 2 < clen
                and reference.base(contig, a) != reference.base(contig, a + 1)
            ):
                ref = reference.seq(contig, a, 2)  # delete the differing next base
                alts = [ref[0]]
            else:
                ref, alts = reference.draw_ref_alt(contig, a, "SNP")

        gts = [draw(genotypes(ploidy, n_alt=len(alts))) for _ in samples]
        b.record(
            contig,
            a + 1,  # 1-based POS
            ref=ref,
            alt=[classify_allele(a) for a in alts],
            gt=gts,
            labels=sorted(labels) if labels else None,
        )
        cursor = a + len(ref) + draw(st.integers(20, 60))

    return b.build()


@st.composite
def documents(
    draw: DrawFn,
    max_samples: int = 3,
    max_records: int = 4,
    max_alt: int = 1,
    *,
    reference: ReferenceSpec | None = None,
    violations: frozenset[str] = frozenset(),
    label_overrides: dict[str, str] | None = None,
    version: VcfVersion = LATEST,
) -> VcfDocument:
    """Draw a small ``VcfDocument`` over a synthetic contig.

    Without ``reference``, records are drawn independently: REF alleles are
    arbitrary single bases and records are not sorted by position.

    When ``reference`` is provided, every REF allele matches the reference
    sequence at its position and records are emitted in position order.
    Violation classes may be opted into to intentionally emit non-canonical
    records alongside the canonical ones.

    Args:
        max_samples: Upper bound on the number of samples drawn.
        max_records: Upper bound on the number of variant records drawn.
        max_alt: Upper bound on the number of ALT alleles per record
            (ignored when ``reference`` is provided).
        reference: When given, draws reference-consistent, position-sorted
            records against this spec.  When ``None``, draws a standalone
            document over a single synthetic ``chr1`` contig.
        violations: Set of violation classes to opt into.  Only meaningful
            when ``reference`` is provided.  Accepted values:

            - ``"multiallelic"`` — emit a record with two ALT alleles;
              labels the record ``multiallelic``.
            - ``"non_atomic"`` — emit a 2-bp MNP that could be decomposed;
              labels the record ``non_atomic``.
            - ``"non_left_aligned"`` — emit a deletion anchored inside a
              planted tandem repeat rather than at its leftmost position;
              labels the record with both ``off_anchor`` and
              ``tandem_repeat``.

            The full emitted label vocabulary is therefore
            ``{multiallelic, non_atomic, off_anchor, tandem_repeat}``.
        label_overrides: Optional mapping from default label strings to
            replacement strings, applied to every label before it is
            attached to a record.  Keys must be drawn from the default
            vocabulary (``multiallelic``, ``non_atomic``, ``off_anchor``,
            ``tandem_repeat``).
    """
    if reference is not None:
        return draw(
            _reference_documents(
                reference,
                violations,
                label_overrides,
                max_samples,
                max_records,
                version,
            )
        )
    # --- existing reference-free body continues unchanged below ---
    n_samples = draw(st.integers(1, max_samples))
    samples = [f"s{i}" for i in range(n_samples)]
    ploidy = draw(st.integers(1, 2))
    b = VcfBuilder(samples=samples, contigs=[("chr1", 100000)], version=version).fmt(
        "GT"
    )

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
        b.record("chr1", pos, ref=ref, alt=[classify_allele(a) for a in alts], gt=gts)
        pos += draw(st.integers(1, 50))
    return b.build()


_SYMBOLIC_SV = sorted(_SV_FIRST_TYPES)


@st.composite
def symbolic_documents(
    draw: DrawFn,
    max_samples: int = 3,
    max_records: int = 4,
    version: VcfVersion = LATEST,
) -> VcfDocument:
    """Draw a ``VcfDocument`` whose records carry symbolic SV alleles.

    Each record holds one symbolic allele drawn from the full set of first-level
    SV types (``DEL``, ``DUP``, ``INS``, ``INV``, ``CNV``, …) or the unspecified
    allele ``<*>``.  Symbolic records include consistent per-allele ``SVLEN`` and
    ``SVCLAIM`` INFO values; unspecified records carry no INFO fields.  Records
    are reference-free (arbitrary single REF base) with positions spaced 2–5 kbp
    apart.

    Args:
        max_samples: Upper bound on the number of samples drawn.
        max_records: Upper bound on the number of variant records drawn.
    """
    n_samples = draw(st.integers(1, max_samples))
    samples = [f"s{i}" for i in range(n_samples)]
    ploidy = draw(st.integers(1, 2))
    has_svclaim = version >= VcfVersion.V4_4
    b = (
        VcfBuilder(samples=samples, contigs=[("chr1", 1_000_000)], version=version)
        .fmt("GT")
        .info("SVLEN")  # version-correct Number/description from the registry
    )
    if has_svclaim:
        b.info("SVCLAIM")
    n_rec = draw(st.integers(1, max_records))
    pos = 1000
    for _ in range(n_rec):
        ref = draw(st.sampled_from(_BASES))
        kind = draw(st.sampled_from([*_SYMBOLIC_SV, "UNSPEC"]))
        if kind == "UNSPEC":
            alts: list[SymbolicAllele | UnspecifiedAllele] = [UnspecifiedAllele()]
            info: dict[str, object] = {}
        else:
            magnitude = draw(st.integers(1, 1000))
            # <= 4.3: SVLEN is the signed REF/ALT length difference (DEL negative).
            # >= 4.4: SVLEN is an unsigned length.
            is_pre_44_del = kind == "DEL" and version < VcfVersion.V4_4
            svlen = -magnitude if is_pre_44_del else magnitude
            alts = [SymbolicAllele(kind)]
            info = {"SVLEN": [svlen]}
            if has_svclaim:
                claim = draw(st.sampled_from(sorted(_SVCLAIM_RULES[kind])))
                info["SVCLAIM"] = [claim]
        gts = [draw(genotypes(ploidy, n_alt=1)) for _ in samples]
        b.record("chr1", pos, ref=ref, alt=alts, gt=gts, info=info)
        pos += draw(st.integers(2000, 5000))
    return b.build()


@st.composite
def reference_and_documents(
    draw: DrawFn,
    *,
    max_samples: int = 3,
    max_records: int = 4,
    violations: frozenset[str] = frozenset(),
    label_overrides: dict[str, str] | None = None,
    max_contigs: int = 2,
    max_contig_len: int = 2000,
    max_repeats: int = 3,
    version: VcfVersion = LATEST,
) -> tuple[ReferenceSpec, VcfDocument, GroundTruth]:
    """Draw a consistent ``(ReferenceSpec, VcfDocument, GroundTruth)`` triple.

    Combines ``references`` and ``documents`` so that the document's REF alleles
    match the drawn reference and records are in position order.  The
    ``GroundTruth`` is derived in the same draw, giving a fully consistent
    triple ready for round-trip testing.

    Args:
        max_samples: Upper bound on the number of samples drawn.
        max_records: Upper bound on the number of variant records drawn.
        violations: Set of violation classes to opt into (see ``documents``).
        label_overrides: Mapping from default label strings to replacement
            strings (see ``documents``).
        max_contigs: Upper bound on the number of reference contigs drawn.
        max_contig_len: Upper bound on per-contig length in base pairs.
        max_repeats: Upper bound on the number of tandem repeats planted.
    """
    spec = draw(
        references(
            max_contigs=max_contigs,
            max_contig_len=max_contig_len,
            max_repeats=max_repeats,
        )
    )
    doc = draw(
        documents(
            max_samples=max_samples,
            max_records=max_records,
            reference=spec,
            violations=violations,
            label_overrides=label_overrides,
            version=version,
        )
    )
    return spec, doc, doc.truth()
