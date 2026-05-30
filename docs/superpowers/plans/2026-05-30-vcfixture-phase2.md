# vcfixture Phase 2 — Matrix & Multiallelic Round-Trip Coverage

> **For agentic workers:** execute task-by-task, TDD, commit per task. Use `.venv/bin/python` (CPython 3.12). NEVER system python. Commit author: `git -c user.name='David' -c user.email='david@standardmodel.bio'`.

**Goal:** Close the gap between the spec's success criteria and what the cyvcf2 round-trip actually validates. Make the strategy generate INFO/FORMAT *values* across the full classic Number×Type matrix and multiallelic records (exercising A/R/G + Number=G end-to-end), round-trip those values through cyvcf2, and add serializer percent-encoding + a reference-aware round-trip.

**Branch:** `impl/vcfixture` (continues the v0.1 foundation; 59 tests passing).

## Empirically-established cyvcf2 decode rules (from a probe of the current library)
These are FACTS to build against — do not re-derive:
- `variant.INFO.get(ID)`: Number=1 scalar → Python `int`/`float`/`str`; Number>1 or A/R/G/`.` → a `tuple`. Float values are **float32-quantized**. Flag → `True` (absent → `None`/`False`).
- `variant.format(ID)`: numeric → `np.ndarray` shape `(n_samples, count)`, dtype `int32` (Integer) or `float32` (Float). **Character → `np.ndarray` shape `(n_samples,)` of strings, each the per-sample values comma-joined** (e.g. `'a,b'`). String FORMAT behaves like Character (joined).
- `variant.genotypes`: `[[a1, a2, phased_bool], ...]`, `-1` for missing — multiallelic indices (e.g. `1/2` → `[1,2,False]`) work.
- **Float exactness:** generate floats with Hypothesis `width=32` so they are float32-representable; then decoded `float32` equals the generated value exactly (no tolerance needed).
- **Strings/chars in the matrix round-trip:** use a SAFE alphabet (alphanumeric only) so no escaping is needed and cyvcf2's percent-decoding is not on the critical path. Escaping is validated separately as a serializer unit test (Task P2).

---

## Task P1: Field-value generation strategy

**Files:** Modify `src/vcfixture/strategies.py`; Test `tests/test_field_value.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_field_value.py
from hypothesis import given, settings, strategies as st
from vcfixture import strategies as S
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture._spec.fielddef import FieldDef

_SAFE = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")

@settings(max_examples=50)
@given(st.data())
def test_field_value_cardinality_and_type(data):
    # Number=R, Float, biallelic diploid -> 2 float32-exact values
    fd = FieldDef("XR", Number.R, Type.FLOAT, "x", "INFO")
    v = data.draw(S.field_value(fd, n_alt=1, ploidy=2))
    assert isinstance(v, list) and len(v) == 2
    import numpy as np
    for x in v:
        assert isinstance(x, float)
        assert np.float32(x) == x  # float32-exact

def test_field_value_flag_is_true():
    fd = FieldDef("XF", Number.FLAG, Type.FLAG, "x", "INFO")
    # Flag carries no value list; strategy returns True deterministically.
    import hypothesis
    @hypothesis.given(hypothesis.strategies.data())
    def inner(data):
        assert data.draw(S.field_value(fd, n_alt=2, ploidy=2)) is True
    inner()

@settings(max_examples=30)
@given(st.data())
def test_field_value_G_count_multiallelic(data):
    fd = FieldDef("PL", Number.G, Type.INTEGER, "x", "FORMAT")
    # triallelic diploid -> C(3+2-1,2) = 6
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
```

- [ ] **Step 2: Run, confirm FAIL** (`S.field_value` missing).

Run: `.venv/bin/python -m pytest tests/test_field_value.py -q`

- [ ] **Step 3: Add to `src/vcfixture/strategies.py`** (append; keep existing code):

```python
# --- field-value generation (for the Number x Type matrix) ------------------
from ._spec.fielddef import FieldDef  # noqa: E402  (grouped with other imports is fine)

_SAFE_ALNUM = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

@st.composite
def _scalar_value(draw, typ: Type):
    if typ is Type.INTEGER:
        return draw(st.integers(min_value=-1000, max_value=1000))
    if typ is Type.FLOAT:
        # width=32 -> float32-exact so cyvcf2's float32 decode matches exactly.
        return draw(st.floats(min_value=-1e6, max_value=1e6,
                              allow_nan=False, allow_infinity=False, width=32))
    if typ is Type.CHARACTER:
        return draw(st.sampled_from(_SAFE_ALNUM))           # single char
    # STRING
    return draw(st.text(alphabet=_SAFE_ALNUM, min_size=1, max_size=6))

@st.composite
def field_value(draw, fielddef: FieldDef, n_alt: int, ploidy: int):
    """A spec-valid value for `fielddef` at a record with n_alt/ploidy.

    Flag -> True. Otherwise a list of `cardinality` scalars (Number=. picks a
    small random count). Values use safe alphabets / float32-exact floats so a
    cyvcf2 round-trip needs no escaping and no float tolerance.
    """
    if fielddef.type is Type.FLAG:
        return True
    card = fielddef.number.cardinality(n_alt, ploidy)
    if card is None:                      # Number=. (variable)
        card = draw(st.integers(min_value=1, max_value=3))
    return [draw(_scalar_value(fielddef.type)) for _ in range(card)]
```

- [ ] **Step 4: Run, confirm PASS.** `.venv/bin/python -m pytest tests/test_field_value.py -q`
- [ ] **Step 5: Commit** `feat: add field-value generation strategy for Number x Type matrix`

---

## Task P2: Serializer percent-encoding for String/Character values

**Files:** Modify `src/vcfixture/serialize.py`; Test add to `tests/test_serialize.py`

VCF requires reserved characters in String/Character INFO/FORMAT values to be percent-encoded. Encode: `%`→`%25`, `:`→`%3A`, `;`→`%3B`, `=`→`%3D`, `,`→`%2C`, `\r`→`%0D`, `\n`→`%0A`, `\t`→`%09`. (`%` MUST be encoded first.) Only str values are encoded; numbers/Genotype are untouched.

- [ ] **Step 1: Add failing test** to `tests/test_serialize.py`:

```python
def test_percent_encoding_of_reserved_chars():
    from vcfixture.serialize import _encode, _fmt_scalar
    assert _encode("a;b") == "a%3Bb"
    assert _encode("a:b,c=d") == "a%3Ab%2Cc%3Dd"
    assert _encode("100%") == "100%25"           # % encoded
    assert _encode("x\ty\n") == "x%09y%0A"
    # _fmt_scalar applies encoding to strings, not to ints/floats
    assert _fmt_scalar("a;b") == "a%3Bb"
    assert _fmt_scalar(5) == "5"
    assert _fmt_scalar(0.5) == "0.5"
```

- [ ] **Step 2: Run, confirm FAIL.** `.venv/bin/python -m pytest tests/test_serialize.py::test_percent_encoding_of_reserved_chars -q`

- [ ] **Step 3: Edit `src/vcfixture/serialize.py`** — add `_encode` and apply it to str values in `_fmt_scalar`:

```python
# Order matters: '%' must be replaced first.
_PERCENT = [("%", "%25"), (":", "%3A"), (";", "%3B"), ("=", "%3D"),
            (",", "%2C"), ("\r", "%0D"), ("\n", "%0A"), ("\t", "%09")]

def _encode(s: str) -> str:
    for ch, rep in _PERCENT:
        s = s.replace(ch, rep)
    return s
```

Then in `_fmt_scalar`, change the final `return str(v)` path so that genuine strings are encoded:

```python
def _fmt_scalar(v: Any) -> str:
    if v is None:
        return "."
    if isinstance(v, bool):           # guard: bool is a subclass of int
        return str(int(v))
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return "."
        return repr(v)
    if isinstance(v, str):
        return _encode(v)
    return str(v)
```

(The `bool` guard prevents a stray boolean from rendering as `True`/`False`; Flag handling in `_render_info` already intercepts real flags before this point, so this is just defensive.)

- [ ] **Step 4: Run, confirm PASS** + full `test_serialize.py` still green.
- [ ] **Step 5: Commit** `feat: percent-encode reserved chars in string values`

---

## Task P3: Multiallelic generation in `documents()`

**Files:** Modify `src/vcfixture/strategies.py`; Test add to `tests/test_strategies.py`

Make `documents()` draw `n_alt` in 1..3, build that many ALTs (mixed classes, but skip a 2nd `*` — at most one spanning deletion and only as the last alt), and draw GT allele indices in `0..n_alt`. Existing single-ALT behavior is the n_alt==1 case.

- [ ] **Step 1: Add failing test** to `tests/test_strategies.py`:

```python
@settings(max_examples=80, suppress_health_check=[HealthCheck.too_slow])
@given(S.documents(max_alt=3))
def test_documents_can_be_multiallelic(doc):
    # Over many examples at least one record should have >1 ALT.
    # (Property: every record's GT indices never exceed its n_alt.)
    from vcfixture.genotype import Genotype
    for rec in doc.records:
        n_alt = len(rec.alts)
        for s in rec.samples:
            gt = s["GT"]
            for a in gt.alleles:
                assert a is None or a <= n_alt
```

(Also keep the existing `test_documents_are_well_formed`.)

- [ ] **Step 2: Run, confirm FAIL** (`documents()` has no `max_alt` kwarg).

- [ ] **Step 3: Edit `documents()` in `src/vcfixture/strategies.py`.** Replace the record loop so it draws `n_alt` and multiple alts. Reference implementation:

```python
@st.composite
def documents(draw, max_samples: int = 3, max_records: int = 4, max_alt: int = 1):
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
            # at most one spanning '*', only as the final alt
            if klass == "SPANNING_DEL" and j != n_alt - 1:
                klass = "SNP"
            _, alt = draw(_ref_alt(klass))
            alts.append(alt)
        ref = draw(_ref_for_alts(alts))
        gts = [draw(genotypes(ploidy, n_alt=n_alt)) for _ in samples]
        b.record("chr1", pos, ref=ref, alt=alts, gt=gts)
        pos += draw(st.integers(1, 50))
    return b.build()
```

Add a small helper that picks a REF compatible with the chosen alts (REF is shared across alts in one record). Simplest correct approach: pick REF independently as a short base string, but ensure GT indices are the only correctness constraint (the builder does not validate REF/ALT shape beyond GT range). A minimal helper:

```python
@st.composite
def _ref_for_alts(draw, alts):
    # A single shared REF base; '*' alts impose no REF constraint, and the
    # builder/serializer/cyvcf2 accept any REF base with multiple ALTs.
    return draw(st.sampled_from(_BASES))
```

NOTE: `_ref_alt` returns `(ref, alt)` but here we only use its `alt`; REF is chosen by `_ref_for_alts`. This keeps multiallelic records simple and valid for the round-trip (which checks genotypes/pos/ref via cyvcf2). If cyvcf2 rejects any generated record, READ the emitted VCF and report the offending line.

- [ ] **Step 4: Run, confirm PASS.** Also re-run `tests/test_roundtrip.py` and `tests/test_strategies.py` — the existing round-trip uses default `max_alt=1` and must stay green.
- [ ] **Step 5: Commit** `feat: generate multiallelic records in documents() strategy`

---

## Task P4: Matrix round-trip through cyvcf2 (the headline)

**Files:** Modify `src/vcfixture/strategies.py` (add `documents_with_fields`); Test `tests/test_matrix_roundtrip.py`

Add a strategy that declares a set of INFO + FORMAT fields spanning the Number×Type matrix, populates each record's values via `field_value`, and a round-trip test that compares cyvcf2's decode to truth using the established decode rules.

- [ ] **Step 1: Add `documents_with_fields` to `src/vcfixture/strategies.py`:**

```python
def _matrix_field_defs():
    """One INFO and one FORMAT FieldDef per classic combo (skip Flag for FORMAT;
    Flag only as INFO). Deterministic, stable IDs."""
    info_defs, format_defs = [], []
    numbers = [("1", Number.ONE), ("2", Number.fixed(2)), ("A", Number.A),
               ("R", Number.R), ("G", Number.G), ("D", Number.DOT)]
    types = [("i", Type.INTEGER), ("f", Type.FLOAT),
             ("c", Type.CHARACTER), ("s", Type.STRING)]
    for nk, num in numbers:
        for tk, typ in types:
            info_defs.append(FieldDef(f"I{nk}{tk}", num, typ, "x", "INFO"))
            format_defs.append(FieldDef(f"F{nk}{tk}", num, typ, "x", "FORMAT"))
    info_defs.append(FieldDef("IFLAG", Number.FLAG, Type.FLAG, "x", "INFO"))
    return info_defs, format_defs

MATRIX_INFO_DEFS, MATRIX_FORMAT_DEFS = _matrix_field_defs()

@st.composite
def documents_with_fields(draw, max_samples: int = 3, max_records: int = 3,
                          max_alt: int = 3):
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
            fmt[fd.id] = [draw(field_value(fd, n_alt=n_alt, ploidy=ploidy))
                          for _ in samples]
        b.record("chr1", pos, ref=ref, alt=alts, gt=gts, info=info, **fmt)
        pos += draw(st.integers(1, 50))
    return b.build()
```

- [ ] **Step 2: Write the round-trip test** `tests/test_matrix_roundtrip.py`:

```python
"""Round-trip the full Number x Type matrix + multiallelic records through the
independent cyvcf2 parser and assert decoded values match the derived truth."""
import tempfile, os
import numpy as np
import pytest
from hypothesis import given, settings, HealthCheck
from vcfixture import strategies as S
from vcfixture._spec.types import Type

cyvcf2 = pytest.importorskip("cyvcf2")

def _info_defs_by_id(doc):
    return {fd.id: fd for fd in doc.info_defs}

def _fmt_defs_by_id(doc):
    return {fd.id: fd for fd in doc.format_defs}

def _norm_info_expected(typ, val):
    # Truth INFO values are lists (or True for Flag).
    if typ is Type.FLAG:
        return val  # True
    return list(val)

def _norm_info_got(typ, got):
    if typ is Type.FLAG:
        return bool(got) is True
    if got is None:
        return None
    if isinstance(got, tuple):
        return list(got)
    return [got]  # Number=1 scalar

@settings(max_examples=60, deadline=None,
          suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large])
@given(S.documents_with_fields())
def test_matrix_round_trips_through_cyvcf2(doc):
    truth = doc.truth()
    info_defs = _info_defs_by_id(doc)
    fmt_defs = _fmt_defs_by_id(doc)
    n_samples = len(doc.samples)

    d = tempfile.mkdtemp()
    path = doc.write(os.path.join(d, "m.vcf.gz"), bgzip=True, index=True)
    vf = cyvcf2.VCF(str(path))

    for ri, variant in enumerate(vf):
        # ---- INFO ----
        for fid, fd in info_defs.items():
            exp = _norm_info_expected(fd.type, truth.info[ri][fid])
            got = _norm_info_got(fd.type, variant.INFO.get(fid))
            if fd.type is Type.FLAG:
                assert got == exp, f"INFO {fid} flag r{ri}: {got} != {exp}"
            elif fd.type is Type.FLOAT:
                np.testing.assert_array_equal(
                    np.float32(got), np.float32(exp),
                    err_msg=f"INFO {fid} r{ri}")
            elif fd.type is Type.INTEGER:
                assert [int(x) for x in got] == [int(x) for x in exp], \
                    f"INFO {fid} r{ri}: {got} != {exp}"
            else:  # CHARACTER / STRING: cyvcf2 may comma-join multi-valued INFO
                got_join = got[0] if len(got) == 1 else ",".join(map(str, got))
                exp_join = ",".join(map(str, exp))
                assert got_join == exp_join, \
                    f"INFO {fid} r{ri}: {got_join!r} != {exp_join!r}"

        # ---- FORMAT ----
        for fid, fd in fmt_defs.items():
            got = variant.format(fid)  # ndarray
            for si in range(n_samples):
                exp = list(truth.format[ri][si][fid])
                if fd.type is Type.FLOAT:
                    np.testing.assert_array_equal(
                        np.float32(list(got[si])), np.float32(exp),
                        err_msg=f"FMT {fid} r{ri} s{si}")
                elif fd.type is Type.INTEGER:
                    assert [int(x) for x in got[si]] == [int(x) for x in exp], \
                        f"FMT {fid} r{ri} s{si}: {list(got[si])} != {exp}"
                else:  # CHARACTER/STRING: cyvcf2 returns one comma-joined str/sample
                    exp_join = ",".join(map(str, exp))
                    got_s = got[si]
                    got_s = got_s if isinstance(got_s, str) else str(got_s)
                    assert got_s == exp_join, \
                        f"FMT {fid} r{ri} s{si}: {got_s!r} != {exp_join!r}"
```

- [ ] **Step 3: Run.** `.venv/bin/python -m pytest tests/test_matrix_roundtrip.py -q --hypothesis-show-statistics`

This is the hardest task. Expect to ITERATE:
- If a Float comparison fails, confirm `width=32` is set in `_scalar_value` (Task P1) — without it float32 decode won't match.
- If a Character/String FORMAT comparison fails, inspect one decoded value (`print(repr(got))`) — adjust the join logic to match cyvcf2's ACTUAL shape (the probe showed `array(['a,b','c,d'])`, i.e. one joined string per sample). Do NOT change truth; adapt the comparison to cyvcf2's real behavior, which is the point of an independent oracle.
- If cyvcf2 reports a malformed file, write one failing `doc.render()` to a file and read it to find the bad line; the bug is in the serializer or an invalid doc the builder accepted — fix the SOURCE, add a focused unit regression test, and keep this property test intact.
- If Hypothesis flags `data_too_large`/too-slow, lower `max_records`/`max_samples` defaults in `documents_with_fields` (keep ≥2 each) rather than weakening assertions.

- [ ] **Step 4: Once green, commit** `test: round-trip full Number x Type matrix + multiallelic through cyvcf2` (include any source fixes in clearly-named commits).

---

## Task P5: Reference-aware round-trip

**Files:** Test `tests/test_reference_roundtrip.py`

Validate that a reference-anchored document (REF drawn from a FASTA) round-trips through cyvcf2 with matching POS/REF — the GVL-facing success criterion.

- [ ] **Step 1: Write the test:**

```python
# tests/test_reference_roundtrip.py
import tempfile, os
import pysam, pytest
from vcfixture.build import VcfBuilder
from vcfixture.reference import Reference

cyvcf2 = pytest.importorskip("cyvcf2")

def _ref(tmp):
    fa = os.path.join(tmp, "ref.fa")
    with open(fa, "w") as f:
        f.write(">chr1\n" + "ACGTACGTAC" * 10 + "\n")
    pysam.faidx(fa)
    return Reference(fa), fa

def test_reference_anchored_doc_round_trips():
    tmp = tempfile.mkdtemp()
    ref, _ = _ref(tmp)
    b = VcfBuilder(samples=["s1"], contigs=[("chr1", 100)]).fmt("GT")
    # Build SNP, DEL, INS, MNP anchored to the reference at distinct positions.
    specs = [(5, "SNP"), (15, "DEL"), (25, "INS"), (35, "MNP")]
    expected = []
    for pos0, klass in specs:
        rref, alts = ref.draw_ref_alt("chr1", pos0, klass=klass)
        b.record("chr1", pos0 + 1, ref=rref, alt=alts, gt=["0|1"])  # 1-based POS
        expected.append((pos0 + 1, rref))
    doc = b.build()
    path = doc.write(os.path.join(tmp, "r.vcf.gz"), bgzip=True, index=True)
    vf = cyvcf2.VCF(str(path))
    got = [(v.POS, v.REF) for v in vf]
    assert got == expected
    # REF must equal the actual reference sequence at each position.
    for (pos1, rref) in expected:
        assert ref.seq("chr1", pos1 - 1, len(rref)) == rref
```

- [ ] **Step 2: Run, confirm PASS.** If REF doesn't match the reference sequence, the bug is in `Reference.draw_ref_alt` (Task 13) — report it.
- [ ] **Step 3: Commit** `test: reference-anchored round-trip through cyvcf2`

---

## After all P-tasks
Run full suite: `.venv/bin/python -m pytest -q --hypothesis-show-statistics`. Report total count and the example counts for the matrix round-trip. Confirm the matrix round-trip exercised multiallelic records (n_alt>1) and all Number kinds.
