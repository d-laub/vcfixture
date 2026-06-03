# Per-version VCF generation (`version=`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `VcfVersion` enum and `version=` argument to `VcfBuilder` and the document-producing strategies so vcfixture generates VCFs that conform to spec versions 4.1–4.5, gating reserved-field availability and the one breaking change (the `SVLEN` definition flip at the 4.3/4.4 boundary).

**Architecture:** One orderable `VcfVersion` enum drives everything. The reserved-field registry (`reserved()`) gains a `version` parameter: it rejects fields not yet introduced and returns the version-correct `SVLEN` `FieldDef`. The model stores `version` instead of a free `fileformat` string; serialization derives the `##fileformat` header from it. `VcfBuilder` threads its version into `reserved()` and into eager validation; strategies thread it into the builder and emit version-correct symbolic-SV metadata. Ground-truth needs no change (it already `abs()`-normalizes `SVLEN`).

**Tech Stack:** Python 3.10+, `uv` (run everything via `uv run`), Hypothesis, pytest, cyvcf2/pysam (independent oracle), pyrefly (strict type-check on `src/`), ruff.

**Reference:** Design spec at `docs/superpowers/specs/2026-06-02-vcf-version-conformance-design.md`. Spec sources tracked at `docs/reference/VCFv4.{1,2,3,4,5}.tex`.

---

## File Structure

- **Create** `src/vcfixture/_spec/version.py` — the `VcfVersion` enum + `LATEST`. One responsibility: the version vocabulary and its ordering.
- **Modify** `src/vcfixture/_spec/reserved.py` — per-field `since` table, `version` param on `reserved()`, the `SVLEN` pre-4.4 variant.
- **Modify** `src/vcfixture/model.py` — `VcfDocument.fileformat: str` → `version: VcfVersion`.
- **Modify** `src/vcfixture/serialize.py` — header line uses `doc.version.value`.
- **Modify** `src/vcfixture/build.py` — `version=` param; thread into `reserved()`, eager validation, and `VcfDocument(...)`.
- **Modify** `src/vcfixture/strategies.py` — `version=` on `documents`, `symbolic_documents`, `reference_and_documents` (+ internal `_reference_documents`); version-correct `SVLEN`/`SVCLAIM`.
- **Modify** `src/vcfixture/__init__.py` — export `VcfVersion`.
- **Create** `tests/test_version.py` — enum behavior.
- **Modify** test files that construct `VcfDocument(...)` directly (rename), and the two cross-parser round-trip tests (parametrize across versions).

**Notes / scope decisions locked here:**
- The design mentioned an `available()` helper for strategies to enumerate reserved IDs. **Omitted as YAGNI** — no strategy enumerates reserved fields generically; version-correctness is achieved by conditional declaration in `symbolic_documents` plus `reserved()` returning the version-correct `SVLEN`. Do not add `available()` unless a later task needs it.
- `VcfVersion` defines all four comparison dunders explicitly (not `functools.total_ordering`) to stay unambiguous under pyrefly strict.

---

## Task 1: `VcfVersion` enum

**Files:**
- Create: `src/vcfixture/_spec/version.py`
- Modify: `src/vcfixture/__init__.py`
- Test: `tests/test_version.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_version.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_version.py -q`
Expected: FAIL — `ImportError: cannot import name 'VcfVersion' from 'vcfixture'`.

- [ ] **Step 3: Create the enum**

Create `src/vcfixture/_spec/version.py`:

```python
from __future__ import annotations

from enum import Enum


class VcfVersion(Enum):
    """A supported VCF spec version.

    ``.value`` is the exact ``##fileformat`` string. Members are orderable in
    chronological (declaration) order, so ``VcfVersion.V4_3 < VcfVersion.V4_4``.
    """

    V4_1 = "VCFv4.1"
    V4_2 = "VCFv4.2"
    V4_3 = "VCFv4.3"
    V4_4 = "VCFv4.4"
    V4_5 = "VCFv4.5"

    def __lt__(self, other: VcfVersion) -> bool:
        return _ORDER[self] < _ORDER[other]

    def __le__(self, other: VcfVersion) -> bool:
        return _ORDER[self] <= _ORDER[other]

    def __gt__(self, other: VcfVersion) -> bool:
        return _ORDER[self] > _ORDER[other]

    def __ge__(self, other: VcfVersion) -> bool:
        return _ORDER[self] >= _ORDER[other]


_ORDER: dict[VcfVersion, int] = {v: i for i, v in enumerate(VcfVersion)}

LATEST = VcfVersion.V4_5
```

- [ ] **Step 4: Export `VcfVersion` from the package root**

In `src/vcfixture/__init__.py`, add the import next to the other `_spec` imports (after line 7, `from ._spec.types import Type`):

```python
from ._spec.version import VcfVersion
```

And add `"VcfVersion"` to `__all__` (keep the list alphabetically ordered — insert it just before `"__version__"` / after `"VcfBuilder"`):

```python
    "VcfBuilder",
    "VcfVersion",
    "strategies",
    "__version__",
```

- [ ] **Step 5: Run tests + type-check**

Run: `uv run pytest tests/test_version.py -q`
Expected: PASS (4 tests).

Run: `uv run pyrefly check`
Expected: PASS (0 errors).

- [ ] **Step 6: Commit**

```bash
git add src/vcfixture/_spec/version.py src/vcfixture/__init__.py tests/test_version.py
git commit -m "feat: add orderable VcfVersion enum"
```

---

## Task 2: Reserved-field ladder + `SVLEN` switch

**Files:**
- Modify: `src/vcfixture/_spec/reserved.py`
- Test: `tests/test_reserved.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reserved.py` (it already imports `reserved`; add the new imports at the top of the file if not present: `import pytest`, `from vcfixture import VcfVersion`, `from vcfixture._spec.number import NumberKind`):

```python
def test_svclaim_rejected_before_4_4():
    with pytest.raises(ValueError, match="introduced in VCFv4.4"):
        reserved("SVCLAIM", "INFO", VcfVersion.V4_3)


def test_svclaim_available_at_4_4():
    assert reserved("SVCLAIM", "INFO", VcfVersion.V4_4).id == "SVCLAIM"


def test_len_rejected_before_4_4():
    with pytest.raises(ValueError, match="introduced in VCFv4.4"):
        reserved("LEN", "FORMAT", VcfVersion.V4_1)


def test_svlen_definition_flips_at_4_4():
    old = reserved("SVLEN", "INFO", VcfVersion.V4_3)
    assert old.number.kind is NumberKind.DOT
    assert "Difference in length" in old.description
    new = reserved("SVLEN", "INFO", VcfVersion.V4_4)
    assert new.number.kind is NumberKind.A
    assert new.description == "Length of structural variant"


def test_unknown_id_still_keyerror():
    with pytest.raises(KeyError):
        reserved("NOPE", "INFO", VcfVersion.V4_5)


def test_default_version_is_latest():
    # existing 2-arg call sites keep working and see the latest definitions
    assert reserved("SVLEN", "INFO").number.kind is NumberKind.A
    assert reserved("SVCLAIM", "INFO").id == "SVCLAIM"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_reserved.py -q`
Expected: FAIL — `reserved()` takes 2 positional args, the new tests pass 3 (`TypeError`).

- [ ] **Step 3: Add the version parameter, `since` table, and SVLEN variant**

In `src/vcfixture/_spec/reserved.py`:

Add the import after line 5 (`from .types import Type`):

```python
from .version import LATEST, VcfVersion
```

Then, after the `_FORMAT = {...}` dict (currently ends at line 56) and before `def reserved`, insert:

```python
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
```

Replace the existing `reserved` function (currently lines 59–61):

```python
def reserved(id: str, kind: str) -> FieldDef:
    table = _INFO if kind == "INFO" else _FORMAT
    return table[id]
```

with:

```python
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
```

- [ ] **Step 4: Run tests + type-check**

Run: `uv run pytest tests/test_reserved.py -q`
Expected: PASS (all, including the pre-existing reserved tests).

Run: `uv run pyrefly check`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/_spec/reserved.py tests/test_reserved.py
git commit -m "feat: gate reserved fields by version; SVLEN definition flips at 4.4"
```

---

## Task 3: Model + serialization (`version` replaces `fileformat`)

**Files:**
- Modify: `src/vcfixture/model.py`
- Modify: `src/vcfixture/serialize.py`
- Test (migrate existing direct constructions): `tests/test_model.py`, `tests/test_truth.py`, `tests/test_serialize.py`, `tests/test_symbolic_truth.py`, `tests/test_repr.py`

> This task renames a model field, which breaks every direct `VcfDocument(...)` construction in tests. `VcfBuilder` is migrated in Task 4, so after this task the builder still references the old name — that is expected; run the targeted test files listed below, not the whole suite, until Task 4.

- [ ] **Step 1: Rename the model field**

In `src/vcfixture/model.py`:

Add the import after line 9 (`from ._spec.fielddef import FieldDef`):

```python
from ._spec.version import VcfVersion
```

Change the field declaration (line 77) from:

```python
    fileformat: str
```

to:

```python
    version: VcfVersion
```

Change the `__repr__` line (line 89) from:

```python
            f"VcfDocument({self.fileformat} samples={len(self.samples)} "
```

to:

```python
            f"VcfDocument({self.version.value} samples={len(self.samples)} "
```

- [ ] **Step 2: Update serialization**

In `src/vcfixture/serialize.py`, change line 95 from:

```python
    lines = [f"##fileformat={doc.fileformat}"]
```

to:

```python
    lines = [f"##fileformat={doc.version.value}"]
```

- [ ] **Step 3: Migrate direct `VcfDocument(...)` constructions in tests**

Each test below constructs `VcfDocument` directly. Add `from vcfixture import VcfVersion` to the imports of each file, then change the version argument:

- `tests/test_model.py`:
  - line 24: `fileformat="VCFv4.5",` → `version=VcfVersion.V4_5,`
  - line ~55 (the positional construction): `"VCFv4.5", (), (), ()` → `VcfVersion.V4_5, (), (), ()`
- `tests/test_truth.py` line ~26 (positional): `"VCFv4.5", (), (), ()` → `VcfVersion.V4_5, (), (), ()`
- `tests/test_serialize.py` line ~31 (positional, first arg on its own line): `"VCFv4.5",` → `VcfVersion.V4_5,`
- `tests/test_symbolic_truth.py` line 24: `fileformat="VCFv4.5",` → `version=VcfVersion.V4_5,`
- `tests/test_repr.py` lines 114 and 157: `fileformat="VCFv4.5",` → `version=VcfVersion.V4_5,`

Do **not** change assertion strings: `"##fileformat=VCFv4.5"` (test_serialize.py:43, test_build.py:30) and `"VcfDocument(VCFv4.5 ...)"` (test_repr.py:122, 166) remain correct because `VcfVersion.V4_5.value == "VCFv4.5"`.

- [ ] **Step 4: Run the affected tests + type-check**

Run: `uv run pytest tests/test_model.py tests/test_truth.py tests/test_serialize.py tests/test_symbolic_truth.py tests/test_repr.py -q`
Expected: PASS.

Run: `uv run pyrefly check`
Expected: FAIL with errors in `build.py` only (it still uses `fileformat=`/`self._fileformat`). This is expected and fixed in Task 4. Confirm the *only* errors are in `build.py`.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/model.py src/vcfixture/serialize.py tests/test_model.py tests/test_truth.py tests/test_serialize.py tests/test_symbolic_truth.py tests/test_repr.py
git commit -m "feat: store VcfVersion on VcfDocument; derive fileformat header from it"
```

---

## Task 4: `VcfBuilder` `version=` parameter

**Files:**
- Modify: `src/vcfixture/build.py`
- Test: `tests/test_build.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_build.py` (add imports at top if absent: `import pytest`, and extend the `vcfixture` import to include `VcfVersion` and `Sym`):

```python
def test_builder_version_sets_header():
    from vcfixture import VcfBuilder, VcfVersion

    doc = VcfBuilder(
        samples=["s1"], contigs=[("chr1", 1000)], version=VcfVersion.V4_2
    ).fmt("GT").build()
    assert doc.render().startswith("##fileformat=VCFv4.2\n")


def test_builder_rejects_svclaim_before_4_4():
    from vcfixture import VcfBuilder, VcfVersion

    b = VcfBuilder(samples=["s1"], contigs=[("chr1", 1000)], version=VcfVersion.V4_3)
    with pytest.raises(ValueError, match="introduced in VCFv4.4"):
        b.info("SVCLAIM")


def test_svlen_number_a_count_enforced_at_4_4():
    from vcfixture import Sym, VcfBuilder, VcfVersion

    b = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 100000)], version=VcfVersion.V4_4)
        .fmt("GT")
        .info("SVLEN")
    )
    with pytest.raises(ValueError, match="cardinality"):
        # one ALT but two SVLEN values: Number=A requires exactly n_alt
        b.record("chr1", 10, ref="A", alt=[Sym("INS")], gt=["0/1"], info={"SVLEN": [50, 60]})


def test_svlen_any_count_allowed_at_4_3():
    from vcfixture import Sym, VcfBuilder, VcfVersion

    b = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 100000)], version=VcfVersion.V4_3)
        .fmt("GT")
        .info("SVLEN")
    )
    # Number=. at 4.3: no cardinality enforcement; single value for single ALT is fine.
    b.record("chr1", 10, ref="A", alt=[Sym("INS")], gt=["0/1"], info={"SVLEN": [30]})
    assert b.build().render().startswith("##fileformat=VCFv4.3\n")


def test_symbolic_del_no_svclaim_required_before_4_4():
    from vcfixture import Sym, VcfBuilder, VcfVersion

    # DEL requires SVCLAIM at >= 4.4, but SVCLAIM does not exist pre-4.4, so the
    # requirement must not apply there.
    b = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 100000)], version=VcfVersion.V4_3)
        .fmt("GT")
        .info("SVLEN")
    )
    b.record("chr1", 10, ref="A", alt=[Sym("DEL")], gt=["0/1"], info={"SVLEN": [-200]})
    assert "DEL" in b.build().render()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_build.py -q`
Expected: FAIL — `VcfBuilder.__init__() got an unexpected keyword argument 'version'`.

- [ ] **Step 3: Thread version through the builder**

In `src/vcfixture/build.py`:

Add the import after line 10 (`from ._spec.types import Type`):

```python
from ._spec.version import VcfVersion
```

Change `__init__` (lines 48–56). Replace the param and the stored attribute:

```python
    def __init__(
        self,
        samples: Iterable[str],
        contigs: Iterable[tuple[str, int | None]],
        version: VcfVersion = VcfVersion.V4_5,
    ):
        self._samples = tuple(samples)
        self._contigs = tuple(ContigDef(c[0], c[1]) for c in contigs)
        self._version = version
```

Convert `_make_def` from a staticmethod to an instance method so it can pass the version. Change lines 91–101 from:

```python
    @staticmethod
    def _make_def(
        id: str,
        number: Number | None,
        type: Type | None,
        description: str | None,
        kind: str,
    ) -> FieldDef:
        if number is None or type is None:
            try:
                return reserved(id, kind)
```

to:

```python
    def _make_def(
        self,
        id: str,
        number: Number | None,
        type: Type | None,
        description: str | None,
        kind: str,
    ) -> FieldDef:
        if number is None or type is None:
            try:
                return reserved(id, kind, self._version)
```

(The two callers `self._make_def(...)` in `info()`/`fmt()` are already instance-method calls — no change needed there.)

Convert `_validate_alleles` from a staticmethod to an instance method and version-gate the SVCLAIM-*required* rule. Change lines 109–114 from:

```python
    @staticmethod
    def _validate_alleles(
        ref: str,
        alts: tuple[Allele, ...],
        info: Mapping[str, object] | None,
    ) -> None:
```

to:

```python
    def _validate_alleles(
        self,
        ref: str,
        alts: tuple[Allele, ...],
        info: Mapping[str, object] | None,
    ) -> None:
```

Then change the SVCLAIM-required branch (lines 138–139) from:

```python
                if a.first_type in _SVCLAIM_REQUIRED and cl is None:
                    raise ValueError(f"SVCLAIM required for {a.render()} (D/J/DJ)")
```

to:

```python
                if (
                    self._version >= VcfVersion.V4_4
                    and a.first_type in _SVCLAIM_REQUIRED
                    and cl is None
                ):
                    raise ValueError(f"SVCLAIM required for {a.render()} (D/J/DJ)")
```

(The call site at line 161, `self._validate_alleles(ref, alts, info)`, is already an instance call — no change. The SVCLAIM *value*-validity branch above it needs no gate: a non-`None` `cl` can only arrive if `SVCLAIM` was declared via `.info()`, which `reserved()` already rejects pre-4.4.)

Finally, in `build()` change line 257 from:

```python
            fileformat=self._fileformat,
```

to:

```python
            version=self._version,
```

- [ ] **Step 4: Run tests + full suite + type-check**

Run: `uv run pytest tests/test_build.py -q`
Expected: PASS.

Run: `uv run pyrefly check`
Expected: PASS (0 errors — the `build.py` errors from Task 3 are now resolved).

Run: `uv run pytest -q`
Expected: PASS (full suite green; strategies still default to V4_5).

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/build.py tests/test_build.py
git commit -m "feat: VcfBuilder version= param gates reserved fields and SVLEN cardinality"
```

---

## Task 5: Strategies `version=` parameter

**Files:**
- Modify: `src/vcfixture/strategies.py`
- Test: `tests/test_strategies.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_strategies.py` (it imports `from vcfixture import strategies as S`; add `from vcfixture import VcfVersion` and ensure `from hypothesis import given, settings` and `from hypothesis import strategies as st` are present — they are used elsewhere in the file):

```python
@settings(max_examples=10, deadline=None)
@given(st.data())
def test_documents_emit_requested_version_header(data):
    for v in VcfVersion:
        doc = data.draw(S.documents(version=v))
        assert doc.render().startswith(f"##fileformat={v.value}\n")
        assert doc.version is v


@settings(max_examples=15, deadline=None)
@given(st.data())
def test_symbolic_documents_version_correct_svlen_sign(data):
    # At <= 4.3, deletions carry a negative SVLEN and no SVCLAIM is declared.
    doc = data.draw(S.symbolic_documents(version=VcfVersion.V4_3))
    assert doc.version is VcfVersion.V4_3
    info_ids = {d.id for d in doc.info_defs}
    assert "SVCLAIM" not in info_ids
    for rec in doc.records:
        svlen = rec.info.get("SVLEN")
        if svlen is None:
            continue
        vals = svlen if isinstance(svlen, (list, tuple)) else [svlen]
        for alt, val in zip(rec.alts, vals):
            if getattr(alt, "first_type", None) == "DEL":
                assert val < 0


@settings(max_examples=10, deadline=None)
@given(st.data())
def test_symbolic_documents_declare_svclaim_at_4_4(data):
    doc = data.draw(S.symbolic_documents(version=VcfVersion.V4_4))
    info_ids = {d.id for d in doc.info_defs}
    assert "SVCLAIM" in info_ids
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_strategies.py -q -k "version"`
Expected: FAIL — `documents() got an unexpected keyword argument 'version'`.

- [ ] **Step 3: Thread version into the strategies**

In `src/vcfixture/strategies.py`:

Add the import after line 8 (`from ._spec.types import Type`):

```python
from ._spec.version import VcfVersion
```

**`documents`** (def at line 339): add a `version` keyword and pass it to both the reference delegation and the reference-free builder. Change the signature (lines 339–348) to add `version: VcfVersion = VcfVersion.V4_5` to the keyword-only block:

```python
def documents(
    draw: DrawFn,
    max_samples: int = 3,
    max_records: int = 4,
    max_alt: int = 1,
    *,
    reference: ReferenceSpec | None = None,
    violations: frozenset[str] = frozenset(),
    label_overrides: dict[str, str] | None = None,
    version: VcfVersion = VcfVersion.V4_5,
) -> VcfDocument:
```

Change the reference delegation (lines 387–392) to forward `version`:

```python
    if reference is not None:
        return draw(
            _reference_documents(
                reference, violations, label_overrides, max_samples, max_records, version
            )
        )
```

Change the reference-free builder line (line 397) from:

```python
    b = VcfBuilder(samples=samples, contigs=[("chr1", 100000)]).fmt("GT")
```

to:

```python
    b = VcfBuilder(samples=samples, contigs=[("chr1", 100000)], version=version).fmt("GT")
```

**`_reference_documents`** (def at line 236): add a trailing `version` parameter and use it in its builder. Change the signature (lines 236–243) to add `version: VcfVersion,` before the return annotation:

```python
def _reference_documents(
    draw: DrawFn,
    reference: ReferenceSpec,
    violations: frozenset[str],
    label_overrides: dict[str, str] | None,
    max_samples: int,
    max_records: int,
    version: VcfVersion,
) -> VcfDocument:
```

Change its builder (lines 261–264) from:

```python
    b = VcfBuilder(
        samples=samples,
        contigs=[(cid, reference.length(cid)) for cid, _ in reference.contigs],
    ).fmt("GT")
```

to:

```python
    b = VcfBuilder(
        samples=samples,
        contigs=[(cid, reference.length(cid)) for cid, _ in reference.contigs],
        version=version,
    ).fmt("GT")
```

**`symbolic_documents`** (def at line 421): add `version`, declare `SVLEN` via the reserved registry (so its definition is version-correct), declare `SVCLAIM` only at >= 4.4, and emit version-correct SVLEN sign + SVCLAIM. Change the signature (lines 421–423) from:

```python
def symbolic_documents(
    draw: DrawFn, max_samples: int = 3, max_records: int = 4
) -> VcfDocument:
```

to:

```python
def symbolic_documents(
    draw: DrawFn,
    max_samples: int = 3,
    max_records: int = 4,
    version: VcfVersion = VcfVersion.V4_5,
) -> VcfDocument:
```

Change the builder construction (lines 429–434) from:

```python
    b = (
        VcfBuilder(samples=samples, contigs=[("chr1", 1_000_000)])
        .fmt("GT")
        .info("SVLEN", Number.A, Type.INTEGER)
        .info("SVCLAIM", Number.A, Type.STRING)
    )
```

to:

```python
    has_svclaim = version >= VcfVersion.V4_4
    b = (
        VcfBuilder(samples=samples, contigs=[("chr1", 1_000_000)], version=version)
        .fmt("GT")
        .info("SVLEN")  # version-correct Number/description from the registry
    )
    if has_svclaim:
        b.info("SVCLAIM")
```

Change the per-record symbolic branch (lines 443–447) from:

```python
        else:
            svlen = draw(st.integers(1, 1000))
            claim = draw(st.sampled_from(sorted(_SVCLAIM_RULES[kind])))
            alts = [SymbolicAllele(kind)]
            info = {"SVLEN": [svlen], "SVCLAIM": [claim]}
```

to:

```python
        else:
            magnitude = draw(st.integers(1, 1000))
            # <= 4.3: SVLEN is the signed REF/ALT length difference (DEL negative).
            # >= 4.4: SVLEN is an unsigned length.
            svlen = -magnitude if (kind == "DEL" and version < VcfVersion.V4_4) else magnitude
            alts = [SymbolicAllele(kind)]
            info = {"SVLEN": [svlen]}
            if has_svclaim:
                claim = draw(st.sampled_from(sorted(_SVCLAIM_RULES[kind])))
                info["SVCLAIM"] = [claim]
```

Note: after this change `Number` / `Type` may become unused imports in `strategies.py` — check and remove them from the import block (lines 7–8) only if no other reference remains (`_build_number_type_combos` and `_matrix_field_defs` still use both, so they will remain used; verify with `uv run ruff check`).

**`reference_and_documents`** (def at line 455): add `version` and forward it to `documents`. Add `version: VcfVersion = VcfVersion.V4_5,` to the keyword-only block (after `max_repeats`), then change the `documents(...)` call (lines 474–482) to pass `version=version`:

```python
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
```

- [ ] **Step 4: Run tests + full suite + lint + type-check**

Run: `uv run pytest tests/test_strategies.py -q`
Expected: PASS.

Run: `uv run ruff check`
Expected: PASS (fix any now-unused import it flags, then re-run).

Run: `uv run pyrefly check`
Expected: PASS.

Run: `uv run pytest -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/strategies.py tests/test_strategies.py
git commit -m "feat: version= on document strategies with version-correct SVLEN/SVCLAIM"
```

---

## Task 6: Cross-parser round-trip across all versions

**Files:**
- Modify: `tests/test_roundtrip.py`
- Modify: `tests/test_symbolic_roundtrip.py`

This keeps the independent-oracle cross-check (serialize → cyvcf2 parse → assert vs `GroundTruth`) intact for every version, using `@parametrize` for the version axis per the project convention.

- [ ] **Step 1: Parametrize the genotype round-trip across versions**

In `tests/test_roundtrip.py`, add `from vcfixture import VcfVersion` to the imports and `from hypothesis import strategies as st` (the file already imports `given, settings`). Replace the decorator + signature (lines 24–26):

```python
@settings(max_examples=75, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(S.documents())
def test_genotypes_round_trip_through_cyvcf2(doc):
    truth = doc.truth()
```

with:

```python
@pytest.mark.parametrize("version", list(VcfVersion))
@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(data=st.data())
def test_genotypes_round_trip_through_cyvcf2(version, data):
    doc = data.draw(S.documents(version=version))
    truth = doc.truth()
```

(`max_examples` is reduced from 75 to 25 because the test now runs once per version — 5× the examples overall. The rest of the test body is unchanged.)

- [ ] **Step 2: Parametrize the symbolic round-trip across versions**

In `tests/test_symbolic_roundtrip.py`, add `from vcfixture import VcfVersion` and `from hypothesis import strategies as st` to the imports. Replace the decorator + signature (lines 14–16):

```python
@settings(max_examples=50, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(S.symbolic_documents())
def test_symbolic_alts_round_trip_through_cyvcf2(doc):
    truth = doc.truth()
```

with:

```python
@pytest.mark.parametrize("version", list(VcfVersion))
@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(data=st.data())
def test_symbolic_alts_round_trip_through_cyvcf2(version, data):
    doc = data.draw(S.symbolic_documents(version=version))
    truth = doc.truth()
```

The existing body already compares `abs(int(got_i)) == at.svlen`, so the signed-vs-unsigned SVLEN difference across versions is handled.

- [ ] **Step 3: Run the round-trip tests**

Run: `uv run pytest tests/test_roundtrip.py tests/test_symbolic_roundtrip.py -q`
Expected: PASS (each test now reports 5 parametrized cases). If cyvcf2 is not installed they are skipped — in that case run `uv run pip list | grep cyvcf2` to confirm, and note the skip rather than treating it as a pass.

- [ ] **Step 4: Full verification**

Run: `uv run pytest -q`
Expected: PASS (full suite).

Run: `uv run ruff check && uv run ruff format --check && uv run pyrefly check`
Expected: all PASS.

Run: `uv run pytest tests/test_benchmark.py -q`
Expected: PASS (benchmark runs at the default version; budget unchanged).

- [ ] **Step 5: Commit**

```bash
git add tests/test_roundtrip.py tests/test_symbolic_roundtrip.py
git commit -m "test: round-trip documents through cyvcf2 across all VCF versions"
```

---

## Final verification (after all tasks)

- [ ] Run `uv run prek run --all-files` — all hooks pass against the whole tree.
- [ ] Run `uv run pytest -q` — full suite green.
- [ ] Confirm the public API exposes `VcfVersion`: `uv run python -c "from vcfixture import VcfVersion; print(list(VcfVersion))"`.
- [ ] Sanity render: `uv run python -c "from vcfixture import VcfBuilder, VcfVersion; print(VcfBuilder(samples=['s1'], contigs=[('chr1',1000)], version=VcfVersion.V4_1).fmt('GT').build().render().splitlines()[0])"` → `##fileformat=VCFv4.1`.

---

## Self-Review

**Spec coverage:**
- `VcfVersion` enum, orderable, exported, `LATEST` → Task 1. ✓
- Model `version` field + serialization header → Task 3. ✓
- `VcfBuilder` `version=` replacing `fileformat=` → Task 4. ✓
- Reserved-field ladder (`since`) + out-of-version rejection → Task 2. ✓
- SVLEN version-variant `FieldDef` driving header + count-validation → Tasks 2 (def) + 4 (validation via builder). ✓
- No `truth.py` change (already `abs()`-normalizes) → confirmed; no task touches it. ✓
- Strategies `version=` on `documents`/`reference_and_documents` (+ `symbolic_documents`), construct-not-reject, version-correct SVLEN sign, `@parametrize` coverage → Tasks 5 + 6. ✓
- Reserved `since` table: SVCLAIM=4.4, LEN=4.4, rest 4.1 → Task 2. ✓
- Out of scope (allele-form gating, percent-encoding, localized alleles, runtime .tex parsing) → not implemented; the SVCLAIM-required gate in Task 4 is the minimal coupling needed so pre-4.4 symbolic generation is valid, not allele-form gating. ✓

**Placeholder scan:** No TBD/TODO; every code step shows complete code. The `available()` helper from the spec is explicitly omitted (YAGNI) with rationale. ✓

**Type/name consistency:** `VcfVersion`, `LATEST`, `reserved(id, kind, version)`, `self._version`, `VcfDocument.version`, `doc.version.value`, `_SVLEN_PRE_4_4`, `_SINCE_INFO`/`_SINCE_FORMAT`, strategy param `version=` — used identically across all tasks. Comparison operators (`<`, `>=`, `<=`) are all defined on `VcfVersion` in Task 1. ✓
