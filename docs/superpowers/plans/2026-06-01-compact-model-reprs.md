# Compact Model Reprs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop the 2,000-line Hypothesis "Falsifying example" blowup and give every model dataclass a compact, human-readable repr that shows in plain `repr()`, pytest diffs, and Hypothesis output.

**Architecture:** First fix the root cause — `Number`'s `ClassVar` singletons leak into `__dataclass_fields__`, so Hypothesis's field-walking pretty-printer recurses into them. Moving those annotations into a `TYPE_CHECKING` block removes the leak. Then add a `CompactRepr` mixin (providing `_repr_pretty_`, which Hypothesis checks *before* field-walking) and a hand-written `__repr__` on each model dataclass.

**Tech Stack:** Python 3.10+, frozen dataclasses, Hypothesis (dev dep), pytest, pyrefly (strict), ruff. Run all tools via `uv run`.

---

## Background for the implementer

- `Number` (`src/vcfixture/_spec/number.py`) is a frozen dataclass whose singletons
  `ONE/A/R/G/DOT/FLAG` are declared as `ClassVar[Number]` *in the class body*. With
  `from __future__ import annotations` active, dataclasses stores them in
  `Number.__dataclass_fields__` as `_FIELD_CLASSVAR` entries that still have
  `init=True`. Hypothesis's vendored pretty-printer
  (`hypothesis/vendor/pretty.py`) field-walks `__dataclass_fields__` filtering only
  on `v.init`, so it recurses into every singleton → exponential blowup.
- Standard `repr()` already excludes those entries; only field-walkers explode.
- Hypothesis's pretty-printer checks for an IPython-style `_repr_pretty_(self, p,
  cycle)` method **before** the dataclass field-walk. Implementing it routes
  Hypothesis through our `__repr__`. It needs **no Hypothesis import** (the clean-room
  wheel smoke test must stay green — nothing in `src/` may import a dev-only dep).
- `@dataclass` does NOT overwrite a `__repr__` defined in the class body, but it
  WOULD shadow one inherited from a base class. Therefore: `__repr__` lives in each
  dataclass body; only `_repr_pretty_` lives on the mixin.

Token rules for `Number.__repr__` (note: this differs from `header_str()`, which
returns `"0"` for FLAG and `str(count)` for FIXED):
- `FIXED` → `str(self.count)` (e.g. `Number(2)`)
- `FLAG` → `"FLAG"` (e.g. `Number(FLAG)`)
- otherwise → `self.kind.value` (`A`, `R`, `G`, `.`)

---

## File Structure

- `src/vcfixture/_repr.py` — **new**. The `CompactRepr` mixin (`_repr_pretty_` only).
- `src/vcfixture/_spec/number.py` — **modify**. TYPE_CHECKING block + mixin + `__repr__`.
- `src/vcfixture/_spec/fielddef.py` — **modify**. Mixin + `__repr__`.
- `src/vcfixture/genotype.py` — **modify**. Mixin + `__repr__`.
- `src/vcfixture/model.py` — **modify**. Mixin + `__repr__` on `ContigDef`, `Record`, `VcfDocument`.
- `src/vcfixture/reference.py` — **modify**. Mixin + `__repr__` on `RepeatFeature`, `ReferenceSpec`.
- `tests/test_repr.py` — **new**. All repr tests (regression guard + per-class formats).

---

## Task 1: Root-cause fix — remove the ClassVar recursion in `Number`

**Files:**
- Modify: `src/vcfixture/_spec/number.py:1-30`
- Test: `tests/test_repr.py` (new)

- [ ] **Step 1: Write the failing regression test**

Create `tests/test_repr.py`:

```python
from hypothesis.vendor.pretty import pretty

from vcfixture._spec.number import Number


def test_number_singletons_not_in_dataclass_fields():
    # ClassVar singletons must not leak into the dataclass field set, or any
    # field-walking pretty-printer (Hypothesis) recurses into them forever.
    assert set(Number.__dataclass_fields__) == {"kind", "count"}


def test_number_pretty_does_not_recurse():
    # Hypothesis's vendored pretty-printer is what produces "Falsifying example"
    # output. It must not explode Number into its singletons.
    out = pretty(Number.G)
    assert "ONE=Number" not in out
    assert len(out) < 200
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest tests/test_repr.py -v`
Expected: both FAIL — `__dataclass_fields__` contains `ONE/A/R/G/DOT/FLAG`, and `pretty(Number.G)` is thousands of chars containing `ONE=Number`.

- [ ] **Step 3: Move the ClassVar annotations into a `TYPE_CHECKING` block**

In `src/vcfixture/_spec/number.py`, change the import on line 6 from:

```python
from typing import ClassVar
```

to:

```python
from typing import TYPE_CHECKING, ClassVar
```

Then replace the singleton declaration block (current lines 23-30):

```python
    # Canonical singletons, assigned after the class body. Declared as ClassVar
    # so type checkers know these attributes exist on Number.
    ONE: ClassVar[Number]
    A: ClassVar[Number]
    R: ClassVar[Number]
    G: ClassVar[Number]
    DOT: ClassVar[Number]
    FLAG: ClassVar[Number]
```

with (wrap them in `if TYPE_CHECKING:` so they are not collected into
`__dataclass_fields__` at runtime, but type checkers still see them):

```python
    # Canonical singletons, assigned after the class body. Declared under
    # TYPE_CHECKING so type checkers know these attributes exist on Number,
    # while keeping them out of __dataclass_fields__ at runtime (otherwise
    # field-walking pretty-printers recurse into them — see tests/test_repr.py).
    if TYPE_CHECKING:
        ONE: ClassVar[Number]
        A: ClassVar[Number]
        R: ClassVar[Number]
        G: ClassVar[Number]
        DOT: ClassVar[Number]
        FLAG: ClassVar[Number]
```

The runtime assignments on lines 61-66 (`Number.ONE = Number(...)` etc.) stay unchanged.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest tests/test_repr.py -v`
Expected: both PASS.

- [ ] **Step 5: Verify type-checking still passes**

Run: `uv run pyrefly check`
Expected: 0 errors (the `TYPE_CHECKING` ClassVar declarations are recognized).

- [ ] **Step 6: Run the existing Number tests (no regressions)**

Run: `uv run pytest tests/test_number.py -v`
Expected: all PASS.

- [ ] **Step 7: Commit**

```bash
git add src/vcfixture/_spec/number.py tests/test_repr.py
git commit -m "fix: stop Number ClassVar singletons leaking into dataclass fields

The ClassVar singletons (ONE/A/R/G/DOT/FLAG) were recorded in
Number.__dataclass_fields__ with init=True, so Hypothesis's pretty-printer
recursed into them and produced ~2k-line Falsifying examples. Moving the
annotations under TYPE_CHECKING keeps type info while dropping them from
the runtime field set.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Add the `CompactRepr` mixin and `Number.__repr__`

**Files:**
- Create: `src/vcfixture/_repr.py`
- Modify: `src/vcfixture/_spec/number.py`
- Test: `tests/test_repr.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_repr.py`:

```python
def test_number_repr_compact():
    assert repr(Number.G) == "Number(G)"
    assert repr(Number.A) == "Number(A)"
    assert repr(Number.R) == "Number(R)"
    assert repr(Number.DOT) == "Number(.)"
    assert repr(Number.ONE) == "Number(1)"
    assert repr(Number.fixed(2)) == "Number(2)"
    assert repr(Number.FLAG) == "Number(FLAG)"


def test_number_pretty_uses_compact_repr():
    # _repr_pretty_ must route Hypothesis's printer through __repr__.
    assert pretty(Number.G) == "Number(G)"
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/test_repr.py -k "compact or pretty_uses" -v`
Expected: FAIL — repr is the default `Number(kind=<NumberKind.G: 'G'>, count=None)`.

- [ ] **Step 3: Create the mixin**

Create `src/vcfixture/_repr.py`:

```python
from __future__ import annotations

from typing import Any


class CompactRepr:
    """Mixin that routes pretty-printers through ``__repr__``.

    Hypothesis's vendored pretty-printer (and IPython's) honor the
    ``_repr_pretty_`` protocol *before* field-walking a dataclass. Implementing
    it here makes our compact ``__repr__`` show up in Falsifying examples too.
    Subclasses must define their own ``__repr__`` in the class body (``@dataclass``
    would shadow an inherited one).
    """

    def _repr_pretty_(self, p: Any, cycle: bool) -> None:
        p.text(repr(self))
```

- [ ] **Step 4: Make `Number` use the mixin and define `__repr__`**

In `src/vcfixture/_spec/number.py`, add the import near the top (after the existing
`from typing import ...` line):

```python
from .._repr import CompactRepr
```

Change the class declaration from:

```python
@dataclass(frozen=True)
class Number:
```

to:

```python
@dataclass(frozen=True)
class Number(CompactRepr):
```

Add this `__repr__` method inside the class (e.g. directly after `header_str`):

```python
    def __repr__(self) -> str:
        if self.kind is NumberKind.FIXED:
            tok = str(self.count)
        elif self.kind is NumberKind.FLAG:
            tok = "FLAG"
        else:
            tok = self.kind.value
        return f"Number({tok})"
```

- [ ] **Step 5: Run to verify pass**

Run: `uv run pytest tests/test_repr.py -v`
Expected: all PASS.

- [ ] **Step 6: Type-check**

Run: `uv run pyrefly check`
Expected: 0 errors.

- [ ] **Step 7: Commit**

```bash
git add src/vcfixture/_repr.py src/vcfixture/_spec/number.py tests/test_repr.py
git commit -m "feat: add CompactRepr mixin and compact Number repr

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Compact `FieldDef.__repr__`

**Files:**
- Modify: `src/vcfixture/_spec/fielddef.py`
- Test: `tests/test_repr.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_repr.py`:

```python
from vcfixture._spec.fielddef import FieldDef
from vcfixture._spec.types import Type


def test_fielddef_repr_compact():
    gt = FieldDef(id="GT", number=Number.ONE, type=Type.STRING,
                  description="Genotype", kind="FORMAT")
    assert repr(gt) == "FieldDef(GT FORMAT Number=1 Type=String)"

    dp = FieldDef(id="DP", number=Number.ONE, type=Type.INTEGER,
                  description="Depth", kind="INFO")
    assert repr(dp) == "FieldDef(DP INFO Number=1 Type=Integer)"
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/test_repr.py -k fielddef -v`
Expected: FAIL — default dataclass repr including the description.

- [ ] **Step 3: Implement**

In `src/vcfixture/_spec/fielddef.py`, add the import after the existing imports
(after line 7, `from .types import Type`):

```python
from .._repr import CompactRepr
```

Change the class declaration from:

```python
@dataclass(frozen=True)
class FieldDef:
```

to:

```python
@dataclass(frozen=True)
class FieldDef(CompactRepr):
```

Add this method inside the class (e.g. after `__post_init__`):

```python
    def __repr__(self) -> str:
        return (
            f"FieldDef({self.id} {self.kind} "
            f"Number={self.number.header_str()} Type={self.type.value})"
        )
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_repr.py -k fielddef -v`
Expected: PASS.

- [ ] **Step 5: Type-check**

Run: `uv run pyrefly check`
Expected: 0 errors.

- [ ] **Step 6: Commit**

```bash
git add src/vcfixture/_spec/fielddef.py tests/test_repr.py
git commit -m "feat: compact FieldDef repr

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Compact `Genotype.__repr__`

**Files:**
- Modify: `src/vcfixture/genotype.py`
- Test: `tests/test_repr.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_repr.py`:

```python
from vcfixture.genotype import Genotype


def test_genotype_repr_compact():
    assert repr(Genotype((0, 1), (True,))) == "Genotype(0|1)"
    assert repr(Genotype((0, 1), (False,))) == "Genotype(0/1)"
    assert repr(Genotype((None, None), (False,))) == "Genotype(./.)"
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/test_repr.py -k genotype -v`
Expected: FAIL — default dataclass repr.

- [ ] **Step 3: Implement**

In `src/vcfixture/genotype.py`, add the import after the existing imports
(after line 6, `_TOKEN_RE = ...` — place the import above that constant, after the
`from dataclasses import dataclass` line):

```python
from ._repr import CompactRepr
```

Change the class declaration from:

```python
@dataclass(frozen=True)
class Genotype:
```

to:

```python
@dataclass(frozen=True)
class Genotype(CompactRepr):
```

Add this method inside the class (e.g. after `render`):

```python
    def __repr__(self) -> str:
        return f"Genotype({self.render()})"
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_repr.py -k genotype -v`
Expected: PASS.

- [ ] **Step 5: Type-check**

Run: `uv run pyrefly check`
Expected: 0 errors.

- [ ] **Step 6: Commit**

```bash
git add src/vcfixture/genotype.py tests/test_repr.py
git commit -m "feat: compact Genotype repr

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Compact reprs for `ContigDef`, `Record`, `VcfDocument`

**Files:**
- Modify: `src/vcfixture/model.py`
- Test: `tests/test_repr.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_repr.py`:

```python
from vcfixture.model import ContigDef, Record, VcfDocument


def test_contigdef_repr_compact():
    assert repr(ContigDef(id="chr1", length=200)) == "ContigDef(chr1:200)"
    assert repr(ContigDef(id="chr1")) == "ContigDef(chr1)"


def _make_record(alts=("T", "G"), labels=frozenset()):
    return Record(
        chrom="chr1", pos=5, ids=None, ref="A", alts=alts, qual=None,
        filters=None, info={}, fmt_keys=("GT",),
        samples=(
            {"GT": Genotype((0, 1), (True,))},
            {"GT": Genotype((1, 1), (True,))},
        ),
        labels=labels,
    )


def test_record_repr_compact():
    rec = _make_record(labels=frozenset({"multiallelic", "snp"}))
    assert repr(rec) == "Record(chr1:5 A>T,G ×2 [multiallelic,snp])"


def test_record_repr_no_labels():
    rec = _make_record(labels=frozenset())
    assert repr(rec) == "Record(chr1:5 A>T,G ×2)"


def test_record_repr_empty_alts():
    rec = _make_record(alts=(), labels=frozenset())
    assert repr(rec) == "Record(chr1:5 A>. ×2)"


def test_vcfdocument_repr_compact():
    gt = FieldDef(id="GT", number=Number.ONE, type=Type.STRING,
                  description="Genotype", kind="FORMAT")
    doc = VcfDocument(
        fileformat="VCFv4.5", info_defs=(), format_defs=(gt,),
        filter_defs=(), contigs=(), samples=("s0", "s1"),
        records=(_make_record(), _make_record()),
    )
    assert repr(doc) == "VcfDocument(VCFv4.5 samples=2 records=2 info=0 format=1)"
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/test_repr.py -k "contigdef or record or vcfdocument" -v`
Expected: FAIL — default dataclass reprs.

- [ ] **Step 3: Implement**

In `src/vcfixture/model.py`, add the import after the existing imports
(after line 10, `from .genotype import Genotype`):

```python
from ._repr import CompactRepr
```

Change `class ContigDef:` (line 17) to `class ContigDef(CompactRepr):` and add inside it:

```python
    def __repr__(self) -> str:
        if self.length is None:
            return f"ContigDef({self.id})"
        return f"ContigDef({self.id}:{self.length})"
```

Change `class Record:` (line 28) to `class Record(CompactRepr):` and add inside it
(e.g. after the `n_alt` property):

```python
    def __repr__(self) -> str:
        alts = ",".join(self.alts) if self.alts else "."
        out = f"Record({self.chrom}:{self.pos} {self.ref}>{alts} ×{len(self.samples)}"
        if self.labels:
            out += f" [{','.join(sorted(self.labels))}]"
        return out + ")"
```

Change `class VcfDocument:` (line 47) to `class VcfDocument(CompactRepr):` and add
inside it (e.g. after `max_ploidy`):

```python
    def __repr__(self) -> str:
        return (
            f"VcfDocument({self.fileformat} samples={len(self.samples)} "
            f"records={len(self.records)} info={len(self.info_defs)} "
            f"format={len(self.format_defs)})"
        )
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_repr.py -k "contigdef or record or vcfdocument" -v`
Expected: all PASS.

- [ ] **Step 5: Type-check**

Run: `uv run pyrefly check`
Expected: 0 errors.

- [ ] **Step 6: Commit**

```bash
git add src/vcfixture/model.py tests/test_repr.py
git commit -m "feat: compact ContigDef/Record/VcfDocument reprs

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Compact reprs for `RepeatFeature` and `ReferenceSpec`

**Files:**
- Modify: `src/vcfixture/reference.py`
- Test: `tests/test_repr.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_repr.py`:

```python
from vcfixture.reference import ReferenceSpec, RepeatFeature


def test_repeatfeature_repr_compact():
    rf = RepeatFeature(contig="chr1", pos0=10, motif="AT", count=3)
    assert repr(rf) == "RepeatFeature(chr1@10 AT×3)"


def test_referencespec_repr_compact():
    spec = ReferenceSpec(
        contigs=(("chr1", "A" * 200),),
        repeats=(RepeatFeature(contig="chr1", pos0=10, motif="AT", count=3),),
    )
    assert repr(spec) == "ReferenceSpec(contigs=[chr1:200bp], repeats=1)"


def test_referencespec_repr_pretty_no_full_sequence():
    spec = ReferenceSpec(contigs=(("chr1", "ACGT" * 50),), repeats=())
    out = pretty(spec)
    assert out == "ReferenceSpec(contigs=[chr1:200bp], repeats=0)"
    assert "ACGTACGT" not in out
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/test_repr.py -k "repeatfeature or referencespec" -v`
Expected: FAIL — default reprs dump the full 200-char sequence.

- [ ] **Step 3: Implement**

In `src/vcfixture/reference.py`, add the import after the existing imports
(after line 10, `from ._typing import StrPath`):

```python
from ._repr import CompactRepr
```

Change `class RepeatFeature:` (line 91) to `class RepeatFeature(CompactRepr):` and
add inside it (e.g. after the `length` property):

```python
    def __repr__(self) -> str:
        return f"RepeatFeature({self.contig}@{self.pos0} {self.motif}×{self.count})"
```

Change `class ReferenceSpec:` (line 105) to `class ReferenceSpec(CompactRepr):` and
add inside it (e.g. after `_seq_for`):

```python
    def __repr__(self) -> str:
        contigs = ", ".join(f"{cid}:{len(seq)}bp" for cid, seq in self.contigs)
        return f"ReferenceSpec(contigs=[{contigs}], repeats={len(self.repeats)})"
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_repr.py -k "repeatfeature or referencespec" -v`
Expected: all PASS.

- [ ] **Step 5: Type-check**

Run: `uv run pyrefly check`
Expected: 0 errors.

- [ ] **Step 6: Commit**

```bash
git add src/vcfixture/reference.py tests/test_repr.py
git commit -m "feat: compact RepeatFeature/ReferenceSpec reprs

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 7: End-to-end regression guard — bounded `VcfDocument` pretty output

**Files:**
- Test: `tests/test_repr.py`

- [ ] **Step 1: Write the test**

Append to `tests/test_repr.py`:

```python
def test_full_document_pretty_is_bounded():
    # The original bug: a VcfDocument printed ~2k lines via Hypothesis's
    # pretty-printer. With compact reprs + _repr_pretty_, the whole document
    # collapses to a single short line.
    gt = FieldDef(id="GT", number=Number.ONE, type=Type.STRING,
                  description="Genotype", kind="FORMAT")
    doc = VcfDocument(
        fileformat="VCFv4.5", info_defs=(), format_defs=(gt,),
        filter_defs=(), contigs=(ContigDef("chr1", 200),),
        samples=("s0", "s1"),
        records=tuple(_make_record() for _ in range(5)),
    )
    out = pretty(doc)
    assert out == "VcfDocument(VCFv4.5 samples=2 records=5 info=0 format=1)"
    assert "ONE=Number" not in out
    assert len(out) < 200
```

- [ ] **Step 2: Run to verify pass**

Run: `uv run pytest tests/test_repr.py::test_full_document_pretty_is_bounded -v`
Expected: PASS (everything it depends on is already implemented).

- [ ] **Step 3: Run the full suite + lint + type-check**

Run:
```bash
uv run pytest
uv run ruff check
uv run ruff format
uv run pyrefly check
```
Expected: all green; `ruff format` reports no changes (or reformats and you re-stage).

- [ ] **Step 4: Commit**

```bash
git add tests/test_repr.py
git commit -m "test: end-to-end guard that VcfDocument pretty output stays bounded

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Self-Review notes

- **Spec coverage:** root-cause fix (Task 1); `CompactRepr` mixin (Task 2);
  compact `__repr__` for Number (T2), FieldDef (T3), Genotype (T4), ContigDef/Record/VcfDocument (T5),
  RepeatFeature/ReferenceSpec (T6); regression + per-class + edge-case tests
  (empty labels, empty alts, `ContigDef` length None) across T1–T7. All spec
  table rows and the test section are covered.
- **No runtime dependency added:** `_repr.py` imports only `typing`; the
  `hypothesis.vendor.pretty` import lives only in `tests/`. Clean-room wheel test unaffected.
- **Type consistency:** mixin method `_repr_pretty_(self, p: Any, cycle: bool)`;
  every dataclass defines `__repr__` in its own body (not inherited) so `@dataclass`
  does not shadow it; `Number.__repr__` token rules differ from `header_str()`
  (FLAG→`FLAG`, FIXED→count) and the tests assert exactly that.
```
