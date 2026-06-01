# Reference-consistent `documents` strategy + variant labels — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reference-consistent variant-generation capability to vcfixture: a mutable `ReferenceBuilder` → frozen `ReferenceSpec`, a general per-variant `labels` mechanism, and composable Hypothesis strategies (`references()`, `documents(reference=…)`, `reference_and_documents()`) that emit deliberately non-canonical but reference-consistent VCFs for downstream tools (GenVarLoader) to test against.

**Architecture:** `ReferenceBuilder` random-fills contigs and supports multi-nucleotide / tandem-repeat writes; `.build()` yields a frozen `ReferenceSpec` (advertising planted repeats) that can draw spec-correct REF/ALT and `write()` a bgzipped+faidx'd FASTA. `Record`/`GroundTruth` gain a domain-agnostic `labels: frozenset[str]`. Strategies draw a `ReferenceSpec`, then draw records whose REFs match it, opting into violation classes (`multiallelic`, `non_atomic`, `non_left_aligned`) that auto-tag descriptive provenance labels. vcfixture never implements normalization and its own tests never invoke bcftools.

**Tech Stack:** Python ≥3.10, `numpy`, `pysam`, `hypothesis`, `pytest`, `uv`, `ruff`, `pyrefly` (strict on `src`), `commitizen`.

**Scope:** vcfixture upstream only. The downstream GenVarLoader Phase 2 property-test module that consumes the released `vcfixture>=0.3.0` is NOT part of this plan.

**Reference design:** `docs/superpowers/specs/2026-05-31-reference-consistent-documents-design.md`.

---

## File Structure

- **Modify** `src/vcfixture/model.py` — add `labels: frozenset[str] = frozenset()` to `Record`.
- **Modify** `src/vcfixture/build.py` — `VcfBuilder.record(..., labels=...)`.
- **Modify** `src/vcfixture/truth.py` — `GroundTruth.labels: list[frozenset[str]]`; populate in `derive_truth`.
- **Modify** `src/vcfixture/reference.py` — add `RepeatFeature`, `ReferenceBuilder`, `ReferenceSpec`; factor `draw_ref_alt` into a shared module helper used by both `Reference` and `ReferenceSpec`.
- **Modify** `src/vcfixture/strategies.py` — add `references()`, a `_reference_documents()` composite, a `reference=`/`violations=`/`label_overrides=` branch on `documents()`, and `reference_and_documents()`.
- **Modify** `src/vcfixture/__init__.py` — export `ReferenceBuilder`, `ReferenceSpec`, `RepeatFeature`.
- **Modify** `pyproject.toml` — version bump 0.2.1 → 0.3.0.
- **Modify** `tests/test_reference.py` — builder/spec/write/draw_ref_alt tests.
- **Modify** `tests/test_truth.py` and **add** `tests/test_labels.py` — labels carried through; absent from `render()`.
- **Modify** `tests/test_strategies.py` — reference-consistency, violation labels, paired output, back-compat.

---

## Task 1: Worktree environment + clean baseline

**Files:** none (environment only).

- [ ] **Step 1: Confirm you are in the worktree**

Run:
```bash
git rev-parse --show-toplevel && git branch --show-current
```
Expected: toplevel ends with `.worktrees/ref-consistent-docs`, branch `feat/reference-consistent-documents`. If not, stop — you are in the wrong checkout.

- [ ] **Step 2: Sync deps and run the baseline suite**

Run:
```bash
uv sync --quiet && uv run pytest -q
```
Expected: PASS — `67 passed` (or more), 0 failures. If anything fails, report and stop before changing code.

- [ ] **Step 3: No commit** (environment-only).

---

## Task 2: General `labels` on the variant model

Add a domain-agnostic per-variant label set carried from `VcfBuilder.record` through `Record` into `GroundTruth`, never serialized into the VCF.

**Files:**
- Test: `tests/test_labels.py` (create)
- Modify: `src/vcfixture/model.py` (Record), `src/vcfixture/build.py` (record), `src/vcfixture/truth.py` (GroundTruth + derive_truth)

- [ ] **Step 1: Write the failing test**

Create `tests/test_labels.py`:
```python
from vcfixture import VcfBuilder


def _builder() -> VcfBuilder:
    return VcfBuilder(samples=["s0"], contigs=[("chr1", 1000)]).fmt("GT")


def test_labels_default_empty():
    b = _builder().record("chr1", 10, ref="A", alt=["C"], gt=["0|1"])
    doc = b.build()
    assert doc.records[0].labels == frozenset()
    assert doc.truth().labels == [frozenset()]


def test_labels_carried_record_to_truth():
    b = _builder().record(
        "chr1", 10, ref="A", alt=["C"], gt=["0|1"], labels=["off_anchor", "x"]
    )
    doc = b.build()
    assert doc.records[0].labels == frozenset({"off_anchor", "x"})
    assert doc.truth().labels == [frozenset({"off_anchor", "x"})]


def test_labels_not_serialized():
    b = _builder().record(
        "chr1", 10, ref="A", alt=["C"], gt=["0|1"], labels=["off_anchor"]
    )
    text = b.build().render()
    assert "off_anchor" not in text
```

- [ ] **Step 2: Run it to verify it fails**

Run:
```bash
uv run pytest tests/test_labels.py -q
```
Expected: FAIL — `TypeError: record() got an unexpected keyword argument 'labels'` (and/or `AttributeError` on `.labels`).

- [ ] **Step 3: Add `labels` to `Record`**

In `src/vcfixture/model.py`, add a field at the END of the `Record` dataclass (after `samples`):
```python
    samples: tuple[Mapping[str, Any], ...]  # per-sample: key -> value(s)/Genotype
    labels: frozenset[str] = frozenset()
```

- [ ] **Step 4: Accept `labels` in `VcfBuilder.record`**

In `src/vcfixture/build.py`, change the `record` signature to add `labels` (keyword), keeping all other params:
```python
    def record(
        self,
        chrom: str,
        pos: int,
        *,
        ref: str,
        alt: Sequence[str],
        ids: Iterable[str] | None = None,
        qual: float | None = None,
        filter: Iterable[str] | None = None,
        gt: Sequence[str] | None = None,
        info: Mapping[str, object] | None = None,
        labels: Iterable[str] | None = None,
        **fmt_fields: Sequence[object],
    ) -> VcfBuilder:
```
Then, in the `Record(...)` construction at the end of `record`, add the field after `samples=...`:
```python
                samples=tuple(samples),
                labels=frozenset(labels) if labels else frozenset(),
```

- [ ] **Step 5: Add `labels` to `GroundTruth` and populate it**

In `src/vcfixture/truth.py`, add a field at the END of the `GroundTruth` dataclass (after `format`):
```python
    format: list[list[dict[str, object]]]  # per record, per sample: id -> value(s)
    labels: list[frozenset[str]]  # per record
```
In `derive_truth`, initialize a list near the other per-record lists:
```python
    fmt: list[list[dict[str, object]]] = []
    labels: list[frozenset[str]] = []
```
Append inside the record loop (after `fmt.append(per_sample)`):
```python
        fmt.append(per_sample)
        labels.append(rec.labels)
```
And add to the returned `GroundTruth(...)` (after `format=fmt,`):
```python
        format=fmt,
        labels=labels,
    )
```

- [ ] **Step 6: Run the test to verify it passes**

Run:
```bash
uv run pytest tests/test_labels.py -q
```
Expected: PASS (3 passed).

- [ ] **Step 7: Run the full suite (no regressions)**

Run:
```bash
uv run pytest -q
```
Expected: all prior tests still pass.

- [ ] **Step 8: Commit**

```bash
git add src/vcfixture/model.py src/vcfixture/build.py src/vcfixture/truth.py tests/test_labels.py
git commit -m "feat: general per-variant labels carried into GroundTruth"
```

---

## Task 3: `RepeatFeature` + `ReferenceBuilder` + `ReferenceSpec`

Add the mutable builder and frozen spec (base/seq/repeats/build). `write()` and `draw_ref_alt` come in Tasks 4–5.

**Files:**
- Modify: `src/vcfixture/reference.py`
- Modify: `tests/test_reference.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_reference.py`:
```python
from vcfixture.reference import ReferenceBuilder, ReferenceSpec, RepeatFeature


def test_builder_random_fill_is_seeded_and_acgt():
    a = ReferenceBuilder(seed=0).add_contig("chr1", 200).build()
    b = ReferenceBuilder(seed=0).add_contig("chr1", 200).build()
    assert a.contigs == b.contigs  # deterministic
    seq = a.seq("chr1", 0, 200)
    assert len(seq) == 200 and set(seq) <= set("ACGT")


def test_set_base_and_set_seq_overwrite():
    spec = (
        ReferenceBuilder(seed=1)
        .add_contig("chr1", 100)
        .set_base("chr1", 10, "A")
        .set_seq("chr1", 20, "GATTACA")
        .build()
    )
    assert spec.base("chr1", 10) == "A"
    assert spec.seq("chr1", 20, 7) == "GATTACA"


def test_tandem_repeat_writes_and_records_feature():
    spec = (
        ReferenceBuilder(seed=2)
        .add_contig("chr1", 100)
        .tandem_repeat("chr1", 30, "AG", 5)
        .build()
    )
    assert spec.seq("chr1", 30, 10) == "AGAGAGAGAG"
    assert spec.repeats == (RepeatFeature("chr1", 30, "AG", 5),)
    assert spec.repeats[0].length == 10


def test_set_seq_out_of_bounds_raises():
    import pytest

    rb = ReferenceBuilder(seed=0).add_contig("chr1", 10)
    with pytest.raises(ValueError):
        rb.set_seq("chr1", 8, "ACGT")  # runs past length 10
```

- [ ] **Step 2: Run it to verify it fails**

Run:
```bash
uv run pytest tests/test_reference.py -q
```
Expected: FAIL — `ImportError: cannot import name 'ReferenceBuilder'`.

- [ ] **Step 3: Implement the builder, spec, and feature**

In `src/vcfixture/reference.py`, update the imports at the top to:
```python
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam

from ._typing import StrPath

_BASES = "ACGT"
_BASES_ARR = np.frombuffer(b"ACGT", dtype="S1")
```
Then append the following to the file (after the existing `Reference` class):
```python
@dataclass(frozen=True)
class RepeatFeature:
    """A tandem repeat planted into a reference, for provenance."""

    contig: str
    pos0: int  # 0-based start of the repeat run
    motif: str
    count: int

    @property
    def length(self) -> int:
        return len(self.motif) * self.count


@dataclass(frozen=True)
class ReferenceSpec:
    """Immutable in-memory reference: contig sequences + planted repeats."""

    contigs: tuple[tuple[str, str], ...]  # (id, sequence)
    repeats: tuple[RepeatFeature, ...] = ()

    def _seq_for(self, contig: str) -> str:
        for cid, seq in self.contigs:
            if cid == contig:
                return seq
        raise KeyError(contig)

    def length(self, contig: str) -> int:
        return len(self._seq_for(contig))

    def base(self, contig: str, pos0: int) -> str:
        return self._seq_for(contig)[pos0]

    def seq(self, contig: str, start0: int, length: int) -> str:
        return self._seq_for(contig)[start0 : start0 + length]


class ReferenceBuilder:
    """Mutable builder for a synthetic reference.

    Random-fills contigs (seeded), supports single-base / multi-nucleotide /
    tandem-repeat overwrites, then ``build()`` freezes a ``ReferenceSpec``.
    """

    def __init__(self, seed: int = 0):
        self._rng = np.random.default_rng(seed)
        self._seqs: dict[str, np.ndarray] = {}
        self._order: list[str] = []
        self._repeats: list[RepeatFeature] = []

    def add_contig(self, id: str, length: int) -> ReferenceBuilder:
        if id in self._seqs:
            raise ValueError(f"contig {id!r} already added")
        self._seqs[id] = self._rng.choice(_BASES_ARR, size=length)
        self._order.append(id)
        return self

    def set_base(self, contig: str, pos0: int, base: str) -> ReferenceBuilder:
        if len(base) != 1:
            raise ValueError(f"set_base expects one base, got {base!r}")
        self._seqs[contig][pos0] = base.encode()
        return self

    def set_seq(self, contig: str, pos0: int, seq: str) -> ReferenceBuilder:
        arr = self._seqs[contig]
        if pos0 < 0 or pos0 + len(seq) > arr.size:
            raise ValueError(
                f"set_seq {contig}:{pos0}+{len(seq)} runs past length {arr.size}"
            )
        arr[pos0 : pos0 + len(seq)] = np.frombuffer(seq.encode(), dtype="S1")
        return self

    def tandem_repeat(
        self, contig: str, pos0: int, motif: str, n: int
    ) -> ReferenceBuilder:
        self.set_seq(contig, pos0, motif * n)
        self._repeats.append(RepeatFeature(contig, pos0, motif, n))
        return self

    def build(self) -> ReferenceSpec:
        contigs = tuple(
            (cid, self._seqs[cid].tobytes().decode()) for cid in self._order
        )
        return ReferenceSpec(contigs=contigs, repeats=tuple(self._repeats))
```
(Note: `field` is imported for forward-compat but unused here; if ruff flags it as unused, drop `field` from the import.)

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run pytest tests/test_reference.py -q
```
Expected: PASS (existing `Reference` tests + 4 new tests).

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/reference.py tests/test_reference.py
git commit -m "feat: ReferenceBuilder/ReferenceSpec with tandem-repeat provenance"
```

---

## Task 4: `ReferenceSpec.write` — bgzipped, faidx'd FASTA roundtrip

**Files:**
- Modify: `src/vcfixture/reference.py` (add `write` to `ReferenceSpec`)
- Modify: `tests/test_reference.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_reference.py`:
```python
import pysam as _pysam


def test_referencespec_write_roundtrips(tmp_path):
    spec = (
        ReferenceBuilder(seed=3)
        .add_contig("chr1", 300)
        .add_contig("chr2", 150)
        .set_seq("chr1", 50, "GATTACA")
        .build()
    )
    out = spec.write(tmp_path / "ref.fa.bgz")
    assert out.exists()
    assert (out.parent / (out.name + ".fai")).exists()
    with _pysam.FastaFile(str(out)) as fa:
        assert fa.fetch("chr1", 50, 57).upper() == "GATTACA"
        assert fa.fetch("chr1", 0, 300).upper() == spec.seq("chr1", 0, 300)
        assert fa.references == ("chr1", "chr2") or set(fa.references) == {
            "chr1",
            "chr2",
        }


def test_referencespec_write_plain(tmp_path):
    spec = ReferenceBuilder(seed=4).add_contig("chr1", 80).build()
    out = spec.write(tmp_path / "ref.fa", bgzip=False)
    assert out.exists()
    with _pysam.FastaFile(str(out)) as fa:
        assert fa.fetch("chr1", 0, 80).upper() == spec.seq("chr1", 0, 80)
```

- [ ] **Step 2: Run it to verify it fails**

Run:
```bash
uv run pytest tests/test_reference.py -k write -q
```
Expected: FAIL — `AttributeError: 'ReferenceSpec' object has no attribute 'write'`.

- [ ] **Step 3: Implement `write`**

Add this method to the `ReferenceSpec` dataclass in `src/vcfixture/reference.py` (after `seq`):
```python
    def write(
        self, path: StrPath, *, bgzip: bool = True, index: bool = True
    ) -> Path:
        """Write a 60-col FASTA; bgzip + faidx it via pysam. Returns the path."""
        path = Path(path)
        text_lines: list[str] = []
        for cid, seq in self.contigs:
            text_lines.append(f">{cid}")
            text_lines.extend(seq[i : i + 60] for i in range(0, len(seq), 60))
        fasta_text = "\n".join(text_lines) + "\n"

        if bgzip:
            plain = path.with_name(path.name + ".tmp.fa")
            plain.write_text(fasta_text)
            pysam.tabix_compress(str(plain), str(path), force=True)
            plain.unlink()
        else:
            path.write_text(fasta_text)

        if index:
            pysam.faidx(str(path))  # writes <path>.fai (+ .gzi when bgzipped)
        return path
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run pytest tests/test_reference.py -k write -q
```
Expected: PASS (2 passed). If `pysam.faidx` errors on the bgzipped file, confirm the file is true bgzip (it is — `tabix_compress` writes BGZF); do not switch to plain gzip.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/reference.py tests/test_reference.py
git commit -m "feat: ReferenceSpec.write bgzipped+faidx FASTA"
```

---

## Task 5: Factor shared `_ref_alt`; add `ReferenceSpec.draw_ref_alt`

Both the file-backed `Reference` and the in-memory `ReferenceSpec` need to draw spec-correct REF/ALT. Factor the existing logic into one module helper over `(base_fn, seq_fn)` accessors; both classes delegate. Existing `Reference.draw_ref_alt` behavior is preserved (its tests must stay green).

**Files:**
- Modify: `src/vcfixture/reference.py`
- Modify: `tests/test_reference.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_reference.py`:
```python
def test_referencespec_draw_ref_alt_matches_sequence():
    spec = ReferenceBuilder(seed=5).add_contig("chr1", 100).set_seq(
        "chr1", 10, "ACGT"
    ).build()
    ref, alts = spec.draw_ref_alt("chr1", pos0=10, klass="SNP", alt_index=1)
    assert ref == "A" and len(alts) == 1 and alts[0] != "A"
    dref, dalts = spec.draw_ref_alt("chr1", pos0=10, klass="DEL", del_len=2)
    assert dref == spec.seq("chr1", 10, 3) and dalts == [dref[0]]
```

- [ ] **Step 2: Run it to verify it fails**

Run:
```bash
uv run pytest tests/test_reference.py -k draw_ref_alt -q
```
Expected: FAIL — `AttributeError: 'ReferenceSpec' object has no attribute 'draw_ref_alt'`.

- [ ] **Step 3: Extract the shared helper and delegate**

In `src/vcfixture/reference.py`, add this module-level function (place it after `_BASES_ARR`, before the `Reference` class). It is the body currently inside `Reference.draw_ref_alt`, parameterized over accessors:
```python
from collections.abc import Callable


def _draw_ref_alt(
    base_fn: Callable[[str, int], str],
    seq_fn: Callable[[str, int, int], str],
    contig: str,
    pos0: int,
    klass: str,
    *,
    alt_index: int = 1,
    del_len: int = 1,
    ins_seq: str = "T",
    mnp_len: int = 2,
) -> tuple[str, list[str]]:
    if klass == "SNP":
        r = base_fn(contig, pos0)
        alt = _BASES[(_BASES.index(r) + alt_index) % 4]
        return r, [alt]
    if klass == "MNP":
        r = seq_fn(contig, pos0, mnp_len)
        alt = "".join(_BASES[(_BASES.index(b) + alt_index) % 4] for b in r)
        return r, [alt]
    if klass == "INS":
        anchor = base_fn(contig, pos0)
        return anchor, [anchor + ins_seq]
    if klass == "DEL":
        r = seq_fn(contig, pos0, del_len + 1)
        return r, [r[0]]
    if klass == "DELINS":
        r = seq_fn(contig, pos0, mnp_len)
        return r, [ins_seq]
    if klass == "SPANNING_DEL":
        return base_fn(contig, pos0), ["*"]
    raise ValueError(f"unknown class {klass!r}")
```
Then REPLACE the body of `Reference.draw_ref_alt` with a delegating call (keep its exact signature):
```python
    def draw_ref_alt(
        self,
        contig: str,
        pos0: int,
        klass: str,
        *,
        alt_index: int = 1,
        del_len: int = 1,
        ins_seq: str = "T",
        mnp_len: int = 2,
    ) -> tuple[str, list[str]]:
        return _draw_ref_alt(
            self.base,
            self.seq,
            contig,
            pos0,
            klass,
            alt_index=alt_index,
            del_len=del_len,
            ins_seq=ins_seq,
            mnp_len=mnp_len,
        )
```
And add the same delegating method to `ReferenceSpec` (after `write`):
```python
    def draw_ref_alt(
        self,
        contig: str,
        pos0: int,
        klass: str,
        *,
        alt_index: int = 1,
        del_len: int = 1,
        ins_seq: str = "T",
        mnp_len: int = 2,
    ) -> tuple[str, list[str]]:
        return _draw_ref_alt(
            self.base,
            self.seq,
            contig,
            pos0,
            klass,
            alt_index=alt_index,
            del_len=del_len,
            ins_seq=ins_seq,
            mnp_len=mnp_len,
        )
```

- [ ] **Step 4: Run the reference tests to verify all pass**

Run:
```bash
uv run pytest tests/test_reference.py -q
```
Expected: PASS — existing `Reference.draw_ref_alt` tests (`test_draw_snp_ref_matches_sequence`, `test_draw_deletion_ref_starts_with_sequence`) AND the new `ReferenceSpec` test.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/reference.py tests/test_reference.py
git commit -m "refactor: share draw_ref_alt between Reference and ReferenceSpec"
```

---

## Task 6: `references()` strategy

A Hypothesis strategy that draws a small `ReferenceSpec` with optional non-overlapping planted tandem repeats.

**Files:**
- Modify: `src/vcfixture/strategies.py`
- Modify: `tests/test_strategies.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_strategies.py`:
```python
from hypothesis import given, settings, HealthCheck

from vcfixture.reference import ReferenceSpec


@settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
@given(S.references())
def test_references_are_well_formed(spec: ReferenceSpec):
    assert isinstance(spec, ReferenceSpec)
    assert len(spec.contigs) >= 1
    for cid, seq in spec.contigs:
        assert len(seq) >= 1 and set(seq) <= set("ACGT")
    # planted repeats actually appear at their advertised loci
    for rf in spec.repeats:
        assert spec.seq(rf.contig, rf.pos0, rf.length) == rf.motif * rf.count
```

- [ ] **Step 2: Run it to verify it fails**

Run:
```bash
uv run pytest tests/test_strategies.py -k references -q
```
Expected: FAIL — `AttributeError: module 'vcfixture.strategies' has no attribute 'references'`.

- [ ] **Step 3: Implement `references()`**

In `src/vcfixture/strategies.py`, add to the imports near the top:
```python
from .reference import ReferenceBuilder, ReferenceSpec
```
Then append:
```python
@st.composite
def references(
    draw: DrawFn,
    *,
    max_contigs: int = 2,
    max_contig_len: int = 2000,
    max_repeats: int = 3,
) -> ReferenceSpec:
    """Draw a small reference-consistent ``ReferenceSpec`` with optional,
    non-overlapping planted tandem repeats (advertised on ``spec.repeats``)."""
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run pytest tests/test_strategies.py -k references -q
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/strategies.py tests/test_strategies.py
git commit -m "feat: references() strategy drawing ReferenceSpec with planted repeats"
```

---

## Task 7: reference-consistent `documents(reference=…, violations=…)`

Add a reference-aware branch to `documents()` that draws records whose REFs match the spec, opting into violation classes that auto-tag provenance labels. The reference-free path is untouched (back-compat).

**Files:**
- Modify: `src/vcfixture/strategies.py`
- Modify: `tests/test_strategies.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_strategies.py` (note: `reference_and_documents` is implemented in Task 8, but these tests only need `documents(reference=…)`; they draw a reference via `references()` and compose inline):
```python
@settings(
    max_examples=40,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much],
)
@given(st.data())
def test_reference_consistent_and_labeled(data):
    spec = data.draw(S.references(max_repeats=3))
    doc = data.draw(
        S.documents(
            reference=spec,
            violations=frozenset({"multiallelic", "non_atomic", "non_left_aligned"}),
        )
    )
    truth = doc.truth()
    # Every REF matches the reference sequence at its position.
    for rec in doc.records:
        assert spec.seq(rec.chrom, rec.pos - 1, len(rec.ref)) == rec.ref
    # Records are position-sorted per contig (norm/consensus/gvl require this).
    last: dict[str, int] = {}
    for rec in doc.records:
        assert rec.pos >= last.get(rec.chrom, 0)
        last[rec.chrom] = rec.pos
    # truth lines up with the document.
    assert truth.genotypes.shape[0] == len(doc.records)
    assert truth.labels == [r.labels for r in doc.records]
    # Provenance labels only use the known vocabulary.
    allowed = {"multiallelic", "non_atomic", "off_anchor", "tandem_repeat"}
    for lbls in truth.labels:
        assert lbls <= allowed


@settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow])
@given(S.documents())
def test_documents_back_compat_unlabeled(doc):
    # Reference-free documents still work and carry no labels.
    assert all(r.labels == frozenset() for r in doc.records)
```
This test also needs `st` in scope — add `from hypothesis import strategies as st` to the imports at the top of `tests/test_strategies.py` if it is not already there.

- [ ] **Step 2: Run it to verify it fails**

Run:
```bash
uv run pytest tests/test_strategies.py -k "reference_consistent or back_compat" -q
```
Expected: FAIL — `documents()` has no `reference`/`violations` kwargs and `reference_and_documents` does not exist yet.

- [ ] **Step 3: Implement the reference-consistent record helper and the `documents` branch**

In `src/vcfixture/strategies.py`, add this composite (it does the reference-aware drawing):
```python
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
    ).fmt("GT")

    enabled = [v for v in ("multiallelic", "non_atomic", "non_left_aligned") if v in violations]

    n_rec = draw(st.integers(1, max_records))
    cursor = 10  # low start so planted repeats (pos0 >= 50) are reachable ahead
    for _ in range(n_rec):
        if cursor + 30 >= clen:
            break
        # Decide what kind of record to emit. Anchor `a` defaults to the cursor;
        # only a forward jump (off_anchor) may move it ahead, never behind, so
        # records stay position-sorted.
        a = cursor
        kind = draw(st.sampled_from(["canonical", *enabled])) if enabled else "canonical"

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
            if want_del and a + 2 < clen and reference.base(
                contig, a
            ) != reference.base(contig, a + 1):
                ref = reference.seq(contig, a, 2)  # delete the differing next base
                alts = [ref[0]]
            else:
                ref, alts = reference.draw_ref_alt(contig, a, "SNP")

        gts = [draw(genotypes(ploidy, n_alt=len(alts))) for _ in samples]
        b.record(
            contig,
            a + 1,  # 1-based POS
            ref=ref,
            alt=alts,
            gt=gts,
            labels=sorted(labels) if labels else None,
        )
        cursor = a + len(ref) + draw(st.integers(20, 60))

    return b.build()
```
Then update the existing `documents` strategy. (1) Replace ONLY its `def documents(...)`/signature lines (currently `def documents(draw, max_samples=3, max_records=4, max_alt=1) -> VcfDocument:`) with the expanded signature below. (2) Insert the `if reference is not None:` block as the FIRST statements of the body, immediately before the existing `n_samples = draw(st.integers(1, max_samples))` line. Leave every existing body line (from `n_samples = ...` onward) exactly as it is.
```python
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
) -> VcfDocument:
    if reference is not None:
        return draw(
            _reference_documents(
                reference, violations, label_overrides, max_samples, max_records
            )
        )
    # --- existing reference-free body continues unchanged below ---
    n_samples = draw(st.integers(1, max_samples))
    samples = [f"s{i}" for i in range(n_samples)]
    ploidy = draw(st.integers(1, 2))
    b = VcfBuilder(samples=samples, contigs=[("chr1", 100000)]).fmt("GT")
    # ... (remaining existing lines unchanged) ...
```

- [ ] **Step 4: Run the tests to verify they pass**

Run:
```bash
uv run pytest tests/test_strategies.py -k "reference_consistent or back_compat" -q
```
Expected: PASS. If `filter_too_much` health-check fires, confirm `max_contig_len` headroom in `references()` and that `cursor` advancement leaves room; the suppress list already includes it.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/strategies.py tests/test_strategies.py
git commit -m "feat: reference-consistent documents() with violation labels"
```

---

## Task 8: `reference_and_documents()` paired strategy

Compose `references()` and `documents(reference=…)` into the `(ReferenceSpec, VcfDocument, GroundTruth)` tuple.

**Files:**
- Modify: `src/vcfixture/strategies.py`
- Modify: `tests/test_strategies.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_strategies.py`:
```python
@settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow])
@given(S.reference_and_documents())
def test_reference_and_documents_tuple(triple):
    spec, doc, truth = triple
    assert isinstance(spec, ReferenceSpec)
    assert truth.genotypes.shape[0] == len(doc.records)
    # Default (no violations) -> reference-consistent, unlabeled.
    for rec in doc.records:
        assert spec.seq(rec.chrom, rec.pos - 1, len(rec.ref)) == rec.ref
        assert rec.labels == frozenset()
```

- [ ] **Step 2: Run it to verify it fails**

Run:
```bash
uv run pytest tests/test_strategies.py -k reference_and_documents -q
```
Expected: FAIL — `AttributeError: ... has no attribute 'reference_and_documents'`.

- [ ] **Step 3: Implement `reference_and_documents()`**

In `src/vcfixture/strategies.py`, add `GroundTruth` to the imports (for the return annotation):
```python
from .truth import GroundTruth
```
Then append:
```python
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
) -> tuple[ReferenceSpec, VcfDocument, GroundTruth]:
    """Draw a consistent ``(ReferenceSpec, VcfDocument, GroundTruth)``."""
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
        )
    )
    return spec, doc, doc.truth()
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run pytest tests/test_strategies.py -k reference_and_documents -q
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/strategies.py tests/test_strategies.py
git commit -m "feat: reference_and_documents() paired strategy"
```

---

## Task 9: Public API exports, version bump, and full verification

**Files:**
- Modify: `src/vcfixture/__init__.py`
- Modify: `pyproject.toml`
- Modify: `tests/test_public_api.py` (if it asserts the export set — check first)

- [ ] **Step 1: Export the new public types**

In `src/vcfixture/__init__.py`, add to the imports:
```python
from .reference import Reference, ReferenceBuilder, ReferenceSpec, RepeatFeature
```
(replacing the existing `from .reference import Reference` line) and add the names to `__all__`:
```python
__all__ = [
    "VcfBuilder",
    "Genotype",
    "Reference",
    "ReferenceBuilder",
    "ReferenceSpec",
    "RepeatFeature",
    "GroundTruth",
    "Number",
    "Type",
    "strategies",
    "__version__",
]
```

- [ ] **Step 2: Check whether the public-API test pins the export set**

Run:
```bash
uv run pytest tests/test_public_api.py -q
```
If it FAILS because it asserts an exact `__all__`/dir set, open `tests/test_public_api.py` and add `ReferenceBuilder`, `ReferenceSpec`, `RepeatFeature` to whatever expected collection it checks. Re-run until green. If it already passes, make no change.

- [ ] **Step 3: Bump the version to 0.3.0**

In `pyproject.toml`, change:
```toml
version = "0.2.1"
```
to:
```toml
version = "0.3.0"
```
and under `[tool.commitizen]`:
```toml
version = "0.3.0"
```
(Both `version` fields — the `[project]` one and the `[tool.commitizen]` one — must read `0.3.0`.)

- [ ] **Step 4: Lint, type-check, and run the full suite**

Run:
```bash
uv run ruff check src tests
uv run ruff format --check src tests
uv run pyrefly check
uv run pytest -q
```
Expected: ruff clean (fix any unused imports it flags, e.g. a leftover `field` import or the illustrative test stubs from Task 7), `pyrefly check` reports no errors on `src`, and the full suite is green. If pyrefly flags the numpy `choice`/buffer assignments in `reference.py`, add a narrowly-scoped `# pyrefly: ignore[<code>]` on that exact line (mirror the existing `# pyrefly: ignore[bad-argument-type]` pattern in `strategies.py`), not a blanket suppression.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/__init__.py pyproject.toml tests/test_public_api.py
git commit -m "feat: export ReferenceBuilder/ReferenceSpec/RepeatFeature; bump 0.3.0"
```

- [ ] **Step 6: Final acceptance check**

Run:
```bash
uv run python -c "from vcfixture import ReferenceBuilder, ReferenceSpec, RepeatFeature, strategies; \
spec, doc, truth = strategies.reference_and_documents().example(); \
print('contigs', len(spec.contigs), 'records', len(doc.records), 'labels', truth.labels[:2])"
uv run pytest -q
```
Expected: prints a small summary with no error; full suite green. This is the acceptance gate for the upstream work.

---

## Notes & follow-ups (not part of this plan)

- **Release:** after merge, tag/release `0.3.0` to PyPI (commitizen `cz bump` / CI). GenVarLoader Phase 2 then bumps its `vcfixture>=0.3.0` dev dep — no path/git pins.
- **Downstream (GenVarLoader Phase 2):** add the property-test module consuming `reference_and_documents()`, running bcftools norm/consensus as the haplotype oracle, asserting gvl haplotypes == consensus and gvl genotypes/AF == `GroundTruth`, and using `labels` to partition canonical vs. deliberately-invalid examples. That work lives in the GenVarLoader repo, not here.
- **Possible later enhancement:** canonical insertions (left-aligned) and richer non-left-aligned constructs (right-shifted MNPs); deferred until a consumer needs them (YAGNI).
