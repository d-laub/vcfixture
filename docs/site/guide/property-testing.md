# Property-based testing

vcfixture ships a `strategies` module built on [Hypothesis](https://hypothesis.readthedocs.io/).
Use it when you want Hypothesis to explore the space of valid VCF documents for
you, rather than writing individual named fixtures by hand.

## Why property-based testing here

Both front-ends — [`VcfBuilder`](vcfbuilder.md) and the Hypothesis strategies — converge on the
same immutable `VcfDocument` hub:

```
 VcfBuilder.record(...)  ──┐
                           ├──▶  VcfDocument  ──▶  render()  →  VCF text
 strategies.documents()  ──┘     (immutable)  ──▶  truth()   →  GroundTruth
```

Because both paths write to the same model, the builder and the fuzzer **cannot
diverge**. [`GroundTruth`](../api/index.md) is derived from the same object that was serialized, so
the expected values are always in sync with the bytes on disk. You do not need to
maintain expected arrays separately — the oracle is free.

## The `strategies` surface

Import the module and call strategies directly:

```python
from vcfixture import strategies as S
```

| Strategy | Returns | Notes |
|---|---|---|
| `S.documents(...)` | `VcfDocument` | The general-purpose entry point. |
| `S.references(...)` | `ReferenceSpec` | A synthetic FASTA with optional planted tandem repeats. |
| `S.reference_and_documents(...)` | `tuple[ReferenceSpec, VcfDocument, GroundTruth]` | Convenience: draws spec + doc + truth together. |
| `S.symbolic_documents(...)` | `VcfDocument` | Documents whose ALTs are symbolic SV alleles (`<DEL>`, `<INS>`, …) with consistent SVLEN/SVCLAIM. |
| `S.documents_with_fields(...)` | `VcfDocument` | Documents carrying one INFO and one FORMAT field for every Number×Type combination. |
| `S.field_value(fielddef, n_alt, ploidy)` | `bool \| list` | A spec-valid value for a single `FieldDef`. |
| `S.genotypes(ploidy, n_alt, ...)` | `str` | A valid GT string (e.g. `"0|1"`, `"./1"`). |

All strategies are `@composite` and compose cleanly with `st.data()`.

### Drawing a document in a test

```python
from hypothesis import given, settings, HealthCheck
from vcfixture import strategies as S

@settings(max_examples=50, deadline=None,
          suppress_health_check=[HealthCheck.too_slow])
@given(S.documents())
def test_document_is_well_formed(doc):
    truth = doc.truth()
    assert truth.genotypes.shape[0] == len(doc.records)
    assert doc.render().startswith("##fileformat=")
```

`documents()` accepts `max_samples`, `max_records`, and `max_alt` to bound the
search space. Pass `reference=spec` (a `ReferenceSpec` drawn from
`S.references()`) to switch to reference-consistent, position-sorted output.

## Construct, never reject

vcfixture strategies **never call `.filter()` or `assume()`**. Instead, they
draw the parameters needed to construct a valid structure, then compute it
directly. For example, `genotypes` draws ploidy and allele count first, then
assembles the GT string; `documents` draws the number of records first, then
loops and draws each one in sequence.

This keeps Hypothesis's shrinking intact (filters break the causal chain between
drawn integers and the example), avoids health-check failures from excessive
rejection, and makes the strategies fast enough for CI.

The rule of thumb: if a value depends on another drawn value, draw the dependency
first and compute the result — do not draw speculatively and filter.

## Matrix coverage via `@parametrize`

`S.ALL_NUMBER_TYPE_COMBOS` is a pre-built list of every valid `(Number, Type,
kind)` triple, where `kind` is `"INFO"` or `"FORMAT"`. Use it to drive a
`@pytest.mark.parametrize` loop when you need exhaustive Number×Type coverage,
reserving Hypothesis for the *values within* a fixed combination:

```python
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from vcfixture import strategies as S
from vcfixture._spec.fielddef import FieldDef

@pytest.mark.parametrize("number,typ,kind", S.ALL_NUMBER_TYPE_COMBOS)
@settings(max_examples=30)
@given(st.data())
def test_field_value_has_correct_cardinality(number, typ, kind, data):
    fd = FieldDef("X", number, typ, "x", kind)
    val = data.draw(S.field_value(fd, n_alt=2, ploidy=2))
    if typ.value == "Flag":
        assert val is True
    else:
        assert isinstance(val, list)
```

This gives you a separate Hypothesis test per combo — deterministic outer loop,
fuzzed inner loop — without blowing up example counts. The full matrix is
currently 49 triples (all Number × all non-Flag Type for INFO and FORMAT, plus
the single Flag/INFO entry).

When you want every combo exercised in a single `@given` test (e.g. for a full
round-trip), use `S.documents_with_fields()` instead: it builds one document
that carries **all** Number×Type combinations simultaneously, so a single
property test covers the entire matrix.

## The oracle / round-trip pattern

The most important property test in vcfixture's own suite serializes an
arbitrary document, parses it with an **independent** third-party parser
(cyvcf2), and asserts that the decoded values match `GroundTruth`. This is what
makes the truth trustworthy rather than merely self-consistent.

The following is adapted from `tests/test_roundtrip.py`:

```python
import tempfile
from pathlib import Path

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings

from vcfixture import strategies as S

cyvcf2 = pytest.importorskip("cyvcf2")


def _genos_from_cyvcf2(variant, n_samples, ploidy):
    out = np.full((n_samples, ploidy), -1, dtype=np.int32)
    for si, g in enumerate(variant.genotypes):
        for ai, allele in enumerate(g[:-1]):
            out[si, ai] = allele
    return out


@settings(max_examples=75, deadline=None,
          suppress_health_check=[HealthCheck.too_slow])
@given(S.documents())
def test_genotypes_round_trip_through_cyvcf2(doc):
    truth = doc.truth()
    d = tempfile.mkdtemp()
    path = doc.write(Path(d) / "x.vcf.gz", bgzip=True, index=True)

    vf = cyvcf2.VCF(str(path))
    n_samples = len(doc.samples)
    ploidy = truth.genotypes.shape[2]
    for ri, variant in enumerate(vf):
        got = _genos_from_cyvcf2(variant, n_samples, ploidy)
        np.testing.assert_array_equal(
            got, truth.genotypes[ri],
            err_msg=f"genotype mismatch at record {ri}",
        )
        assert variant.POS == int(truth.pos[ri])
        assert variant.REF == truth.ref[ri]
```

Key points:

- `doc.truth()` returns a `GroundTruth` with `genotypes` (shape
  `[n_records, n_samples, ploidy]`, `-1` = missing), `pos`, `ref`, and
  per-record `info` / per-sample `format` dicts.
- `doc.write(path, bgzip=True, index=True)` serializes to `.vcf.gz` and writes
  a `.csi` index, ready for any parser that reads bgzipped VCF.
- The round-trip test uses `bgzip=True, index=True` because cyvcf2 requires the
  index to random-access records; if you only iterate sequentially, plain text
  suffices: `doc.write(path)` or `doc.render()` for an in-memory string.
- Missing alleles serialize as `.` and are decoded back to `-1` by both cyvcf2
  and pysam, matching vcfixture's convention.

`tests/test_matrix_roundtrip.py` extends this pattern to the full Number×Type
matrix via `S.documents_with_fields()`, verifying that every INFO and FORMAT
field type round-trips correctly through cyvcf2.
