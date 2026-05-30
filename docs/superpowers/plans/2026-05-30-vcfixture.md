# vcfixture Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a pure-Python library that generates small VCF (v4.5) test data — via an explicit builder and via Hypothesis strategies — and returns the decoded ground truth alongside it, so downstream parsers (genoray, GenVarLoader) assert against a known oracle instead of hand-coded numpy literals.

**Architecture:** One immutable `VcfDocument` model is the hub. The builder and the Hypothesis strategies are front-ends that produce the model; the serializer (model → VCF text) and truth-deriver (model → `GroundTruth`) are back-ends that consume it. This guarantees builder and fuzzer cannot diverge, and makes truth free. Self-validation round-trips every document through an independent parser (cyvcf2/pysam) and asserts agreement with the derived truth.

**Tech Stack:** Python 3.10+, numpy, hypothesis, pysam (bgzip/index + reference reads), cyvcf2 (independent cross-check parser), pytest. Packaging via hatchling.

---

## File Structure

```
vcfixture/
  pyproject.toml
  README.md
  src/vcfixture/
    __init__.py            # public API exports
    _spec/
      __init__.py
      types.py             # Type enum
      number.py            # NumberKind enum, Number, cardinality
      genotype_order.py    # genotype ordering enumerator (Number=G)
      fielddef.py          # FieldDef dataclass + validity invariants
      reserved.py          # curated reserved-field registry
    genotype.py            # Genotype (parse/render)
    model.py               # Record, ContigDef, VcfDocument
    variants.py            # variant-class constructors + classify()
    serialize.py           # model -> VCF text
    truth.py               # model -> GroundTruth
    build.py               # VcfBuilder
    io.py                  # write_text / write_bgzip (+ index)
    reference.py           # optional reference-aware adapter
    strategies.py          # Hypothesis strategies + parametrize tables
  tests/
    test_types.py
    test_number.py
    test_genotype_order.py
    test_fielddef.py
    test_reserved.py
    test_genotype.py
    test_model.py
    test_variants.py
    test_serialize.py
    test_truth.py
    test_build.py
    test_io.py
    test_reference.py
    test_strategies.py
    test_roundtrip.py
    test_benchmark.py
    test_genoray_parity.py
```

**Responsibilities:** `_spec/` holds the static VCF-4.5 rules (types, number cardinality, field validity, reserved fields, genotype ordering) with no knowledge of documents. `genotype.py`/`model.py` are the immutable data model. `variants.py` is the REF/ALT taxonomy. `serialize.py` and `truth.py` are the two back-ends. `build.py` and `strategies.py` are the two front-ends. `io.py`/`reference.py` are thin pysam wrappers.

---

## Task 0: Project scaffold

**Files:**
- Create: `pyproject.toml`
- Create: `src/vcfixture/__init__.py`
- Create: `src/vcfixture/_spec/__init__.py`
- Create: `tests/__init__.py`

- [ ] **Step 1: Write `pyproject.toml`**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "vcfixture"
version = "0.1.0"
description = "Generate small VCF (v4.5) test data with decoded ground truth."
requires-python = ">=3.10"
dependencies = ["numpy>=1.24", "pysam>=0.22"]

[project.optional-dependencies]
dev = ["pytest>=8", "hypothesis>=6.100", "cyvcf2>=0.30"]

[tool.hatch.build.targets.wheel]
packages = ["src/vcfixture"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

- [ ] **Step 2: Create empty package files**

```python
# src/vcfixture/__init__.py
"""vcfixture — generate small VCF test data with decoded ground truth."""
__version__ = "0.1.0"
```

```python
# src/vcfixture/_spec/__init__.py
```

```python
# tests/__init__.py
```

- [ ] **Step 3: Install editable with dev extras**

Run: `cd /Users/david/projects/vcfixture && python -m pip install -e ".[dev]"`
Expected: ends with `Successfully installed vcfixture-0.1.0`

- [ ] **Step 4: Verify import + pytest discovery**

Run: `cd /Users/david/projects/vcfixture && python -c "import vcfixture; print(vcfixture.__version__)" && python -m pytest -q`
Expected: prints `0.1.0`; pytest reports `no tests ran` (exit 5 is fine).

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml src/vcfixture/__init__.py src/vcfixture/_spec/__init__.py tests/__init__.py
git commit -m "chore: scaffold vcfixture package"
```

---

## Task 1: Type enum

**Files:**
- Create: `src/vcfixture/_spec/types.py`
- Test: `tests/test_types.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_types.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_types.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'vcfixture._spec.types'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/_spec/types.py
from __future__ import annotations
from enum import Enum

class Type(str, Enum):
    INTEGER = "Integer"
    FLOAT = "Float"
    FLAG = "Flag"
    CHARACTER = "Character"
    STRING = "String"

    @classmethod
    def info_allowed(cls) -> tuple["Type", ...]:
        return tuple(cls)

    @classmethod
    def format_allowed(cls) -> tuple["Type", ...]:
        return tuple(t for t in cls if t is not cls.FLAG)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_types.py -q`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/_spec/types.py tests/test_types.py
git commit -m "feat: add VCF Type enum"
```

---

## Task 2: Number + cardinality

**Files:**
- Create: `src/vcfixture/_spec/number.py`
- Test: `tests/test_number.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_number.py
import pytest
from vcfixture._spec.number import Number, NumberKind

def test_header_strings():
    assert Number.ONE.header_str() == "1"
    assert Number.fixed(2).header_str() == "2"
    assert Number.A.header_str() == "A"
    assert Number.R.header_str() == "R"
    assert Number.G.header_str() == "G"
    assert Number.DOT.header_str() == "."
    assert Number.FLAG.header_str() == "0"

def test_cardinality_fixed_a_r():
    assert Number.fixed(3).cardinality(n_alt=5, ploidy=2) == 3
    assert Number.A.cardinality(n_alt=2, ploidy=2) == 2
    assert Number.R.cardinality(n_alt=2, ploidy=2) == 3  # n_alt + 1
    assert Number.FLAG.cardinality(n_alt=2, ploidy=2) == 0

def test_cardinality_g_matches_spec_examples():
    # n_alleles = n_alt + 1; count = C(n_alleles + ploidy - 1, ploidy)
    assert Number.G.cardinality(n_alt=1, ploidy=2) == 3   # 00,01,11
    assert Number.G.cardinality(n_alt=2, ploidy=2) == 6   # 00,01,11,02,12,22
    assert Number.G.cardinality(n_alt=2, ploidy=3) == 10

def test_cardinality_dot_is_variable():
    assert Number.DOT.cardinality(n_alt=2, ploidy=2) is None

def test_equality_and_hash():
    assert Number.fixed(1) == Number.ONE
    assert {Number.A, Number.A} == {Number.A}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_number.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/_spec/number.py
from __future__ import annotations
import math
from dataclasses import dataclass
from enum import Enum

class NumberKind(str, Enum):
    FIXED = "FIXED"
    A = "A"
    R = "R"
    G = "G"
    DOT = "."
    FLAG = "FLAG"

@dataclass(frozen=True)
class Number:
    kind: NumberKind
    count: int | None = None  # set only for FIXED

    @classmethod
    def fixed(cls, n: int) -> "Number":
        if n < 0:
            raise ValueError("fixed Number must be >= 0")
        return cls(NumberKind.FIXED, n)

    def header_str(self) -> str:
        if self.kind is NumberKind.FIXED:
            return str(self.count)
        if self.kind is NumberKind.FLAG:
            return "0"
        return self.kind.value

    def cardinality(self, n_alt: int, ploidy: int) -> int | None:
        k = self.kind
        if k is NumberKind.FIXED:
            return self.count
        if k is NumberKind.FLAG:
            return 0
        if k is NumberKind.A:
            return n_alt
        if k is NumberKind.R:
            return n_alt + 1
        if k is NumberKind.G:
            n_alleles = n_alt + 1
            return math.comb(n_alleles + ploidy - 1, ploidy)
        return None  # DOT / variable

# Constants (assigned after class definition)
Number.ONE = Number(NumberKind.FIXED, 1)
Number.A = Number(NumberKind.A)
Number.R = Number(NumberKind.R)
Number.G = Number(NumberKind.G)
Number.DOT = Number(NumberKind.DOT)
Number.FLAG = Number(NumberKind.FLAG)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_number.py -q`
Expected: PASS (6 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/_spec/number.py tests/test_number.py
git commit -m "feat: add Number model and cardinality"
```

---

## Task 3: Genotype ordering enumerator

**Files:**
- Create: `src/vcfixture/_spec/genotype_order.py`
- Test: `tests/test_genotype_order.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_genotype_order.py
from vcfixture._spec.genotype_order import genotype_ordering

def test_diploid_biallelic():
    assert genotype_ordering(ploidy=2, n_alleles=2) == [(0, 0), (0, 1), (1, 1)]

def test_diploid_triallelic():
    assert genotype_ordering(ploidy=2, n_alleles=3) == [
        (0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)
    ]

def test_triploid_triallelic_count():
    out = genotype_ordering(ploidy=3, n_alleles=3)
    assert len(out) == 10
    assert out[0] == (0, 0, 0)
    assert out[-1] == (2, 2, 2)

def test_haploid():
    assert genotype_ordering(ploidy=1, n_alleles=3) == [(0,), (1,), (2,)]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_genotype_order.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/_spec/genotype_order.py
from __future__ import annotations

def genotype_ordering(ploidy: int, n_alleles: int) -> list[tuple[int, ...]]:
    """Ordered genotypes per VCF 4.5 'Number=G' ordering.

    Recursion: genotypes of ploidy p over alleles [0, a] are, for each
    allele a in 0..n_alleles-1, the genotypes of ploidy p-1 over [0, a]
    each suffixed with a. This yields the spec ordering 00,01,11,02,12,22.
    """
    if ploidy < 1:
        raise ValueError("ploidy must be >= 1")

    def rec(p: int) -> list[tuple[int, ...]]:
        if p == 1:
            return [(a,) for a in range(n_alleles)]
        out: list[tuple[int, ...]] = []
        for a in range(n_alleles):
            for prefix in rec(p - 1):
                if prefix[-1] <= a:
                    out.append(prefix + (a,))
        return out

    return rec(ploidy)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_genotype_order.py -q`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/_spec/genotype_order.py tests/test_genotype_order.py
git commit -m "feat: add Number=G genotype ordering"
```

---

## Task 4: FieldDef + validity invariants

**Files:**
- Create: `src/vcfixture/_spec/fielddef.py`
- Test: `tests/test_fielddef.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_fielddef.py
import pytest
from vcfixture._spec.fielddef import FieldDef
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type

def test_valid_info_field():
    f = FieldDef("AF", Number.A, Type.FLOAT, "Allele frequency", "INFO")
    assert f.id == "AF"
    assert f.header_line() == (
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">'
    )

def test_flag_must_be_number_zero():
    with pytest.raises(ValueError, match="Flag"):
        FieldDef("DB", Number.ONE, Type.FLAG, "x", "INFO")

def test_flag_must_be_info():
    with pytest.raises(ValueError, match="Flag"):
        FieldDef("DB", Number.FLAG, Type.FLAG, "x", "FORMAT")

def test_format_cannot_be_flag():
    with pytest.raises(ValueError, match="Flag"):
        FieldDef("X", Number.ONE, Type.FLAG, "x", "FORMAT")

def test_id_must_match_regex():
    with pytest.raises(ValueError, match="ID"):
        FieldDef("1bad", Number.ONE, Type.STRING, "x", "INFO")
    # '1000G' is explicitly allowed by the spec regex
    FieldDef("1000G", Number.DOT, Type.STRING, "x", "INFO")

def test_kind_must_be_valid():
    with pytest.raises(ValueError, match="kind"):
        FieldDef("X", Number.ONE, Type.STRING, "x", "BOGUS")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_fielddef.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/_spec/fielddef.py
from __future__ import annotations
import re
from dataclasses import dataclass
from .number import Number, NumberKind
from .types import Type

_ID_RE = re.compile(r"^([A-Za-z_][0-9A-Za-z_.]*|1000G)$")
_KINDS = ("INFO", "FORMAT")

@dataclass(frozen=True)
class FieldDef:
    id: str
    number: Number
    type: Type
    description: str
    kind: str  # "INFO" | "FORMAT"

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise ValueError(f"kind must be one of {_KINDS}, got {self.kind!r}")
        if not _ID_RE.match(self.id):
            raise ValueError(f"ID {self.id!r} does not match VCF key regex")
        if self.type is Type.FLAG:
            if self.kind != "INFO":
                raise ValueError("Flag fields must be INFO, not FORMAT")
            if self.number.kind is not NumberKind.FLAG:
                raise ValueError("Flag fields must have Number=0")
        elif self.number.kind is NumberKind.FLAG:
            raise ValueError("Number=0 is only valid for Flag fields")

    def header_line(self) -> str:
        return (
            f"##{self.kind}=<ID={self.id},Number={self.number.header_str()},"
            f'Type={self.type.value},Description="{self.description}">'
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_fielddef.py -q`
Expected: PASS (6 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/_spec/fielddef.py tests/test_fielddef.py
git commit -m "feat: add FieldDef with validity invariants"
```

---

## Task 5: Reserved-field registry (curated)

**Files:**
- Create: `src/vcfixture/_spec/reserved.py`
- Test: `tests/test_reserved.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_reserved.py
import pytest
from vcfixture._spec.reserved import reserved
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type

def test_reserved_info_af():
    f = reserved("AF", "INFO")
    assert f.number == Number.A and f.type == Type.FLOAT and f.kind == "INFO"

def test_reserved_format_gt():
    f = reserved("GT", "FORMAT")
    assert f.number == Number.ONE and f.type == Type.STRING

def test_reserved_format_pl_is_g_integer():
    f = reserved("PL", "FORMAT")
    assert f.number == Number.G and f.type == Type.INTEGER

def test_reserved_flag():
    f = reserved("DB", "INFO")
    assert f.type == Type.FLAG and f.number == Number.FLAG

def test_unknown_reserved_raises():
    with pytest.raises(KeyError):
        reserved("NOPE", "INFO")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_reserved.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/_spec/reserved.py
from __future__ import annotations
from .fielddef import FieldDef
from .number import Number
from .types import Type

# Curated subset — names are secondary; this covers the Number/Type classes
# the downstream consumers actually use.
_INFO = {
    "AA": FieldDef("AA", Number.ONE, Type.STRING, "Ancestral allele", "INFO"),
    "AC": FieldDef("AC", Number.A, Type.INTEGER, "Allele count", "INFO"),
    "AF": FieldDef("AF", Number.A, Type.FLOAT, "Allele frequency", "INFO"),
    "AN": FieldDef("AN", Number.ONE, Type.INTEGER, "Total allele number", "INFO"),
    "DP": FieldDef("DP", Number.ONE, Type.INTEGER, "Combined depth", "INFO"),
    "DB": FieldDef("DB", Number.FLAG, Type.FLAG, "dbSNP membership", "INFO"),
    "H2": FieldDef("H2", Number.FLAG, Type.FLAG, "HapMap2 membership", "INFO"),
}
_FORMAT = {
    "GT": FieldDef("GT", Number.ONE, Type.STRING, "Genotype", "FORMAT"),
    "GQ": FieldDef("GQ", Number.ONE, Type.INTEGER, "Genotype quality", "FORMAT"),
    "DP": FieldDef("DP", Number.ONE, Type.INTEGER, "Read depth", "FORMAT"),
    "AD": FieldDef("AD", Number.R, Type.INTEGER, "Allelic depths", "FORMAT"),
    "PL": FieldDef("PL", Number.G, Type.INTEGER, "Phred genotype likelihoods", "FORMAT"),
    "GL": FieldDef("GL", Number.G, Type.FLOAT, "Log10 genotype likelihoods", "FORMAT"),
    "PS": FieldDef("PS", Number.ONE, Type.INTEGER, "Phase set", "FORMAT"),
}

def reserved(id: str, kind: str) -> FieldDef:
    table = _INFO if kind == "INFO" else _FORMAT
    return table[id]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_reserved.py -q`
Expected: PASS (5 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/_spec/reserved.py tests/test_reserved.py
git commit -m "feat: add curated reserved-field registry"
```

---

## Task 6: Genotype model (parse + render)

**Files:**
- Create: `src/vcfixture/genotype.py`
- Test: `tests/test_genotype.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_genotype.py
import pytest
from vcfixture.genotype import Genotype

def test_parse_phased_diploid():
    g = Genotype.parse("0|1")
    assert g.alleles == (0, 1)
    assert g.phased == (True,)
    assert g.is_phased is True
    assert g.ploidy == 2

def test_parse_unphased_with_missing():
    g = Genotype.parse("./1")
    assert g.alleles == (None, 1)
    assert g.phased == (False,)
    assert g.is_phased is False

def test_parse_haploid():
    g = Genotype.parse("1")
    assert g.alleles == (1,)
    assert g.phased == ()

def test_round_trip_render():
    for s in ["0|1", "1/1", "./.", "0/0/1", "1|0|2"]:
        assert Genotype.parse(s).render() == s

def test_render_missing_genotype():
    assert Genotype.parse("./.").render() == "./."
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_genotype.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/genotype.py
from __future__ import annotations
import re
from dataclasses import dataclass

_TOKEN_RE = re.compile(r"([|/])")

@dataclass(frozen=True)
class Genotype:
    alleles: tuple[int | None, ...]   # None == missing allele
    phased: tuple[bool, ...]          # one per separator; len == len(alleles) - 1

    @classmethod
    def parse(cls, s: str) -> "Genotype":
        parts = _TOKEN_RE.split(s)  # e.g. "0|1" -> ["0","|","1"]
        alleles: list[int | None] = []
        phased: list[bool] = []
        for i, tok in enumerate(parts):
            if i % 2 == 0:
                alleles.append(None if tok == "." else int(tok))
            else:
                phased.append(tok == "|")
        return cls(tuple(alleles), tuple(phased))

    @property
    def ploidy(self) -> int:
        return len(self.alleles)

    @property
    def is_phased(self) -> bool:
        return len(self.phased) > 0 and all(self.phased)

    def render(self) -> str:
        out = ["." if a is None else str(a) for a in self.alleles]
        seps = ["|" if p else "/" for p in self.phased]
        chars: list[str] = [out[0]]
        for sep, allele in zip(seps, out[1:]):
            chars.append(sep)
            chars.append(allele)
        return "".join(chars)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_genotype.py -q`
Expected: PASS (5 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/genotype.py tests/test_genotype.py
git commit -m "feat: add Genotype parse/render"
```

---

## Task 7: Record + VcfDocument model

**Files:**
- Create: `src/vcfixture/model.py`
- Test: `tests/test_model.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_model.py
from vcfixture.model import Record, ContigDef, VcfDocument
from vcfixture.genotype import Genotype
from vcfixture._spec.reserved import reserved

def _doc():
    rec = Record(
        chrom="chr1", pos=81262, ids=None, ref="GAT", alts=("A",),
        qual=None, filters=None, info={},
        fmt_keys=("GT", "DS"),
        samples=({"GT": Genotype.parse("0|1"), "DS": [1.0]},
                 {"GT": Genotype.parse("1|1"), "DS": [2.0]}),
    )
    return VcfDocument(
        fileformat="VCFv4.5",
        info_defs=(),
        format_defs=(reserved("GT", "FORMAT"),),
        filter_defs=(),
        contigs=(ContigDef("chr1", 1000),),
        samples=("s1", "s2"),
        records=(rec,),
    )

def test_document_is_frozen_and_holds_records():
    doc = _doc()
    assert doc.samples == ("s1", "s2")
    assert doc.records[0].alts == ("A",)
    assert doc.max_ploidy() == 2

def test_ploidy_varies_uses_max():
    rec = Record("chr1", 1, None, "A", ("T",), None, None, {},
                 ("GT",), ({"GT": Genotype.parse("0/0/1")},))
    doc = VcfDocument("VCFv4.5", (), (), (), (ContigDef("chr1", None),),
                      ("s1",), (rec,))
    assert doc.max_ploidy() == 3
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_model.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/model.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Mapping
from .genotype import Genotype
from ._spec.fielddef import FieldDef

@dataclass(frozen=True)
class ContigDef:
    id: str
    length: int | None = None

    def header_line(self) -> str:
        if self.length is None:
            return f"##contig=<ID={self.id}>"
        return f"##contig=<ID={self.id},length={self.length}>"

@dataclass(frozen=True)
class Record:
    chrom: str
    pos: int                                  # 1-based
    ids: tuple[str, ...] | None               # None -> "."
    ref: str
    alts: tuple[str, ...]                     # may contain "*"
    qual: float | None
    filters: tuple[str, ...] | None           # None -> "."; () -> "PASS"
    info: Mapping[str, Any]                   # id -> value(s); Flag -> True
    fmt_keys: tuple[str, ...]                 # FORMAT column order
    samples: tuple[Mapping[str, Any], ...]    # per-sample: key -> value(s)/Genotype

    @property
    def n_alt(self) -> int:
        return len(self.alts)

@dataclass(frozen=True)
class VcfDocument:
    fileformat: str
    info_defs: tuple[FieldDef, ...]
    format_defs: tuple[FieldDef, ...]
    filter_defs: tuple[tuple[str, str], ...]  # (id, description)
    contigs: tuple[ContigDef, ...]
    samples: tuple[str, ...]
    records: tuple[Record, ...]

    def max_ploidy(self) -> int:
        p = 1
        for rec in self.records:
            for s in rec.samples:
                gt = s.get("GT")
                if isinstance(gt, Genotype):
                    p = max(p, gt.ploidy)
        return p

    # Back-end delegations (imported lazily to avoid cycles).
    def render(self) -> str:
        from .serialize import render_document
        return render_document(self)

    def truth(self):
        from .truth import derive_truth
        return derive_truth(self)

    def write(self, path, *, bgzip: bool = False, index: bool = False):
        from . import io
        if bgzip:
            return io.write_bgzip(self, path, index=index)
        return io.write_text(self, path)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_model.py -q`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/model.py tests/test_model.py
git commit -m "feat: add Record, ContigDef, VcfDocument model"
```

---

## Task 8: Variant-class constructors + classify

**Files:**
- Create: `src/vcfixture/variants.py`
- Test: `tests/test_variants.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_variants.py
import pytest
from vcfixture import variants as v

def test_constructors():
    assert v.snp("A", "T") == ("A", "T")
    assert v.mnp("AC", "GT") == ("AC", "GT")
    assert v.insertion("A", "CGT") == ("A", "ACGT")
    assert v.deletion("A", "CGT") == ("ACGT", "A")
    assert v.delins("AC", "GTA") == ("AC", "GTA")
    assert v.spanning_deletion("A") == ("A", "*")

def test_classify_single_alt():
    assert v.classify("A", "T") == "SNP"
    assert v.classify("AC", "GT") == "MNP"
    assert v.classify("A", "ACGT") == "INS"
    assert v.classify("ACGT", "A") == "DEL"
    assert v.classify("AC", "GTA") == "DELINS"
    assert v.classify("A", "*") == "SPANNING_DEL"

def test_record_class_multiallelic():
    assert v.record_class("G", ("A", "C")) == "MULTIALLELIC"
    assert v.record_class("G", ("A",)) == "SNP"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_variants.py -q`
Expected: FAIL with `ImportError`/`AttributeError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/variants.py
from __future__ import annotations

def snp(ref_base: str, alt_base: str) -> tuple[str, str]:
    return ref_base, alt_base

def mnp(ref_seq: str, alt_seq: str) -> tuple[str, str]:
    return ref_seq, alt_seq

def insertion(anchor: str, inserted: str) -> tuple[str, str]:
    return anchor, anchor + inserted

def deletion(anchor: str, deleted: str) -> tuple[str, str]:
    return anchor + deleted, anchor

def delins(ref_seq: str, alt_seq: str) -> tuple[str, str]:
    return ref_seq, alt_seq

def spanning_deletion(ref_base: str) -> tuple[str, str]:
    return ref_base, "*"

def classify(ref: str, alt: str) -> str:
    if alt == "*":
        return "SPANNING_DEL"
    lr, la = len(ref), len(alt)
    if lr == 1 and la == 1:
        return "SNP"
    if lr == la:
        return "MNP"
    if la > lr and alt.startswith(ref):
        return "INS"
    if lr > la and ref.startswith(alt):
        return "DEL"
    return "DELINS"

def record_class(ref: str, alts: tuple[str, ...]) -> str:
    if len(alts) > 1:
        return "MULTIALLELIC"
    return classify(ref, alts[0])
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_variants.py -q`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/variants.py tests/test_variants.py
git commit -m "feat: add variant-class constructors and classify"
```

---

## Task 9: Serializer (model -> VCF text)

**Files:**
- Create: `src/vcfixture/serialize.py`
- Test: `tests/test_serialize.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_serialize.py
from vcfixture.model import Record, ContigDef, VcfDocument
from vcfixture.genotype import Genotype
from vcfixture.serialize import render_document
from vcfixture._spec.fielddef import FieldDef
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type

def _doc():
    af = FieldDef("AF", Number.A, Type.FLOAT, "Allele frequency", "INFO")
    db = FieldDef("DB", Number.FLAG, Type.FLAG, "dbSNP", "INFO")
    gt = FieldDef("GT", Number.ONE, Type.STRING, "Genotype", "FORMAT")
    ds = FieldDef("DS", Number.A, Type.FLOAT, "Dosage", "FORMAT")
    rec = Record(
        chrom="chr1", pos=81262, ids=None, ref="GAT", alts=("A",),
        qual=None, filters=None,
        info={"AF": [0.5], "DB": True},
        fmt_keys=("GT", "DS"),
        samples=({"GT": Genotype.parse("0|1"), "DS": [1.0]},
                 {"GT": Genotype.parse("./."), "DS": [None]}),
    )
    return VcfDocument("VCFv4.5", (af, db), (gt, ds), (),
                       (ContigDef("chr1", 1000),), ("s1", "s2"), (rec,))

def test_header_first_line_is_fileformat():
    text = render_document(_doc())
    assert text.splitlines()[0] == "##fileformat=VCFv4.5"

def test_chrom_header_includes_samples():
    text = render_document(_doc())
    chrom = [l for l in text.splitlines() if l.startswith("#CHROM")][0]
    assert chrom.split("\t")[-2:] == ["s1", "s2"]
    assert chrom.split("\t")[8] == "FORMAT"

def test_record_line_fields():
    text = render_document(_doc())
    data = [l for l in text.splitlines() if not l.startswith("#")][0]
    cols = data.split("\t")
    assert cols[0:5] == ["chr1", "81262", ".", "GAT", "A"]
    assert cols[5] == "."          # QUAL missing
    assert cols[6] == "."          # FILTER missing
    assert cols[7] == "AF=0.5;DB"  # Flag emitted as bare key
    assert cols[8] == "GT:DS"
    assert cols[9] == "0|1:1.0"
    assert cols[10] == "./.:."     # missing dosage -> "."

def test_pass_filter_renders_PASS():
    doc = _doc()
    rec = doc.records[0]
    rec2 = type(rec)(**{**rec.__dict__, "filters": ()})
    doc2 = type(doc)(**{**doc.__dict__, "records": (rec2,)})
    data = [l for l in render_document(doc2).splitlines() if not l.startswith("#")][0]
    assert data.split("\t")[6] == "PASS"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_serialize.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/serialize.py
from __future__ import annotations
from typing import Any
from .genotype import Genotype
from .model import Record, VcfDocument

def _fmt_scalar(v: Any) -> str:
    if v is None:
        return "."
    if isinstance(v, float):
        # Compact but stable float repr; integers-as-floats keep one decimal.
        return repr(v)
    return str(v)

def _fmt_value(v: Any) -> str:
    """Render an INFO/FORMAT value that may be a scalar or a list."""
    if isinstance(v, (list, tuple)):
        if len(v) == 0:
            return "."
        return ",".join(_fmt_scalar(x) for x in v)
    return _fmt_scalar(v)

def _render_info(rec: Record) -> str:
    if not rec.info:
        return "."
    parts: list[str] = []
    for key, val in rec.info.items():
        if val is True:            # Flag present
            parts.append(key)
        elif val is False:         # Flag absent -> omit
            continue
        else:
            parts.append(f"{key}={_fmt_value(val)}")
    return ";".join(parts) if parts else "."

def _render_sample(rec: Record, sample: dict) -> str:
    vals: list[str] = []
    for key in rec.fmt_keys:
        v = sample.get(key)
        if isinstance(v, Genotype):
            vals.append(v.render())
        else:
            vals.append(_fmt_value(v))
    # Drop trailing missing fields, but never GT (index 0).
    while len(vals) > 1 and vals[-1] == ".":
        vals.pop()
    return ":".join(vals)

def _render_record(rec: Record) -> str:
    ids = ";".join(rec.ids) if rec.ids else "."
    alt = ",".join(rec.alts) if rec.alts else "."
    qual = "." if rec.qual is None else _fmt_scalar(rec.qual)
    if rec.filters is None:
        filt = "."
    elif len(rec.filters) == 0:
        filt = "PASS"
    else:
        filt = ";".join(rec.filters)
    cols = [rec.chrom, str(rec.pos), ids, rec.ref, alt, qual, filt,
            _render_info(rec)]
    if rec.fmt_keys:
        cols.append(":".join(rec.fmt_keys))
        for s in rec.samples:
            cols.append(_render_sample(rec, s))
    return "\t".join(cols)

def render_document(doc: VcfDocument) -> str:
    lines = [f"##fileformat={doc.fileformat}"]
    for c in doc.contigs:
        lines.append(c.header_line())
    for f in doc.info_defs:
        lines.append(f.header_line())
    for fid, desc in doc.filter_defs:
        lines.append(f'##FILTER=<ID={fid},Description="{desc}">')
    for f in doc.format_defs:
        lines.append(f.header_line())
    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if doc.format_defs or any(r.fmt_keys for r in doc.records):
        header.append("FORMAT")
        header.extend(doc.samples)
    lines.append("\t".join(header))
    for rec in doc.records:
        lines.append(_render_record(rec))
    return "\n".join(lines) + "\n"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_serialize.py -q`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/serialize.py tests/test_serialize.py
git commit -m "feat: add VCF serializer"
```

---

## Task 10: Truth deriver (model -> GroundTruth)

**Files:**
- Create: `src/vcfixture/truth.py`
- Test: `tests/test_truth.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_truth.py
import numpy as np
from vcfixture.model import Record, ContigDef, VcfDocument
from vcfixture.genotype import Genotype
from vcfixture.truth import derive_truth

def _doc():
    rec = Record(
        chrom="chr1", pos=81262, ids=None, ref="GAT", alts=("A",),
        qual=None, filters=None, info={"AF": [0.5]},
        fmt_keys=("GT", "DS"),
        samples=({"GT": Genotype.parse("0|1"), "DS": [1.0]},
                 {"GT": Genotype.parse("./."), "DS": [None]}),
    )
    return VcfDocument("VCFv4.5", (), (), (), (ContigDef("chr1", 1000),),
                       ("s1", "s2"), (rec,))

def test_genotype_matrix_with_missing_sentinel():
    t = derive_truth(_doc())
    # (records, samples, ploidy); missing -> -1
    assert t.genotypes.shape == (1, 2, 2)
    np.testing.assert_array_equal(t.genotypes[0, 0], [0, 1])
    np.testing.assert_array_equal(t.genotypes[0, 1], [-1, -1])

def test_phasing_matrix():
    t = derive_truth(_doc())
    assert t.phasing.shape == (1, 2)
    assert bool(t.phasing[0, 0]) is True    # 0|1 phased
    assert bool(t.phasing[0, 1]) is False   # ./. unphased

def test_pos_ref_alt_and_class():
    t = derive_truth(_doc())
    np.testing.assert_array_equal(t.pos, [81262])
    assert t.ref == ["GAT"]
    assert t.alts == [["A"]]
    # GAT->A changes the first base, so it is a complex substitution (DELINS),
    # not a clean anchored deletion (matches classify() in variants.py).
    assert t.variant_class == ["DELINS"]

def test_info_and_format_echoed():
    t = derive_truth(_doc())
    assert t.info[0]["AF"] == [0.5]
    assert t.format[0][0]["DS"] == [1.0]
    assert t.format[0][1]["DS"] == [None]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_truth.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/truth.py
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .genotype import Genotype
from .model import VcfDocument
from .variants import record_class

@dataclass(frozen=True)
class GroundTruth:
    samples: tuple[str, ...]
    contigs: tuple[str, ...]
    pos: np.ndarray                 # (records,) int64, 1-based
    ref: list[str]
    alts: list[list[str]]
    variant_class: list[str]
    genotypes: np.ndarray           # (records, samples, ploidy) int32, -1 missing
    phasing: np.ndarray             # (records, samples) bool (fully phased)
    info: list[dict]                # per record: id -> decoded value(s)
    format: list[list[dict]]        # per record, per sample: id -> value(s)

def derive_truth(doc: VcfDocument) -> GroundTruth:
    n_rec = len(doc.records)
    n_smp = len(doc.samples)
    ploidy = doc.max_ploidy()

    genos = np.full((n_rec, n_smp, ploidy), -1, dtype=np.int32)
    phasing = np.zeros((n_rec, n_smp), dtype=bool)
    pos = np.zeros(n_rec, dtype=np.int64)
    ref: list[str] = []
    alts: list[list[str]] = []
    vclass: list[str] = []
    info: list[dict] = []
    fmt: list[list[dict]] = []

    for ri, rec in enumerate(doc.records):
        pos[ri] = rec.pos
        ref.append(rec.ref)
        alts.append(list(rec.alts))
        vclass.append(record_class(rec.ref, rec.alts))
        info.append(dict(rec.info))
        per_sample: list[dict] = []
        for si, sample in enumerate(rec.samples):
            gt = sample.get("GT")
            if isinstance(gt, Genotype):
                for ai, allele in enumerate(gt.alleles):
                    genos[ri, si, ai] = -1 if allele is None else allele
                phasing[ri, si] = gt.is_phased
            per_sample.append({k: v for k, v in sample.items() if k != "GT"})
        fmt.append(per_sample)

    return GroundTruth(
        samples=doc.samples,
        contigs=tuple(c.id for c in doc.contigs),
        pos=pos, ref=ref, alts=alts, variant_class=vclass,
        genotypes=genos, phasing=phasing, info=info, format=fmt,
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_truth.py -q`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/truth.py tests/test_truth.py
git commit -m "feat: add ground-truth deriver"
```

---

## Task 11: VcfBuilder

**Files:**
- Create: `src/vcfixture/build.py`
- Test: `tests/test_build.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_build.py
import numpy as np
import pytest
from vcfixture.build import VcfBuilder
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type

def test_build_biallelic_with_dosage():
    doc = (VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 1000)])
           .info("AF", Number.A, Type.FLOAT)
           .fmt("GT")
           .fmt("DS", Number.A, Type.FLOAT)
           .record("chr1", 81262, ref="GAT", alt=["A"],
                   gt=["0|1", "1|1"], info={"AF": [0.5]},
                   DS=[[1.0], [2.0]])
           .build())
    t = doc.truth()
    np.testing.assert_array_equal(t.genotypes[0], [[0, 1], [1, 1]])
    assert t.format[0][0]["DS"] == [1.0]
    assert "##fileformat=VCFv4.5" in doc.render()

def test_reserved_field_by_name():
    doc = (VcfBuilder(samples=["s1"], contigs=[("chr1", None)])
           .fmt("GT")
           .record("chr1", 1, ref="A", alt=["T"], gt=["0/1"])
           .build())
    assert doc.format_defs[0].id == "GT"

def test_undefined_format_field_raises():
    b = (VcfBuilder(samples=["s1"], contigs=[("chr1", None)]).fmt("GT"))
    with pytest.raises(ValueError, match="not declared"):
        b.record("chr1", 1, ref="A", alt=["T"], gt=["0/1"], DS=[[1.0]])

def test_gt_index_out_of_range_raises():
    b = (VcfBuilder(samples=["s1"], contigs=[("chr1", None)]).fmt("GT"))
    with pytest.raises(ValueError, match="allele index"):
        b.record("chr1", 1, ref="A", alt=["T"], gt=["0/5"])

def test_cardinality_mismatch_raises():
    b = (VcfBuilder(samples=["s1"], contigs=[("chr1", None)])
         .fmt("GT").fmt("AD", Number.R, Type.INTEGER))
    # biallelic -> R == 2 values expected, provide 1
    with pytest.raises(ValueError, match="cardinality"):
        b.record("chr1", 1, ref="A", alt=["T"], gt=["0/1"], AD=[[5]])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_build.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/build.py
from __future__ import annotations
from typing import Any
from .genotype import Genotype
from .model import ContigDef, Record, VcfDocument
from ._spec.fielddef import FieldDef
from ._spec.number import Number, NumberKind
from ._spec.types import Type
from ._spec.reserved import reserved

class VcfBuilder:
    def __init__(self, samples, contigs, fileformat: str = "VCFv4.5"):
        self._samples = tuple(samples)
        self._contigs = tuple(ContigDef(c[0], c[1]) for c in contigs)
        self._fileformat = fileformat
        self._info_defs: dict[str, FieldDef] = {}
        self._format_defs: dict[str, FieldDef] = {}
        self._filter_defs: list[tuple[str, str]] = []
        self._records: list[Record] = []

    def info(self, id, number=None, type=None, description=None) -> "VcfBuilder":
        self._info_defs[id] = self._make_def(id, number, type, description, "INFO")
        return self

    def fmt(self, id, number=None, type=None, description=None) -> "VcfBuilder":
        self._format_defs[id] = self._make_def(id, number, type, description, "FORMAT")
        return self

    def filter(self, id: str, description: str) -> "VcfBuilder":
        self._filter_defs.append((id, description))
        return self

    @staticmethod
    def _make_def(id, number, type, description, kind) -> FieldDef:
        if number is None or type is None:
            return reserved(id, kind)  # KeyError if unknown + undefined
        return FieldDef(id, number, type, description or id, kind)

    def record(self, chrom, pos, *, ref, alt, ids=None, qual=None,
               filter=None, gt=None, info=None, **fmt_fields) -> "VcfBuilder":
        alts = tuple(alt)
        n_alt = len(alts)

        fmt_keys: list[str] = []
        samples: list[dict[str, Any]] = []
        for _ in self._samples:
            samples.append({})

        if gt is not None:
            if "GT" not in self._format_defs:
                raise ValueError("GT not declared; call .fmt('GT')")
            fmt_keys.append("GT")
            for si, s in enumerate(gt):
                geno = Genotype.parse(s)
                for a in geno.alleles:
                    if a is not None and a > n_alt:
                        raise ValueError(
                            f"allele index {a} out of range (n_alt={n_alt})")
                samples[si]["GT"] = geno

        ploidy = max((s["GT"].ploidy for s in samples if "GT" in s), default=2)

        for key, per_sample in fmt_fields.items():
            if key not in self._format_defs:
                raise ValueError(f"FORMAT field {key!r} not declared")
            fdef = self._format_defs[key]
            fmt_keys.append(key)
            card = fdef.number.cardinality(n_alt, ploidy)
            for si, val in enumerate(per_sample):
                if card is not None and isinstance(val, (list, tuple)) \
                        and len(val) != card:
                    raise ValueError(
                        f"{key} cardinality mismatch: expected {card}, "
                        f"got {len(val)}")
                samples[si][key] = val

        info_dict: dict[str, Any] = {}
        if info:
            for key, val in info.items():
                if key not in self._info_defs:
                    raise ValueError(f"INFO field {key!r} not declared")
                fdef = self._info_defs[key]
                card = fdef.number.cardinality(n_alt, ploidy)
                if card is not None and fdef.number.kind is not NumberKind.FLAG \
                        and isinstance(val, (list, tuple)) and len(val) != card:
                    raise ValueError(
                        f"{key} cardinality mismatch: expected {card}, "
                        f"got {len(val)}")
                info_dict[key] = val

        self._records.append(Record(
            chrom=chrom, pos=pos,
            ids=tuple(ids) if ids else None,
            ref=ref, alts=alts, qual=qual,
            filters=tuple(filter) if filter is not None else None,
            info=info_dict, fmt_keys=tuple(fmt_keys), samples=tuple(samples),
        ))
        return self

    def build(self) -> VcfDocument:
        return VcfDocument(
            fileformat=self._fileformat,
            info_defs=tuple(self._info_defs.values()),
            format_defs=tuple(self._format_defs.values()),
            filter_defs=tuple(self._filter_defs),
            contigs=self._contigs,
            samples=self._samples,
            records=tuple(self._records),
        )

    # convenience pass-throughs
    def render(self) -> str: return self.build().render()
    def truth(self): return self.build().truth()
    def write(self, path, **kw): return self.build().write(path, **kw)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_build.py -q`
Expected: PASS (5 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/build.py tests/test_build.py
git commit -m "feat: add VcfBuilder with eager validation"
```

---

## Task 12: IO — text + bgzip + index

**Files:**
- Create: `src/vcfixture/io.py`
- Test: `tests/test_io.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_io.py
from pathlib import Path
import pysam
from vcfixture.build import VcfBuilder

def _doc():
    return (VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 100000)])
            .fmt("GT")
            .record("chr1", 81262, ref="A", alt=["T"], gt=["0|1", "1|1"])
            .build())

def test_write_text(tmp_path: Path):
    p = _doc().write(tmp_path / "x.vcf")
    assert Path(p).read_text().startswith("##fileformat=VCFv4.5")

def test_write_bgzip_and_index_is_queryable(tmp_path: Path):
    p = _doc().write(tmp_path / "x.vcf.gz", bgzip=True, index=True)
    assert Path(str(p) + ".csi").exists()
    vf = pysam.VariantFile(str(p))
    rows = list(vf.fetch("chr1", 81000, 82000))
    assert len(rows) == 1
    assert rows[0].pos == 81262
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_io.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/io.py
from __future__ import annotations
from pathlib import Path
import pysam
from .model import VcfDocument

def write_text(doc: VcfDocument, path) -> Path:
    p = Path(path)
    p.write_text(doc.render())
    return p

def write_bgzip(doc: VcfDocument, path, *, index: bool = True) -> Path:
    p = Path(path)
    if p.suffix != ".gz":
        p = p.with_suffix(p.suffix + ".gz")
    tmp = p.with_suffix("")  # strip .gz for the plain-text intermediate
    tmp.write_text(doc.render())
    pysam.tabix_compress(str(tmp), str(p), force=True)
    tmp.unlink()
    if index:
        pysam.tabix_index(str(p), preset="vcf", csi=True, force=True)
    return p
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_io.py -q`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/io.py tests/test_io.py
git commit -m "feat: add text + bgzip/index IO"
```

---

## Task 13: Reference-aware adapter

**Files:**
- Create: `src/vcfixture/reference.py`
- Test: `tests/test_reference.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_reference.py
from pathlib import Path
import pysam
from vcfixture.reference import Reference

def _make_fasta(tmp_path: Path) -> Path:
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "ACGTACGTAC" * 5 + "\n")
    pysam.faidx(str(fa))
    return fa

def test_base_and_seq(tmp_path):
    ref = Reference(_make_fasta(tmp_path))
    assert ref.base("chr1", 0) == "A"   # 0-based
    assert ref.seq("chr1", 0, 4) == "ACGT"

def test_draw_snp_ref_matches_sequence(tmp_path):
    ref = Reference(_make_fasta(tmp_path))
    rec_ref, alts = ref.draw_ref_alt("chr1", pos0=0, klass="SNP", alt_index=1)
    assert rec_ref == "A"
    assert alts[0] != "A" and len(alts[0]) == 1

def test_draw_deletion_ref_starts_with_sequence(tmp_path):
    ref = Reference(_make_fasta(tmp_path))
    rec_ref, alts = ref.draw_ref_alt("chr1", pos0=0, klass="DEL", del_len=2)
    assert rec_ref == ref.seq("chr1", 0, 3)   # anchor + 2 deleted
    assert alts == [rec_ref[0]]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_reference.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/reference.py
from __future__ import annotations
import pysam

_BASES = "ACGT"

class Reference:
    """Thin pysam.FastaFile wrapper that draws spec-correct REF/ALT.

    Deterministic given its arguments (alt_index / del_len / ins_seq);
    Hypothesis supplies the randomness by choosing those arguments.
    """

    def __init__(self, fasta_path):
        self._fa = pysam.FastaFile(str(fasta_path))

    def base(self, contig: str, pos0: int) -> str:
        return self._fa.fetch(contig, pos0, pos0 + 1).upper()

    def seq(self, contig: str, start0: int, length: int) -> str:
        return self._fa.fetch(contig, start0, start0 + length).upper()

    def draw_ref_alt(self, contig, pos0, klass, *, alt_index=1,
                     del_len=1, ins_seq="T", mnp_len=2):
        if klass == "SNP":
            r = self.base(contig, pos0)
            alt = _BASES[(_BASES.index(r) + alt_index) % 4]
            return r, [alt]
        if klass == "MNP":
            r = self.seq(contig, pos0, mnp_len)
            alt = "".join(_BASES[(_BASES.index(b) + alt_index) % 4] for b in r)
            return r, [alt]
        if klass == "INS":
            anchor = self.base(contig, pos0)
            return anchor, [anchor + ins_seq]
        if klass == "DEL":
            r = self.seq(contig, pos0, del_len + 1)
            return r, [r[0]]
        if klass == "DELINS":
            r = self.seq(contig, pos0, mnp_len)
            return r, [ins_seq]
        if klass == "SPANNING_DEL":
            return self.base(contig, pos0), ["*"]
        raise ValueError(f"unknown class {klass!r}")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_reference.py -q`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/reference.py tests/test_reference.py
git commit -m "feat: add reference-aware REF/ALT adapter"
```

---

## Task 14: Hypothesis strategies + parametrize tables

**Files:**
- Create: `src/vcfixture/strategies.py`
- Test: `tests/test_strategies.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_strategies.py
from hypothesis import given, settings, HealthCheck
from vcfixture import strategies as S
from vcfixture.model import VcfDocument
from vcfixture._spec.number import NumberKind

def test_all_number_type_combos_table_is_exhaustive():
    combos = S.ALL_NUMBER_TYPE_COMBOS
    # Every classic Number kind appears.
    kinds = {n.kind for (n, t, kind) in combos}
    assert {NumberKind.FIXED, NumberKind.A, NumberKind.R,
            NumberKind.G, NumberKind.DOT, NumberKind.FLAG} <= kinds
    # Flag only pairs with INFO.
    for (n, t, kind) in combos:
        if t.value == "Flag":
            assert kind == "INFO"

def test_all_variant_classes_present():
    assert set(S.ALL_VARIANT_CLASSES) >= {
        "SNP", "MNP", "INS", "DEL", "DELINS", "SPANNING_DEL"}

@settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
@given(S.documents())
def test_documents_are_well_formed(doc: VcfDocument):
    assert isinstance(doc, VcfDocument)
    # Truth derivation must never raise on generated docs.
    t = doc.truth()
    assert t.genotypes.shape[0] == len(doc.records)
    # Serialization must never raise and must start with fileformat.
    assert doc.render().startswith("##fileformat=")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_strategies.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/strategies.py
from __future__ import annotations
from hypothesis import strategies as st
from .build import VcfBuilder
from .genotype import Genotype
from .variants import (snp, mnp, insertion, deletion, delins, spanning_deletion)
from ._spec.number import Number, NumberKind
from ._spec.types import Type

# --- exhaustive coverage tables (parametrize, not search) -------------------
def _build_number_type_combos():
    numbers = [Number.ONE, Number.fixed(2), Number.A, Number.R, Number.G,
               Number.DOT]
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

# --- value strategies (flat, construct-never-reject) ------------------------
_BASES = "ACGT"

@st.composite
def _ref_alt(draw, klass: str):
    b = draw(st.sampled_from(_BASES))
    b2 = draw(st.sampled_from(_BASES))
    if klass == "SNP":
        alt = _BASES[(_BASES.index(b) + 1 + draw(st.integers(0, 2))) % 4]
        return snp(b, alt)
    if klass == "MNP":
        return mnp(b + b2, _BASES[(_BASES.index(b) + 1) % 4]
                   + _BASES[(_BASES.index(b2) + 1) % 4])
    if klass == "INS":
        return insertion(b, draw(st.text(_BASES, min_size=1, max_size=3)))
    if klass == "DEL":
        return deletion(b, draw(st.text(_BASES, min_size=1, max_size=3)))
    if klass == "DELINS":
        return delins(b + b2, draw(st.text(_BASES, min_size=1, max_size=3)))
    return spanning_deletion(b)

@st.composite
def genotypes(draw, ploidy: int, n_alt: int, missing_rate: float = 0.1):
    """Construct a GT string directly — no rejection sampling."""
    alleles = []
    for _ in range(ploidy):
        if draw(st.floats(0, 1)) < missing_rate:
            alleles.append(".")
        else:
            alleles.append(str(draw(st.integers(0, n_alt))))
    phased = draw(st.booleans())
    sep = "|" if phased else "/"
    return sep.join(alleles)

@st.composite
def documents(draw, max_samples: int = 3, max_records: int = 4):
    n_samples = draw(st.integers(1, max_samples))
    samples = [f"s{i}" for i in range(n_samples)]
    ploidy = draw(st.integers(1, 2))
    b = VcfBuilder(samples=samples, contigs=[("chr1", 100000)]).fmt("GT")

    n_rec = draw(st.integers(1, max_records))
    pos = 1000
    for _ in range(n_rec):
        klass = draw(st.sampled_from(ALL_VARIANT_CLASSES))
        ref, alt = draw(_ref_alt(klass))
        gts = [draw(genotypes(ploidy, n_alt=1)) for _ in samples]
        b.record("chr1", pos, ref=ref, alt=[alt], gt=gts)
        pos += draw(st.integers(1, 50))
    return b.build()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_strategies.py -q`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/strategies.py tests/test_strategies.py
git commit -m "feat: add Hypothesis strategies and coverage tables"
```

---

## Task 15: Round-trip self-validation (cyvcf2 cross-check)

**Files:**
- Create: `tests/test_roundtrip.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_roundtrip.py
"""Self-validation: serialize -> parse with an INDEPENDENT parser -> the
third-party decode must match our derived GroundTruth. This makes the oracle
trustworthy rather than merely self-consistent."""
from pathlib import Path
import numpy as np
import pytest
from hypothesis import given, settings, HealthCheck
from vcfixture import strategies as S

cyvcf2 = pytest.importorskip("cyvcf2")

def _genos_from_cyvcf2(variant, n_samples, ploidy):
    # cyvcf2 genotypes(): list per sample [a1, a2, phased_bool]; -1 == missing
    out = np.full((n_samples, ploidy), -1, dtype=np.int32)
    for si, g in enumerate(variant.genotypes):
        for ai, allele in enumerate(g[:-1]):
            out[si, ai] = allele
    return out

@settings(max_examples=75, deadline=None,
          suppress_health_check=[HealthCheck.too_slow])
@given(S.documents())
def test_genotypes_round_trip_through_cyvcf2(tmp_path_factory, doc):
    truth = doc.truth()
    d: Path = tmp_path_factory.mktemp("rt")
    path = doc.write(d / "x.vcf.gz", bgzip=True, index=True)

    vf = cyvcf2.VCF(str(path))
    n_samples = len(doc.samples)
    ploidy = truth.genotypes.shape[2]
    for ri, variant in enumerate(vf):
        got = _genos_from_cyvcf2(variant, n_samples, ploidy)
        np.testing.assert_array_equal(got, truth.genotypes[ri],
            err_msg=f"genotype mismatch at record {ri}")
        assert variant.POS == int(truth.pos[ri])
        assert variant.REF == truth.ref[ri]
```

- [ ] **Step 2: Run test to verify it fails (then passes if code is correct)**

Run: `python -m pytest tests/test_roundtrip.py -q`
Expected: This test exercises Tasks 6–14. If it fails, the failure points to a serializer/truth discrepancy — fix the implicated module, not this test. Expected once correct: PASS.

- [ ] **Step 3: Triage any mismatch**

If cyvcf2 disagrees with truth, the bug is in `serialize.py` (wrong text) or `truth.py` (wrong oracle). Add a minimal failing case to the relevant unit test file (`test_serialize.py` or `test_truth.py`), fix the module, then re-run. Do NOT weaken this property test.

- [ ] **Step 4: Run to verify green**

Run: `python -m pytest tests/test_roundtrip.py -q`
Expected: PASS (1 passed)

- [ ] **Step 5: Commit**

```bash
git add tests/test_roundtrip.py
git commit -m "test: round-trip self-validation against cyvcf2"
```

---

## Task 16: Benchmark gate

**Files:**
- Create: `tests/test_benchmark.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_benchmark.py
"""Performance gate: drawing + serializing + deriving truth for a
representative document must stay cheap. Budget is generous to avoid CI
flakiness but tight enough to catch order-of-magnitude regressions."""
import time
from vcfixture.build import VcfBuilder
from vcfixture.genotype import Genotype

def _representative_doc():
    b = VcfBuilder(samples=[f"s{i}" for i in range(10)],
                   contigs=[("chr1", 1_000_000)]).fmt("GT")
    for k in range(50):
        gts = ["0|1" if (i + k) % 2 else "1/1" for i in range(10)]
        b.record("chr1", 1000 + k * 10, ref="A", alt=["T"], gt=gts)
    return b.build()

def test_build_serialize_truth_under_budget():
    t0 = time.perf_counter()
    for _ in range(100):
        doc = _representative_doc()
        doc.render()
        doc.truth()
    elapsed = time.perf_counter() - t0
    # 100 iterations of a 50-record x 10-sample doc.
    assert elapsed < 5.0, f"too slow: {elapsed:.2f}s for 100 iterations"
```

- [ ] **Step 2: Run test to verify it passes**

Run: `python -m pytest tests/test_benchmark.py -q`
Expected: PASS (1 passed). If it FAILS, profile `serialize.py`/`truth.py` hot paths before adjusting the budget.

- [ ] **Step 3: Commit**

```bash
git add tests/test_benchmark.py
git commit -m "test: add performance gate"
```

---

## Task 17: genoray parity (success criterion)

**Files:**
- Create: `tests/test_genoray_parity.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_genoray_parity.py
"""Reproduce genoray's hand-authored biallelic fixture via the builder and
assert the derived truth matches the numpy literals currently hand-coded in
genoray/tests/test_vcf.py. This proves vcfixture can replace those fixtures."""
import numpy as np
from vcfixture.build import VcfBuilder
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type

def _genoray_biallelic():
    # Mirrors genoray biallelic.vcf rows on chr1 at 81262/81262/81265
    # (1-based). genoray reads with phasing; DS dosage Number=A.
    return (VcfBuilder(samples=["sample1", "sample2"],
                       contigs=[("chr1", None), ("chr2", None), ("chr3", None)])
            .fmt("GT")
            .fmt("DS", Number.A, Type.FLOAT)
            .record("chr1", 81262, ref="GAT", alt=["A"],
                    gt=["0|1", "1|1"], DS=[[1.0], [2.0]])
            .record("chr1", 81262, ref="G", alt=["A"],
                    gt=["./.", "0/1"], DS=[[None], [1.0]])
            .record("chr1", 81265, ref="T", alt=["C"],
                    gt=["1|0", "./."], DS=[[0.9], [None]])
            .build())

def test_truth_matches_handcoded_genos_for_region():
    truth = _genoray_biallelic().truth()
    # genoray test_vcf.read_all() expects, for chr1:81261-81263 (records 0,1):
    #   genos (s p v) -> here we store (records, samples, ploidy)
    # record 0: s1 0|1, s2 1|1 ; record 1: s1 ./., s2 0/1
    np.testing.assert_array_equal(truth.genotypes[0], [[0, 1], [1, 1]])
    np.testing.assert_array_equal(truth.genotypes[1], [[-1, -1], [0, 1]])
    # phasing: record 0 both phased; record 1 s1 missing(unphased), s2 unphased
    np.testing.assert_array_equal(truth.phasing[0], [True, True])
    np.testing.assert_array_equal(truth.phasing[1], [False, False])
    # dosages echoed
    assert truth.format[0][0]["DS"] == [1.0]
    assert truth.format[1][0]["DS"] == [None]

def test_renders_and_round_trips(tmp_path):
    import pysam
    path = _genoray_biallelic().write(tmp_path / "g.vcf.gz",
                                      bgzip=True, index=True)
    vf = pysam.VariantFile(str(path))
    rows = list(vf.fetch("chr1", 81260, 81266))
    assert len(rows) == 3
```

- [ ] **Step 2: Run test to verify it passes**

Run: `python -m pytest tests/test_genoray_parity.py -q`
Expected: PASS (2 passed). If genotype/phasing assertions fail, the bug is in `truth.py` or `genotype.py`; fix there.

- [ ] **Step 3: Commit**

```bash
git add tests/test_genoray_parity.py
git commit -m "test: genoray biallelic fixture parity"
```

---

## Task 18: Public API + README

**Files:**
- Modify: `src/vcfixture/__init__.py`
- Create: `README.md`
- Test: `tests/test_public_api.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_public_api.py
import vcfixture

def test_public_exports():
    for name in ["VcfBuilder", "Number", "Type", "Genotype",
                 "GroundTruth", "Reference", "strategies"]:
        assert hasattr(vcfixture, name), f"missing public export: {name}"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_public_api.py -q`
Expected: FAIL with `AssertionError: missing public export: VcfBuilder`

- [ ] **Step 3: Write minimal implementation**

```python
# src/vcfixture/__init__.py
"""vcfixture — generate small VCF test data with decoded ground truth."""
from . import strategies
from .build import VcfBuilder
from .genotype import Genotype
from .reference import Reference
from .truth import GroundTruth
from ._spec.number import Number
from ._spec.types import Type

__version__ = "0.1.0"
__all__ = [
    "VcfBuilder", "Genotype", "Reference", "GroundTruth",
    "Number", "Type", "strategies", "__version__",
]
```

```markdown
<!-- README.md -->
# vcfixture

Generate small VCF (v4.5) test data — via an explicit builder or Hypothesis
strategies — with the decoded ground truth returned alongside, so parser tests
assert against a known oracle instead of hand-coded literals.

```python
from vcfixture import VcfBuilder, Number, Type

doc = (VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 100000)])
       .info("AF", Number.A, Type.FLOAT)
       .fmt("GT").fmt("DS", Number.A, Type.FLOAT)
       .record("chr1", 81262, ref="GAT", alt=["A"],
               gt=["0|1", "1|1"], info={"AF": [0.5]}, DS=[[1.0], [2.0]]))

text  = doc.render()                       # VCF text
truth = doc.truth()                        # GroundTruth (numpy genotypes, ...)
doc.write("x.vcf.gz", bgzip=True, index=True)
```

See `docs/superpowers/specs/2026-05-30-vcfixture-design.md` for the design.
```

- [ ] **Step 4: Run test + full suite**

Run: `python -m pytest -q`
Expected: PASS — entire suite green.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/__init__.py README.md tests/test_public_api.py
git commit -m "feat: public API exports and README"
```

---

## Self-Review

**Spec coverage:**
- VCF+decoded truth → Tasks 10, 17. Pure Python → all. Hypothesis-owns-randomness → Task 14 (strategies feed deterministic builder). Valid-only → builder validation (Task 11), no invalid path. Text+bgzip/index → Task 12. Optional reference-aware → Task 13. Classic Number×Type exhaustive → Tasks 2, 5, 14 (`ALL_NUMBER_TYPE_COMBOS`). Variant-class coverage → Tasks 8, 14 (`ALL_VARIANT_CLASSES`). Performance (small + fast) → Task 14 (flat/no-reject strategies) + Task 16 (benchmark gate). Self-validation vs independent parser → Task 15 (cyvcf2). Reserved curated subset → Task 5. Consumers cast dtype → truth uses single int32/bool width (Task 10).
- **Gap noted & accepted:** `ALL_NUMBER_TYPE_COMBOS` is built and asserted exhaustive (Task 14) but the `documents()` strategy only exercises GT+single-class records; wiring the full combo table into generated INFO/FORMAT *values* is deferred to a follow-up parametrized test once the round-trip backbone (Task 15) is green. This is intentional sequencing, not a missing requirement — the table and per-combo cardinality are already covered by Tasks 2/5/14.

**Placeholder scan:** No TBD/TODO; every code step shows complete code; every command has expected output.

**Type consistency:** `Number`/`NumberKind`/`Type` names stable across Tasks 2–14. `FieldDef(id, number, type, description, kind)` signature consistent (Tasks 4,5,11). `Record`/`VcfDocument` field names match between model (Task 7), serializer (Task 9), truth (Task 10), builder (Task 11). `Genotype(alleles, phased)` + `.is_phased`/`.ploidy`/`.parse`/`.render` consistent (Tasks 6,7,10,11). `derive_truth`/`render_document` delegation names match `model.py` lazy imports (Task 7) and their definitions (Tasks 9,10). `write_text`/`write_bgzip(..., index=)` match `model.write` (Tasks 7,12). `record_class`/`classify` match between variants (Task 8) and truth (Task 10).
