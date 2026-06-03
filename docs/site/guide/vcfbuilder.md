# VcfBuilder guide

`VcfBuilder` is the explicit-construction front-end for vcfixture. Use it when
you want a specific, named fixture — a particular SNP at a known position, a
deletion with a known genotype matrix, a symbolic SV with exact SVLEN — and you
want to assert against decoded ground truth rather than hand-coded literals.

## Mental model

Everything in vcfixture converges on one immutable hub: `VcfDocument`. Two
front-ends produce it (the builder and Hypothesis strategies); two back-ends
consume it (the serializer and the truth deriver).

```
 VcfBuilder.record(...)  ──┐
                           ├──▶  VcfDocument  ──▶  render()  →  VCF text
 strategies.documents()  ──┘     (immutable)  ──▶  truth()   →  GroundTruth
```

Because `GroundTruth` is derived from the same model object that was serialized,
the expected values are never out of sync with the bytes on disk. You get the
oracle for free.

## Construct a builder

Pass your sample names and contigs. Contigs take `(name, length)` pairs; use
`None` for the length when it is not known.

```python
from vcfixture import VcfBuilder

b = VcfBuilder(
    samples=["s1", "s2"],
    contigs=[("chr1", 248956422), ("chr2", 242193529)],
)
```

The builder is mutable. Every declaration and record call returns `self`, so you
can chain them.

## Declare fields

Declare every INFO and FORMAT field you intend to use before adding records.
Reserved fields (`AC`, `AF`, `AN`, `DP`, `AD`, `PL`, `GL`, `GT`, `GQ`, `PS`,
`CN`, `LEN`, `DB`, `H2`, `AA`, `END`, `SVLEN`, `SVCLAIM`, `SVTYPE`, `CIPOS`,
`CIEND`, `CILEN`, `MATEID`, `PARID`, `IMPRECISE`) resolve their `Number` and
`Type` automatically — just pass the ID.

```python
from vcfixture import VcfBuilder, Number, Type

b = (
    VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 100_000)])
    # Reserved fields — Number/Type auto-resolved:
    .info("AF")         # Number=A, Type=FLOAT
    .info("AN")         # Number=1, Type=INTEGER
    .fmt("GT")          # Number=1, Type=STRING
    .fmt("AD")          # Number=R, Type=INTEGER
    .fmt("GQ")          # Number=1, Type=INTEGER
    # Custom field — Number/Type required:
    .fmt("DS", Number.A, Type.FLOAT)
    # FILTER and ALT metadata:
    .filter("LowQual", "Low-quality call")
    .alt("DEL", "Deletion relative to the reference")
)
```

Attempting to use an undeclared field ID in a record raises `ValueError`
immediately.

## Add records

Call `.record(chrom, pos, ...)` to append a variant. ALT alleles must be typed
`Allele` instances — plain strings are not accepted. Genotypes are passed as
`"0|1"`-style strings (phased with `|`, unphased with `/`, missing as `./.`).
Extra keyword arguments become per-sample FORMAT values.

```python
from vcfixture import VcfBuilder, Number, Type, Seq

doc = (
    VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 100_000)])
    .info("AF", Number.A, Type.FLOAT)
    .fmt("GT")
    .fmt("DS", Number.A, Type.FLOAT)
    .record(
        "chr1", 81262,
        ref="GAT",
        alt=[Seq("A")],          # sequence ALT — typed instance required
        gt=["0|1", "1|1"],       # one string per sample
        info={"AF": [0.5]},
        DS=[[1.0], [2.0]],       # per-sample, inner list is per-ALT
    )
)
text  = doc.render()
truth = doc.truth()
```

(`doc` is still the builder here — `render()`, `truth()`, and `write()` are builder methods that act on the document it has accumulated.)

### Allele types

| Class | Alias | Example | Serializes as |
|-------|-------|---------|---------------|
| `SequenceAllele` | `Seq` | `Seq("A")` | `A` |
| `SymbolicAllele` | `Sym` | `Sym.deletion()` | `<DEL>` |
| `SpanningDeletion` | `Star` | `Star()` | `*` |
| `UnspecifiedAllele` | `Unspecified` | `Unspecified()` | `<*>` |
| `BreakendAllele` | `Bnd` | `Bnd.parse("T[chr2:5[")` | `T[chr2:5[` |

`SymbolicAllele` has class methods for the five SV first-types:
`Sym.deletion()`, `Sym.insertion()`, `Sym.duplication()`, `Sym.inversion()`,
`Sym.cnv()`. Each accepts optional subtypes, e.g. `Sym.duplication("TANDEM")`.

## Outputs

### `render()` — VCF text

```python
from vcfixture import VcfBuilder, Seq

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 1_000)])
    .fmt("GT")
    .record("chr1", 10, ref="A", alt=[Seq("T")], gt=["0/1"])
)
text = doc.render()
assert "##fileformat=VCFv4.5" in text
assert "chr1\t10" in text
```

### `truth()` — `GroundTruth`

`GroundTruth` is a frozen dataclass. Key fields:

- `genotypes` — `(records, samples, ploidy)` int32 ndarray; **`-1` = missing**.
- `phasing` — `(records, samples)` bool ndarray; `True` iff fully phased.
- `pos` — `(records,)` int64 ndarray; 1-based.
- `ref` / `alts` — lists of strings.
- `variant_class` — one label per record (see [Variant classes](#variant-classes)).
- `info` / `format` — lists of `{id: value}` dicts mirroring the VCF fields.
- `alts_truth` — per-ALT `AlleleTruth` with `.kind`, `.is_sequence`, `.sv_type`,
  `.svlen`, `.sv_end`.

```python
import numpy as np
from vcfixture import VcfBuilder, Number, Type, Seq

doc = (
    VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 100_000)])
    .info("AF", Number.A, Type.FLOAT)
    .fmt("GT")
    .fmt("DS", Number.A, Type.FLOAT)
    .record(
        "chr1", 81262,
        ref="GAT", alt=[Seq("A")],
        gt=["0|1", "1|1"],
        info={"AF": [0.5]},
        DS=[[1.0], [2.0]],
    )
)
t = doc.truth()
np.testing.assert_array_equal(t.genotypes[0], [[0, 1], [1, 1]])
assert t.phasing[0, 0] == True   # phased
assert t.info[0]["AF"] == [0.5]
assert t.format[0][0]["DS"] == [1.0]
```

### `write(path, *, bgzip=False, index=False)`

Write to disk. Returns the `Path` that was written.

```python
import tempfile, pathlib
from vcfixture import VcfBuilder, Seq

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 1_000)])
    .fmt("GT")
    .record("chr1", 10, ref="A", alt=[Seq("T")], gt=["0/1"])
)

with tempfile.TemporaryDirectory() as d:
    # Plain text VCF:
    p = doc.write(pathlib.Path(d) / "out.vcf")
    assert p.exists()

    # bgzip-compressed + tabix-indexed:
    p2 = doc.write(pathlib.Path(d) / "out.vcf.gz", bgzip=True, index=True)
    assert p2.exists()
```

## Eager validation

The builder validates at `.record()` time, not at build time. Errors surface
as `ValueError` immediately so the source of the mistake is obvious.

**Undeclared field used in a record:**
```python
# Raises ValueError: FORMAT field 'DS' not declared
from vcfixture import VcfBuilder, Seq
b = VcfBuilder(samples=["s1"], contigs=[("chr1", None)]).fmt("GT")
# b.record("chr1", 1, ref="A", alt=[Seq("T")], gt=["0/1"], DS=[[1.0]])
```

**GT allele index out of range:**
```python
# Raises ValueError: allele index 5 out of range (n_alt=1)
# b.record("chr1", 1, ref="A", alt=[Seq("T")], gt=["0/5"])
```

**Value count does not match resolved cardinality:**
```python
# AD is Number=R (one per allele including REF), so biallelic needs 2 values.
# Raises ValueError: AD cardinality mismatch: expected 2, got 1
from vcfixture import VcfBuilder, Number, Type, Seq
b = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", None)])
    .fmt("GT").fmt("AD", Number.R, Type.INTEGER)
)
# b.record("chr1", 1, ref="A", alt=[Seq("T")], gt=["0/1"], AD=[[5]])
```

**Symbolic SV validation:**

- Symbolic and breakend ALTs require a single-base REF padding base.
- Every symbolic allele requires `SVLEN` in INFO.
- `<DEL>` and `<DUP>` require `SVCLAIM` (values `"D"`, `"J"`, or `"DJ"`).
- `SVLEN` must be absent for breakend, spanning-deletion, and unspecified alleles.

## Variant classes

The `GroundTruth.variant_class` list assigns one label per record. For
single-ALT records the label reflects the allele type; for multi-ALT records it
is always `"MULTIALLELIC"`.

### SNP

Single-base substitution.

```python
from vcfixture import VcfBuilder, Seq

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 1_000)])
    .fmt("GT")
    .record("chr1", 100, ref="A", alt=[Seq("T")], gt=["0/1"])
)
assert doc.truth().variant_class == ["SNP"]
```

### MNP

Multi-nucleotide polymorphism — `len(ref) == len(alt) > 1`.

```python
from vcfixture import VcfBuilder, Seq

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 1_000)])
    .fmt("GT")
    .record("chr1", 100, ref="AC", alt=[Seq("TG")], gt=["0/1"])
)
assert doc.truth().variant_class == ["MNP"]
```

### Insertion

`alt` starts with `ref` and is longer — left-anchored convention.

```python
from vcfixture import VcfBuilder, Seq

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 1_000)])
    .fmt("GT")
    .record("chr1", 100, ref="A", alt=[Seq("ATG")], gt=["0/1"])
)
assert doc.truth().variant_class == ["INS"]
```

### Deletion

`ref` starts with `alt` and is longer — left-anchored convention.

```python
from vcfixture import VcfBuilder, Seq

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 1_000)])
    .fmt("GT")
    .record("chr1", 100, ref="ATG", alt=[Seq("A")], gt=["0/1"])
)
assert doc.truth().variant_class == ["DEL"]
```

### Complex substitution (DELINS)

Different lengths, no shared prefix — the catch-all for complex alleles.

```python
from vcfixture import VcfBuilder, Seq

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 1_000)])
    .fmt("GT")
    .record("chr1", 100, ref="GAT", alt=[Seq("A")], gt=["0/1"])
)
assert doc.truth().variant_class == ["DELINS"]
```

### Spanning deletion

The `*` allele: overlaps a deletion from an upstream record.

```python
from vcfixture import VcfBuilder, Star

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 1_000)])
    .fmt("GT")
    .record("chr1", 100, ref="A", alt=[Star()], gt=["0/1"])
)
assert doc.truth().variant_class == ["SPANNING_DEL"]
```

### Multiallelic

Any record with more than one ALT allele.

```python
from vcfixture import VcfBuilder, Seq

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 1_000)])
    .fmt("GT")
    .record("chr1", 100, ref="A", alt=[Seq("T"), Seq("G")], gt=["0/1"])
)
assert doc.truth().variant_class == ["MULTIALLELIC"]
```

### Symbolic SV

Symbolic alleles (`<DEL>`, `<DUP>`, `<INS>`, `<INV>`, `<CNV>`) need a
single-base REF anchor, `SVLEN` in INFO, and — for `<DEL>` and `<DUP>` —
`SVCLAIM` (`"D"`, `"J"`, or `"DJ"`).

```python
from vcfixture import VcfBuilder, Sym

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 100_000)])
    .fmt("GT")
    .info("SVLEN")    # reserved: Number=A, Type=INTEGER
    .info("SVCLAIM")  # reserved: Number=A, Type=STRING
    .record(
        "chr1", 1000,
        ref="G",
        alt=[Sym.deletion()],
        gt=["0/1"],
        info={"SVLEN": [500], "SVCLAIM": ["DJ"]},
    )
)
t = doc.truth()
assert t.variant_class == ["SV_DEL"]
at = t.alts_truth[0][0]
assert at.kind == "SYMBOLIC"
assert at.sv_type == "DEL"
assert at.svlen == 500
assert at.sv_end == 1500   # POS + SVLEN
```

The five class methods on `SymbolicAllele` / `Sym`:

```python
from vcfixture import Sym

Sym.deletion()             # <DEL>
Sym.insertion()            # <INS>
Sym.duplication()          # <DUP>
Sym.inversion()            # <INV>
Sym.cnv()                  # <CNV>
Sym.duplication("TANDEM")  # <DUP:TANDEM>
```

### Breakend

Breakend alleles use `Bnd.parse(s)` with any valid VCF 4 replacement string.

```python
from vcfixture import VcfBuilder, Number, Type, Bnd

doc = (
    VcfBuilder(samples=["s1"], contigs=[("chr1", 1_000), ("chr2", 1_000)])
    .fmt("GT")
    .info("MATEID", Number.A, Type.STRING)
    .record(
        "chr1", 100,
        ref="T",
        alt=[Bnd.parse("T[chr2:5[")],
        gt=["0/1"],
        info={"MATEID": ["bnd_mate"]},
    )
)
t = doc.truth()
assert t.variant_class == ["BND"]
assert t.alts_truth[0][0].kind == "BND"
assert t.alts_truth[0][0].is_sequence is False
```

## Why the truth is trustworthy

`GroundTruth` is derived from the same `VcfDocument` that is serialized to text,
so the oracle and the bytes are always in sync. That alone makes it
self-consistent — but self-consistency is not enough: a bug that serializes
incorrectly *and* derives truth incorrectly would pass silently.

The test suite breaks that symmetry with an independent cross-check: every
generated document is serialized to a temporary file, re-parsed by a third-party
library (cyvcf2 or pysam), and the decoded values are asserted against
`GroundTruth`. Because the parser has no knowledge of vcfixture internals, any
discrepancy between the serialized bytes and the oracle is caught immediately.
There is also a Hypothesis property test that runs this round-trip for arbitrary
documents.

This is what makes the oracle trustworthy rather than merely self-consistent.
