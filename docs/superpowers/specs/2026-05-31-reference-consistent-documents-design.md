# Reference-consistent `documents` strategy + variant labels — design

**Date:** 2026-05-31
**Repo:** vcfixture
**Branch / worktree:** `feat/reference-consistent-documents`
(`.worktrees/ref-consistent-docs`, based on `main` @ 0.2.1)
**Status:** approved design; implementation plan to follow.

## Problem & intent

vcfixture's current `strategies.documents` draws REF bases at random, so the
generated VCFs are **not reference-consistent**: their REF alleles do not match
any FASTA. Downstream tools that validate REF against a reference (`bcftools
norm`, `bcftools consensus`, gvl) reject them. GenVarLoader's Phase 2 test work
needs reference-consistent fixtures to drive property tests.

**The consuming intent (from GenVarLoader):** gvl *requires* every input
variant to already be **left-aligned, atomized, and bi-allelic** before
consumption, and assumes those invariants everywhere downstream. The point of
this capability is the **opposite** of feeding gvl clean input — it is to
generate *deliberately invalid* (non-canonical) but still **reference-consistent
and spec-valid** VCFs, so gvl's test suite can:

1. relax those invariants experimentally and measure throughput degradation, and
2. assert gvl fails — and characterize *how* it fails — on invalid input.

**Out of scope for both vcfixture and gvl:** implementing variant normalization.
vcfixture never re-implements `bcftools norm`. The haplotype oracle stays
`bcftools consensus` (which applies variants by POS/REF regardless of
left-alignment, as long as REF matches the reference — guaranteed by
reference-consistency), and lives downstream in gvl. The genotype/AF oracle
stays vcfixture's as-authored `GroundTruth`.

## Key decisions (resolved during brainstorming)

- **The strategy generates the reference** and emits paired output
  `(ReferenceSpec, VcfDocument, GroundTruth)`. Owning the reference bases lets
  the strategy construct the exact contexts (e.g. tandem repeats) that make a
  variant non-canonical.
- **Reference construction is a mutable `ReferenceBuilder` → frozen
  `ReferenceSpec`.** Multi-nucleotide writes and tandem repeats are first-class,
  because tandem repeats are the contexts where left-alignment matters.
- **A general per-variant `labels` mechanism**, domain-agnostic: vcfixture
  attaches descriptive provenance labels for what it generated; gvl assigns its
  own invariant-semantics labels. Labels are **never serialized** into the VCF.
- **Composable strategy layer:** `references()` + `documents(reference=…)` +
  a thin `reference_and_documents()` wrapper, so a reference can be drawn once
  and reused across many drawn documents.
- **vcfixture's own tests do not run bcftools** (not a dev dep); they assert
  reference-consistency in-process.

## Components

### 1. Reference construction — `reference.py`

Two new public types alongside the existing pysam file-backed `Reference`
(unchanged):

- **`ReferenceBuilder`** (mutable):
  - `__init__(seed: int = 0)`
  - `add_contig(id: str, length: int)` — random-fills ACGT from a seeded numpy
    `default_rng`.
  - `set_base(contig, pos0, base)` — single-base overwrite (0-based).
  - `set_seq(contig, pos0, seq)` — multi-nucleotide run overwrite.
  - `tandem_repeat(contig, pos0, motif, n)` — writes `motif * n` at `pos0` and
    records a `RepeatFeature(contig, pos0, motif, n)` so the built spec can
    advertise it.
  - `build() -> ReferenceSpec`.
- **`ReferenceSpec`** (frozen dataclass):
  - `contigs: tuple[tuple[str, str], ...]` — `(id, sequence)`.
  - `repeats: tuple[RepeatFeature, ...]` — planted-repeat provenance.
  - `base(contig, pos0) -> str`, `seq(contig, start0, length) -> str`.
  - `write(path, *, bgzip: bool = True, index: bool = True) -> Path` — writes a
    60-col FASTA, then bgzip + faidx via pysam (reuse `io.py` patterns). Returns
    the written path.
- **`RepeatFeature`** (frozen): `contig: str, pos0: int, motif: str, count: int`.

The reference-consistent REF/ALT logic currently in `Reference.draw_ref_alt` is
factored into a module-level helper parameterized over a `(base, seq)` accessor,
so both `Reference` (file-backed) and `ReferenceSpec` (in-memory) reuse it with
no duplication. `Reference.draw_ref_alt`'s public behavior is preserved.

### 2. General variant labels — `model.py` / `build.py` / `truth.py`

- `Record` gains `labels: frozenset[str] = frozenset()`.
- `VcfBuilder.record(..., labels: Iterable[str] | None = None)` — stores the
  frozenset on the record.
- `GroundTruth` gains `labels: list[frozenset[str]]` (one entry per record);
  `derive_truth` populates it from `rec.labels`.
- **Labels are out-of-band test metadata and are never written to the VCF.**
  `VcfDocument.render()` output is unchanged. vcfixture attaches no semantics;
  gvl interprets the strings.

### 3. Strategies — `strategies.py`

- **`references(*, max_contigs=2, max_contig_len=2000, max_repeats=3)`** → draws
  a `ReferenceSpec`. Draws a seed integer + small contigs (length drawn up to
  `max_contig_len`), optionally plants
  tandem repeats at drawn loci (motif and count drawn within bounds), advertised
  via `spec.repeats`. Contigs are kept small to keep Hypothesis examples small.
- **`documents(max_samples=3, max_records=4, max_alt=1, *, reference=None,
  violations=frozenset(), label_overrides=None)`** (existing positional params
  unchanged; new params are keyword-only and additive):
  - `reference=None` → existing reference-free behavior preserved exactly
    (back-compat for current tests/consumers).
  - `reference=<ReferenceSpec>` → REFs are drawn from the spec
    (reference-consistent), reusing the shared REF/ALT helper.
  - `violations ⊆ {"multiallelic", "non_atomic", "non_left_aligned"}` — opt into
    each non-canonical class. Enabled violations auto-tag descriptive provenance
    labels on the affected records: `"multiallelic"`, `"non_atomic"`,
    `"off_anchor"`, `"tandem_repeat"`. `"non_left_aligned"` places a deletion off
    the leftmost copy of a planted repeat drawn from `spec.repeats`.
  - `label_overrides: Mapping[str, str] | None` — remap the default provenance
    label strings; gvl may instead/also add its own labels downstream.
- **`reference_and_documents(…)`** → composes the two and returns
  `(ReferenceSpec, VcfDocument, GroundTruth)`.

**Coupling resolution:** a non-left-aligned indel requires a repeat at its
locus. `references()` plants repeats and `ReferenceSpec` *advertises* them via
`spec.repeats`, so the composable `documents(reference=spec)` can place
violations against a spec it did not itself build. (Rejected alternative: fuse
planting + placement into only the paired wrapper — loses standalone
`documents(reference=spec)` violation support.)

### 4. Public API + release — `__init__.py`, `pyproject.toml`

- Export `ReferenceBuilder`, `ReferenceSpec`, `RepeatFeature` from `__init__.py`.
  `Reference`, `strategies`, `VcfBuilder`, `Genotype`, `GroundTruth`, `Number`,
  `Type` stay exported.
- Conventional-commit feature work → commitizen **minor bump 0.2.1 → 0.3.0** +
  changelog. Actual PyPI release is a separate manual step after merge; gvl
  consumes via a version bump (no path/git pins).

## Testing (TDD)

vcfixture's own suite stays bcftools-free.

- `test_reference.py` (extend): `ReferenceBuilder.set_base/set_seq/tandem_repeat`
  produce the expected bases; `ReferenceSpec.base/seq`; `ReferenceSpec.write`
  roundtrips (read back via `pysam.FastaFile`); `repeats` provenance is correct;
  the shared REF/ALT helper still satisfies the existing `draw_ref_alt` tests.
- `test_model.py` / `test_build.py` / `test_truth.py`: `labels` carried
  `record → truth`; `labels` absent from `render()` output; default empty set.
- `test_strategies.py` (extend): `references()` specs are self-consistent
  (planted repeats present at advertised loci); `documents(reference=spec)` REFs
  equal `spec.seq(...)` at each record (assert in-process, no bcftools); each
  violation kind yields its provenance label; `reference_and_documents()` tuple
  is well-formed and `truth().genotypes.shape[0] == len(doc.records)`;
  reference-free `documents()` behavior is unchanged.

## Architecture notes

- **Determinism under Hypothesis:** `references()` draws a seed integer and
  structural parameters via Hypothesis (not a fixed seed), so shrinking and
  reproducibility work through Hypothesis' own machinery.
- **No I/O in strategies:** strategies are pure value producers; FASTA
  materialization happens only when a consumer calls `ReferenceSpec.write(...)`.
- **Back-compat:** `documents()`/`documents_with_fields()` keep their current
  signatures and reference-free behavior; `reference=` and `violations=` are
  additive keyword-only parameters with inert defaults.

## Downstream (GenVarLoader Phase 2 — not part of this work)

gvl consumes the released `vcfixture>=0.3.0` to write a property-test module:
per Hypothesis example, draw `(reference, document, truth)`, write the FASTA +
VCF, run `bcftools norm`/`consensus` as the haplotype oracle, `gvl.write`, then
assert gvl haplotypes == consensus and gvl genotypes/AF == `GroundTruth`; use
`labels` to partition canonical vs. deliberately-invalid examples and assert
gvl's failure/degradation behavior on the latter.
