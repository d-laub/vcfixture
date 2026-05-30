# vcfixture — Design

**Date:** 2026-05-30
**Status:** Approved design; ready for implementation planning.

## Purpose

A pure-Python library for generating VCF test data — both for **specific named
test cases** (an explicit builder API) and for **property-based / fuzz testing**
(Hypothesis strategies). Every generated VCF comes with its **decoded ground
truth**, so consumers assert parser output against a known oracle instead of
hand-deriving expected values.

### Problem it solves

Two downstream repos — `GenVarLoader` and `genoray` — currently maintain small
VCF files as on-disk fixtures and **hand-code the expected decoded values** as
numpy literals in their test files (e.g. `genos = np.array([[[0,-1],[1,-1]], ...])`,
`dosages = [[1.0, nan],[2.0,1.0]]`). The VCF text and the expected values drift
independently, coverage of the INFO/FORMAT `Number`×`Type` matrix is ad hoc, and
maintaining many tiny VCFs is a smell. `genoray` derives expected arrays by hand
from VCF text; `GenVarLoader` computes an oracle via `bcftools norm` + `bcftools
consensus`.

vcfixture cures this by making the generator the source of truth: because it
constructs the VCF, it knows the decoded semantics by construction.

## Key decisions (from brainstorming)

1. **Output = VCF + decoded truth.** The generator returns the VCF *and* a
   structured ground-truth object (genotype matrices, phasing, resolved
   INFO/FORMAT values, missing sentinels). This is what removes the drift smell.
2. **Pure Python, no Rust/PyO3.** Rust was the original premise but is YAGNI
   here: Hypothesis provides shrinking for free, Python's type system models the
   `Number`×`Type` matrix adequately, and a Rust generator would only pay off if
   the downstream consumers were also Rust (they are not). Generating a few-KB
   test VCF is not performance-sensitive. Revisit only if profiling proves
   generation is a real bottleneck.
3. **Hypothesis owns randomness; the model is deterministic.** Strategies draw
   parameters and feed deterministic builders, so shrinking works end-to-end.
4. **Valid-only generation.** Only spec-conformant VCFs in v1. Negative testing
   becomes a future mutation layer over valid output.
5. **Text + bgzip/index artifacts.** Emit `.vcf` text always; optionally write
   bgzipped `.vcf.gz` with a `.csi`/`.tbi` index via pysam. Downstream
   PGEN/SVAR/consensus transcoding stays in each consumer's own scripts.
6. **Optional reference-aware mode.** Reference-free by default (arbitrary REF —
   fine for genoray, which only reads genotypes). When given a FASTA + region,
   draw POS and spec-correct REF/ALT from the actual sequence so GVL's
   `bcftools norm`/`consensus` oracle accepts the output.
7. **Classic Number vocabulary, exhaustive.** Cover `1`, fixed-N, `A`, `R`, `G`,
   `.`, and `0`(Flag) across all Types, with correct per-record cardinality.
   Localized-allele codes (`LA`/`LR`/`LG`) are out of scope. **Field *names* are
   secondary** — a curated subset of reserved fields suffices; the priority is
   exhaustive coverage of the Number×Type matrix.
8. **Variant-class coverage is first-class** (see below), not a toggle
   afterthought: SNP, MNP, all indel classes, complex/non-atomic substitutions,
   and spanning deletions are explicit, individually-addressable generation
   targets with their own truth.
9. **Small outputs, optimized hot path** (see Performance): all generated data is
   tiny, and the Hypothesis strategy path is engineered to avoid the slowness
   that complex-type fuzzing usually incurs.

## Spec grounding (VCF 4.5)

Modeled against
`https://github.com/samtools/hts-specs/blob/master/VCFv4.5.tex`:

- **Number values:** integer (fixed), `A` (per alt), `R` (per allele incl. ref),
  `G` (per genotype), `.` (unknown/variable), `0` (Flag only).
- **Cardinality:** `A = n_alt`; `R = n_alt + 1`;
  `G = C(n_alt + 1 + ploidy − 1, ploidy)` (binomial over alleles-incl-ref and
  ploidy). Examples: diploid biallelic → 3; diploid triallelic → 6; triploid
  triallelic → 10. The genotype ordering enumerator is required to place
  `PL`/`GL` values correctly.
- **Types:** INFO = Integer, Float, Flag, Character, String; FORMAT = same minus
  Flag. **Flag must be Number=0 and INFO-only.**
- **Missing:** single `.`; arrays may use a single `.` or a list of dots
  (`.,.,.`). Zero-length fields are illegal — a dot is required. Trailing FORMAT
  fields may be dropped except GT.
- **GT:** allele indices separated by `/` (unphased) or `|` (phased); `0`=REF,
  `1..n`=ALTs in order, `.`=missing allele, `*`=spanning deletion ALT. The first
  phasing indicator may be omitted (implicitly `/` if any separator is `/`, else
  `|`).
- **Reserved fields** with canonical Number/Type: INFO `AA`(1,String),
  `AC`(A,Integer), `AF`(A,Float), `AN`(1,Integer), `DP`(1,Integer), `DB`(0,Flag),
  `H2`(0,Flag); FORMAT `GT`, `GQ`(1,Integer), `DP`(1,Integer), `AD`(R,Integer),
  `PL`(G,Integer), `GL`(G,Float), `PS`(1,Integer).
- **Header:** `##fileformat=VCFv4.5` must be the first line; other meta lines in
  any order; `#CHROM ... INFO [FORMAT samples...]` header line; UTF-8, no BOM;
  percent-encoding for reserved characters.

## Variant-class coverage (first-class)

`REF`/`ALT` generation is organized around an explicit, individually-addressable
taxonomy. Each class is its own strategy and builder shortcut, and each carries
unambiguous truth (REF length, ALT length, and the implied edit). Reference-aware
mode draws REF from the sequence; reference-free mode constructs abstract bases.

| Class | Shape (REF→ALT) | Example | Notes |
|-------|-----------------|---------|-------|
| **SNP** | 1→1, differ | `A→T` | |
| **MNP** | k→k, k>1, every base may differ | `AC→GT` | non-atomic substitution |
| **Insertion** | anchor base + inserted seq | `A→ACGT` | left-anchored per spec |
| **Deletion** | anchor + deleted seq → anchor | `ACGT→A` | |
| **Complex / delins (non-atomic)** | m→n, m≠n, m,n≥1, not pure ins/del | `AC→GTA` | combined substitution+indel |
| **Spanning deletion** | ALT = `*` | `A→*` | overlapping upstream deletion |
| **Multiallelic** | multiple ALTs of mixed classes at one site | `G→A,C` / `GAT→A,GA` | drives `A`/`R`/`G` cardinality |

"Non-atomic" = variants that are not minimal/left-aligned single edits (MNPs and
delins), which exercise the consumers' normalization and span logic. The
generator can emit them pre- or post-normalization; reference-aware output stays
acceptable to `bcftools norm`. Strategies expose per-class weights so a test can
say "only indels" or "mix all classes."

## Performance

All generated data is intentionally **small** — bounded sample counts, ploidy,
record counts, and allele/sequence lengths. There is no use case for large
outputs; size bounds are tight defaults, not just caps.

Hypothesis fuzzing over complex types is historically slow, so the strategy path
is engineered to avoid the usual culprits:

- **Construct, never reject.** No `.filter()` / `assume()` rejection sampling on
  the hot path — every draw is mapped to a valid construction directly (e.g. draw
  ploidy and allele count, then *compute* genotypes, rather than draw-and-reject).
- **Flat, non-recursive strategies.** Avoid recursive/deeply-nested strategies;
  prefer flat integer/index draws composed into structures, which shrink and
  generate far faster.
- **Matrix coverage via parametrization, not search.** Exhaustive Number×Type and
  variant-class coverage is driven by `@pytest.mark.parametrize` over enumerated
  combos, reserving Hypothesis for the *values within* a fixed combo — this keeps
  each Hypothesis search space small and fast instead of asking it to explore the
  whole matrix.
- **Validation-skipping fast constructor.** The eager builder validation runs for
  the user-facing builder; the strategy path uses an internal constructor that
  skips re-validation, since strategies are valid by construction. The model's hot
  path stays allocation-light.
- **Benchmark gate.** The suite includes a micro-benchmark asserting that drawing
  + serializing + deriving truth for a representative document stays within a
  target budget, so regressions are caught.

One immutable **document model** is the hub. The builder and the Hypothesis
strategies are two *front-ends* producing the model; the serializer and the
truth-deriver are two *back-ends* consuming it. This guarantees the builder and
fuzzer cannot diverge, and makes ground truth free.

```
 Builder kwargs ─┐                          ┌─→ Serializer ─→ .vcf text ─→ [bgzip + .csi/.tbi]
                 ├─→  VcfDocument (model) ───┤
 Hypothesis draws┘     (immutable)           └─→ Truth deriver ─→ GroundTruth (numpy + dicts)
```

### Components

1. **Spec layer** (`_spec/`)
   - `Number`: `ONE | FIXED(n) | A | R | G | DOT | FLAG`.
   - `Type`: `Integer | Float | Flag | Character | String`.
   - `FieldDef`: `id`, `number`, `type`, `description`, `kind` (INFO|FORMAT).
   - `cardinality(number, n_alt, ploidy) -> int | None` (None for `.`/variable).
   - Genotype-ordering enumerator (for `Number=G` placement of PL/GL).
   - Validity invariants: Flag ⇒ Number=0 & INFO-only; no zero-length fields;
     INFO-key regex `^([A-Za-z_][0-9A-Za-z_.]*|1000G)$`; no duplicate IDs/samples.
   - **Reserved-field registry** so builders/strategies can reference fields by
     name with canonical Number/Type.

2. **Document model** (`model.py`) — immutable dataclasses:
   - `Genotype`: allele indices (`int | None` for missing), per-separator phasing.
   - `Record`: chrom, pos, id, ref, alts, qual, filter, info dict, format key
     order, per-sample format values.
   - `VcfDocument`: fileformat, declared field defs, contigs, sample names,
     records.

3. **Serializer** (`serialize.py`) — model → spec-correct text. Handles missing
   `.`, multi-dot arrays, GT first-phasing-omission, trailing-FORMAT dropping
   (never GT), percent-encoding.

4. **Truth deriver** (`truth.py`) — model → `GroundTruth`:
   - Genotype matrix as `np.int*`, **`-1` = missing** (matches both consumers).
   - Phasing matrix as `bool`.
   - Per-record INFO dict and per-sample FORMAT dict, cardinality resolved.

5. **Builder API** (`build.py`) — explicit construction for named tests:
   ```python
   doc = (VcfBuilder(samples=["s1", "s2"], contigs=["chr1", "chr2"])
          .info("AF", Number.A, Type.Float)
          .fmt("GT").fmt("DS", Number.A, Type.Float)
          .record("chr1", 81262, ref="GAT", alt=["A"],
                  gt=["0|1", "1|1"], DS=[[1.0], [2.0]]))
   text  = doc.render()                       # str
   truth = doc.truth()                        # GroundTruth
   doc.write("x.vcf.gz", bgzip=True, index=True)
   ```
   Eager validation: Flag/Number mismatch, undefined ID used in a record, GT
   index out of range, value count ≠ resolved cardinality.

6. **Hypothesis strategies** (`strategies.py`) — `@composite` strategies feeding
   the *same* builders, so shrinking works end-to-end. Knobs: ploidy, n_samples,
   contigs, chosen Number×Type combos, phasing, missingness rate,
   multiallelic/indel/spanning-deletion toggles, optional reference. Plus
   `all_number_type_combos()` for exhaustive `@parametrize` coverage of the
   matrix.

7. **Reference adapter** (`reference.py`, optional) — given FASTA + region, draws
   POS + matching REF + realistic ALTs (SNP/indel/MNP/spanning-del) so GVL's
   `bcftools norm`/`consensus` accepts the output. Reference-free by default.

8. **IO** — text always; bgzip + `.csi`/`.tbi` via pysam on request.

No CLI in v1 (YAGNI).

## Data flow

Parameters (builder kwargs *or* Hypothesis draws) → `VcfDocument` →
{Serializer → VCF text → optional bgzip/index} **and** {Truth deriver →
`GroundTruth`}. Consumer tests compare parser output against `GroundTruth`.

## Error handling

- Builder validates **eagerly** and raises on invariant violations (Flag with
  Number≠0, undefined INFO/FORMAT ID in a record, GT allele index out of range,
  value count ≠ resolved cardinality).
- Strategies produce valid documents by construction.
- Reference mode raises if a drawn REF cannot match the sequence (should not
  occur, since REF is read from the reference).

## Self-validation (how vcfixture earns trust)

vcfixture's own test suite **round-trips every document through an independent
parser** (`pysam`/`cyvcf2`): serialize → parse → assert the third-party decode
matches our `GroundTruth`. A Hypothesis property test asserts this for arbitrary
documents. This makes the oracle trustworthy: truth is cross-checked against a
real VCF parser, not merely self-consistent. The bgzip/index path is validated
by opening the indexed file and performing a region query.

## Scope boundaries (explicitly NOT building in v1)

- Invalid/malformed generation (future mutation layer over valid output).
- Localized-allele codes `LA`/`LR`/`LG` and the `LAA` field.
- PGEN/SVAR/consensus transcoding (stays in consumer repos).
- A CLI.

## Dependencies

`numpy`, `hypothesis`, `pysam` (bgzip/index + reference reading and round-trip
validation). `cyvcf2` optionally as a second independent parser for
cross-checking.

## Success criteria

- A genoray-style fixture (biallelic with `DS` dosages, phasing, missing/half
  calls, multiallelic, sample reorder) can be expressed via the builder, and its
  `GroundTruth` matches what the current hand-coded numpy literals assert.
- A GVL-style reference-anchored VCF generates cleanly through `bcftools
  norm`/`consensus`.
- The classic Number×Type matrix and every variant class (SNP, MNP, insertion,
  deletion, complex/non-atomic, spanning deletion, multiallelic) are each covered
  by parametrized cases that round-trip through pysam with truth agreement.
- The benchmark gate passes: draw + serialize + derive-truth for a representative
  document stays within the target budget.
