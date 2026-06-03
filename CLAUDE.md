# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

vcfixture generates small VCF test data conforming to spec versions 4.1–4.5
(selectable via `version=`, default latest) — via an explicit builder or
Hypothesis strategies — and returns the **decoded ground truth** alongside, so
parser tests assert against a known oracle instead of hand-coded numpy literals.
It exists to replace the hand-maintained VCF fixtures + hand-derived expected
arrays in downstream repos `genoray` and `GenVarLoader`.

## Commands

`uv` is **required** (build backend + locked dev toolchain); bare `pip`/`venv` is
not supported. Run tools through `uv run` so they use the locked environment.

```bash
uv sync                       # create .venv, install project + dev deps
uv run prek install           # install pre-commit + commit-msg git hooks

uv run pytest                 # full test suite
uv run pytest tests/test_truth.py::test_name   # a single test
uv run pytest -q -k roundtrip # filter by name

uv run ruff check             # lint
uv run ruff format            # format
uv run pyrefly check          # strict type-check (src/ only)
uv run prek run --all-files   # all hooks against the whole tree
```

Commits must follow Conventional Commits (`feat:`, `fix:`, `docs:`, `test:`,
`chore:`) — the `commit-msg` hook (commitizen) rejects non-conforming messages.

**Do not bump the version or edit `CHANGELOG.md` by hand.** Versioning and the
changelog are handled in CI: the manual `release.yml` workflow
(`workflow_dispatch`) runs commitizen, which reads conventional commits since the
last tag to compute the next version, updates the changelog, tags, pushes to
`main`, and publishes to PyPI via OIDC. (Early-history bump commits that touched
`pyproject.toml`/`CHANGELOG.md` manually predate this — don't follow them.)

## Architecture

One immutable **document model** (`model.py`: `VcfDocument`/`Record`, plus
`Genotype`) is the hub. Two front-ends produce it; two back-ends consume it.
This guarantees the builder and the fuzzer cannot diverge, and makes ground
truth free.

```
 Builder kwargs ─┐                          ┌─→ serialize.py ─→ .vcf text ─→ io.py [bgzip + .csi/.tbi]
                 ├─→  VcfDocument (model) ───┤
 Hypothesis draws┘     (immutable)           └─→ truth.py ─→ GroundTruth (numpy + dicts)
```

- **`_spec/`** — the spec layer, covering VCF 4.1–4.5. `version.py` (`VcfVersion`:
  orderable enum `V4_1..V4_5` + `LATEST`), `number.py` (`Number`:
  `ONE | FIXED(n) | A | R | G | DOT | FLAG` + `cardinality()`), `types.py`
  (`Type`), `fielddef.py` (`FieldDef` + validity invariants), `genotype_order.py`
  (`Number=G` enumerator for placing PL/GL), `reserved.py` (curated reserved-field
  registry: AC/AF/AN/DP/AD/PL/GL/GT/...; `reserved(id, kind, version)` gates fields
  by introduction version and returns the version-correct `SVLEN` — its definition
  flips at the 4.3/4.4 boundary). Everything else builds on these.
- **`build.py`** (`VcfBuilder`) — explicit construction for named tests. Validates
  **eagerly**: Flag⇒Number=0 & INFO-only, undefined INFO/FORMAT ID used in a
  record, GT index out of range, value count ≠ resolved cardinality.
- **`strategies.py`** — Hypothesis `@composite` strategies feeding the *same*
  builders so shrinking works end-to-end (`documents`, `references`,
  `reference_and_documents`, `field_value`, `genotypes`). Plus
  `all_number_type_combos()`-style helpers for exhaustive `@parametrize` of the
  Number×Type matrix.
- **`serialize.py`** — model → spec-correct text (missing `.`, multi-dot arrays,
  GT first-phasing omission, trailing-FORMAT dropping except GT, percent-encoding).
- **`truth.py`** (`GroundTruth`) — model → numpy: genotype matrix with
  **`-1` = missing** (matches both consumers), bool phasing matrix, resolved
  per-record INFO and per-sample FORMAT.
- **`reference.py`** (`Reference`/`ReferenceBuilder`/`ReferenceSpec`/`RepeatFeature`)
  — optional reference-aware mode: given a FASTA + region, draws POS + matching
  REF + realistic ALTs so GVL's `bcftools norm`/`consensus` accepts the output.
  Reference-free by default (arbitrary REF — fine for genoray).
- **`io.py`** — text always; bgzip `.vcf.gz` + `.csi`/`.tbi` index via pysam on request.

The public API (`__init__.py`) re-exports `VcfBuilder`, `VcfVersion`, `Genotype`,
`Reference` family, `GroundTruth`, `Number`, `Type`, and the `strategies` module.
`_spec/` and `_typing.py` are private (`VcfVersion` is the one `_spec/` type
re-exported; `LATEST` stays internal).

## Conventions that matter here

- **Type-safety is a shipped guarantee.** The library has `py.typed` and `src/`
  must pass `pyrefly check` in **strict** mode. Tests are deliberately out of the
  type-check scope (they interrogate untyped cyvcf2/numpy output). New code in
  `src/` needs accurate annotations.
- **The test suite is the oracle's self-validation.** Tests round-trip every
  generated document through an *independent* parser (pysam/cyvcf2): serialize →
  parse → assert the third-party decode matches our `GroundTruth`. Keep that
  cross-check intact — it is what makes the truth trustworthy rather than merely
  self-consistent. There is a Hypothesis property test asserting this for
  arbitrary documents.
- **Valid-only generation.** Only spec-conformant VCFs (v1 scope). No
  invalid/malformed generation, no localized-allele `LA`/`LR`/`LG`, no CLI.
- **Strategy hot path: construct, never reject.** No `.filter()`/`assume()`
  rejection sampling — draw parameters then *compute* valid structures. Prefer
  flat, non-recursive strategies. Drive matrix coverage via `@parametrize`,
  reserving Hypothesis for values within a fixed combo. `test_benchmark.py` gates
  draw+serialize+derive-truth time against a budget.
- **Variant classes are first-class** (SNP, MNP, insertion, deletion,
  complex/delins, spanning deletion `*`, multiallelic) — each is an
  individually-addressable generation target carrying its own truth.

## CI

`.github/workflows/test.yml` runs lint+type-check, pytest across Python
3.10–3.13 (+ one macOS run), and a **clean-room wheel smoke test** that imports
the built wheel in a runtime-deps-only venv — this catches any `src/` module
reaching for a dev-only dependency. If you add an import in `src/`, ensure it is
declared in `[project].dependencies`, not just the dev group.

See `docs/superpowers/specs/2026-05-30-vcfixture-design.md` for full design
rationale and `CONTRIBUTING.md` for setup.
