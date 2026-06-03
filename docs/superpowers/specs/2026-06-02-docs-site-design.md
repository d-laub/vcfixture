# vcfixture documentation site — design

**Date:** 2026-06-02
**Status:** approved (brainstorming) — pending implementation plan
**Branch:** `docs/site`

## Goal

Publish a documentation site for vcfixture, hosted on GitHub Pages, built with
[Zensical](https://zensical.org/) (the Material for MkDocs team's static site
generator). The site must cover **both** front-ends of the library:

1. **`VcfBuilder`** — explicit, named-test construction.
2. **Property-based testing** — the Hypothesis `strategies` module and the
   testing philosophy around it.

It also publishes an auto-generated **API reference** for the public surface.

## Non-goals

- No publishing of internal design specs/plans (`docs/superpowers`, `docs/plans`,
  `docs/reference`) — those stay private.
- No versioned-docs (mike) setup in v1.
- No CHANGELOG/version automation changes (handled by `release.yml`).

## Decisions (from brainstorming)

| Question | Decision |
|----------|----------|
| Content style | Hand-written guides **plus** auto-generated API reference. |
| API reference source | Write/expand docstrings across the public API **now**, then auto-render. |
| Site source location | Dedicated `docs/site/` source dir; internal specs stay private. |
| Tooling | Zensical + `mkdocstrings-python`, driven through `uv`. |
| Hosting | GitHub Pages via a new `.github/workflows/docs.yml`. |

## Architecture & tooling

- **Generator:** Zensical, configured via `zensical.toml` at the repo root.
- **API docs:** `mkdocstrings` with the Python handler
  (`mkdocstrings-python`), `docstring_style = "google"`, `paths = ["src"]`.
  Griffe analyzes `src/` statically — rendering does not require importing the
  package, but the package is installed in the docs env anyway for accurate
  annotation resolution.
- **Dependency management:** docs deps live in a **`docs` dependency group** in
  `pyproject.toml` (`[dependency-groups]`), never in `[project].dependencies`.
  This keeps the runtime wheel clean (the clean-room wheel smoke test in
  `test.yml` must stay green). All doc commands run through `uv run`.

## Repository layout

```
zensical.toml                 # repo root; docs_dir = "docs/site", site_dir = "site"
.github/workflows/docs.yml    # build (PRs) + deploy (push to main)
docs/
  site/                       # docs_dir — the published site source
    index.md                  # what vcfixture is, install, 30-second example
    guide/
      vcfbuilder.md           # explicit construction guide
      property-testing.md     # Hypothesis strategies + testing philosophy
    api/
      index.md                # mkdocstrings render of the public API
    stylesheets/              # (only if custom CSS is needed)
  superpowers/                # internal — NOT published (outside docs_dir)
  plans/                      # internal — NOT published
  reference/                  # internal (VCF 4.5 .tex) — NOT published
site/                         # build output — gitignored
```

Because `docs_dir = "docs/site"`, the internal trees are simply outside the
build root. Nothing about them needs an exclude rule — the build is hermetic by
construction and cannot leak internal design docs.

## `zensical.toml` (shape)

```toml
[project]
site_name = "vcfixture"
site_description = "Generate small VCF test data with decoded ground truth."
site_url = "https://d-laub.github.io/vcfixture/"
repo_url = "https://github.com/d-laub/vcfixture"
repo_name = "d-laub/vcfixture"
edit_uri = "edit/main/docs/site/"
docs_dir = "docs/site"
site_dir = "site"
nav = [
  { "Home" = "index.md" },
  { "Guide" = [
      "guide/vcfbuilder.md",
      "guide/property-testing.md",
  ]},
  { "API Reference" = "api/index.md" },
]

[project.theme]
variant = "classic"            # Material for MkDocs look
features = [
  "navigation.sections",
  "navigation.top",
  "navigation.footer",
  "content.action.view",
  "search.highlight",
]

[[project.theme.palette]]
scheme = "default"
toggle.icon = "lucide/sun"
toggle.name = "Switch to dark mode"

[[project.theme.palette]]
scheme = "slate"
toggle.icon = "lucide/moon"
toggle.name = "Switch to light mode"

[project.plugins.mkdocstrings.handlers.python]
paths = ["src"]
inventories = ["https://docs.python.org/3/objects.inv"]

[project.plugins.mkdocstrings.handlers.python.options]
docstring_style = "google"
show_source = false
show_root_heading = true
members_order = "source"
```

(Exact feature/option set may be tuned during implementation; this captures
intent.)

## Content plan

### `index.md` (Home)
- One-paragraph "what is this / why" (lifted/adapted from README + CLAUDE.md).
- Install (`uv add --dev vcfixture` / `pip install vcfixture`).
- The 30-second example from the README (builder → `render()` / `truth()` /
  `write()`).
- Links into the two guides and the API reference.

### `guide/vcfbuilder.md`
- The mental model: one immutable `VcfDocument` is the hub; the builder is one
  front-end that produces it.
- Walkthrough: construct a builder (samples/contigs), declare INFO/FORMAT/FILTER
  fields, add records, then `render()` / `truth()` / `write()`.
- Eager validation guarantees (Flag ⇒ Number=0 & INFO-only, undefined IDs, GT
  index range, count vs cardinality).
- Variant classes as first-class targets (SNP, MNP, ins/del, complex/delins,
  spanning deletion `*`, multiallelic) and symbolic/breakend alleles via the
  `Allele` family.
- What `GroundTruth` gives you (numpy genotype matrix with `-1` = missing,
  phasing matrix, resolved INFO/FORMAT) and why it is trustworthy (round-tripped
  through an independent parser in the test suite).

### `guide/property-testing.md`
- Why property-based testing here: builder and fuzzer share the same model, so
  they cannot diverge; ground truth is free.
- The `strategies` module surface: `documents`, `references`,
  `reference_and_documents`, `field_value`, `genotypes`, and the
  `all_number_type_combos()`-style matrix helpers.
- The hot-path rule: **construct, never reject** — no `.filter()`/`assume()`;
  draw parameters then compute valid structures. Prefer flat, non-recursive
  strategies.
- Driving the Number×Type matrix via `@parametrize`, reserving Hypothesis for
  values within a fixed combo.
- The oracle/round-trip pattern: serialize → parse with cyvcf2/pysam → assert
  the third-party decode matches `GroundTruth`. Include the existing
  property-test example.

### `api/index.md`
- mkdocstrings render of the public API re-exported from `vcfixture/__init__.py`:
  `VcfBuilder`, `Genotype`, the `Allele` family (`Allele`, `Seq`/`SequenceAllele`,
  `Sym`/`SymbolicAllele`, `Bnd`/`BreakendAllele`, `Star`/`SpanningDeletion`,
  `Unspecified`/`UnspecifiedAllele`), `Reference`/`ReferenceBuilder`/
  `ReferenceSpec`/`RepeatFeature`, `GroundTruth`, `Number`, `Type`, and the
  `strategies` module.
- Grouped logically (Builder · Alleles · Reference · Truth · Spec types ·
  Strategies) rather than one flat dump.

## Docstring work

The public symbols are currently thinly documented (`VcfBuilder`/`__init__` have
no docstrings; only `strategies.py` uses a `Args:` block). To make the API
reference useful, add/expand **Google-style** docstrings across the public API:

- `VcfBuilder` (class + `__init__`, `info`, `fmt`, `filter`, `alt`, `record`,
  `render`, `truth`, `write`, and any other public methods).
- The `Allele` family classes and their construction helpers.
- `strategies` public entry points (normalize the existing ones to a consistent
  Google style).
- `GroundTruth`, `Number`, `Type`, and the `Reference` family.

Constraints:
- **Docstrings only — no behavior changes.**
- `src/` must continue to pass `pyrefly check` (strict) and `ruff check`/`format`.
- Use examples in docstrings where they aid the reference; keep them accurate.

## CI & deployment

New `.github/workflows/docs.yml`:

- **On pull request:** `uv sync --group docs` + `uv run zensical build --clean`
  as a **build-only check** (no deploy) so broken docs fail CI.
- **On push to `main`:** build, then `upload-pages-artifact` →
  `deploy-pages`, publishing to `https://d-laub.github.io/vcfixture/`.
- Permissions: `contents: read`, `pages: write`, `id-token: write`;
  `github-pages` environment.
- Uses `astral-sh/setup-uv` + `uv run zensical build` to honor the
  uv-required toolchain (not the skill's generic `pip install zensical`).

`.gitignore`: add `site/` (build output).

GitHub Pages must be set to "GitHub Actions" as the source (one-time repo
setting — flagged as a manual follow-up for the maintainer).

## Verification

- `uv run zensical build --clean` succeeds with **no warnings**.
- `uv run zensical serve` renders Home, both guides, and the API reference;
  API pages show real descriptions (not bare signatures).
- `uv run pyrefly check` and `uv run ruff check` still pass after docstring
  additions; `uv run pytest` stays green (docstrings shouldn't affect tests, but
  doctest-like examples must be correct if any are executed).
- Clean-room wheel smoke test logic unaffected (docs deps are in the `docs`
  group, not runtime deps).

## Open follow-ups (post-v1, out of scope)

- Versioned docs (mike) if/when releases warrant it.
- "Was this page helpful?" feedback / analytics.
- Expanding examples into a cookbook of variant-class recipes.
