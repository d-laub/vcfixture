# Per-version VCF generation (`version=`)

**Date:** 2026-06-02
**Status:** Approved design, ready for implementation plan

## Problem

`vcfixture` generates only VCF 4.5 documents. The `fileformat` field on
`VcfDocument` (and the `fileformat=` arg on `VcfBuilder`) is a free string that
is echoed verbatim into the `##fileformat` header and gates nothing — the entire
`_spec/` layer is modeled against 4.5.

Downstream consumers (`genoray`, `GenVarLoader`) must parse real-world VCFs
written against older spec versions. To test that, we need **faithful per-version
fixtures**: a document declared as `VCFv4.1` must not contain features that did
not exist in 4.1, and version-sensitive constructs must serialize the way a real
file of that version would.

## Scope of version differences

The four added reference specs (`docs/reference/VCFv4.{1,2,3,4}.tex`, alongside
the already-tracked `VCFv4.5.tex`) were consulted to determine what actually
changes across 4.1 → 4.5 within this library's modeling scope.

**There is exactly one truly backward-incompatible change.** Everything else is
*additive* — a feature that did not exist yet — which only "breaks" in one
direction (a newer-feature file fails an older parser; older files stay valid).

### The one breaking change: `SVLEN` (4.3 → 4.4)

| | 4.1 / 4.2 / 4.3 | 4.4 / 4.5 |
|---|---|---|
| `Number` | `.` | `A` |
| Description | "Difference in length between REF and ALT alleles" | "Length of structural variant" |
| Sign for a deletion | negative (`SVLEN=-205`) | positive (`SVLEN=205`) |

Same INFO key, opposite sign for deletions. VCF 4.4 itself adds a back-compat
clause ("the absolute value of SVLEN should be taken and a negative SVLEN should
be treated as positive"), which is the tell-tale of a breaking reinterpretation.

### Additive ladder (sourced from the specs)

- **Percent-encoding** of special characters: introduced 4.3 (absent in 4.1/4.2).
- **`Number=R`** terminology: 4.3+.
- **Spanning-deletion `*` allele**: 4.3+.
- **Breakend (BND) notation**: present since 4.1 (`MATEID`/`PARID` defined in 4.1).
- **`SVCLAIM`** reserved INFO: genuinely new in 4.4 (no mention in ≤4.3).
- **`<*>` reference block** (and the `LEN` reserved FORMAT field that describes it):
  introduced in 4.4.
- **Localized alleles** (LAA/LR/LAD…): 4.5-only — already out of scope per CLAUDE.md.

Note (corrected during planning against the real `.tex` source): `CILEN` and `PS`
are **not** 4.4-only. `CILEN` exists since 4.1 (`##INFO=<ID=CILEN,Number=2,…>`,
a different historical definition); `PS` is defined in 4.1 prose. They are
therefore `since V4_1` and not gated.

## What v1 gates (decided)

**SVLEN switch + reserved-field ladder.** The version arg gates:

1. The **SVLEN representation** (the `FieldDef` it resolves to — `Number`/description).
2. **Reserved-field availability** — which reserved INFO/FORMAT IDs exist at a
   given version.

Out of scope for v1 (the version arg does **not** enforce these): allele-form
gating (`*`, BND, symbolic-SV revisions), version-specific percent-encoding, and
localized alleles. A 4.1 document may still use a newer allele form; only
reserved-field availability and the SVLEN representation are version-gated.

## Design

### 1. `VcfVersion` enum — `_spec/version.py` (new), exported publicly

An **orderable** enum so ladder comparisons (`field.since <= version`) work:

```python
class VcfVersion(Enum):
    V4_1 = "VCFv4.1"
    V4_2 = "VCFv4.2"
    V4_3 = "VCFv4.3"
    V4_4 = "VCFv4.4"
    V4_5 = "VCFv4.5"   # latest / default
```

- Ordering implemented via declaration order (an internal index used by
  `__lt__`/`__le__`), so `VcfVersion.V4_3 < VcfVersion.V4_4` holds. `Enum` +
  explicit comparison methods, *not* `IntEnum` (the public `.value` must be the
  exact `##fileformat` string, not an int).
- `.value` is the exact `##fileformat` string.
- A module-level `LATEST = VcfVersion.V4_5` constant is the single source of the
  default.
- Re-exported from `__init__.py` alongside `Number`/`Type`. `_spec/` stays
  private; the enum is surfaced through the package root.

### 2. Model & serialization

- **`model.py`**: `VcfDocument.fileformat: str` → `version: VcfVersion`. The
  `__repr__` line updated to print `self.version.value`.
- **`serialize.py:95`**: `##fileformat={doc.version.value}`.
- **`build.py`**: `VcfBuilder.__init__` param `fileformat: str = "VCFv4.5"` →
  `version: VcfVersion = VcfVersion.V4_5`; stored as `self._version`; passed at
  the single `VcfDocument(...)` construction site (`build.py:256`). The old
  `fileformat=` param is **removed** (acceptable pre-1.0 break).

### 3. Reserved-field ladder + SVLEN switch — `_spec/reserved.py`

The chokepoint `reserved(id, kind)` is called from the builder and ~33 sites
(mostly tests). It gains an optional version param with a default so existing
callers are unaffected:

```python
def reserved(id: str, kind: str, version: VcfVersion = LATEST) -> FieldDef: ...
```

- Each reserved entry is annotated with the **version it was introduced**
  (`since`). Looking up a field whose `since > version` raises the same
  `KeyError` the builder already converts to a `ValueError` ("… is not a known
  reserved field …"), with the message extended to name the introducing version.
- **SVLEN** carries **version-specific variants**: resolves to
  `Number=.`, *"Difference in length between REF and ALT alleles"* for ≤ V4_3, and
  `Number=A`, *"Length of structural variant"* for ≥ V4_4. This single resolved
  `FieldDef` drives both the emitted `##INFO` header and the builder's eager
  count-validation (`Number=.` allows any count; `Number=A` requires one per ALT).
- **`since` table (verified against the tracked `.tex` source during planning):**
  - `SVCLAIM` (INFO) → `V4_4` (genuinely new in 4.4).
  - `LEN` (FORMAT) → `V4_4` (describes the `<*>` reference block, a 4.4 feature).
  - **All other** currently registered reserved fields → `V4_1` (each is defined
    in the 4.1 spec, whether as a `##INFO`/`##FORMAT` line or in prose).
  - **Known fidelity limitation (accepted for v1):** a few fields' *definitions*
    (not availability) shifted across versions — notably `CILEN`
    (`Number=2`/breakend description in ≤4.3 → `Number=.`/SVLEN description in
    4.4+). v1 keeps these at their latest definition with `since=V4_1`; only
    `SVLEN` gets a true version-variant definition. This mirrors the
    additive-ladder scope (we gate availability + the one SVLEN break, not every
    historical definition drift).
  - The reference `.tex` files are the authority for this table; tracked in
    `docs/reference/`, consulted at author time, never parsed at runtime.

No `truth.py` changes: `_allele_truth` already `abs()`-normalizes SVLEN when
computing the SV span/END (`truth.py:45`), and the verbatim per-record INFO
decode records whatever sign was serialized. Ground-truth is therefore already
version-independent for SVLEN.

The builder does **not** auto-flip a caller-supplied SVLEN sign — the caller owns
the sign on values they pass, and the verbatim INFO truth reflects it.

### 4. Strategies — `strategies.py`

- `documents(...)` and `reference_and_documents(...)` (and any other
  doc-producing strategy) gain `version: VcfVersion = VcfVersion.V4_5`. The
  version is a **fixed parameter per draw**, not drawn or filtered.
- **Construct, never reject** (per CLAUDE.md): a helper exposes the reserved IDs
  *available* at a given version (those with `since <= version`); strategies draw
  reserved fields only from that set. No `.filter()`/`assume()`.
- When a strategy generates a symbolic SV carrying `SVLEN`, it emits the
  **version-correct sign** (negative for spanning types like DEL at ≤ 4.3,
  positive at ≥ 4.4) so fixtures look like real files of that version.
- Exhaustive version coverage is driven by `@parametrize` over `VcfVersion` in
  tests; Hypothesis stays within a fixed version (matches the existing
  parametrize-for-matrix-axes convention).

## Testing

- The cross-parser round-trip property test (serialize → cyvcf2/pysam parse →
  assert against `GroundTruth`) runs **parametrized across all five versions**,
  keeping the independent-parser cross-check intact per version.
- Unit tests:
  - `##fileformat` header string matches the version.
  - `reserved()` rejects out-of-version IDs (e.g. `SVCLAIM` at V4_1) and accepts
    them at/after their `since`.
  - SVLEN `FieldDef` Number/description flips at the 4.3/4.4 boundary.
  - Builder count-validation honors the version-specific SVLEN `Number`
    (`Number=A` rejects a count ≠ n_alt at ≥ 4.4; `Number=.` allows any count at
    ≤ 4.3).
  - Existing `reserved(...)` call sites still pass with the defaulted version.
- `test_benchmark.py` runs against the default version; its budget is unchanged.

## Components affected

| File | Change |
|---|---|
| `_spec/version.py` | **new** — `VcfVersion` enum + `LATEST` |
| `_spec/reserved.py` | per-entry `since`; version param; SVLEN variant; available-IDs helper |
| `model.py` | `fileformat: str` → `version: VcfVersion`; repr |
| `serialize.py` | header line uses `doc.version.value` |
| `build.py` | `version=` param replaces `fileformat=`; thread into `reserved()` and `VcfDocument(...)` |
| `strategies.py` | `version=` param; draw reserved fields from available set; version-correct SVLEN sign |
| `__init__.py` | export `VcfVersion` |
| `docs/reference/VCFv4.{1,2,3,4}.tex` | tracked reference specs (LaTeX source) |
| tests | parametrize round-trip + new version unit tests |

## Out of scope (v1)

- Allele-form gating (`*`, BND, symbolic-SV revision differences).
- Version-specific percent-encoding (strict 4.3+, looser before).
- Localized alleles (LAA/LR/LAD) — 4.5-only, already excluded project-wide.
- Runtime parsing of the `.tex` specs (they are author-time reference only).
