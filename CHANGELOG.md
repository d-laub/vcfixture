## v0.6.0 (2026-06-03)

### Feat

- add vcfixture skill for consuming-lib maintainers
- version= on document strategies with version-correct SVLEN/SVCLAIM
- VcfBuilder version= param gates reserved fields and SVLEN cardinality
- store VcfVersion on VcfDocument; derive fileformat header from it
- gate reserved fields by version; SVLEN definition flips at 4.4
- add orderable VcfVersion enum
- export Allele construction vocabulary; benchmark symbolic path
- Hypothesis strategy for symbolic SV documents + cyvcf2 oracle
- emit ##ALT header lines for symbolic alleles
- eager builder validation for symbolic/breakend alleles
- per-allele AlleleTruth with is_sequence flag and SV geometry
- add reserved structural-variant INFO/FORMAT fields
- typed Allele union with classifier and smart constructors

### Refactor

- default VcfBuilder/strategies version to LATEST; document 4.1-4.5 support
- derive symbolic SV type list from canonical set; pin SVCLAIM keys
- store typed Allele objects in Record.alts

## v0.5.0 (2026-06-01)

### Feat

- compact RepeatFeature/ReferenceSpec reprs
- compact ContigDef/Record/VcfDocument reprs
- compact Genotype repr
- compact FieldDef repr
- add CompactRepr mixin and compact Number repr

### Fix

- stop Number ClassVar singletons leaking into dataclass fields

## v0.4.1 (2026-06-01)

### Fix

- __version__

## v0.4.0 (2026-06-01)

### Feat

- export ReferenceBuilder/ReferenceSpec/RepeatFeature; bump 0.3.0
- reference_and_documents() paired strategy
- reference-consistent documents() with violation labels
- references() strategy drawing ReferenceSpec with planted repeats
- ReferenceSpec.write bgzipped+faidx FASTA
- ReferenceBuilder/ReferenceSpec with tandem-repeat provenance
- general per-variant labels carried into GroundTruth

### Refactor

- share draw_ref_alt between Reference and ReferenceSpec

## v0.2.1 (2026-05-31)

### Fix

- declare hypothesis as a runtime dependency

## v0.2.0 (2026-05-31)

### Feat

- generate multiallelic records in documents() strategy
- percent-encode reserved chars in string values
- add field-value generation strategy for Number x Type matrix
- public API exports and README
- add Hypothesis strategies and coverage tables
- add reference-aware REF/ALT adapter
- add text + bgzip/index IO
- add VcfBuilder with eager validation
- add ground-truth deriver
- add VCF serializer
- add variant-class constructors and classify
- add Record, ContigDef, VcfDocument model
- add Genotype parse/render
- add curated reserved-field registry
- add FieldDef with validity invariants
- add Number=G genotype ordering
- add Number model and cardinality
- add VCF Type enum

### Fix

- serialize non-finite floats as missing; clearer error for unknown reserved field
