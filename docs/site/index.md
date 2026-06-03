# vcfixture

Generate small VCF (v4.5) test data — via an explicit builder or Hypothesis
strategies — with the **decoded ground truth** returned alongside, so parser
tests assert against a known oracle instead of hand-coded literals.

## Install

```bash
uv add --dev vcfixture   # or: pip install vcfixture
```

## 30-second example

```python
from vcfixture import VcfBuilder, Number, Type, Seq

doc = (VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 100000)])
       .info("AF", Number.A, Type.FLOAT)
       .fmt("GT").fmt("DS", Number.A, Type.FLOAT)
       .record("chr1", 81262, ref="GAT", alt=[Seq("A")],
               gt=["0|1", "1|1"], info={"AF": [0.5]}, DS=[[1.0], [2.0]]))

text  = doc.render()                       # VCF text
truth = doc.truth()                        # GroundTruth (numpy genotypes, ...)
doc.write("x.vcf.gz", bgzip=True, index=True)
```

`Seq(...)` builds a sequence ALT allele; symbolic SV and breakend allele types are also available — see the [API reference](api/index.md).

## Where to go next

- [VcfBuilder guide](guide/vcfbuilder.md) — explicit construction for named tests.
- [Property-based testing](guide/property-testing.md) — Hypothesis strategies.
- [API reference](api/index.md) — the full public surface.
