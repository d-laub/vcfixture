# vcforge

Generate small VCF (v4.5) test data — via an explicit builder or Hypothesis
strategies — with the decoded ground truth returned alongside, so parser tests
assert against a known oracle instead of hand-coded literals.

```python
from vcforge import VcfBuilder, Number, Type

doc = (VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 100000)])
       .info("AF", Number.A, Type.FLOAT)
       .fmt("GT").fmt("DS", Number.A, Type.FLOAT)
       .record("chr1", 81262, ref="GAT", alt=["A"],
               gt=["0|1", "1|1"], info={"AF": [0.5]}, DS=[[1.0], [2.0]]))

text  = doc.render()                       # VCF text
truth = doc.truth()                        # GroundTruth (numpy genotypes, ...)
doc.write("x.vcf.gz", bgzip=True, index=True)
```

See `docs/superpowers/specs/2026-05-30-vcforge-design.md` for the design, and
[CONTRIBUTING.md](CONTRIBUTING.md) for development setup (uv is required).
