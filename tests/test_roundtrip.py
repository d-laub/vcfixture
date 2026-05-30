"""Self-validation: serialize -> parse with an INDEPENDENT parser -> the
third-party decode must match our derived GroundTruth."""
from pathlib import Path
import tempfile
import numpy as np
import pytest
from hypothesis import given, settings, HealthCheck
from vcforge import strategies as S

cyvcf2 = pytest.importorskip("cyvcf2")

def _genos_from_cyvcf2(variant, n_samples, ploidy):
    out = np.full((n_samples, ploidy), -1, dtype=np.int32)
    for si, g in enumerate(variant.genotypes):
        for ai, allele in enumerate(g[:-1]):
            out[si, ai] = allele
    return out

@settings(max_examples=75, deadline=None,
          suppress_health_check=[HealthCheck.too_slow])
@given(S.documents())
def test_genotypes_round_trip_through_cyvcf2(doc):
    truth = doc.truth()
    d = tempfile.mkdtemp()
    path = doc.write(Path(d) / "x.vcf.gz", bgzip=True, index=True)

    vf = cyvcf2.VCF(str(path))
    n_samples = len(doc.samples)
    ploidy = truth.genotypes.shape[2]
    for ri, variant in enumerate(vf):
        got = _genos_from_cyvcf2(variant, n_samples, ploidy)
        np.testing.assert_array_equal(got, truth.genotypes[ri],
            err_msg=f"genotype mismatch at record {ri}")
        assert variant.POS == int(truth.pos[ri])
        assert variant.REF == truth.ref[ri]
