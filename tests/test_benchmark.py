"""Performance gate: drawing + serializing + deriving truth for a
representative document must stay cheap."""

import time

from vcfixture.build import VcfBuilder


def _representative_doc():
    b = VcfBuilder(
        samples=[f"s{i}" for i in range(10)], contigs=[("chr1", 1_000_000)]
    ).fmt("GT")
    for k in range(50):
        gts = ["0|1" if (i + k) % 2 else "1/1" for i in range(10)]
        b.record("chr1", 1000 + k * 10, ref="A", alt=["T"], gt=gts)
    return b.build()


def test_build_serialize_truth_under_budget():
    t0 = time.perf_counter()
    for _ in range(100):
        doc = _representative_doc()
        doc.render()
        doc.truth()
    elapsed = time.perf_counter() - t0
    assert elapsed < 5.0, f"too slow: {elapsed:.2f}s for 100 iterations"
