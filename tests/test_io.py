from pathlib import Path
import pysam
from vcforge.build import VcfBuilder

def _doc():
    return (VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 100000)])
            .fmt("GT")
            .record("chr1", 81262, ref="A", alt=["T"], gt=["0|1", "1|1"])
            .build())

def test_write_text(tmp_path: Path):
    p = _doc().write(tmp_path / "x.vcf")
    assert Path(p).read_text().startswith("##fileformat=VCFv4.5")

def test_write_bgzip_and_index_is_queryable(tmp_path: Path):
    p = _doc().write(tmp_path / "x.vcf.gz", bgzip=True, index=True)
    assert Path(str(p) + ".csi").exists()
    vf = pysam.VariantFile(str(p))
    rows = list(vf.fetch("chr1", 81000, 82000))
    assert len(rows) == 1
    assert rows[0].pos == 81262
