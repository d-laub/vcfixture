from __future__ import annotations
from pathlib import Path
import pysam
from .model import VcfDocument

def write_text(doc: VcfDocument, path) -> Path:
    p = Path(path)
    p.write_text(doc.render())
    return p

def write_bgzip(doc: VcfDocument, path, *, index: bool = True) -> Path:
    p = Path(path)
    if p.suffix != ".gz":
        p = p.with_suffix(p.suffix + ".gz")
    tmp = p.with_suffix("")  # strip .gz for the plain-text intermediate
    tmp.write_text(doc.render())
    pysam.tabix_compress(str(tmp), str(p), force=True)
    tmp.unlink()
    if index:
        pysam.tabix_index(str(p), preset="vcf", csi=True, force=True)
    return p
