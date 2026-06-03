from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from .genotype import Genotype
from .model import Record, VcfDocument

# Order matters: '%' must be replaced first.
_PERCENT = [
    ("%", "%25"),
    (":", "%3A"),
    (";", "%3B"),
    ("=", "%3D"),
    (",", "%2C"),
    ("\r", "%0D"),
    ("\n", "%0A"),
    ("\t", "%09"),
]


def _encode(s: str) -> str:
    for ch, rep in _PERCENT:
        s = s.replace(ch, rep)
    return s


def _fmt_scalar(v: Any) -> str:
    if v is None:
        return "."
    if isinstance(v, bool):
        return str(int(v))
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return "."
        return repr(v)
    if isinstance(v, str):
        return _encode(v)
    return str(v)


def _fmt_value(v: Any) -> str:
    if isinstance(v, (list, tuple)):
        if len(v) == 0:
            return "."
        return ",".join(_fmt_scalar(x) for x in v)
    return _fmt_scalar(v)


def _render_info(rec: Record) -> str:
    if not rec.info:
        return "."
    parts: list[str] = []
    for key, val in rec.info.items():
        if val is True:
            parts.append(key)
        elif val is False:
            continue
        else:
            parts.append(f"{key}={_fmt_value(val)}")
    return ";".join(parts) if parts else "."


def _render_sample(rec: Record, sample: Mapping[str, Any]) -> str:
    vals: list[str] = []
    for key in rec.fmt_keys:
        v = sample.get(key)
        if isinstance(v, Genotype):
            vals.append(v.render())
        else:
            vals.append(_fmt_value(v))
    return ":".join(vals)


def _render_record(rec: Record) -> str:
    ids = ";".join(rec.ids) if rec.ids else "."
    alt = ",".join(a.render() for a in rec.alts) if rec.alts else "."
    qual = "." if rec.qual is None else _fmt_scalar(rec.qual)
    if rec.filters is None:
        filt = "."
    elif len(rec.filters) == 0:
        filt = "PASS"
    else:
        filt = ";".join(rec.filters)
    cols = [rec.chrom, str(rec.pos), ids, rec.ref, alt, qual, filt, _render_info(rec)]
    if rec.fmt_keys:
        cols.append(":".join(rec.fmt_keys))
        for s in rec.samples:
            cols.append(_render_sample(rec, s))
    return "\t".join(cols)


def render_document(doc: VcfDocument) -> str:
    lines = [f"##fileformat={doc.version.value}"]
    for c in doc.contigs:
        lines.append(c.header_line())
    for f in doc.info_defs:
        lines.append(f.header_line())
    for fid, desc in doc.filter_defs:
        lines.append(f'##FILTER=<ID={fid},Description="{desc}">')
    for f in doc.format_defs:
        lines.append(f.header_line())
    for ad in doc.alt_defs:
        lines.append(ad.header_line())
    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if doc.format_defs or any(r.fmt_keys for r in doc.records):
        header.append("FORMAT")
        header.extend(doc.samples)
    lines.append("\t".join(header))
    for rec in doc.records:
        lines.append(_render_record(rec))
    return "\n".join(lines) + "\n"
