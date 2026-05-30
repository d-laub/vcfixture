"""Round-trip the full Number x Type matrix + multiallelic records through the
independent cyvcf2 parser and assert decoded values match the derived truth."""

import os
import tempfile

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings

from vcfixture import strategies as S
from vcfixture._spec.number import NumberKind
from vcfixture._spec.types import Type

cyvcf2 = pytest.importorskip("cyvcf2")

# cyvcf2 pads variable-length (Number=.) FORMAT arrays across samples with
# these sentinel values:  INT32_MIN+1 for integer, NaN for float.
_INT_MISSING = np.int32(-2147483647)
_INT_END_OF_VECTOR = np.int32(-2147483646)


def _info_defs_by_id(doc):
    return {fd.id: fd for fd in doc.info_defs}


def _fmt_defs_by_id(doc):
    return {fd.id: fd for fd in doc.format_defs if fd.id != "GT"}


def _norm_info_expected(typ, val):
    if typ is Type.FLAG:
        return val
    return list(val)


def _norm_info_got(typ, got):
    if typ is Type.FLAG:
        return bool(got) is True
    if got is None:
        return None
    if isinstance(got, tuple):
        return list(got)
    return [got]


def _strip_fmt_padding_int(arr):
    """Strip trailing INT_MISSING / INT_END_OF_VECTOR sentinels cyvcf2 appends
    when padding variable-length FORMAT fields across samples."""
    vals = list(arr)
    while vals and vals[-1] in (_INT_MISSING, _INT_END_OF_VECTOR):
        vals.pop()
    return vals


def _strip_fmt_padding_float(arr):
    """Strip trailing NaN sentinels from variable-length FORMAT float fields."""
    vals = list(arr)
    while vals and np.isnan(vals[-1]):
        vals.pop()
    return vals


@settings(
    max_examples=60,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
@given(S.documents_with_fields())
def test_matrix_round_trips_through_cyvcf2(doc):
    truth = doc.truth()
    info_defs = _info_defs_by_id(doc)
    fmt_defs = _fmt_defs_by_id(doc)
    n_samples = len(doc.samples)

    d = tempfile.mkdtemp()
    path = doc.write(os.path.join(d, "m.vcf.gz"), bgzip=True, index=True)
    vf = cyvcf2.VCF(str(path))

    for ri, variant in enumerate(vf):
        for fid, fd in info_defs.items():
            exp = _norm_info_expected(fd.type, truth.info[ri][fid])
            got = _norm_info_got(fd.type, variant.INFO.get(fid))
            if fd.type is Type.FLAG:
                assert got == exp, f"INFO {fid} flag r{ri}: {got} != {exp}"
            elif fd.type is Type.FLOAT:
                np.testing.assert_array_equal(
                    np.float32(got), np.float32(exp), err_msg=f"INFO {fid} r{ri}"
                )
            elif fd.type is Type.INTEGER:
                assert [int(x) for x in got] == [int(x) for x in exp], (
                    f"INFO {fid} r{ri}: {got} != {exp}"
                )
            else:
                # Character/String INFO: cyvcf2 returns a single comma-joined
                # string (or a bare scalar for Number=1).
                got_join = got[0] if len(got) == 1 else ",".join(map(str, got))
                exp_join = ",".join(map(str, exp))
                assert got_join == exp_join, (
                    f"INFO {fid} r{ri}: {got_join!r} != {exp_join!r}"
                )

        for fid, fd in fmt_defs.items():
            got = variant.format(fid)
            is_dot = fd.number.kind is NumberKind.DOT
            for si in range(n_samples):
                exp = list(truth.format[ri][si][fid])
                if fd.type is Type.FLOAT:
                    # For Number=. cyvcf2 pads with NaN to the longest sample
                    got_vals = (
                        _strip_fmt_padding_float(got[si]) if is_dot else list(got[si])
                    )
                    np.testing.assert_array_equal(
                        np.float32(got_vals),
                        np.float32(exp),
                        err_msg=f"FMT {fid} r{ri} s{si}",
                    )
                elif fd.type is Type.INTEGER:
                    # For Number=. cyvcf2 pads with INT_MISSING sentinel
                    got_vals = (
                        _strip_fmt_padding_int(got[si])
                        if is_dot
                        else [int(x) for x in got[si]]
                    )
                    assert [int(x) for x in got_vals] == [int(x) for x in exp], (
                        f"FMT {fid} r{ri} s{si}: {got_vals} != {exp}"
                    )
                else:
                    # Character/String FORMAT: cyvcf2 returns one per-sample
                    # string with values comma-joined.
                    exp_join = ",".join(map(str, exp))
                    got_s = got[si]
                    got_s = got_s if isinstance(got_s, str) else str(got_s)
                    assert got_s == exp_join, (
                        f"FMT {fid} r{ri} s{si}: {got_s!r} != {exp_join!r}"
                    )
