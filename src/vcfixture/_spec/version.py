from __future__ import annotations

from enum import Enum


class VcfVersion(Enum):
    """A supported VCF spec version.

    ``.value`` is the exact ``##fileformat`` string. Members are orderable in
    chronological (declaration) order, so ``VcfVersion.V4_3 < VcfVersion.V4_4``.
    """

    V4_1 = "VCFv4.1"
    V4_2 = "VCFv4.2"
    V4_3 = "VCFv4.3"
    V4_4 = "VCFv4.4"
    V4_5 = "VCFv4.5"

    def __lt__(self, other: VcfVersion) -> bool:
        return _ORDER[self] < _ORDER[other]

    def __le__(self, other: VcfVersion) -> bool:
        return _ORDER[self] <= _ORDER[other]

    def __gt__(self, other: VcfVersion) -> bool:
        return _ORDER[self] > _ORDER[other]

    def __ge__(self, other: VcfVersion) -> bool:
        return _ORDER[self] >= _ORDER[other]


_ORDER: dict[VcfVersion, int] = {v: i for i, v in enumerate(VcfVersion)}

LATEST = VcfVersion.V4_5
