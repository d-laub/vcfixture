from __future__ import annotations

from enum import Enum


class Type(str, Enum):
    """VCF value type for an INFO or FORMAT field.

    Each member corresponds to the ``Type=`` token in a VCF
    meta-information line.  ``Flag`` is only valid in INFO fields; all
    other types are valid in both INFO and FORMAT.

    Use ``Type.info_allowed()`` or ``Type.format_allowed()`` to obtain
    the permitted subsets for header validation.
    """

    INTEGER = "Integer"
    """Whole-number values (``Type=Integer``)."""

    FLOAT = "Float"
    """Floating-point values (``Type=Float``)."""

    FLAG = "Flag"
    """Presence/absence flag (``Type=Flag``); INFO-only, no value written."""

    CHARACTER = "Character"
    """Single-character values (``Type=Character``)."""

    STRING = "String"
    """Arbitrary string values (``Type=String``)."""

    @classmethod
    def info_allowed(cls) -> tuple[Type, ...]:
        """Return all Types that are valid in an INFO field.

        Returns:
            All five ``Type`` members in definition order.
        """
        return tuple(cls)

    @classmethod
    def format_allowed(cls) -> tuple[Type, ...]:
        """Return all Types that are valid in a FORMAT field.

        ``Flag`` is excluded because the VCF spec does not permit
        ``Type=Flag`` in FORMAT meta-information lines.

        Returns:
            All ``Type`` members except ``Flag``.
        """
        return tuple(t for t in cls if t is not cls.FLAG)
