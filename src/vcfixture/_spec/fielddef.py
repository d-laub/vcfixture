from __future__ import annotations

import re
from dataclasses import dataclass

from .number import Number, NumberKind
from .types import Type

_ID_RE = re.compile(r"^([A-Za-z_][0-9A-Za-z_.]*|1000G)$")
_KINDS = ("INFO", "FORMAT")


@dataclass(frozen=True)
class FieldDef:
    id: str
    number: Number
    type: Type
    description: str
    kind: str  # "INFO" | "FORMAT"

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise ValueError(f"kind must be one of {_KINDS}, got {self.kind!r}")
        if not _ID_RE.match(self.id):
            raise ValueError(f"ID {self.id!r} does not match VCF key regex")
        if self.type is Type.FLAG:
            if self.kind != "INFO":
                raise ValueError("Flag fields must be INFO, not FORMAT")
            if self.number.kind is not NumberKind.FLAG:
                raise ValueError("Flag fields must have Number=0")
        elif self.number.kind is NumberKind.FLAG:
            raise ValueError("Number=0 is only valid for Flag fields")

    def header_line(self) -> str:
        return (
            f"##{self.kind}=<ID={self.id},Number={self.number.header_str()},"
            f'Type={self.type.value},Description="{self.description}">'
        )
