from __future__ import annotations

from enum import Enum


class Type(str, Enum):
    INTEGER = "Integer"
    FLOAT = "Float"
    FLAG = "Flag"
    CHARACTER = "Character"
    STRING = "String"

    @classmethod
    def info_allowed(cls) -> tuple[Type, ...]:
        return tuple(cls)

    @classmethod
    def format_allowed(cls) -> tuple[Type, ...]:
        return tuple(t for t in cls if t is not cls.FLAG)
