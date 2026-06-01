from __future__ import annotations

import sys
from typing import Any

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing_extensions import override

__all__ = ["CompactRepr", "override"]


class CompactRepr:
    """Mixin that routes pretty-printers through ``__repr__``.

    Hypothesis's vendored pretty-printer (and IPython's) honor the
    ``_repr_pretty_`` protocol *before* field-walking a dataclass. Implementing
    it here makes our compact ``__repr__`` show up in Falsifying examples too.
    Subclasses must define their own ``__repr__`` in the class body (``@dataclass``
    would shadow an inherited one).
    """

    def _repr_pretty_(self, p: Any, cycle: bool) -> None:
        p.text(repr(self))
