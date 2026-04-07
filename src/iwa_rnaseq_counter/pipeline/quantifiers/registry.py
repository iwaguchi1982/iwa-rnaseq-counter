from __future__ import annotations

from .base import BaseQuantifier
from .salmon import SalmonQuantifier
from .star import StarQuantifier
from .hisat2 import Hisat2Quantifier


_QUANTIFIERS: dict[str, type[BaseQuantifier]] = {
    "salmon": SalmonQuantifier,
    "star": StarQuantifier,
    "hisat2": Hisat2Quantifier,
    # v0.7.3+
    # "kallisto": KallistoQuantifier,
}


def get_quantifier(name: str) -> BaseQuantifier:
    key = str(name).strip().lower()
    try:
        return _QUANTIFIERS[key]()
    except KeyError as e:
        supported = ", ".join(sorted(_QUANTIFIERS))
        raise NotImplementedError(
            f"Unsupported quantifier: {name!r}. Supported: {supported}"
        ) from e


def list_quantifiers() -> list[str]:
    return sorted(_QUANTIFIERS.keys())
