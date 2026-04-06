from __future__ import annotations

from .base import BaseQuantifier
from .salmon import SalmonQuantifier


_QUANTIFIERS: dict[str, type[BaseQuantifier]] = {
    "salmon": SalmonQuantifier,
    # v0.7.1+
    # "star": StarQuantifier,
    # "hisat2": Hisat2Quantifier,
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
