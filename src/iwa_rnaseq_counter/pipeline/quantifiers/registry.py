from .salmon import SalmonQuantifier
from .base import BaseQuantifier

_QUANTIFIERS = {
    "salmon": SalmonQuantifier,
}


def get_quantifier(name: str) -> BaseQuantifier:
    try:
        return _QUANTIFIERS[name]()
    except KeyError as e:
        raise NotImplementedError(f"Unsupported quantifier: {name}") from e