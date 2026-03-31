from .salmon import SalmonQuantifier

_QUANTIFIERS = {
    "salmon": SalmonQuantifier,
}


def get_quantifier(name: str):
    try:
        return _QUANTIFIERS[name]()
    except KeyError as e:
        raise NotImplementedError(f"Unsupported quantifier: {name}") from e