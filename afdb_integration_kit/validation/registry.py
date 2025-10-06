from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional


ValidationCallable = Callable[..., tuple[bool, Dict[str, Any]]]
FormatCallable = Callable[..., str]


@dataclass(frozen=True)
class ValidationHook:
    name: str
    run: ValidationCallable
    formatter: Optional[FormatCallable] = None
    default_kwargs: Dict[str, Any] = field(default_factory=dict)
    description: Optional[str] = None

    def build_kwargs(self, overrides: Optional[Mapping[str, Any]] = None) -> Dict[str, Any]:
        if not overrides:
            return dict(self.default_kwargs)
        merged: Dict[str, Any] = dict(self.default_kwargs)
        merged.update(overrides)
        return merged


_REGISTRY: Dict[str, ValidationHook] = {}


def register_validator(hook: ValidationHook) -> None:
    if hook.name in _REGISTRY:
        raise ValueError(f"Validator '{hook.name}' is already registered")
    _REGISTRY[hook.name] = hook


def get_validator(name: str) -> ValidationHook:
    try:
        return _REGISTRY[name]
    except KeyError as exc:
        available = ", ".join(sorted(_REGISTRY)) or "<none>"
        raise KeyError(f"Unknown validator '{name}'. Available: {available}") from exc


def iter_validators() -> Iterable[ValidationHook]:
    return _REGISTRY.values()


def list_validators() -> List[ValidationHook]:
    return list(_REGISTRY.values())


def list_validator_names() -> List[str]:
    return list(_REGISTRY.keys())


__all__ = [
    "FormatCallable",
    "ValidationCallable",
    "ValidationHook",
    "get_validator",
    "iter_validators",
    "list_validator_names",
    "list_validators",
    "register_validator",
]
