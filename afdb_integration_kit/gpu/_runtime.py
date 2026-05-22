# Copyright 2026 Maciej Majewski
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# SPDX-License-Identifier: Apache-2.0

"""Runtime helpers for optional GPU analysis dependencies."""
from __future__ import annotations

import importlib


_INSTALL_HINT = "Install production dependencies with `uv pip install '.[production]'`."
_TORCH_CLUSTER_HINT = (
    "Install `torch_cluster` separately with a wheel that matches your "
    "PyTorch/CUDA build if you want the accelerator."
)


def import_optional_dependency(module_name: str, feature: str):
    """Import an optional dependency with an actionable error message."""
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        missing = exc.name or module_name
        hint = _TORCH_CLUSTER_HINT if missing == "torch_cluster" else _INSTALL_HINT
        raise ModuleNotFoundError(
            f"{feature} requires optional dependency `{missing}`. {hint}"
        ) from exc


def require_torch(feature: str):
    """Import torch with a production-analysis-specific error message."""
    return import_optional_dependency("torch", feature)


def resolve_device(
    device: str,
    *,
    torch_module=None,
    feature: str = "GPU clash/interface analysis",
) -> str:
    """Resolve ``cpu``/``cuda``/``auto`` into an executable torch device."""
    torch = torch_module if torch_module is not None else require_torch(feature)

    normalized = device.strip().lower()
    if normalized == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    if normalized == "cpu":
        return "cpu"
    if normalized.startswith("cuda"):
        if not torch.cuda.is_available():
            raise RuntimeError(
                f"device={device!r} was requested for {feature}, but CUDA is "
                "unavailable. Install a CUDA-enabled PyTorch build or use "
                "device='cpu' or device='auto'."
            )
        return normalized
    raise ValueError(
        f"Unsupported device {device!r}. Expected 'cpu', 'cuda', or 'auto'."
    )
