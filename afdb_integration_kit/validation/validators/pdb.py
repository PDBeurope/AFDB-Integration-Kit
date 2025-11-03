from __future__ import annotations

from pathlib import Path
from typing import List

from ..context import ValidationContext
from ..registry import register_check
from ..results import ValidationResult


@register_check("pdb")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    return []


__all__ = ["run"]
