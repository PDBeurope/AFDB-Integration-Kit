from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Set


@dataclass(slots=True)
class ValidationContext:
    root: Path
    config: Dict[str, Any]
    selected_checks: Optional[Set[str]]
    logger: logging.Logger


__all__ = ["ValidationContext"]
