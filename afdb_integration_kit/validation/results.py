from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Optional


class Level(str, Enum):
    ERROR = "ERROR"
    WARN = "WARN"
    INFO = "INFO"


@dataclass(frozen=True, slots=True)
class ValidationResult:
    check: str
    file: Path
    level: Level
    code: str
    message: str
    location: Optional[str] = None
    suggested_fix: Optional[str] = None
    metrics: Optional[Dict[str, float]] = None

    def to_dict(self) -> Dict[str, Any]:
        data: Dict[str, Any] = {
            "check": self.check,
            "file": str(self.file),
            "level": self.level.value,
            "code": self.code,
            "message": self.message,
        }
        if self.location is not None:
            data["location"] = self.location
        if self.suggested_fix is not None:
            data["suggested_fix"] = self.suggested_fix
        if self.metrics is not None:
            data["metrics"] = dict(self.metrics)
        return data


__all__ = ["Level", "ValidationResult"]
