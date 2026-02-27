# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Protocol


@dataclass
class AnalysisContext:
    """Single-folder runtime context passed to each analysis method."""

    folder_name: str
    folder_path: str
    species_record: Optional[Dict[str, Any]]
    complex_object: str
    pymol_cmd: Any
    chain_map: Dict[str, str]


@dataclass
class MethodExecutionResult:
    """Standard method output, serialized into method_results."""

    success: bool
    message: str
    metrics: Dict[str, Any] = field(default_factory=dict)
    artifacts: Dict[str, Any] = field(default_factory=dict)
    legacy_interface_analysis: Optional[Dict[str, Any]] = None
    legacy_local_params: Optional[Dict[str, Any]] = None

    def to_dict(self, method_id: str, display_name: str, version: str) -> Dict[str, Any]:
        return {
            "method_id": method_id,
            "display_name": display_name,
            "version": version,
            "success": self.success,
            "message": self.message,
            "metrics": self.metrics,
            "artifacts": self.artifacts,
        }


class AnalysisMethod(Protocol):
    """Plug-in contract for a new analysis method."""

    method_id: str
    display_name: str
    version: str

    def run(self, context: AnalysisContext) -> MethodExecutionResult:
        ...
