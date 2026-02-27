# -*- coding: utf-8 -*-

from typing import Dict, List, Optional, Type

from .base import AnalysisMethod
from .local_interface_strict import LocalInterfaceStrictMethod


_METHOD_REGISTRY: Dict[str, Type[AnalysisMethod]] = {
    "local_interface_strict": LocalInterfaceStrictMethod,
}


def list_available_methods() -> List[str]:
    return sorted(_METHOD_REGISTRY.keys())


def build_methods(method_ids: Optional[List[str]] = None) -> List[AnalysisMethod]:
    """
    Build method instances from ids.

    If method_ids is None, all registered methods are enabled.
    """
    if method_ids is None:
        method_ids = list_available_methods()

    methods: List[AnalysisMethod] = []
    for method_id in method_ids:
        method_cls = _METHOD_REGISTRY.get(method_id)
        if method_cls is None:
            print(f"[警告] 未注册的方法：{method_id}，已跳过。")
            continue
        methods.append(method_cls())
    return methods
