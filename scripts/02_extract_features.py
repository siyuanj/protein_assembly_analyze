# -*- coding: utf-8 -*-
"""
可扩展批处理特征提取脚本

目标：
1. 统一处理目录扫描、最佳模型选择、结果增量保存。
2. 通过 analysis_methods 插件机制执行分析方法，便于后续扩展新方法。
3. 保留兼容字段（interface_analysis / local_params）以兼容旧可视化脚本。
"""

import json
import os
from typing import Any, Dict, List, Optional, Tuple

try:
    import pymol
    from pymol import cmd
except Exception:
    pymol = None
    cmd = None

from analysis_methods.base import AnalysisContext, MethodExecutionResult
from analysis_methods.registry import build_methods, list_available_methods


# ======================================================================
# Config
# ======================================================================
CHAIN_MAP = {"A": "APH1", "B": "NCT", "C": "PEN2", "D": "PS1"}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
BASE_DIR = os.path.join(PROJECT_ROOT, "data", "raw_structures", "af_samllsample_1021")
OUTPUT_JSON = os.path.join(PROJECT_ROOT, "results", "parsed_data", "batch_analysis_results.json")
ALL_DATA_JSON = os.path.join(PROJECT_ROOT, "data", "metadata", "all_data_1101.json")
EXCLUDE_DIRS = {".git", ".vscode", "__pycache__", "venv", "node_modules", ".idea"}

# 方法选择（可选）
# - 默认（空）= 启用全部注册方法
# - 例如：ANALYSIS_METHODS=local_interface_strict
METHODS_ENV = os.getenv("ANALYSIS_METHODS", "").strip()

# 是否强制重跑已存在结果
FORCE_REANALYZE = os.getenv("FORCE_REANALYZE", "0").strip().lower() in {"1", "true", "yes"}


# ======================================================================
# JSON helpers
# ======================================================================
def _load_existing_results(path: str) -> Dict[str, Any]:
    if not os.path.isfile(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _save_results_atomic(data: Dict[str, Any], path: str) -> None:
    tmp = f"{path}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)
    os.replace(tmp, path)


# ======================================================================
# Metadata helpers
# ======================================================================
def _load_all_data_records(path: str) -> List[Dict[str, Any]]:
    if not os.path.isfile(path):
        print(f"[警告] 未找到 metadata 文件：{path}")
        return []
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception as e:
        print(f"[警告] 读取 metadata 失败：{e}")
        return []
    if isinstance(data, list):
        return [x for x in data if isinstance(x, dict)]
    if isinstance(data, dict):
        return [data]
    return []


def _build_species_record_index(records: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    index: Dict[str, Dict[str, Any]] = {}
    for rec in records:
        for key in ("name_final", "name_cleaned", "name_original"):
            name = rec.get(key)
            if isinstance(name, str) and name and name not in index:
                index[name] = rec
    return index


def _resolve_species_record(record_index: Dict[str, Dict[str, Any]], folder_name: str) -> Optional[Dict[str, Any]]:
    return record_index.get(folder_name)


# ======================================================================
# Model selection
# ======================================================================
def find_best_model_and_cif(folder_path: str) -> Tuple[Optional[Dict[str, Any]], Optional[str], float, str]:
    """
    在目录中扫描 summary_confidences_{0..4}.json，选择 iptm 最高的模型，
    并返回对应 cif 路径（兼容 fold_* 前缀）。
    """
    best_iptm = -1.0
    best_summary = None
    best_idx = -1
    best_prefix = ""

    print(f"   - 1. 正在扫描最佳模型：{folder_path}")
    if not os.path.isdir(folder_path):
        print(f"     [错误] 目录不存在：{folder_path}")
        return None, None, -1.0, ""

    files = os.listdir(folder_path)
    for i in range(5):
        summary_name = None
        for filename in files:
            if filename.endswith(f"summary_confidences_{i}.json"):
                summary_name = filename
                break
        if not summary_name:
            continue

        try:
            with open(os.path.join(folder_path, summary_name), "r", encoding="utf-8") as f:
                data = json.load(f)
            iptm = float(data.get("iptm", 0.0))
            if iptm > best_iptm:
                best_iptm = iptm
                best_summary = data
                best_idx = i
                best_prefix = summary_name[: -len(f"summary_confidences_{i}.json")]
        except Exception as e:
            print(f"     [警告] 读取 {summary_name} 失败：{e}")

    if not best_summary:
        print("     [错误] 未找到可用 summary_confidences 文件。")
        return None, None, -1.0, ""

    plain_cif = os.path.join(folder_path, f"model_{best_idx}.cif")
    prefixed_cif = os.path.join(folder_path, f"{best_prefix}model_{best_idx}.cif")
    if os.path.isfile(plain_cif):
        best_cif = plain_cif
    elif os.path.isfile(prefixed_cif):
        best_cif = prefixed_cif
    else:
        print(f"     [错误] 找到最佳 iptm 但缺少 model_{best_idx}.cif")
        return best_summary, None, best_iptm, best_prefix

    print(f"   - 完成：最佳模型 {os.path.basename(best_cif)}（iptm={best_iptm:.4f}）")
    return best_summary, best_cif, best_iptm, best_prefix


def _prepare_complex_for_analysis(cif_path: str, complex_obj: str = "my_complex") -> None:
    cmd.reinitialize()
    cmd.load(os.path.abspath(cif_path), complex_obj)
    cmd.alter("all", "segi=''")

    # 统一 SASA 参数
    cmd.set("dot_solvent", 1)
    cmd.set("solvent_radius", 1.4)
    cmd.set("surface_cavity_mode", 0)
    cmd.set("dot_density", 8)
    cmd.set("surface_quality", 2)
    cmd.set("dot_hydrogens", 0)


# ======================================================================
# Method selection and execution
# ======================================================================
def _parse_method_ids_from_env() -> Optional[List[str]]:
    if not METHODS_ENV:
        return None
    ids = [x.strip() for x in METHODS_ENV.split(",") if x.strip()]
    return ids if ids else None


def _calc_overall_status(method_payload: Dict[str, Dict[str, Any]]) -> str:
    if not method_payload:
        return "Success"
    has_success = any(bool(v.get("success")) for v in method_payload.values())
    return "Success" if has_success else "Failed (all methods)"


def analyze_single_folder(
    folder_name: str,
    folder_path: str,
    species_record: Optional[Dict[str, Any]],
    methods: List[Any],
) -> Dict[str, Any]:
    summary_data, cif_path, iptm_score, _ = find_best_model_and_cif(folder_path)
    if not cif_path:
        return {
            "folder_name": folder_name,
            "status": "Failed (no model)",
            "interface_analysis": {},
            "method_results": {},
        }

    _prepare_complex_for_analysis(cif_path, complex_obj="my_complex")
    context = AnalysisContext(
        folder_name=folder_name,
        folder_path=folder_path,
        species_record=species_record,
        complex_object="my_complex",
        pymol_cmd=cmd,
        chain_map=CHAIN_MAP,
    )

    method_results: Dict[str, Dict[str, Any]] = {}
    merged_interface_analysis: Dict[str, Any] = {}
    merged_local_params: Dict[str, Any] = {}

    for method in methods:
        print(f"   - 3. 执行方法：{method.method_id}（{method.display_name}）")
        try:
            run_result = method.run(context)
        except Exception as e:
            run_result = MethodExecutionResult(success=False, message=f"方法执行异常：{e}")
            print(f"     [错误] 方法 {method.method_id} 执行失败：{e}")

        method_results[method.method_id] = run_result.to_dict(
            method_id=method.method_id,
            display_name=method.display_name,
            version=method.version,
        )

        if run_result.legacy_interface_analysis:
            merged_interface_analysis.update(run_result.legacy_interface_analysis)
        if run_result.legacy_local_params and not merged_local_params:
            merged_local_params = run_result.legacy_local_params

    result = {
        "folder_name": folder_name,
        "best_model_file": os.path.basename(cif_path),
        "total_iptm": iptm_score,
        "metrics_from_json": {
            "chain_pair_iptm": summary_data.get("chain_pair_iptm") if summary_data else None,
            "chain_pair_pae_min": summary_data.get("chain_pair_pae_min") if summary_data else None,
        },
        "interface_analysis": merged_interface_analysis,  # 兼容旧脚本
        "method_results": method_results,  # 新增：可扩展方法输出
        "analysis_method_versions": {m.method_id: m.version for m in methods},
        "status": _calc_overall_status(method_results),
    }
    if merged_local_params:
        result["local_params"] = merged_local_params  # 兼容旧脚本
    return result


# ======================================================================
# Main
# ======================================================================
def main_batch_analysis() -> None:
    if pymol is None or cmd is None:
        print("[错误] 当前环境未安装 pymol，无法执行结构分析。")
        print("请先安装可用的 PyMOL Python 环境后重试。")
        return

    all_results = _load_existing_results(OUTPUT_JSON)

    enabled_ids = _parse_method_ids_from_env()
    methods = build_methods(enabled_ids)
    if not methods:
        print("[错误] 没有可用分析方法，程序终止。")
        return

    print("=== 开始批处理特征提取（可扩展方法架构）===")
    print(f"可用方法：{', '.join(list_available_methods())}")
    print(f"启用方法：{', '.join([m.method_id for m in methods])}")
    print(f"强制重跑：{'是' if FORCE_REANALYZE else '否'}")

    records = _load_all_data_records(ALL_DATA_JSON)
    record_index = _build_species_record_index(records)
    print(f"已加载 metadata 记录：{len(record_index)}")

    try:
        print("正在以无界面模式初始化 PyMOL（-c -q）...")
        pymol.finish_launching(["pymol", "-c", "-q"])

        target_folders = []
        for name in os.listdir(BASE_DIR):
            full_path = os.path.join(BASE_DIR, name)
            if os.path.isdir(full_path) and name not in EXCLUDE_DIRS:
                target_folders.append(name)
        target_folders.sort()

        if FORCE_REANALYZE:
            pending = target_folders
        else:
            pending = [x for x in target_folders if x not in all_results]

        print(
            f"共发现 {len(target_folders)} 个目标目录，"
            f"已完成 {len(all_results)} 个，待处理 {len(pending)} 个。"
        )

        for i, folder_name in enumerate(pending, 1):
            print(f"\n--- [ {i} / {len(pending)} ] 正在处理：{folder_name} ---")
            folder_path = os.path.join(BASE_DIR, folder_name)
            species_record = _resolve_species_record(record_index, folder_name)

            folder_result = analyze_single_folder(
                folder_name=folder_name,
                folder_path=folder_path,
                species_record=species_record,
                methods=methods,
            )
            all_results[folder_name] = folder_result
            print(f"   -> 结果：{folder_result.get('status')}")

            try:
                _save_results_atomic(all_results, OUTPUT_JSON)
                print(f"   -> 已增量保存：{OUTPUT_JSON}")
            except Exception as e:
                print(f"   -> 增量保存失败：{e}")

    except Exception as e:
        print(f"\n!!! 发生严重错误：{e} !!!")
        print("批处理提前终止。")
    finally:
        print("\n正在关闭 PyMOL...")
        try:
            cmd.quit()
        except Exception:
            pass
        print("PyMOL 已关闭。")

    print(f"\n正在进行最终保存：{OUTPUT_JSON} ...")
    try:
        _save_results_atomic(all_results, OUTPUT_JSON)
        print("结果保存成功。")
    except Exception as e:
        print(f"保存 JSON 失败：{e}")

    print("\n=== 批处理特征提取完成 ===")


if __name__ == "__main__":
    main_batch_analysis()
