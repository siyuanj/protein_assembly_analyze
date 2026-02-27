# -*- coding: utf-8 -*-

import json
import os
import re

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

MAPPING_FILE_PATH = os.path.join(PROJECT_ROOT, 'data', 'metadata', 'job_name_mapping.json')
RESULTS_PARENT_FOLDER = os.path.join(PROJECT_ROOT, 'data', 'raw_structures', 'af_samllsample_1021')


def sanitize_af_name_final(job_name_from_mapping: str) -> str:
    """Apply AF-style normalization to folder names."""
    sanitized_name = job_name_from_mapping.lower()
    sanitized_name = re.sub(r'[^a-z0-9]+', '_', sanitized_name)
    return sanitized_name.strip('_')


def verify_downloaded_results() -> None:
    print('--- 开始核对 AlphaFold 结果目录 ---')

    try:
        with open(MAPPING_FILE_PATH, 'r', encoding='utf-8') as f:
            name_mapping = json.load(f)
    except FileNotFoundError:
        print(f"错误：找不到映射文件：{MAPPING_FILE_PATH}")
        return
    except Exception as e:
        print(f"错误：加载映射文件失败：{e}")
        return

    if not isinstance(name_mapping, dict) or not name_mapping:
        print(f"错误：映射文件为空或格式无效：{MAPPING_FILE_PATH}")
        return

    print(f"已加载 {len(name_mapping)} 条映射记录。")

    expected_folders_set = set()
    reverse_name_map = {}
    for original_name, job_name in name_mapping.items():
        final_name = sanitize_af_name_final(str(job_name))
        expected_folders_set.add(final_name)
        reverse_name_map[final_name] = original_name

    try:
        found_folders_set = {
            entry
            for entry in os.listdir(RESULTS_PARENT_FOLDER)
            if os.path.isdir(os.path.join(RESULTS_PARENT_FOLDER, entry))
        }
    except Exception as e:
        print(f"错误：扫描结果目录失败：{e}")
        return

    print(f"在 {RESULTS_PARENT_FOLDER} 下发现 {len(found_folders_set)} 个目录。")

    missing_folders = expected_folders_set - found_folders_set
    extra_folders = found_folders_set - expected_folders_set

    print('\n--- 核对报告 ---')

    if not missing_folders and not extra_folders:
        print('通过：目录集合与映射文件一致。')

    if missing_folders:
        print(f"缺失目录（{len(missing_folders)}）：")
        for i, san_name in enumerate(sorted(missing_folders), 1):
            original_file = reverse_name_map.get(san_name, '[未知来源]')
            print(f"  {i}. {san_name}（来源：{original_file}）")

    if extra_folders:
        print(f"多余目录（{len(extra_folders)}）：")
        for i, folder_name in enumerate(sorted(extra_folders), 1):
            print(f"  {i}. {folder_name}")

    print('--- 核对结束 ---')


if __name__ == '__main__':
    verify_downloaded_results()
