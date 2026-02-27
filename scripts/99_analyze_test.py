# -*- coding: utf-8 -*-
import os
import json
import pymol
from pymol import cmd
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
TARGET_FOLDER_NAME = '6_10_3_7_tupaia_chinensis'
TARGET_FOLDER = os.path.join(PROJECT_ROOT, 'data', 'raw_structures', 'af_samllsample_1021', TARGET_FOLDER_NAME)
ALL_DATA_JSON = os.path.join(PROJECT_ROOT, 'data', 'metadata', 'all_data_1101.json')
LOCAL_RADIUS  = 5.0
HBOND_CUTOFF  = 3.2
HBOND_ANGLE   = 140
CHAIN_NAMES   = ['APH1', 'NCT', 'PEN2', 'PS1']


def load_anchor_from_all_data(all_json_path, target_name):
    fallback = (121, 122)
    if not os.path.isfile(all_json_path):
        print(f"[警告] 未找到 '{all_json_path}'，使用默认锚点 {fallback}")
        return fallback
    try:
        with open(all_json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception as e:
        print(f"[警告] 读取 '{all_json_path}' 失败：{e}；使用默认锚点 {fallback}")
        return fallback

    items = data if isinstance(data, list) else [data]
    for item in items:
        if item.get('name_final') == target_name \
           or item.get('name_cleaned') == target_name \
           or item.get('name_original') == target_name:
            pa = item.get('pen2_alignment', {})
            try:
                i94 = int(pa['pos94']['ungapped_index'])
                i95 = int(pa['pos95']['ungapped_index'])
                print(f"[信息] 从比对数据读取 '{target_name}' 的锚点：({i94}, {i95})")
                return (i94, i95)
            except Exception:
                print(f"[警告] pen2_alignment 不完整，使用默认锚点 {fallback}")
                return fallback
    print(f"[警告] 在比对数据中未找到 '{target_name}'，使用默认锚点 {fallback}")
    return fallback


def print_matrix(title, matrix, names):
    print(f"\n--- {title} ---")
    header = " " * 6 + "".join([f"{n:>10}" for n in names])
    print(header)
    print(" " * 5 + "-" * (len(names) * 10))
    for i, row in enumerate(matrix):
        row_str = "".join([f"{v:10.2f}" for v in row])
        print(f" {names[i]:>5} |{row_str}")


def find_best_model_and_cif(folder_path):
    """
    兼容两种命名：
    1) summary_confidences_0.json / model_0.cif
    2) fold_xxx_summary_confidences_0.json / fold_xxx_model_0.cif
    """
    best_iptm = -1.0
    best_summary = None
    best_idx = -1
    best_prefix = ""   # "" 或 "fold_6_10_3_7_tupaia_chinensis_" 这种

    print(f"--- 1. 扫描 '{folder_path}' 的最佳模型 ---")
    if not os.path.isdir(folder_path):
        print(f"[错误] 目录不存在：{folder_path}")
        return None, None, -1.0

    files = os.listdir(folder_path)

    for i in range(5):
        # 形式1：summary_confidences_i.json
        cand1 = f"summary_confidences_{i}.json"
        # 形式2：*summary_confidences_i.json
        cand2 = None
        for fn in files:
            if fn.endswith(f"summary_confidences_{i}.json"):
                cand2 = fn
                break

        target_json = None
        prefix_used = ""
        if cand1 in files:
            target_json = cand1
            prefix_used = ""
        elif cand2:
            target_json = cand2
            # 从 cand2 里把前缀取出来，比如 fold_6_10_3_7_tupaia_chinensis_
            prefix_used = cand2[:-len(f"summary_confidences_{i}.json")]
        else:
            continue

        json_path = os.path.join(folder_path, target_json)
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            iptm = data.get('iptm', 0.0)
            print(f"  已发现 {target_json}：iptm={iptm:.4f}")
            if iptm > best_iptm:
                best_iptm = iptm
                best_summary = data
                best_idx = i
                best_prefix = prefix_used
        except Exception as e:
            print(f"  [警告] 读取 {target_json} 失败：{e}")

    if best_summary is None:
        print("[错误] 未找到有效的 summary_confidences_*.json（含或不含 fold_ 前缀）")
        return None, None, -1.0

    # 找对应的 cif，也要兼容两种名字
    # 形式1：model_i.cif
    plain_cif = os.path.join(folder_path, f"model_{best_idx}.cif")
    # 形式2：<prefix>model_i.cif
    pref_cif = os.path.join(folder_path, f"{best_prefix}model_{best_idx}.cif")

    if os.path.isfile(plain_cif):
        cif_path = plain_cif
    elif os.path.isfile(pref_cif):
        cif_path = pref_cif
    else:
        print(f"[错误] 未找到最佳模型对应 cif：尝试过 '{plain_cif}' 和 '{pref_cif}'")
        return best_summary, None, best_iptm

    print(f"[完成] 最佳模型：{os.path.basename(cif_path)}（iptm={best_iptm:.4f}）")

    # 打印矩阵
    iptm_mat = best_summary.get('chain_pair_iptm')
    if iptm_mat:
        print_matrix("界面 iPTM（越高越好）", iptm_mat, CHAIN_NAMES)
    else:
        print("[警告] summary 中缺少 chain_pair_iptm")

    pae_min_mat = best_summary.get('chain_pair_pae_min')
    if pae_min_mat:
        print_matrix("界面最小 PAE（Å，越低越好）", pae_min_mat, CHAIN_NAMES)
    else:
        print("[警告] summary 中缺少 chain_pair_pae_min")

    return best_summary, cif_path, best_iptm


def _sum_residue_sasa(obj, sele):
    res_sasa = defaultdict(float)
    cmd.iterate(f"({obj}) and ({sele})", "res_sasa[(chain, resi)] += b",
                space={'res_sasa': res_sasa})
    return res_sasa


def _per_residue_delta_sasa(iso_obj, complex_obj, sele):
    cmd.get_area(iso_obj, 1, 1)
    iso_res = _sum_residue_sasa(iso_obj, sele)
    cmd.get_area(complex_obj, 1, 1)
    cmp_res = _sum_residue_sasa(complex_obj, sele)

    total = 0.0
    keys = set(list(iso_res.keys()) + list(cmp_res.keys()))
    for k in keys:
        d = iso_res.get(k, 0.0) - cmp_res.get(k, 0.0)
        if d > 0:
            total += d
    return total


def calc_local_cd(complex_obj,
                  anchors=(121, 122),
                  r=5.0, cutoff=3.2, angle=140.0):
    anchor_chain = 'C'
    partner_chain = 'D'

    resistr = '+'.join(map(str, anchors))
    anchors_sel = f"{complex_obj} and chain {anchor_chain} and resi {resistr}"
    pocket_sel  = (f"{complex_obj} and chain {partner_chain} within {r} of "
                   f"({complex_obj} and chain {anchor_chain} and resi {resistr})")

    obj_c_iso = "__iso_C_tmp"
    obj_d_iso = "__iso_D_tmp"
    for tmp in (obj_c_iso, obj_d_iso):
        if tmp in cmd.get_names('all'):
            cmd.delete(tmp)
    cmd.create(obj_c_iso, f"{complex_obj} and chain {anchor_chain}")
    cmd.create(obj_d_iso, f"{complex_obj} and chain {partner_chain}")

    anchors_on_c = f"{obj_c_iso} and chain {anchor_chain} and resi {resistr}"
    pocket_on_d  = f"{obj_d_iso} and byres ({pocket_sel})"

    total_c = _per_residue_delta_sasa(obj_c_iso, complex_obj, anchors_on_c)
    total_d = _per_residue_delta_sasa(obj_d_iso, complex_obj, pocket_on_d)
    local_bsa_two = total_c + total_d

    # 氢键（局部 + 严格链间）
    cmd.select("__anchors", anchors_sel)
    cmd.select("__pocketD", pocket_sel)
    cmd.select("__don_A", "__anchors and donor")
    cmd.select("__acc_B", "__pocketD and acceptor")
    cmd.select("__don_B", "__pocketD and donor")
    cmd.select("__acc_A", "__anchors and acceptor")

    if cmd.count_atoms(f"{complex_obj} and elem H") > 0:
        cmd.select("__don_A", "__don_A and neighbor elem H")
        cmd.select("__don_B", "__don_B and neighbor elem H")

    pairs_AB = cmd.find_pairs("__don_A", "__acc_B", mode=1, cutoff=cutoff, angle=angle)
    pairs_BA = cmd.find_pairs("__don_B", "__acc_A", mode=1, cutoff=cutoff, angle=angle)

    seen = set()
    hbond_count = 0
    for pset in (pairs_AB, pairs_BA):
        for p in pset:
            i, j = p[0][1], p[1][1]
            d = cmd.get_distance(f"{complex_obj} and index {i}", f"{complex_obj} and index {j}")
            if d is None or d > (cutoff + 0.1):
                continue
            key = tuple(sorted((i, j)))
            if key in seen:
                continue
            seen.add(key)
            hbond_count += 1

    # 清理
    for s in ("__anchors", "__pocketD", "__don_A", "__acc_B", "__don_B", "__acc_A"):
        if s in cmd.get_names('selections'):
            cmd.delete(s)
    for tmp in (obj_c_iso, obj_d_iso):
        if tmp in cmd.get_names('all'):
            cmd.delete(tmp)

    return local_bsa_two, hbond_count


def main():
    pymol.finish_launching(['pymol', '-c', '-q'])

    anchors = load_anchor_from_all_data(ALL_DATA_JSON, TARGET_FOLDER_NAME)

    summary_data, cif_path, iptm = find_best_model_and_cif(TARGET_FOLDER)
    if not cif_path:
        print("[严重错误] 未找到可用 cif，程序退出。")
        cmd.quit()
        return

    cmd.reinitialize()
    cmd.load(os.path.abspath(cif_path), 'my_complex')
    cmd.alter('all', "segi=''")

    cmd.set('dot_solvent', 1)
    cmd.set('solvent_radius', 1.4)
    cmd.set('surface_cavity_mode', 0)
    cmd.set('dot_density', 8)
    cmd.set('surface_quality', 2)
    cmd.set('dot_hydrogens', 0)

    print("\n--- 2. 计算局部界面 C-D ---")
    local_bsa_two, local_hb = calc_local_cd(
        'my_complex',
        anchors=anchors,
        r=LOCAL_RADIUS,
        cutoff=HBOND_CUTOFF,
        angle=HBOND_ANGLE
    )

    print(f" C 链（PEN2）锚点：{anchors}")
    print(f" 局部 BSA（双侧）= {local_bsa_two:.2f} Å^2")
    print(f" 界面面积（单侧）= {local_bsa_two/2.0:.2f} Å^2（BSA/2）")
    print(f" 局部链间氢键数 = {local_hb}")

    print("\n=== 运行完成 ===")
    cmd.quit()


if __name__ == "__main__":
    main()
