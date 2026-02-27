# -*- coding: utf-8 -*-

from collections import defaultdict
from typing import Any, Dict, List, Tuple

from .base import AnalysisContext, MethodExecutionResult


class LocalInterfaceStrictMethod:
    """
    现有方法迁移版：
    - 3 对界面局部 BSA + H-Bond
    - 严格锚点完整性
    """

    method_id = "local_interface_strict"
    display_name = "局部界面严格锚点分析"
    version = "1.0.0"

    TARGET_PAIRS: List[Tuple[str, str]] = [
        ("C", "D"),  # PEN2-PS1: PEN2 94/95
        ("D", "A"),  # PS1-APH1: PS1 465/466/467
        ("C", "B"),  # PEN2-NCT: PEN2 100
    ]
    REQUIRED_LEN: Dict[Tuple[str, str], int] = {
        ("C", "D"): 2,
        ("D", "A"): 3,
        ("C", "B"): 1,
    }
    LOCAL_RADIUS = 5.0
    HBOND_CUTOFF = 3.2
    HBOND_ANGLE = 140

    @staticmethod
    def _safe_idx(node: Any) -> int:
        if not isinstance(node, dict):
            return None
        aa = node.get("aa")
        if aa in ("-", None):
            return None
        try:
            idx = int(node.get("ungapped_index"))
            if idx <= 0:
                return None
            return idx
        except Exception:
            return None

    def _pair_name(self, chain_map: Dict[str, str], chain_a: str, chain_b: str) -> str:
        return f"{chain_map[chain_a]}-{chain_map[chain_b]}"

    def _load_anchor_sets(self, species_record: Dict[str, Any], target_name: str) -> Dict[Tuple[str, str], List[int]]:
        anchors = {("C", "D"): [], ("D", "A"): [], ("C", "B"): []}
        if not species_record:
            print(f"     [警告] 缺少 '{target_name}' 的 metadata 记录；本方法将全部跳过。")
            return anchors

        # C-D: PEN2 94/95
        try:
            site = species_record.get("pen2_alignment", {}).get("site_94_95", {})
            c94 = self._safe_idx(site.get("pos94"))
            c95 = self._safe_idx(site.get("pos95"))
            vals = [v for v in (c94, c95) if v is not None]
            if len(vals) == self.REQUIRED_LEN[("C", "D")]:
                anchors[("C", "D")] = vals
            else:
                print(f"     [警告] 严格模式：'{target_name}' 的 PEN2 94/95 不完整（{vals}）；C-D 跳过。")
        except Exception:
            print(f"     [警告] '{target_name}' 的 PEN2 site_94_95 解析失败（C-D 跳过）。")

        # D-A: PS1 465/466/467
        try:
            ps1a = species_record.get("ps1_alignment", {})
            p465 = self._safe_idx(ps1a.get("pos465"))
            p466 = self._safe_idx(ps1a.get("pos466"))
            p467 = self._safe_idx(ps1a.get("pos467"))
            vals = [v for v in (p465, p466, p467) if v is not None]
            if len(vals) == self.REQUIRED_LEN[("D", "A")]:
                anchors[("D", "A")] = vals
            else:
                print(f"     [警告] 严格模式：'{target_name}' 的 PS1 465/466/467 不完整（{vals}）；D-A 跳过。")
        except Exception:
            print(f"     [警告] '{target_name}' 的 PS1 465/466/467 解析失败（D-A 跳过）。")

        # C-B: PEN2 100
        try:
            p100 = species_record.get("pen2_alignment", {}).get("site_100", {}).get("pos100", {})
            i100 = self._safe_idx(p100)
            if i100 is not None:
                anchors[("C", "B")] = [i100]
            else:
                print(f"     [警告] 严格模式：'{target_name}' 的 PEN2 100 缺失或无效；C-B 跳过。")
        except Exception:
            print(f"     [警告] '{target_name}' 的 PEN2 site_100 解析失败（C-B 跳过）。")

        print(f"     [信息] '{target_name}' 的严格锚点：{anchors}")
        return anchors

    @staticmethod
    def _sum_residue_sasa(cmd_api: Any, obj: str, sele: str) -> Dict[Tuple[str, str], float]:
        res_sasa = defaultdict(float)
        cmd_api.iterate(f"({obj}) and ({sele})", "res_sasa[(chain, resi)] += b", space={"res_sasa": res_sasa})
        return res_sasa

    def _per_residue_delta_sasa(self, cmd_api: Any, iso_obj: str, complex_obj: str, sele: str) -> float:
        cmd_api.get_area(iso_obj, 1, 1)
        iso_res = self._sum_residue_sasa(cmd_api, iso_obj, sele)
        cmd_api.get_area(complex_obj, 1, 1)
        cmp_res = self._sum_residue_sasa(cmd_api, complex_obj, sele)
        total = 0.0
        for key in set(list(iso_res.keys()) + list(cmp_res.keys())):
            delta = iso_res.get(key, 0.0) - cmp_res.get(key, 0.0)
            if delta > 0:
                total += delta
        return total

    def _calculate_local_interface(
        self,
        cmd_api: Any,
        complex_obj: str,
        anchor_chain: str,
        partner_chain: str,
        anchor_resi: List[int],
    ) -> Tuple[float, int]:
        if not anchor_resi:
            return None, None

        resistr = "+".join(map(str, anchor_resi))
        anchors_sel = f"{complex_obj} and chain {anchor_chain} and resi {resistr}"
        pocket_sel = (
            f"{complex_obj} and chain {partner_chain} within {self.LOCAL_RADIUS} of "
            f"({complex_obj} and chain {anchor_chain} and resi {resistr})"
        )

        obj_anchor_iso = "__iso_anchor_tmp"
        obj_partner_iso = "__iso_partner_tmp"
        for tmp in (obj_anchor_iso, obj_partner_iso):
            if tmp in cmd_api.get_names("all"):
                cmd_api.delete(tmp)

        cmd_api.create(obj_anchor_iso, f"{complex_obj} and chain {anchor_chain}")
        cmd_api.create(obj_partner_iso, f"{complex_obj} and chain {partner_chain}")

        anchors_on_a = f"{obj_anchor_iso} and chain {anchor_chain} and resi {resistr}"
        pocket_on_p = f"{obj_partner_iso} and byres ({pocket_sel})"

        total_a = self._per_residue_delta_sasa(cmd_api, obj_anchor_iso, complex_obj, anchors_on_a)
        total_p = self._per_residue_delta_sasa(cmd_api, obj_partner_iso, complex_obj, pocket_on_p)
        local_bsa_two = total_a + total_p

        cmd_api.select("__anchors", anchors_sel)
        cmd_api.select("__pocket", pocket_sel)
        cmd_api.select("__don_A", "__anchors and donor")
        cmd_api.select("__acc_P", "__pocket and acceptor")
        cmd_api.select("__don_P", "__pocket and donor")
        cmd_api.select("__acc_A", "__anchors and acceptor")

        if cmd_api.count_atoms(f"{complex_obj} and elem H") > 0:
            cmd_api.select("__don_A", "__don_A and neighbor elem H")
            cmd_api.select("__don_P", "__don_P and neighbor elem H")

        pairs_ap = cmd_api.find_pairs("__don_A", "__acc_P", mode=1, cutoff=self.HBOND_CUTOFF, angle=self.HBOND_ANGLE)
        pairs_pa = cmd_api.find_pairs("__don_P", "__acc_A", mode=1, cutoff=self.HBOND_CUTOFF, angle=self.HBOND_ANGLE)

        seen = set()
        hbond_count = 0
        for pset in (pairs_ap, pairs_pa):
            for pair in pset:
                atom_i, atom_j = pair[0][1], pair[1][1]
                dist = cmd_api.get_distance(f"{complex_obj} and index {atom_i}", f"{complex_obj} and index {atom_j}")
                if dist is None or dist > (self.HBOND_CUTOFF + 0.1):
                    continue
                key = tuple(sorted((atom_i, atom_j)))
                if key in seen:
                    continue
                seen.add(key)
                hbond_count += 1

        for sel in ("__anchors", "__pocket", "__don_A", "__acc_P", "__don_P", "__acc_A"):
            if sel in cmd_api.get_names("selections"):
                cmd_api.delete(sel)
        for tmp in (obj_anchor_iso, obj_partner_iso):
            if tmp in cmd_api.get_names("all"):
                cmd_api.delete(tmp)

        return local_bsa_two, hbond_count

    def run(self, context: AnalysisContext) -> MethodExecutionResult:
        cmd_api = context.pymol_cmd
        chain_map = context.chain_map
        chain_ids = sorted(chain_map.keys())
        anchors = self._load_anchor_sets(context.species_record, context.folder_name)

        print("   - 2. 执行方法：局部界面严格锚点分析")
        results: Dict[str, Dict[str, Any]] = {}

        for i, chain_a in enumerate(chain_ids):
            for chain_b in chain_ids[i + 1 :]:
                results[self._pair_name(chain_map, chain_a, chain_b)] = {
                    "BSA": None,
                    "H_Bonds": None,
                    "Method": "未分析",
                }

        for anchor_chain, partner_chain in self.TARGET_PAIRS:
            label = self._pair_name(chain_map, *sorted((anchor_chain, partner_chain)))
            anchor_resi = anchors.get((anchor_chain, partner_chain), [])
            required = self.REQUIRED_LEN[(anchor_chain, partner_chain)]
            if len(anchor_resi) != required:
                print(
                    f"     {label}: [跳过] 严格模式需要 {required} 个锚点，当前仅 {len(anchor_resi)} 个（{anchor_resi}）。"
                )
                results[label] = {
                    "BSA": None,
                    "Interface_one_sided": None,
                    "H_Bonds": None,
                    "Method": f"局部严格模式（锚点在 {chain_map[anchor_chain]}_{anchor_chain}）",
                    "Note": f"已跳过：需要 {required} 个锚点，但实际为 {len(anchor_resi)}",
                }
                continue

            bsa_two, h_bonds = self._calculate_local_interface(
                cmd_api,
                context.complex_object,
                anchor_chain=anchor_chain,
                partner_chain=partner_chain,
                anchor_resi=anchor_resi,
            )
            if bsa_two is None:
                print(f"     {label}: [跳过] 计算失败。")
                results[label] = {
                    "BSA": None,
                    "Interface_one_sided": None,
                    "H_Bonds": None,
                    "Method": f"局部严格模式（锚点在 {chain_map[anchor_chain]}_{anchor_chain}）",
                    "Note": "计算失败",
                }
                continue

            print(
                f"     {label}: 锚点链 {anchor_chain} = {tuple(anchor_resi)} | "
                f"局部 BSA（双侧）={bsa_two:.2f} Å^2，界面面积（单侧）={bsa_two/2.0:.2f} Å^2，氢键数={h_bonds}"
            )
            results[label] = {
                "BSA": bsa_two,
                "Interface_one_sided": bsa_two / 2.0,
                "H_Bonds": h_bonds,
                "Method": (
                    f"局部严格模式（锚点在 {chain_map[anchor_chain]}_{anchor_chain}，"
                    f"对比 {chain_map[partner_chain]}_{partner_chain}）"
                ),
                "Anchors_used": list(anchor_resi),
            }

        local_params = {
            "radius_r_A": self.LOCAL_RADIUS,
            "hbond_cutoff_A": self.HBOND_CUTOFF,
            "hbond_angle_deg": self.HBOND_ANGLE,
            "strict_required_len": {f"{a}-{b}": n for (a, b), n in self.REQUIRED_LEN.items()},
            "anchors_used": {
                "PEN2-PS1 (C-D)": list(anchors.get(("C", "D"), ())),
                "APH1-PS1 (A-D)": list(anchors.get(("D", "A"), ())),
                "NCT-PEN2 (B-C)": list(anchors.get(("C", "B"), ())),
            },
        }

        return MethodExecutionResult(
            success=True,
            message="方法执行完成",
            metrics={
                "interface_analysis": results,
                "local_params": local_params,
            },
            legacy_interface_analysis=results,
            legacy_local_params=local_params,
        )
