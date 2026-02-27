import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from matplotlib.ticker import FixedLocator

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# ==============================================================================
# --- 閰嶇疆鍖?(CONFIGURATION AREA) ---
# ==============================================================================

# 1. Input JSON file path
JSON_FILE_PATH = os.path.join(PROJECT_ROOT, 'results', 'parsed_data', 'batch_analysis_results.json')

# 2. Chart output settings
FIGURES_ROOT = os.path.join(PROJECT_ROOT, 'results', 'figures')
COMBINED_LINE_DIR = os.path.join(FIGURES_ROOT, 'combined_lineplots')
COMBINED_BOX_DIR = os.path.join(FIGURES_ROOT, 'combined_boxplots')
INDIVIDUAL_BOX_DIR = os.path.join(FIGURES_ROOT, 'individual_boxplots')
INDIVIDUAL_LINE_DIR = os.path.join(FIGURES_ROOT, 'individual_lineplots')
FIGURE_FORMAT = 'png'   # 'png' | 'svg' | 'pdf' ...
FIGURE_DPI = 300

# 3. Species classification information
CLASSIFICATION_MAP = {
    '3.1':      'Fungi: Chytridiomycota',
    '3.4':      'Fungi: Mucoromycota',
   # '4.3':      'Protists: SAR Supergroup',
   # '5.6':      'Plants: Angiosperms',
    '6.2':      'Animalia: Cnidaria',
    '6.4':      'Animalia: Nematoda',
    '6.6':      'Animalia: Mollusca',
    '6.7':      'Animalia: Arthropoda',
    '6.8':      'Animalia: Echinodermata',
    '6.10.3.3': 'Animalia: Osteichthyes',
    '6.10.3.4': 'Animalia: Amphibia',
    '6.10.3.5': 'Animalia: Reptilia',
    '6.10.3.6': 'Animalia: Aves',
    '6.10.3.7': 'Animalia: Mammalia'
}

# 4. Key protein interfaces to analyze
KEY_INTERFACES = [
    "APH1-PS1",
    "APH1-NCT",
    "PEN2-PS1",
    "NCT-PEN2"
]

# 5. Chain order from the AlphaFold prediction
CHAIN_ORDER = {'APH1': 0, 'NCT': 1, 'PEN2': 2, 'PS1': 3}

COLOR_PALETTE = 'tab10'

# ==============================================================================
# --- Data Processing and Visualization Functions (DO NOT EDIT BELOW) ---
# ==============================================================================

def set_matplotlib_font():
    """Set a default font to avoid CJK issues."""
    try:
        plt.rcParams['font.family'] = 'Arial'
    except Exception as e:
        print(f"无法将字体设置为 Arial，将使用默认字体。错误：{e}")

def get_classification_key(folder_name: str) -> str:
    """Parse the classification key from the folder name like '6_10_3_7_species' -> '6.10.3.7'."""
    parts = folder_name.split('_')
    key_parts = []
    for part in parts:
        if part.isdigit():
            key_parts.append(part)
        else:
            break
    return '.'.join(key_parts)

def is_success_status(status_value) -> bool:
    status = str(status_value).strip().lower()
    return status == "success" or status.startswith("success") or status == "成功"

def load_and_parse_data(json_path, classification_map, key_interfaces, chain_order):
    """Load JSON and convert to a tidy DataFrame."""
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"错误：未找到文件 '{json_path}'。")
        return None
    except json.JSONDecodeError:
        print(f"错误：文件 '{json_path}' 不是有效的 JSON。")
        return None

    processed = []
    classification_order = list(classification_map.values())

    for folder_name, details in data.items():
        if not isinstance(details, dict) or not is_success_status(details.get("status")):
            continue

        class_key = get_classification_key(folder_name)
        if class_key not in classification_map:
            continue

        classification = classification_map[class_key]
        record = {
            'folder_name': folder_name,
            'classification': classification,
            'total_iptm': details.get('total_iptm')  # may be None
        }

        pae_matrix = details.get('metrics_from_json', {}).get('chain_pair_pae_min', [])
        interface_analysis = details.get('interface_analysis', {}) or {}
        if not interface_analysis:
            for method_payload in (details.get("method_results", {}) or {}).values():
                metrics = method_payload.get("metrics", {}) if isinstance(method_payload, dict) else {}
                if isinstance(metrics.get("interface_analysis"), dict):
                    interface_analysis.update(metrics["interface_analysis"])

        for interface in key_interfaces:
            chain1, chain2 = interface.split('-')
            # PAE
            try:
                idx1, idx2 = chain_order[chain1], chain_order[chain2]
                record[f'PAE_{interface}'] = pae_matrix[idx1][idx2]
            except (IndexError, KeyError, TypeError):
                record[f'PAE_{interface}'] = None

            # Interface metrics (BSA / H_Bonds)
            ida = interface_analysis.get(interface) or interface_analysis.get(f"{chain2}-{chain1}")
            if isinstance(ida, dict):
                record[f'BSA_{interface}'] = ida.get('BSA')
                record[f'H_Bonds_{interface}'] = ida.get('H_Bonds')
            else:
                record[f'BSA_{interface}'] = None
                record[f'H_Bonds_{interface}'] = None

        processed.append(record)

    if not processed:
        print("警告：未从 JSON 中解析到有效数据。")
        return pd.DataFrame()

    df = pd.DataFrame(processed)
    df['classification'] = pd.Categorical(df['classification'], categories=classification_order, ordered=True)
    return df.sort_values('classification').reset_index(drop=True)

# ---------- helpers: annotate sample counts on x-axis ----------
def _apply_xtick_counts(ax, counts_by_class: pd.Series):
    """
    Append sample count to x tick labels: "Class\n(n=XX)".
    Use FixedLocator to avoid set_ticklabels warnings.
    """
    try:
        labels = [tick.get_text() for tick in ax.get_xticklabels()]
        locs = ax.get_xticks()
        new_labels = []
        for lab in labels:
            n = int(counts_by_class.get(lab, 0))
            new_labels.append(f"{lab}\n(n={n})")
        ax.xaxis.set_major_locator(FixedLocator(locs))
        ax.set_xticklabels(new_labels)
    except Exception:
        # best-effort; ignore if fails
        pass

def plot_combined_visuals(df, metric_prefix, title_suffix, y_label, combined_line_dir, combined_box_dir, file_format, dpi, y_limit=None):
    """Generate plots comparing all key interfaces together."""
    metric_cols = [c for c in df.columns if c.startswith(metric_prefix)]
    if not metric_cols:
        print(f"跳过 {title_suffix} 的组合图：不存在以 '{metric_prefix}' 开头的列。")
        return

    df_melted = df.melt(
        id_vars=['classification'],
        value_vars=metric_cols,
        var_name='interface',
        value_name='value'
    )
    # 灏嗗€艰浆涓烘暟鍊煎苟涓㈠純 NaN
    df_melted['value'] = pd.to_numeric(df_melted['value'], errors='coerce')
    df_melted = df_melted.dropna(subset=['value'])
    if df_melted.empty:
        print(f"跳过 {title_suffix} 的组合图：所有值均为空。")
        return

    df_melted['interface'] = df_melted['interface'].str.replace(f'{metric_prefix}_', '', regex=False)

    plt.style.use('seaborn-v0_8-whitegrid')

    # Box Plot
    plt.figure(figsize=(20, 10))
    ax = sns.boxplot(data=df_melted, x='classification', y='value', hue='interface', palette=COLOR_PALETTE)
    plt.title(f'{title_suffix} Distribution (Combined)', fontsize=20, pad=20)
    plt.xlabel('Species Classification', fontsize=14)
    plt.ylabel(y_label, fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.legend(title='Protein Interface', bbox_to_anchor=(1.02, 1), loc='upper left')
    if y_limit is not None:
        plt.ylim(0, y_limit)
    # 鏍囨敞鏍锋湰鏁扮洰锛堟寜褰撳墠鍥句腑瀹為檯浣跨敤鐨勬暟鎹偣璁℃暟锛?
    counts = df_melted.groupby('classification', observed=True)['value'].count()
    _apply_xtick_counts(ax, counts)
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    plt.savefig(os.path.join(combined_box_dir, f'combined_boxplot_{metric_prefix.lower()}.{file_format}'), dpi=dpi)
    plt.close()

    # Line Plot (mean with sd)
    plt.figure(figsize=(20, 10))
    ax = sns.lineplot(data=df_melted, x='classification', y='value', hue='interface', marker='o', errorbar='sd', palette=COLOR_PALETTE)
    plt.title(f'Mean {title_suffix} Trend (Combined)', fontsize=20, pad=20)
    plt.xlabel('Species Classification', fontsize=14)
    plt.ylabel(f'Mean {y_label} (with Std. Dev.)', fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.legend(title='Protein Interface', bbox_to_anchor=(1.02, 1), loc='upper left')
    if y_limit is not None:
        plt.ylim(0, y_limit)
    counts = df_melted.groupby('classification', observed=True)['value'].count()
    _apply_xtick_counts(ax, counts)
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    plt.savefig(os.path.join(combined_line_dir, f'combined_lineplot_{metric_prefix.lower()}.{file_format}'), dpi=dpi)
    plt.close()

    print(f"已生成 {title_suffix} 的组合图。")

    # --- BSA 涓撶敤锛氭寜鈥滄瘡涓晫闈㈢殑鏈€灏忓垎绫诲潎鍊尖€濆綊涓€鍖栵紙鍊煎潎 鈮?1锛夛紝鍐嶆眹鎬诲埌涓€寮犲浘 ---
    if metric_prefix.upper() == 'BSA':
        # 1) 璁＄畻姣忎釜 interface 鍦ㄥ悇鍒嗙被涓婄殑鍧囧€?
        class_means = (
            df_melted
            .groupby(['interface', 'classification'], observed=True)['value']
            .mean()
            .reset_index()
        )
        # 2) 姣忎釜 interface 鐨勨€滄鍧囧€兼渶灏忓€尖€濅綔涓鸿 interface 鐨勫熀绾?
        pos_means = class_means[class_means['value'] > 0]
        if pos_means.empty:
            print("跳过归一化 BSA 图：没有任何界面存在正的分类均值。")
            return
        baseline = (
            pos_means
            .groupby('interface', observed=True)['value']
            .min()
            .rename('baseline')
            .reset_index()
        )
        # 3) 鍚堝苟鍩虹嚎骞惰繘琛屽綊涓€鍖栵紙value / baseline_ifc锛夛紝涓㈠純鏃犲熀绾跨殑鐣岄潰
        df_norm = df_melted.merge(baseline, on='interface', how='left').dropna(subset=['baseline'])
        if df_norm.empty:
            print("跳过归一化 BSA 图：所有界面都缺少正的基线值（按分类均值计算）。")
            return
        df_norm['value'] = df_norm['value'] / df_norm['baseline']  # >= 1

        # Normalized Box Plot锛堟瘡涓晫闈㈢嫭绔嬪綊涓€鍖栧悗鍚堝苟锛?
        plt.figure(figsize=(20, 10))
        ax = sns.boxplot(data=df_norm, x='classification', y='value', hue='interface', palette=COLOR_PALETTE)
        plt.title('Interface BSA Distribution (Combined, normalized by min class mean per interface, 鈮?)', fontsize=20, pad=20)
        plt.xlabel('Species Classification', fontsize=14)
        plt.ylabel('Buried Surface Area (normalized, 鈮?)', fontsize=14)
        plt.xticks(rotation=45, ha='right', fontsize=12)
        plt.legend(title='Protein Interface', bbox_to_anchor=(1.02, 1), loc='upper left')
        counts = df_norm.groupby('classification', observed=True)['value'].count()
        _apply_xtick_counts(ax, counts)
        plt.tight_layout(rect=[0, 0, 0.88, 1])
        plt.savefig(os.path.join(combined_box_dir, f'combined_boxplot_{metric_prefix.lower()}_normalized_by_min_class_mean.{file_format}'), dpi=dpi)
        plt.close()

        # Normalized Line Plot锛堝潎鍊间笌鏍囧噯宸熀浜庡綊涓€鍖栧悗鐨勬牱鏈噸鏂拌绠楋級
        plt.figure(figsize=(20, 10))
        ax = sns.lineplot(data=df_norm, x='classification', y='value', hue='interface', marker='o', errorbar='sd', palette=COLOR_PALETTE)
        plt.title('Mean Interface BSA Trend (Combined, normalized by min class mean per interface, 鈮?)', fontsize=20, pad=20)
        plt.xlabel('Species Classification', fontsize=14)
        plt.ylabel('Mean Buried Surface Area (normalized, 鈮?)', fontsize=14)
        plt.xticks(rotation=45, ha='right', fontsize=12)
        plt.legend(title='Protein Interface', bbox_to_anchor=(1.02, 1), loc='upper left')
        counts = df_norm.groupby('classification', observed=True)['value'].count()
        _apply_xtick_counts(ax, counts)
        plt.tight_layout(rect=[0, 0, 0.88, 1])
        plt.savefig(os.path.join(combined_line_dir, f'combined_lineplot_{metric_prefix.lower()}_normalized_by_min_class_mean.{file_format}'), dpi=dpi)
        plt.close()

        print("已生成按“每个界面最小分类均值”归一化的组合 BSA 图。")

def plot_individual_visuals(df, metric_prefix, title_suffix, y_label, individual_line_dir, individual_box_dir, file_format, dpi, y_limit=None):
    """Generate a separate plot for each interface column."""
    metric_cols = [c for c in df.columns if c.startswith(metric_prefix)]
    if not metric_cols:
        print(f"跳过 {title_suffix} 的单界面图：不存在以 '{metric_prefix}' 开头的列。")
        return

    for col in metric_cols:
        # 鏇寸ǔ濡ュ啓娉曪細浠呭綋鍓嶇紑瀛樺湪鏃跺幓鎺夊畠
        if col.startswith(f'{metric_prefix}_'):
            interface_name = col[len(metric_prefix) + 1:]
        else:
            interface_name = col
        # 鍙敤璇ュ垪鏈夊€肩殑鏁版嵁
        df_col = df[['classification', col]].copy()
        df_col[col] = pd.to_numeric(df_col[col], errors='coerce')
        df_col = df_col.dropna(subset=[col])
        if df_col.empty:
            print(f"跳过 {interface_name}：{metric_prefix} 全部为空值。")
            continue
        # 鍘绘帀鏈娇鐢ㄧ殑鍒嗙被姘村钩
        if hasattr(df_col['classification'], 'cat'):
            df_col['classification'] = df_col['classification'].cat.remove_unused_categories()

        plt.style.use('seaborn-v0_8-whitegrid')

        # Box Plot
        plt.figure(figsize=(18, 9))
        ax = sns.boxplot(data=df_col, x='classification', y=col, color='skyblue')
        plt.title(f'{interface_name}: {title_suffix} Distribution', fontsize=20, pad=20)
        plt.xlabel('Species Classification', fontsize=14)
        plt.ylabel(y_label, fontsize=14)
        plt.xticks(rotation=45, ha='right', fontsize=12)
        if y_limit is not None:
            plt.ylim(0, y_limit)
        # 鏍锋湰鏁帮細鎸夊綋鍓嶅垪鍦ㄥ悇鍒嗙被涓殑闈炵┖璁℃暟
        counts = df_col.groupby('classification', observed=True)[col].count()
        _apply_xtick_counts(ax, counts)
        plt.tight_layout()
        plt.savefig(os.path.join(individual_box_dir, f'individual_boxplot_{metric_prefix.lower()}_{interface_name}.{file_format}'), dpi=dpi)
        plt.close()

        # Line Plot
        plt.figure(figsize=(18, 9))
        ax = sns.lineplot(data=df_col, x='classification', y=col, marker='o', color='royalblue', errorbar='sd')
        plt.title(f'{interface_name}: Mean {title_suffix} Trend', fontsize=20, pad=20)
        plt.xlabel('Species Classification', fontsize=14)
        plt.ylabel(f'Mean {y_label} (with Std. Dev.)', fontsize=14)
        plt.xticks(rotation=45, ha='right', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        if y_limit is not None:
            plt.ylim(0, y_limit)
        counts = df_col.groupby('classification', observed=True)[col].count()
        _apply_xtick_counts(ax, counts)
        plt.tight_layout()
        plt.savefig(os.path.join(individual_line_dir, f'individual_lineplot_{metric_prefix.lower()}_{interface_name}.{file_format}'), dpi=dpi)
        plt.close()

    print(f"已完成 {title_suffix} 的单界面图生成。")

def main():
    """Main entry."""
    set_matplotlib_font()
    for out_dir in (COMBINED_LINE_DIR, COMBINED_BOX_DIR, INDIVIDUAL_BOX_DIR, INDIVIDUAL_LINE_DIR):
        os.makedirs(out_dir, exist_ok=True)

    df = load_and_parse_data(JSON_FILE_PATH, CLASSIFICATION_MAP, KEY_INTERFACES, CHAIN_ORDER)
    if df is None or df.empty:
        print("未加载到可用数据，程序终止。")
        sys.exit(0)

    print(f"成功加载并解析 {len(df)} 条有效记录。")
    print("-" * 40)

    # Combined
    print("正在生成组合图（多界面对比）...")
    plot_combined_visuals(df, 'total_iptm', 'Overall iPTM', 'iPTM Score', COMBINED_LINE_DIR, COMBINED_BOX_DIR, FIGURE_FORMAT, FIGURE_DPI)
    plot_combined_visuals(df, 'PAE', 'Interface PAE', 'PAE Score', COMBINED_LINE_DIR, COMBINED_BOX_DIR, FIGURE_FORMAT, FIGURE_DPI)
    plot_combined_visuals(df, 'BSA', 'Interface BSA', 'Buried Surface Area (脜虏)', COMBINED_LINE_DIR, COMBINED_BOX_DIR, FIGURE_FORMAT, FIGURE_DPI)
    plot_combined_visuals(df, 'H_Bonds', 'Interface H-Bonds', 'Number of H-Bonds', COMBINED_LINE_DIR, COMBINED_BOX_DIR, FIGURE_FORMAT, FIGURE_DPI)

    # Individual
    print("\n正在生成单界面图（每个界面一组图）...")
    plot_individual_visuals(df, 'PAE', 'Interface PAE', 'PAE Score', INDIVIDUAL_LINE_DIR, INDIVIDUAL_BOX_DIR, FIGURE_FORMAT, FIGURE_DPI)
    plot_individual_visuals(df, 'BSA', 'Interface BSA', 'Buried Surface Area (脜虏)', INDIVIDUAL_LINE_DIR, INDIVIDUAL_BOX_DIR, FIGURE_FORMAT, FIGURE_DPI)
    plot_individual_visuals(df, 'H_Bonds', 'Interface H-Bonds', 'Number of H-Bonds', INDIVIDUAL_LINE_DIR, INDIVIDUAL_BOX_DIR, FIGURE_FORMAT, FIGURE_DPI)

    print("-" * 40)
    print(f"所有图表已保存到：{FIGURES_ROOT}")

if __name__ == '__main__':
    main()


