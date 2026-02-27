# gamma_sec_evolution 项目说明

## 目录结构

```text
.
├── data/
│   ├── raw_structures/af_samllsample_1021/   # AlphaFold 原始结构目录（只读）
│   └── metadata/                              # 元数据
├── scripts/                                   # 执行脚本（按流程编号）
│   └── analysis_methods/                      # 可扩展分析方法插件
├── results/                                   # 结果输出
│   ├── parsed_data/
│   └── figures/
└── logs/                                      # 运行日志
```

## 运行顺序（01 -> 02 -> 03）

请在项目根目录执行（即本 README 所在目录）。

### 1. 数据完整性核对

```powershell
python scripts/01_verify_data.py
```

写入日志：

```powershell
python scripts/01_verify_data.py *> logs/01_verify.log
```

### 2. 结构特征提取（批处理）

```powershell
python scripts/02_extract_features.py
```

写入日志（推荐长任务使用）：

```powershell
python scripts/02_extract_features.py *> logs/02_extract_run.log
```

仅记录错误：

```powershell
python scripts/02_extract_features.py 2> logs/02_extract_error.log
```

输出结果文件：

- `results/parsed_data/batch_analysis_results.json`

### 3. 趋势可视化

```powershell
python scripts/03_visualize_trends.py
```

输出图表目录：

- `results/figures/combined_lineplots`
- `results/figures/combined_boxplots`
- `results/figures/individual_lineplots`
- `results/figures/individual_boxplots`

## 可扩展分析方法（重点）

`scripts/02_extract_features.py` 已改为“调度器 + 插件”架构：

- 调度器：负责目录遍历、最佳模型选择、结果保存。
- 插件方法：放在 `scripts/analysis_methods/`，每个方法一个文件。
- 注册入口：`scripts/analysis_methods/registry.py`

当前内置方法：

- `local_interface_strict`（局部界面严格锚点分析）

### 新增方法步骤

1. 在 `scripts/analysis_methods/` 新建一个方法文件，例如 `global_geometry.py`。  
2. 参考 `local_interface_strict.py`，实现：
   - `method_id`
   - `display_name`
   - `version`
   - `run(context) -> MethodExecutionResult`
3. 在 `scripts/analysis_methods/registry.py` 的 `_METHOD_REGISTRY` 中注册新方法。  
4. 重新运行 `python scripts/02_extract_features.py`。

### 方法启用控制（环境变量）

- 默认不设置：运行全部已注册方法。
- 仅运行指定方法：

```powershell
$env:ANALYSIS_METHODS="local_interface_strict"
python scripts/02_extract_features.py
```

- 强制重跑（忽略已有结果）：

```powershell
$env:FORCE_REANALYZE="1"
python scripts/02_extract_features.py
```

## 环境依赖

- Python 3.9+
- `pymol`（用于 `02_extract_features.py`）
- `pandas`、`matplotlib`、`seaborn`（用于 `03_visualize_trends.py`）

## 说明

- `data/raw_structures/` 视为原始数据层，建议只读，不直接修改。
- 历史图表批次保存在 `results/figures/archive_runs/`。
